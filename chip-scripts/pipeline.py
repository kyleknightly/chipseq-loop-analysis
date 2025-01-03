"""
This is the culmination of all of my matrix-producing scripts.

Each step is denumerated with one of these comments.

You can start this pipeline at different points by commenting out 
preceeding steps and uncommenting the file importer for the desired step.

Inputs: ChIP-seq tracks (with data), loop list
Output: Protein contact matrices in cis, trans, with various filters.
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import ast
import pickle
import os
import seaborn as sns
from collections import defaultdict

# MARK: Step 1
"""
ChIP-seq tracks + loop list -> paired-anchor-TFs file
"""


# MARK: Step 2
"""
Paired-anchor-TFs file -> proportions
"""
# File importer
paired_anchor_tfs_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/flag-tag-catchseq/beds/paired-anchor-TFs.bed'

paired_anchor_tfs = pd.read_csv(paired_anchor_tfs_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
paired_anchor_tfs['prots1'] = paired_anchor_tfs['prots1'].apply(ast.literal_eval)
paired_anchor_tfs['prots2'] = paired_anchor_tfs['prots2'].apply(ast.literal_eval)
print('paired-anchor-TFs created')

loop_end_counts = defaultdict(int)

for _, row in paired_anchor_tfs.iterrows():
    for prot in row['prots1']:
        loop_end_counts[prot] += 1
    for prot in row['prots2']:
        loop_end_counts[prot] += 1

loops = len(paired_anchor_tfs)
loop_ends = loops*2
loop_end_proportions = {protein: float(count) / loop_ends for protein, count in loop_end_counts.items()}

# Unpairing loop anchors
end1 = paired_anchor_tfs[['ch1', 'start1', 'end1', 'prots1']]
end2 = paired_anchor_tfs[['ch2', 'start2', 'end2', 'prots2']]
end1.columns = ['ch', 'start', 'end', 'prots']
end2.columns = ['ch', 'start', 'end', 'prots']
unpaired_anchor_tfs = pd.concat([end1, end2])
unpaired_anchor_tfs['prots'] = unpaired_anchor_tfs['prots'].apply(tuple)
unpaired_anchor_tfs = unpaired_anchor_tfs.drop_duplicates()
unpaired_anchor_tfs['prots'] = unpaired_anchor_tfs['prots'].apply(list)

anchor_counts = defaultdict(int)

for _, row in unpaired_anchor_tfs.iterrows():
    for protein in row['prots']:
        anchor_counts[protein] += 1

loop_anchors = len(unpaired_anchor_tfs)
anchor_proportions = {protein: float(count) / loop_anchors for protein, count in anchor_counts.items()}

print('proportions created')

# MARK: Step 3
"""
Paired-anchor-TFs file -> matrices
"""

proteins = set(paired_anchor_tfs['prots1'].explode()).union(set(paired_anchor_tfs['prots2'].explode()))
proteins = {x for x in proteins if pd.notna(x)}
proteins=sorted(proteins)

cis_matrix = pd.DataFrame(1, index=proteins, columns=proteins)

for prots in unpaired_anchor_tfs['prots']:
    for i in range(len(prots)):
        for j in range(i, len(prots)):
            protein1 = prots[i]
            protein2 = prots[j]
            cis_matrix.at[protein1, protein2] += 1
            cis_matrix.at[protein2, protein1] += 1

print('cis_matrix created')

proteins = set(paired_anchor_tfs['prots1'].explode()).union(set(paired_anchor_tfs['prots2'].explode()))
trans_matrix = pd.DataFrame(1, index=proteins, columns=proteins)
for _, row in paired_anchor_tfs.iterrows():
    proteins1 = set(row['prots1'])
    proteins2 = set(row['prots2'])
    for protein in proteins1:
        for other_protein in proteins2:
            trans_matrix.at[protein, other_protein] += 1
    for protein in proteins2:
        for other_protein in proteins1:
            trans_matrix.at[protein, other_protein] += 1
    # Remove rows and columns with no name
    trans_matrix = trans_matrix[trans_matrix.index.notna()]
    trans_matrix = trans_matrix.loc[:, trans_matrix.columns.notna()]

print('trans_matrix created')

# MARK: Step 4
"""
Proportions, matrices -> Enrichment matrices
"""

cis_enrichment_matrix = cis_matrix.copy()
trans_enrichment_matrix = trans_matrix.copy()

for row in cis_enrichment_matrix.index:
    for col in cis_enrichment_matrix.columns:
        if row in proportions and col in proportions:
            cis_enrichment_matrix.loc[row, col] /= (loop_anchors * float(proportions[row]) * float(proportions[col])+1) 
        else:
            print('no proportion for:' + str(row)+ ', ' + str(col))

for row in trans_enrichment_matrix.index:
    for col in trans_enrichment_matrix.columns:
        if row in proportions and col in proportions:
            trans_enrichment_matrix.loc[row, col] /= (loops * float(proportions[row]) * float(proportions[col])+1) 
        else:
            print('no proportion for:' + str(row)+ ', ' + str(col))


# MARK: Step 5
"""
Enrichment matrices -> filtered enrichment matrices
"""
geq10k_proteins = [protein for protein, count in loop_end_counts.items() if count >= 10000]
geq10k_cis_enrichment_matrix = cis_enrichment_matrix.loc[geq10k_proteins, geq10k_proteins]
geq10k_trans_enrichment_matrix = trans_enrichment_matrix.loc[geq10k_proteins, geq10k_proteins]

# MARK: Step 6
"""
All matrices -> plotted matrices
"""

dfs = {'cis_enrichment_matrix': cis_enrichment_matrix,
        'trans_enrichment_matrix': trans_enrichment_matrix,
        'geq10k_cis_enrichment_matrix': geq10k_cis_enrichment_matrix,
        'geq10k_trans_enrichment_matrix': geq10k_trans_enrichment_matrix}
for name, enrichment_df in dfs:
    max = enrichment_df.max().max()
    #log, replace 0 with a 1x-10
    log10_enrichment_df = np.log10(enrichment_df.replace(0, 10 ** (-np.log10(max))))

    #calc vmin vmax values
    vmin = log10_enrichment_df.min().min()
    vmax = log10_enrichment_df.max().max()

    linkage_matrix = linkage(enrichment_df, method='ward')

    # Apply Optimal Leaf Ordering to the linkage matrix
    linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, enrichment_df)

    # Get the ordered indices after optimal leaf ordering
    ordered_index = leaves_list(linkage_matrix_olo)

    ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
    print(ordered_protein_names)

    ordered_log10_enrichment_df = log10_enrichment_df.iloc[ordered_index, ordered_index]

    fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), gridspec_kw={'width_ratios': [1, 15]})
    ax_heatmap.set_aspect('equal')
    dendro = dendrogram(linkage_matrix_olo, labels=ordered_log10_enrichment_df.index, orientation='left', ax=ax_dendro)
    ax_dendro.invert_yaxis()  # Reverse the y-axis to make the dendrogram go from top to bottom
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])

    heatmap = sns.heatmap(
        ordered_log10_enrichment_df, 
        cmap='RdBu_r', # 'rdbu_r' "Spectral_r",    
        annot=False,      
        linewidths=0.05,   
        # cbar_kws={'label': 'Log10 Enrichment Score'}, 
        center=0,  
        ax=ax_heatmap,
        xticklabels=True,
        yticklabels=True,
        vmin=vmin,  
        vmax=vmax   
    )

    controls = ['CTCF', 'STAG1', 'SMC3', 'RAD21']
    highlight_labels = ['RNF219', 'ZNF615', 'ZSCAN30', 'SSRP1', 'AFF4', 'ZC3H4', 'ZBTB2', 'ZBTB25', 'ZNF639', 'ZNF816', 'ZNF670', 'ZFP82', 'ZNF608', 'EED', 'ZNF747', 'ZNF772', 'ZSCAN31', 'E2F1', 'TFDP1', 'TFDP2', 'ZNF430', 'ZNF221', 'ZNF605', 'PBX2', 'SP2', 'NFYA', 'NFYB', 'NFYC', 'ZNF646', 'ZNF865', 'ZMYM2', 'KDM1A', 'ZMYM3', 'ZC3H13', 'AKAP8', 'AKAP8L', 'JUN', 'FOSL1', 'JUNB', 'ATF3', 'JUND', 'CEBPD', 'CEBPA', 'CEBPG', 'HLF', 'NFIL3']
    for label in ax_heatmap.get_yticklabels():
        #print(label)
        if label.get_text() in highlight_labels:
            label.set_weight('bold')
            label.set_color('red')
            label.set_fontsize(12)
        if label.get_text() in controls:
            label.set_weight('bold')
            label.set_color('blue')
            label.set_fontsize(12)
    # Stagger the y-tick labels
    yticks = np.arange(len(ordered_log10_enrichment_df.index)) + 0.5
    yticklabels = ordered_log10_enrichment_df.index
    for i, label in enumerate(ax_heatmap.get_yticklabels()):
        if i % 2 == 0:
            label.set_x(-0.0005)  # Shift to the left
        else:
            label.set_x(-0.018)  # Shift further to the left

    # Adjust tick lengths for staggered labels
    ax_heatmap.yaxis.set_tick_params(which='both', length=0)
    for i, tick in enumerate(ax_heatmap.yaxis.get_major_ticks()):
        if i % 2 == 0:
            tick.tick1line.set_markersize(2.75)  # Set longer tick length
        else:
            tick.tick1line.set_markersize(47)  # Set shorter tick length

    #ax_heatmap.yaxis.tick_right()
    #ax_heatmap.yaxis.set_label_position('right')

    heatmap.set_xticks(np.arange(len(ordered_log10_enrichment_df.columns)) + 0.5)
    heatmap.set_yticks(np.arange(len(ordered_log10_enrichment_df.index)) + 0.5)
    heatmap.set_xticklabels(ordered_log10_enrichment_df.columns, rotation=90, fontsize=6)
    heatmap.set_yticklabels(ordered_log10_enrichment_df.index, rotation=0, fontsize=9)

    # Access the colorbar
    colorbar = heatmap.collections[0].colorbar
    # Set the font size for the colorbar label
    colorbar.ax.yaxis.label.set_size(30)  # Set the desired font size
    # Optional: set the font size for the colorbar ticks
    colorbar.ax.tick_params(labelsize=30)
    log_ticks = colorbar.get_ticks()

    # Convert log10 ticks back to the original scale
    true_ticks = [10 ** tick for tick in log_ticks]

    # Format the tick labels
    # For values >= 1, display as integers; for small values, use scientific notation
    formatted_ticks = [
        f"{t:.2e}" if t < 1 else f"{int(t):,}" for t in true_ticks
    ]

    # Update colorbar with formatted labels
    colorbar.set_ticks(log_ticks)
    colorbar.set_ticklabels(formatted_ticks)# Format as integers with commas

    plt.subplots_adjust(wspace=0.07)
    # plt.suptitle('Enrichment Heatmap of TF Overlaps with Clustering (Log10 Scale)')
    # plt.xlabel('Transcription Factors')
    # plt.ylabel('Transcription Factors')
    plt.savefig('olo-log10' + name[:-3] + 'png', dpi=300, bbox_inches='tight')