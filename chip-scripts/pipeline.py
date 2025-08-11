"""
This is the culmination of all of my matrix-producing scripts.

Each step is denumerated with one of these comments.

You can start this pipeline at different points by commenting out 
preceeding steps and uncommenting the file importer for the desired step.

Inputs: ChIP-seq peak calls (with data), loop list
Output: Protein contact matrices in cis, trans, with various filters.
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import ast
import pickle
import os
import csv
import seaborn as sns
from collections import defaultdict
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering


# MARK: Paired-Anchor-TFs
"""
ChIP-seq peak calls + loop list -> paired-anchor-TFs file
"""
chipdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/K562/beds/merged-filtered'
paired_anchors_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe'

def lookup_dic(chipdir):
    interval_origins = defaultdict(list) #the dic
    for bed_filename in os.listdir(chipdir):
        if bed_filename.endswith('.bed'): #verify filetype
            bed_filepath = os.path.join(chipdir, bed_filename) #get path
            with open(bed_filepath, 'r') as bed_file: #open
                reader = csv.reader(bed_file, delimiter='\t') #make csv
                for row in reader: #parse
                    interval = (row[0], row[1], row[2])  # (ch2, start2, end2) in origin's inbed var
                    interval_origins[interval].append(bed_filename.split('-')[0])
    print("made lookup_dic")
    return interval_origins

def find_overlaps(interval, lookup_dict):
    overlap_values = []
    chrom, start, end = interval
    for (lookup_chrom, lookup_start, lookup_end), values in lookup_dict.items():
        if chrom == lookup_chrom and int(start) <= int(lookup_end) and int(end) >= int(lookup_start):
            overlap_values.extend(values)

    return list(set(overlap_values))

paired_anchors = pd.read_csv(paired_anchors_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'ch2', 'start2', 'end2'])
interval_prots = lookup_dic(chipdir)

seen = {}

prots1_list = []
prots2_list = []
n = len(paired_anchors)
    
for i, (_, row) in enumerate(paired_anchors.iterrows(), start=1):
    print(f"Processing row {i} of {n}")

    # Interval 1
    interval1 = (row['ch1'], row['start1'], row['end1'])
    if interval1 not in seen:
        seen[interval1] = find_overlaps(interval1, interval_prots)
    prots1_list.append(seen[interval1])
    
    # Interval 2
    interval2 = (row['ch2'], row['start2'], row['end2'])
    if interval2 not in seen:
        seen[interval2] = find_overlaps(interval2, interval_prots)
    prots2_list.append(seen[interval2])

# Add the new columns
paired_anchors['prots1'] = prots1_list
paired_anchors['prots2'] = prots2_list

paired_anchor_tfs=paired_anchors

print('paired-anchor-TFs created')
extrusion = ['CTCF', 'RAD21', 'SMC3', 'STAG1']
extrusion_set = set(extrusion)  # set for O(1) membership checks
mask = (
    paired_anchor_tfs['prots1'].apply(lambda prots: set(prots).isdisjoint(extrusion_set))
    &
    paired_anchor_tfs['prots2'].apply(lambda prots: set(prots).isdisjoint(extrusion_set))
)

# Keep only the rows that pass the mask
noctcf_paired_anchor_tfs = paired_anchor_tfs[mask]

# File Exporter
noctcf_paired_anchor_tfs.to_csv(
    'no-ctcf-paired-anchor-tfs.tsv',
    sep='\t',
    index=False,
    header=False,
    columns=['ch1', 'start1', 'end2', 'prots1', 'ch2', 'end1', 'start2', 'prots2']
)

paired_anchor_tfs.to_csv(
    'paired-anchor-tfs.tsv',
    sep='\t',
    index=False,
    header=False,
    columns=['ch1', 'start1', 'end2', 'prots1', 'ch2', 'end1', 'start2', 'prots2']
)

# MARK: Proportions
"""
Paired-anchor-TFs file -> proportions
"""

# File importer
# paired_anchor_tfs_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/pipeline/paired-anchor-tfs.tsv'

# paired_anchor_tfs = pd.read_csv(paired_anchor_tfs_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'ch2', 'start2', 'end2', 'prots1', 'prots2'])
# paired_anchor_tfs['prots1'] = paired_anchor_tfs['prots1'].apply(ast.literal_eval)
# paired_anchor_tfs['prots2'] = paired_anchor_tfs['prots2'].apply(ast.literal_eval)


extrusion = ['CTCF'] # 'RAD21', 'SMC3', 'STAG1'
extrusion_set = set(extrusion)  # set for O(1) membership checks
mask = (
    paired_anchor_tfs['prots1'].apply(lambda prots: set(prots).isdisjoint(extrusion_set))
    &
    paired_anchor_tfs['prots2'].apply(lambda prots: set(prots).isdisjoint(extrusion_set))
)

# Keep only the rows that pass the mask
noctcf_paired_anchor_tfs = paired_anchor_tfs[mask]
# print(paired_anchor_tfs)
# print(noctcf_paired_anchor_tfs)
print('filtered paired-anchor-TFs created')

def proportion_maker(paired_anchor_tfs):
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
    return anchor_proportions, loop_end_proportions, loop_end_counts, unpaired_anchor_tfs, loop_anchors, loop_ends

anchor_proportions, loop_end_proportions, loop_end_counts, unpaired_anchor_tfs, loop_anchors, loop_ends = proportion_maker(paired_anchor_tfs)
noctcf_anchor_proportions, noctcf_loop_end_proportions, noctcf_loop_end_counts, noctcf_unpaired_anchor_tfs, noctcf_loop_anchors, noctcf_loop_ends = proportion_maker(noctcf_paired_anchor_tfs)


# MARK: Matrices
"""
Paired-anchor-TFs file -> matrices
"""
def matrix_maker(paired_anchor_tfs, unpaired_anchor_tfs):
    proteins = set(paired_anchor_tfs['prots1'].explode()).union(set(paired_anchor_tfs['prots2'].explode()))
    proteins = {x for x in proteins if pd.notna(x)}
    proteins=sorted(proteins) # Converts set to list

    cis_matrix = pd.DataFrame(1, index=proteins, columns=proteins)

    for prots in unpaired_anchor_tfs['prots']:
        cis_matrix.loc[prots, prots] +=1

    print('cis_matrix created')

    trans_matrix = pd.DataFrame(1, index=proteins, columns=proteins)
    for _, row in paired_anchor_tfs.iterrows():
        proteins1 = row['prots1']
        proteins2 = row['prots2']
        trans_matrix.loc[proteins1, proteins2] += 1
        trans_matrix.loc[proteins2, proteins1] += 1
    # Remove rows and columns with no name
    trans_matrix = trans_matrix[trans_matrix.index.notna()]
    trans_matrix = trans_matrix.loc[:, trans_matrix.columns.notna()]
    return cis_matrix, trans_matrix

# # Make Files
cis_matrix, trans_matrix = matrix_maker(paired_anchor_tfs, unpaired_anchor_tfs)
noctcf_cis_matrix, noctcf_trans_matrix = matrix_maker(noctcf_paired_anchor_tfs, noctcf_unpaired_anchor_tfs)

# # File exporter
cis_matrix.to_csv("cis_matrix.tsv",sep="\t",index=True, header=True)
trans_matrix.to_csv("trans_matrix.tsv",sep="\t",index=True, header=True)
noctcf_cis_matrix.to_csv("noctcf_cis_matrix.tsv",sep="\t",index=True, header=True)
noctcf_trans_matrix.to_csv("noctcf_trans_matrix.tsv",sep="\t",index=True, header=True)

# print('trans_matrix created')

# MARK: Enrichments
"""
Proportions, matrices -> Enrichment matrices
"""

# File importer
# cis_matrix   = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/basic-matrices/cis_matrix.tsv", sep="\t",index_col=0)
# trans_matrix   = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/basic-matrices/trans_matrix.tsv", sep="\t",index_col=0)
# noctcf_cis_matrix   = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/basic-matrices/noctcf_cis_matrix.tsv", sep="\t",index_col=0)
# noctcf_trans_matrix   = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/basic-matrices/noctcf_trans_matrix.tsv", sep="\t",index_col=0)


def enricher(anchor_proportions, loop_end_proportions, cis_matrix, trans_matrix, loop_anchors, loop_ends):
    cis_enrichment_matrix = cis_matrix.copy()
    cis_enrichment_matrix = cis_enrichment_matrix.astype(float)
    trans_enrichment_matrix = trans_matrix.copy()
    trans_enrichment_matrix = trans_enrichment_matrix.astype(float)

    for row in cis_enrichment_matrix.index:
        for col in cis_enrichment_matrix.columns:
            cis_enrichment_matrix.loc[row, col] /= (float(loop_anchors) * float(anchor_proportions[row]) * float(anchor_proportions[col])+1) 

    for row in trans_enrichment_matrix.index:
        for col in trans_enrichment_matrix.columns:
            trans_enrichment_matrix.loc[row, col] /= (float(loop_ends) * float(loop_end_proportions[row]) * float(loop_end_proportions[col])+1) 
    return cis_enrichment_matrix, trans_enrichment_matrix

cis_enrichment_matrix, trans_enrichment_matrix = enricher(anchor_proportions, loop_end_proportions, cis_matrix, trans_matrix, loop_anchors, loop_ends)
noctcf_cis_enrichment_matrix, noctcf_trans_enrichment_matrix = enricher(noctcf_anchor_proportions, noctcf_loop_end_proportions, noctcf_cis_matrix, noctcf_trans_matrix, noctcf_loop_anchors, noctcf_loop_ends)

print('enriched')

# MARK: Filters
"""
Enrichment matrices -> filtered enrichment matrices
"""
cutoffs = [5, 10, 20, 50]

dfs = {"cis_enrichment_matrix": cis_enrichment_matrix,
        "trans_enrichment_matrix": trans_enrichment_matrix,
        "noctcf_cis_enrichment_matrix": noctcf_cis_enrichment_matrix,
        "noctcf_trans_enrichment_matrix": noctcf_trans_enrichment_matrix}

for cutoff in cutoffs:
    label = f"geq{cutoff}k"

    # Select proteins that meet the cutoff
    filtered_proteins = [
        protein for protein, count in loop_end_counts.items() 
        if count >= (cutoff * 1000)
    ]

    # Slice the matrices
    filtered_cis_matrix = cis_enrichment_matrix.loc[filtered_proteins, filtered_proteins]
    filtered_trans_matrix = trans_enrichment_matrix.loc[filtered_proteins, filtered_proteins]
    filtered_proteins = [
        protein for protein, count in noctcf_loop_end_counts.items() 
        if loop_end_counts[protein] >= (cutoff * 1000)
    ]
    filtered_noctcf_cis_matrix = noctcf_cis_enrichment_matrix.loc[filtered_proteins, filtered_proteins]
    filtered_noctcf_trans_matrix = noctcf_trans_enrichment_matrix.loc[filtered_proteins, filtered_proteins]
    dfs[f"geq{cutoff}k_cis_enrichment_matrix"] = filtered_cis_matrix
    dfs[f"geq{cutoff}k_trans_enrichment_matrix"] = filtered_trans_matrix
    dfs[f"geq{cutoff}k_noctcf_cis_enrichment_matrix"] = filtered_noctcf_cis_matrix
    dfs[f"geq{cutoff}k_noctcf_trans_enrichment_matrix"] = filtered_noctcf_trans_matrix



# MARK: Plots
"""
All matrices -> plotted matrices
"""
def plot(name, enrichment_df):
    maxval = enrichment_df.max().max()
    log10_enrichment_df = np.log10(enrichment_df.replace(0, 10 ** (-np.log10(maxval))))

    vmin = log10_enrichment_df.min().min()
    vmax = log10_enrichment_df.max().max()

    linkage_matrix = linkage(enrichment_df, method='ward')

    linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, enrichment_df)

    ordered_index = leaves_list(linkage_matrix_olo)

    ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
    print(name)
    print(ordered_protein_names)

    ordered_log10_enrichment_df = log10_enrichment_df.iloc[ordered_index, ordered_index]

    fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), gridspec_kw={'width_ratios': [1, 15]})
    ax_heatmap.set_aspect('equal')
    dendro = dendrogram(linkage_matrix_olo, labels=ordered_log10_enrichment_df.index, orientation='left', ax=ax_dendro)

    ax_dendro.invert_yaxis()  # Reverses the y-axis to make the dendrogram go from top to bottom
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
    """
    Option A: Highlight certain labels
    """
    # controls = ['CTCF', 'STAG1', 'SMC3', 'RAD21']
    # highlight_labels = ['RNF219', 'ZNF615', 'ZSCAN30', 'SSRP1', 'AFF4', 'ZC3H4', 'ZBTB2', 'ZBTB25', 'ZNF639', 'ZNF816', 'ZNF670', 'ZFP82', 'ZNF608', 'EED', 'ZNF747', 'ZNF772', 'ZSCAN31', 'E2F1', 'TFDP1', 'TFDP2', 'ZNF430', 'ZNF221', 'ZNF605', 'PBX2', 'SP2', 'NFYA', 'NFYB', 'NFYC', 'ZNF646', 'ZNF865', 'ZMYM2', 'KDM1A', 'ZMYM3', 'ZC3H13', 'AKAP8', 'AKAP8L', 'JUN', 'FOSL1', 'JUNB', 'ATF3', 'JUND', 'CEBPD', 'CEBPA', 'CEBPG', 'HLF', 'NFIL3']
    # for label in ax_heatmap.get_yticklabels():
    #     #print(label)
    #     if label.get_text() in highlight_labels:
    #         label.set_weight('bold')
    #         label.set_color('red')
    #         label.set_fontsize(12)
    #     if label.get_text() in controls:
    #         label.set_weight('bold')
    #         label.set_color('blue')
    #         label.set_fontsize(12)
    """
    Option B: Colormap from average distance to TSS
    """
    distance_df = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/enhancer-promoter/avg_tss_dist.csv")  # or pass this as an argument
    distance_dict = dict(zip(distance_df["protein"], distance_df["average_distance_to_tss"]))

    # Normalize distances to 0–1 for colormap
    distances = [distance_dict[prot] for prot in enrichment_df.index if prot in distance_dict and not pd.isna(distance_dict[prot])]
    vmin_dist = min(distances)
    vmax_dist = max(distances)
    norm = mcolors.Normalize(vmin=vmin_dist, vmax=vmax_dist)
    cmap = mpl.colormaps["viridis"]
    for label in ax_heatmap.get_yticklabels():
        protein = label.get_text()
        distance = distance_dict.get(protein)

        if distance is not None and not pd.isna(distance):
            color = cmap(norm(distance))
            label.set_color(color)
            label.set_fontsize(9)
        else:
            label.set_color("black")
            label.set_fontsize(9)


    
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
    plt.savefig(name, dpi=300, bbox_inches='tight')
    return(ordered_protein_names)

orders = {'tss': ['ZFP62', 'POLR2AphosphoS2', 'PAF1', 'H3K79me2', 'ZNF597', 'PPARGC1A', 'HDAC6', 'ZNF326', 'MLLT10', 'ELK4', 'MATR3', 'AKAP8L', 'AFF4', 'ZC3H13', 'ZNF713', 'NCOA5', 'ZNF513', 'H3K9ac', 'SREBF2', 'KDM5A', 'ZNF839', 'HMGA1', 'AKAP8', 'ZNF674', 'KIAA2018', 'CREBL2', 'ZNF138', 'CAMTA2', 'PAWR', 'GMEB2', 'ZFP41', 'ZNF571', 'ZNF569', 'ZC3H4', 'SAFB2', 'REL', 'H4K20me1', 'IRF3', 'DMTF1', 'ZNF30', 'KLF13', 'NR0B2', 'ZNF446', 'EEA1', 'SMAD1', 'DZIP1', 'SOX18', 'ZNF827', 'SP110', 'ATF5', 'ZNF235', 'FUBP1', 'SIN3A', 'E2F2', 'ZNF567', 'ZNF879', 'ZNF778', 'ZBTB2', 'DR1', 'ZFP82', 'ZNF639', 'IRF5', 'ZIK1', 'ZNF776', 'ZNF790', 'ZNF468', 'PAX8', 'THYN1', 'ZNF546', 'ZNF451', 'KLF15', 'THAP4', 'ARHGAP35', 'THAP7', 'TIGD3', 'KLF9', 'CENPT', 'NFATC3', 'PRMT3', 'ZNF234', 'ZNF558', 'ZNF333', 'ZNF678', 'SMAD9', 'ZNF367', 'ZXDC', 'SNAI1', 'ZSCAN30', 'ZNF275', 'REXO4', 'PREB', 'MXD1', 'GLI4', 'ZBTB8A', 'ZZZ3', 'ZNF383', 'ZNF883', 'ZNF577', 'TRAFD1', 'ZNF550', 'THAP8', 'CHCHD3', 'MTERF4', 'KAT7', 'EED', 'HCFC1', 'ZNF704', 'TSC22D1', 'MZF1', 'SIN3B', 'TMF1', 'ZNF552', 'POLR2AphosphoS5', 'MXD3', 'ZNF749', 'ZNF34', 'SRY', 'ZNF180', 'ZNF354B', 'TIGD6', 'ZNF619', 'IKZF4', 'CSRNP1', 'CSRNP2', 'BHLHA15', 'ZNF83', 'PLSCR1', 'ZNF44', 'ZNF570', 'JARID2', 'GLYR1', 'ZNF572', 'NFE2L1', 'ZNF672', 'FUBP3', 'HIVEP1', 'IRF9', 'ZNF782', 'SREBF1', 'KMT2A', 'NR1H2', 'MTERF2', 'KLF12', 'DNMT1', 'CERS6', 'ZNF318', 'GMEB1', 'THAP9', 'ZNF547', 'ZNF527', 'CENPBD1', 'ISX', 'BATF2', 'ZNF747', 'ZNF33B', 'NRL', 'ZNF510', 'H3K4me3', 'ZNF891', 'NR3C1', 'MXD4', 'E2F1', 'GTF3A', 'ZNF407', 'ZNF224', 'SP4', 'ZBTB46', 'ZNF485', 'ELK1', 'ELF4', 'PHF20', 'ZNF589', 'IRX3', 'DBP', 'ARID2', 'ZFYVE20', 'ZBTB1', 'ZNF607', 'ZNF253', 'ZNF598', 'ZNF225', 'ZHX2', 'LBX2', 'BCL3', 'PIN1', 'POLR2A', 'ZNF75D', 'ZNF296', 'MXI1', 'NAIF1', 'ZNF136', 'SSRP1', 'ZNF484', 'E2F5', 'HES4', 'AKNA', 'FOXC1', 'NFKB2', 'AHR', 'ARNT2', 'ZNF737', 'HOXA9', 'ZNF660', 'ZNF784', 'ZNF781', 'DMAP1', 'MAF1', 'ZNF608', 'H3K4me2', 'TAF1', 'SALL2', 'CC2D1A', 'ZNF25', 'ZMAT3', 'ZNF786', 'ZNF256', 'ZNF251', 'ZNF483', 'SNAPC4', 'BRD4', 'ZNF280B', 'H3K36me3', 'ZNF10', 'ZC3H8', 'SP2', 'ZSCAN31', 'ZNF773', 'SUZ12', 'ZNF343', 'ZBTB38', 'ZNF709', 'FLYWCH1', 'ZNF230', 'GRHL1', 'YEATS4', 'ZNF530', 'ZHX1', 'TOPORS', 'NPAS2', 'ETV6', 'NFAT5', 'RNF219', 'SPEN', 'YEATS2', 'ZNF181', 'MAZ', 'KDM5B', 'SMAD7', 'KMT2B', 'FOXO4', 'HINFP', 'ZNF563', 'ZNF276', 'JDP2', 'ZNF274', 'ARID4A', 'CCDC6', 'ZSCAN12', 'CREB3', 'FOXM1', 'BRCA1', 'MYNN', 'ZNF772', 'ZNF605', 'DRAP1', 'ZNF496', 'ZNF724P', 'ZNF746', 'MTF1', 'ZBTB7A', 'ZSCAN21', 'BAZ2A', 'ZNF460', 'HBP1', 'MYC', 'CEBPZ', 'TBP', 'ZNF335', 'MED8', 'ZNF350', 'PHF8', 'CHD2', 'ZNF615', 'ZNF766', 'ZNF788', 'KLF6', 'BRF2', 'ZNF430', 'ZFP90', 'ZNF221', 'BORCS8', 'ZBTB4', 'ZNF414', 'ZBTB3', 'ZNF543', 'SFPQ', 'RARG', 'XBP1', 'ZNF501', 'MTA1', 'NFIB', 'ZNF142', 'PBX3', 'ZNF18', 'RFXAP', 'ADNP', 'ZNF761', 'GATAD1', 'ZBTB10', 'PHF5A', 'ARID4B', 'ZBTB42', 'ZNF691', 'ZNF697', 'TFDP2', 'TEAD2', 'ZNF33A', 'ZNF548', 'YY1', 'ZFP36L1', 'FBXL19', 'ZNF580', 'ZNF670', 'MYRF', 'HMGXB4', 'ZFP37', 'TSC22D2', 'ZNF770', 'ATRX', 'ZSCAN20', 'NCOA1', 'POGK', 'ZMAT5', 'RFXANK', 'ZHX3', 'ZNF740', 'ZBED4', 'TGIF2', 'RCOR1', 'GPN1', 'KLF11', 'RELA', 'ZSCAN25', 'GLMP', 'ZNF780A', 'ZNF556', 'ZFP91', 'SNAPC2', 'MNX1', 'NFYA', 'ZNF850', 'ZFP14', 'ZUFSP', 'E2F4', 'ZFX', 'ZNF232', 'UBTF', 'ZBTB25', 'ZNF101', 'KDM3A', 'ZFY', 'ZNF121', 'JRK', 'HOXA7', 'SMYD3', 'STAT6', 'SNAPC5', 'ZNF816', 'E2F8', 'MEIS1', 'ZSCAN29', 'ZNF720', 'ZNF329', 'ZNF48', 'ATAD3A', 'PRDM10', 'SATB2', 'HMGXB3', 'ZCCHC11', 'ZNF511', 'HIC2', 'ZFP36L2', 'H3K27ac', 'RNF2', 'TFDP1', 'ZBTB14', 'ZNF557', 'ATF6', 'POU2F1', 'HOXA10', 'ZFAT', 'KDM2A', 'MYPOP', 'ZBTB24', 'ZNF260', 'HOXA3', 'CBX1', 'HOXA5', 'ZBTB39', 'ZBTB7B', 'PRDM4', 'ARNTL', 'TERF1', 'ZNF703', 'CBFB', 'REPIN1', 'ZSCAN22', 'STAT5B', 'ERF', 'HDAC1', 'KLF16', 'ZMYM2', 'DNMT3B', 'ZNF20', 'KDM4B', 'MTA3', 'ZNF17', 'KAT8', 'ZNF775', 'IKZF5', 'ZNF777', 'ZBTB49', 'RBPJ', 'ZNF337', 'MBD1', 'ZNF143', 'ZNF431', 'MED13', 'ZBED5', 'PATZ1', 'ZMYM3', 'SAP130', 'EGR1', 'HNF1B', 'ZNF574', 'SP5', 'ZNF878', 'ZNF280D', 'NKX3', 'DDIT3', 'HOXD1', 'THAP11', 'ZNF441', 'ZNF562', 'HSF1', 'ESRRA', 'MLXIP', 'ZNF564', 'ELF1', 'CHD4', 'LIN54', 'TCF3', 'NR2C2', 'TCF12', 'SKIL', 'ZSCAN5A', 'PHF21A', 'MEIS2', 'ZBTB33', 'ZNF616', 'SP140L', 'ZBTB44', 'ZNF687', 'NFYC', 'ZNF792', 'U2AF1', 'ATF1', 'MGA', 'ZSCAN9', 'SRF', 'IRF2', 'GATAD2B', 'ZNF205', 'ETV4', 'ZGPAT', 'ZNF576', 'MEF2D', 'TBX2', 'H2AFZ', 'ZEB1', 'CBX5', 'ZBTB34', 'ZKSCAN8', 'MIER3', 'ETV5', 'ZNF263', 'GTF2F1', 'FIP1L1', 'ZNF768', 'ZKSCAN5', 'CENPX', 'ZBTB26', 'SRSF4', 'NFYB', 'MAX', 'LRRFIP1', 'SMAD3', 'ZNF512B', 'ZNF614', 'ZNF841', 'MED1', 'ZBTB37', 'ZNF526', 'HSF2', 'ZNF549', 'ZBTB43', 'ZNF124', 'MYBL2', 'AGO1', 'IRF1', 'ZNF264', 'NFKBIZ', 'ZBTB21', 'HNRNPH1', 'NACC2', 'TARDBP', 'GABPB1', 'SRSF1', 'RXRA', 'GPBP1L1', 'ZNF710', 'ZFP64', 'ZNF800', 'RFX5', 'GZF1', 'TFAP4', 'ZNF432', 'RREB1', 'SALL1', 'TFE3', 'PITX1', 'FOXK1', 'THRA', 'AGO2', 'MTF2', 'ATF3', 'ELF3', 'NRF1', 'THRB', 'WIZ', 'ETS1', 'RCOR2', 'MLX', 'ZKSCAN1', 'FOXK2', 'ZFP1', 'ZNF3', 'SKI', 'RFX3', 'HDAC2', 'ZNF292', 'SP1', 'ZBTB20', 'GABPA', 'SMAD4', 'CREM', 'BHLHE40', 'FOXQ1', 'PRDM15', 'FUS', 'KDM1A', 'ZNF362', 'ZNF331', 'ZNF512', 'ZNF217', 'H3K4me1', 'RXRB', 'MBNL3', 'ZNF382', 'PBX2', 'GATA4', 'PTBP1', 'ZNF219', 'HNRNPL', 'LCORL', 'GFI1', 'TP53', 'RORA', 'RBM22', 'TEAD3', 'PPARG', 'FOXO1', 'ZNF503', 'NR2F6', 'ZNF160', 'ZNF629', 'TBL1XR1', 'TAF15', 'TEF', 'LCOR', 'RBFOX2', 'CCAR2', 'PRPF4', 'ZNF281', 'DLX6', 'POLR2G', 'POGZ', 'ZNF609', 'ZNF646', 'MEF2A', 'NCOR1', 'PROX1', 'USF1', 'HNF4G', 'HNF1A', 'NONO', 'GATAD2A', 'CEBPD', 'TEAD4', 'RARA', 'FOXP4', 'ZNF317', 'HMG20B', 'RBM39', 'ZMYM4', 'ATF4', 'NR2F1', 'ARID3A', 'DPF2', 'NR5A1', 'ISL2', 'JUNB', 'PAXIP1', 'RBAK', 'NFIL3', 'KDM6A', 'SIX4', 'NCOA2', 'MIXL1', 'MAFG', 'NFRKB', 'ZNF384', 'CEBPG', 'NFIA', 'BCL6', 'ZNF652', 'ZNF707', 'FOSL2', 'SOX6', 'NR2F2', 'FOXJ3', 'ZNF865', 'HOMEZ', 'XRCC5', 'PCBP1', 'TEAD1', 'RERE', 'SIX1', 'HNRNPLL', 'EP300', 'ARID5B', 'ZNF12', 'TCF7L2', 'SOX5', 'RFX1', 'MBD4', 'ASH2L', 'USF2', 'FOXA3', 'NFE2', 'ZFP36', 'CBX2', 'ZFHX3', 'SOX13', 'FOSL1', 'TBX3', 'ONECUT2', 'HNF4A', 'ZNF7', 'HMG20A', 'ZBTB40', 'ZNF24', 'TRIM24', 'GATA2', 'ZNF644', 'ZNF207', 'CEBPA', 'NFIC', 'TCF7', 'FOXP1', 'AHDC1', 'SMARCC2', 'SETDB1', 'MIER2', 'REST', 'HHEX', 'U2AF2', 'HNRNPUL1', 'FOXA1', 'SCMH1', 'H3K27me3', 'ONECUT1', 'ZNF490', 'EHMT2', 'EZH2', 'SMC3', 'RAD21', 'SMARCE1', 'HLF', 'KAT2B', 'PCBP2', 'STAG1', 'FOXA2', 'CTCF', 'JUN', 'CEBPB', 'SNRNP70', 'JUND', 'HNRNPK', 'CREB1', 'TOE1', 'ZNF146', 'MNT', 'NFE2L2', 'RAD51', 'YBX1', 'ATF7', 'ATM', 'MAFF', 'ATF2', 'MAFK', 'NBN', 'IKZF1', 'TRIM22', 'SYNCRIP', 'SRSF9', 'ZNF282', 'PLRG1', 'ARNT', 'PHB2', 'BCLAF1', 'H3K9me3']}
for name, enrichment_df in dfs.items():
    print(name, enrichment_df)
    orders[name]=plot(name, enrichment_df)
    # print(orders)
    enrichment_df.to_csv(name,sep="\t",index=True, header=True)
    

order_lock_dfs = {"locked_noctcf_cis_enrichment_matrix": (dfs['noctcf_cis_enrichment_matrix'], orders['cis_enrichment_matrix']),
                "locked_noctcf_trans_enrichment_matrix": (dfs['noctcf_trans_enrichment_matrix'], orders['trans_enrichment_matrix']),
                "locked_geq25k_noctcf_cis_enrichment_matrix": (dfs['geq25k_noctcf_cis_enrichment_matrix'], orders['geq25k_cis_enrichment_matrix']),
                "locked_geq25k_noctcf_trans_enrichment_matrix": (dfs['geq25k_noctcf_trans_enrichment_matrix'], orders['geq25k_trans_enrichment_matrix']),
                "locked_geq50k_noctcf_cis_enrichment_matrix": (dfs['geq50k_noctcf_cis_enrichment_matrix'], orders['geq50k_cis_enrichment_matrix']),
                "locked_geq50k_noctcf_trans_enrichment_matrix": (dfs['geq50k_noctcf_trans_enrichment_matrix'], orders['geq50k_trans_enrichment_matrix']),
                "tss_cis_enrichment_matrix": (dfs['cis_enrichment_matrix'], orders['tss']),
                "tss_trans_enrichment_matrix": (dfs['trans_enrichment_matrix'], orders['tss']),
                "tss_noctcf_cis_enrichment_matrix": (dfs['noctcf_cis_enrichment_matrix'], orders['tss']),
                "tss_noctcf_trans_enrichment_matrix": (dfs['noctcf_trans_enrichment_matrix'], orders['tss']),
                "tss_geq50k_cis_enrichment_matrix": (dfs['geq50k_cis_enrichment_matrix'], orders['tss']),
                "tss_geq50k_trans_enrichment_matrix": (dfs['geq50k_trans_enrichment_matrix'], orders['tss']),
                "tss_geq50k_noctcf_cis_enrichment_matrix": (dfs['geq50k_noctcf_cis_enrichment_matrix'], orders['tss']),
                "tss_geq50k_noctcf_trans_enrichment_matrix": (dfs['geq50k_noctcf_trans_enrichment_matrix'], orders['tss']),
                }

for name, (enrichment_df, order) in order_lock_dfs.items():
    maxval = enrichment_df.max().max()
    log10_enrichment_df = np.log10(enrichment_df.replace(0, 10 ** (-np.log10(maxval))))

    vmin = log10_enrichment_df.min().min()
    vmax = log10_enrichment_df.max().max()

    linkage_matrix = linkage(enrichment_df, method='ward')

    linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, enrichment_df)

    # order = ['JARID2', 'MYNN', 'ZBTB1', 'ZBED5', 'SP140L', 'ZCCHC11', 'ATF6', 'PRDM4', 'STAT5B', 'ZNF720', 'MYRF', 'ZBTB49', 'ZBTB44', 'ATAD3A', 'LRRFIP1', 'ZFP36L2', 'GPN1', 'CERS6', 'TEAD2', 'DBP', 'AKNA', 'ISX', 'ZNF367', 'PREB', 'ZNF180', 'ELK4', 'XBP1', 'PBX3', 'MTF1', 'ZNF773', 'ZNF256', 'ZSCAN12', 'KLF15', 'E2F2', 'FOXM1', 'TP53', 'MAF1', 'ZNF557', 'ZFP37', 'MLXIP', 'CREB3', 'NFAT5', 'SATB2', 'ZSCAN20', 'ZFAT', 'ZFP14', 'ZNF251', 'SMAD9', 'ZNF571', 'ZNF660', 'ZNF484', 'SP110', 'ZNF577', 'ZNF570', 'ZNF101', 'IRF1', 'DMTF1', 'ATF5', 'PAX8', 'ZNF678', 'ZNF318', 'PAWR', 'ZBTB8A', 'RARG', 'THYN1', 'TSC22D1', 'ZNF354B', 'SNAPC5', 'ZNF146', 'CC2D1A', 'CENPBD1', 'ZNF737', 'TIGD3', 'ZNF25', 'SOX18', 'NFATC3', 'MZF1', 'KIAA2018', 'ZNF569', 'FLYWCH1', 'ZNF562', 'ZNF841', 'PIN1', 'HOXA7', 'CSRNP2', 'JDP2', 'KLF13', 'ZFP41', 'ZNF83', 'ZNF30', 'CREBL2', 'LBX2', 'CHCHD3', 'ZNF253', 'CENPT', 'NR1H2', 'ZFYVE20', 'ZNF468', 'TERF1', 'REXO4', 'PLSCR1', 'GLMP', 'ZNF530', 'ZNF75D', 'ZNF17', 'ZNF672', 'HOXA9', 'BATF2', 'REL', 'ZNF589', 'ZNF780A', 'THAP8', 'ZNF234', 'ZNF485', 'THAP4', 'IRX3', 'ZNF616', 'IKZF4', 'ZNF552', 'ARHGAP35', 'ZNF138', 'EEA1', 'DZIP1', 'ZNF827', 'ZNF879', 'PRMT3', 'ZNF674', 'NFE2L1', 'ZNF778', 'ZNF20', 'ZNF839', 'ZNF707', 'ZNF276', 'ZNF782', 'ZNF510', 'FUBP3', 'ZNF44', 'ZNF786', 'ZNF33B', 'GLI4', 'ZNF337', 'ZNF136', 'ZBTB42', 'ZNF704', 'ZNF275', 'ZNF160', 'ZIK1', 'ZNF451', 'ZNF724P', 'ZNF749', 'MTERF2', 'ZNF333', 'ZNF235', 'CAMTA2', 'ZNF567', 'SMAD1', 'ZNF619', 'ZNF527', 'ZNF550', 'ZNF713', 'TRAFD1', 'FUBP1', 'RELA', 'NR3C1', 'CBX1', 'PHF5A', 'KDM4B', 'ZNF124', 'RFXANK', 'TOE1', 'RCOR1', 'DDIT3', 'RNF219', 'ZNF615', 'ZSCAN30', 'ZHX3', 'KLF11', 'TSC22D2', 'ETV6', 'AHR', 'KLF9', 'DNMT3B', 'TCF3', 'CBX2', 'ELK1', 'BRCA1', 'ZNF766', 'NKX3', 'ZNF232', 'ZNF697', 'WIZ', 'DNMT1', 'ZNF225', 'EED', 'SAFB2', 'IRF5', 'ZNF546', 'ZNF329', 'ZNF790', 'ZNF446', 'TMF1', 'ZNF181', 'BRF2', 'ZNF431', 'ZXDC', 'ZNF639', 'E2F5', 'ZNF608', 'ZFP82', 'ZNF776', 'NR0B2', 'TIGD6', 'MTERF4', 'HES4', 'ZNF558', 'GZF1', 'ZNF781', 'THAP7', 'ZNF224', 'SNAI1', 'ZNF703', 'ZNF576', 'ZBTB24', 'NFIB', 'ZNF296', 'ZNF432', 'ZMAT5', 'MTA3', 'ATRX', 'ZUFSP', 'GRHL1', 'ZNF496', 'NPAS2', 'AKAP8L', 'ZNF326', 'NCOA5', 'MLLT10', 'ZNF549', 'ZNF490', 'ZNF775', 'ZNF383', 'GTF3A', 'ZNF34', 'ZC3H8', 'ZZZ3', 'ZBTB46', 'CCAR2', 'FIP1L1', 'U2AF1', 'PTBP1', 'ZFP62', 'ZNF597', 'ZNF513', 'SETDB1', 'GTF2F1', 'HCFC1', 'SIN3B', 'GMEB2', 'ZNF526', 'HOXA5', 'RFXAP', 'CCDC6', 'ZSCAN21', 'HOXA10', 'HIVEP1', 'ZNF772', 'ZSCAN31', 'MXD4', 'BORCS8', 'ARNT2', 'NFKB2', 'GLYR1', 'ZNF441', 'ZNF816', 'ZNF670', 'ZNF18', 'NAIF1', 'ELF4', 'ZC3H13', 'ZBTB2', 'SALL2', 'KLF12', 'E2F1', 'DR1', 'ZNF891', 'NRL', 'ZNF407', 'MXD1', 'ZNF563', 'MTF2', 'ZNF260', 'ZNF512B', 'STAT6', 'NCOA1', 'ZBTB39', 'FOXQ1', 'HOXD1', 'MED13', 'ZBTB10', 'HBP1', 'KLF6', 'ZHX1', 'ZHX2', 'MED8', 'ZNF414', 'ZNF691', 'BAZ2A', 'ZNF747', 'BRD4', 'KAT7', 'SNAPC4', 'JRK', 'ZBTB4', 'AKAP8', 'ZNF572', 'FOXO4', 'TEF', 'HSF2', 'ZKSCAN5', 'RBAK', 'BHLHA15', 'ZBTB3', 'ZNF483', 'ZNF761', 'ZNF543', 'HINFP', 'ZNF230', 'ZNF343', 'ZNF548', 'ZSCAN22', 'ARID2', 'CSRNP1', 'TOPORS', 'SMYD3', 'SNAPC2', 'FBXL19', 'ZSCAN5A', 'ZFP36L1', 'ZNF784', 'ZNF883', 'ZNF274', 'BCL3', 'ZNF547', 'ZNF598', 'ZMAT3', 'ZNF709', 'ZNF607', 'FOXC1', 'IRF9', 'MXD3', 'ZNF280B', 'SRY', 'ZNF564', 'ZFP90', 'ZNF142', 'RNF2', 'ZNF850', 'ZNF33A', 'SRSF4', 'POLR2AphosphoS2', 'AFF4', 'ZC3H4', 'PAF1', 'SSRP1', 'MEF2A', 'ZFP1', 'NFKBIZ', 'ZNF511', 'ZSCAN9', 'ZNF335', 'HMGXB3', 'ZBTB26', 'TGIF2', 'ZBED4', 'ZBTB25', 'ZNF792', 'ZNF205', 'CBX5', 'CBFB', 'KDM3A', 'ZNF580', 'MBD1', 'ZBTB34', 'ZNF878', 'ZNF746', 'TCF12', 'MEIS2', 'ZBTB43', 'ZNF264', 'ZBTB20', 'ZNF280D', 'ZNF317', 'REPIN1', 'ARNTL', 'ESRRA', 'ZSCAN29', 'MEIS1', 'ZNF556', 'SFPQ', 'SRF', 'ZNF221', 'ZNF430', 'POGK', 'ZBTB40', 'ZNF10', 'ZNF512', 'CHD2', 'SP2', 'SIN3A', 'ZNF143', 'THAP9', 'PHF20', 'KMT2A', 'KMT2B', 'YEATS2', 'ATF7', 'E2F8', 'MNT', 'TFDP1', 'TFDP2', 'SP4', 'GMEB1', 'ZNF350', 'DMAP1', 'SPEN', 'ZNF788', 'MYPOP', 'ZNF605', 'ZSCAN25', 'MTA1', 'ZNF777', 'ZBTB7A', 'ZBTB38', 'ZNF768', 'KDM5B', 'SMAD7', 'TAF15', 'PCBP1', 'SRSF1', 'AGO1', 'FUS', 'ZBTB14', 'ZNF800', 'PRDM15', 'ZNF740', 'ZFP64', 'SREBF1', 'U2AF2', 'HNRNPH1', 'NFRKB', 'IRF3', 'CEBPZ', 'SUZ12', 'HNRNPK', 'AGO2', 'NRF1', 'PHF8', 'TAF1', 'ZNF501', 'YY1', 'POLR2AphosphoS5', 'KDM5A', 'E2F4', 'YEATS4', 'MAZ', 'ZNF48', 'MNX1', 'MYC', 'MXI1', 'GATAD1', 'DRAP1', 'TBP', 'NFYA', 'NFYC', 'NFYB', 'RFX5', 'RFX3', 'HDAC6', 'RBM22', 'HNRNPL', 'ZNF770', 'HMGXB4', 'GABPB1', 'ZFY', 'ZFP91', 'UBTF', 'KDM2A', 'RBM39', 'HNRNPLL', 'ZBTB33', 'ZMYM4', 'ARID4B', 'PATZ1', 'XRCC5', 'ARID4A', 'LIN54', 'ASH2L', 'MGA', 'ZFX', 'RBFOX2', 'POLR2G', 'POLR2A', 'H3K9ac', 'MATR3', 'USF1', 'BHLHE40', 'POGZ', 'ZGPAT', 'TBX2', 'HDAC1', 'GABPA', 'ELF1', 'THAP11', 'HOXA3', 'POU2F1', 'ZNF574', 'EGR1', 'NONO', 'SP5', 'ETV4', 'MYBL2', 'NR2C2', 'KAT8', 'HNF1B', 'MBNL3', 'ZNF207', 'RAD51', 'PRPF4', 'PCBP2', 'MAFG', 'NFE2', 'ZNF7', 'ZNF3', 'ZNF382', 'CENPX', 'SCMH1', 'HNRNPUL1', 'PRDM10', 'IKZF5', 'ERF', 'KLF16', 'ZBTB7B', 'TARDBP', 'RBPJ', 'IRF2', 'ZNF710', 'ATF4', 'ZNF281', 'ZNF331', 'ATF2', 'FOXK1', 'TFE3', 'ATF3', 'SKIL', 'ZNF629', 'ZMYM2', 'HIC2', 'ZNF865', 'ZNF263', 'ZNF121', 'USF2', 'ZKSCAN1', 'ZNF24', 'ZBTB21', 'REST', 'MIER2', 'ZNF12', 'SMAD3', 'ZKSCAN8', 'MIER3', 'ZNF614', 'MEF2D', 'NACC2', 'GATAD2B', 'TBL1XR1', 'GATA4', 'RERE', 'ONECUT2', 'JUNB', 'ETS1', 'DLX6', 'ZNF503', 'MLX', 'SIX4', 'HMG20B', 'ZBTB37', 'RORA', 'H3K9me3', 'NFE2L2', 'MAFK', 'MAFF', 'H3K4me1', 'ZNF384', 'HNF4A', 'CEBPA', 'FOXA1', 'FOXA2', 'CEBPB', 'HLF', 'H3K79me2', 'GPBP1L1', 'HMGA1', 'NBN', 'H3K27me3', 'SRSF9', 'PLRG1', 'SMC3', 'CTCF', 'STAG1', 'RAD21', 'SYNCRIP', 'YBX1', 'IKZF1', 'PHB2', 'TRIM22', 'H4K20me1', 'H3K36me3', 'EZH2', 'ZFP36', 'CHD4', 'NR2F2', 'NCOR1', 'NFIC', 'TCF7', 'PROX1', 'FOXP1', 'HMG20A', 'ZNF644', 'GATA2', 'ONECUT1', 'PAXIP1', 'ARID3A', 'SOX6', 'TEAD3', 'FOSL2', 'CEBPG', 'SOX13', 'ARID5B', 'GATAD2A', 'FOXA3', 'TEAD1', 'NFIL3', 'EP300', 'KDM1A', 'AHDC1', 'NR2F1', 'NR5A1', 'HNF4G', 'TCF7L2', 'SMAD4', 'NR2F6', 'RARA', 'SOX5', 'BCL6', 'TEAD4', 'ZNF362', 'JUND', 'ZNF219', 'RXRB', 'LCOR', 'FOXP4', 'NCOA2', 'NFIA', 'CEBPD', 'SNRNP70', 'CREM', 'SP1', 'CREB1', 'SALL1', 'TFAP4', 'SAP130', 'MAX', 'H3K27ac', 'ZNF687', 'H3K4me3', 'H3K4me2', 'H2AFZ', 'ZNF282', 'ARNT', 'PPARGC1A', 'RFX1', 'SREBF2', 'HSF1', 'DPF2', 'KDM6A', 'HHEX', 'ZFHX3', 'SIX1', 'RCOR2', 'FOXJ3', 'TBX3', 'THRB', 'RREB1', 'HOMEZ', 'PPARG', 'MIXL1', 'ZNF652', 'FOSL1', 'JUN', 'EHMT2', 'PBX2', 'FOXK2', 'ATF1', 'ZNF609', 'ZNF217', 'FOXO1', 'THRA', 'HNF1A', 'ISL2', 'PITX1', 'GFI1', 'ZEB1', 'SKI', 'PHF21A', 'RXRA', 'LCORL', 'HDAC2', 'MED1', 'ELF3', 'ETV5', 'ZNF292', 'ZMYM3', 'ZNF460', 'ZNF646', 'ADNP', 'ATM', 'TRIM24', 'MBD4', 'KAT2B', 'SMARCE1', 'SMARCC2']
    # 
    # ['HMGA1', 'ATRX', 'ZNF10', 'ZNF660', 'ATAD3A', 'GPN1', 'ZNF251', 'FUBP1', 'SP110', 'KIAA2018', 'THYN1', 'FIP1L1', 'PTBP1', 'ZNF550', 'LRRFIP1', 'ZNF720', 'ZFP36L2', 'ATF6', 'TRAFD1', 'MZF1', 'ZNF180', 'PAWR', 'ZNF577', 'ZFP62', 'ZNF597', 'NCOA5', 'ZNF260', 'ZNF484', 'ZNF326', 'MTA3', 'MLLT10', 'SRSF4', 'AKAP8L', 'ZNF549', 'CERS6', 'PREB', 'NR3C1', 'ZCCHC11', 'ZBTB44', 'MYNN', 'ZNF713', 'TIGD3', 'ZNF571', 'ZNF383', 'RARG', 'ZNF513', 'ELK4', 'SMAD9', 'JARID2', 'ZNF674', 'ZNF569', 'ZFP14', 'TEAD2', 'E2F2', 'ISX', 'XBP1', 'AKNA', 'DBP', 'ZSCAN12', 'CC2D1A', 'ZNF570', 'ZNF354B', 'CBX2', 'ZBTB49', 'PRDM4', 'SP140L', 'HNRNPH1', 'DMTF1', 'SOX18', 'STAT5B', 'ZC3H4', 'KLF13', 'ZBTB1', 'MATR3', 'SNRNP70', 'NBN', 'ZNF496', 'PPARGC1A', 'TRIM22', 'NPAS2', 'ZNF17', 'CSRNP2', 'FLYWCH1', 'THAP8', 'NR1H2', 'BHLHA15', 'CEBPZ', 'IRF3', 'ZNF207', 'MBNL3', 'POLR2AphosphoS2', 'NFRKB', 'ZNF572', 'THAP7', 'PCBP1', 'MTERF4', 'ZZZ3', 'SRSF1', 'ZNF749', 'TIGD6', 'ZNF451', 'IRX3', 'ZNF615', 'ZSCAN22', 'ZNF343', 'ZNF483', 'ZNF548', 'ZNF781', 'SNAPC4', 'ZNF225', 'ZNF598', 'ZNF230', 'ZNF224', 'GTF2F1', 'TCF3', 'KDM4B', 'ZNF124', 'SMYD3', 'ARID2', 'SNAPC2', 'NCOA1', 'ZNF780A', 'PIN1', 'E2F5', 'ZNF761', 'HBP1', 'TOPORS', 'ZNF883', 'ZNF766', 'BRCA1', 'ELK1', 'RCOR1', 'ZNF576', 'ZNF512B', 'RBAK', 'ZNF512', 'SSRP1', 'RFXANK', 'THAP4', 'ZNF724P', 'AHR', 'ZNF160', 'WIZ', 'ZNF276', 'ZHX3', 'SATB2', 'ZNF333', 'ZNF827', 'MTERF2', 'ZIK1', 'ZNF704', 'ZNF275', 'ZNF616', 'ZNF253', 'ZC3H8', 'ZNF552', 'ZNF546', 'AKAP8', 'AGO1', 'SAFB2', 'ZNF782', 'EED', 'ZNF527', 'ZNF296', 'ZNF138', 'SNAI1', 'ZNF446', 'ZNF790', 'TMF1', 'MLXIP', 'MTF1', 'ZNF30', 'ZNF707', 'ZBTB3', 'KLF11', 'TP53', 'CBX1', 'ZNF181', 'GZF1', 'STAT6', 'TOE1', 'ZNF146', 'ZNF490', 'ZBTB24', 'SREBF2', 'PCBP2', 'PAF1', 'AFF4', 'ZNF775', 'ZKSCAN5', 'MAF1', 'CREB3', 'NFAT5', 'ZFP37', 'DZIP1', 'ZNF485', 'IKZF4', 'ARHGAP35', 'TSC22D2', 'ZNF136', 'GLI4', 'ZNF619', 'ZNF567', 'ETV6', 'ZNF510', 'ZNF33B', 'SMAD1', 'IRF5', 'ZBTB42', 'RELA', 'PHF5A', 'KLF9', 'ZNF235', 'ZNF34', 'DNMT1', 'ZNF678', 'ZNF879', 'FUBP3', 'CHCHD3', 'U2AF2', 'ZBTB46', 'ZFYVE20', 'ZNF367', 'ZNF20', 'ZNF778', 'NFE2L1', 'IRF1', 'ZBED5', 'MYRF', 'EEA1', 'LBX2', 'PBX3', 'KLF15', 'ZNF256', 'ZNF773', 'ZNF83', 'HOXA7', 'DNMT3B', 'U2AF1', 'PAX8', 'ZNF839', 'ZBTB8A', 'ZNF101', 'ATF5', 'CAMTA2', 'ZNF337', 'ZNF318', 'CREBL2', 'ZNF786', 'ZNF44', 'ZFP41', 'ZNF25', 'NFATC3', 'GTF3A', 'CCAR2', 'DDIT3', 'CENPBD1', 'FOXM1', 'REXO4', 'ZNF468', 'JDP2', 'PLSCR1', 'TERF1', 'CENPT', 'ZNF75D', 'ZFAT', 'ZNF557', 'ZSCAN20', 'ZSCAN30', 'SNAPC5', 'ZNF737', 'ZNF234', 'TSC22D1', 'ZNF850', 'ZNF530', 'GLMP', 'ZNF432', 'ZMAT5', 'ZNF382', 'CENPX', 'ZNF282', 'HDAC6', 'ZUFSP', 'GRHL1', 'ARNT', 'MBD1', 'HSF2', 'ZFP64', 'RBM22', 'ZKSCAN1', 'ZNF24', 'ZNF317', 'ZNF205', 'ZNF280D', 'NFKBIZ', 'TCF12', 'ZSCAN29', 'ZSCAN5A', 'ZNF746', 'MED8', 'ESRRA', 'ZNF792', 'ZNF264', 'ZBTB21', 'PRDM10', 'SFPQ', 'SRF', 'MTA1', 'ZNF772', 'ZNF430', 'ZNF221', 'HOXA10', 'ZBTB25', 'NAIF1', 'MED13', 'ZSCAN9', 'RNF219', 'ZBTB39', 'ZNF703', 'RFXAP', 'ZNF607', 'ARNT2', 'MXD4', 'SALL2', 'ZSCAN21', 'ZNF563', 'ZNF407', 'NRL', 'ZNF441', 'RNF2', 'ZNF800', 'ZSCAN25', 'KDM3A', 'ZNF580', 'ZBED4', 'PRDM15', 'ZNF740', 'ZBTB38', 'GLYR1', 'BORCS8', 'NFKB2', 'TAF15', 'ZNF526', 'ZNF33A', 'ZNF841', 'ZC3H13', 'BATF2', 'ZNF608', 'ZNF558', 'ZNF776', 'REL', 'ZNF589', 'HOXA9', 'ZFP82', 'SREBF1', 'KAT7', 'ZNF329', 'CSRNP1', 'ZNF639', 'BRF2', 'RFX5', 'ZNF335', 'ELF4', 'BRD4', 'ZNF414', 'ZFP36L1', 'ZBTB4', 'ZNF784', 'JRK', 'ZNF697', 'ZNF543', 'ZFP90', 'ZNF232', 'KLF6', 'ZNF431', 'ZNF274', 'ZNF280B', 'ZNF547', 'SRY', 'MXD3', 'NR0B2', 'HES4', 'IRF9', 'BCL3', 'KLF12', 'ZNF891', 'FOXC1', 'HSF1', 'PRPF4', 'MTF2', 'ZBTB26', 'ZNF691', 'BAZ2A', 'ZNF142', 'ZMAT3', 'ZNF564', 'ZNF747', 'FBXL19', 'ZNF816', 'ZNF556', 'HOXA5', 'HMGXB3', 'HOXD1', 'MEIS1', 'ZNF511', 'NKX3', 'TEF', 'ZFP1', 'FOXO4', 'ZBTB10', 'ZHX1', 'ZNF670', 'ZNF18', 'ZHX2', 'ZNF709', 'ZSCAN31', 'HIVEP1', 'HINFP', 'SIN3B', 'GMEB2', 'ZNF672', 'FUS', 'ZXDC', 'PRMT3', 'NFIB', 'ZNF562', 'TGIF2', 'CCDC6', 'ZBTB2', 'DR1', 'ZBTB14', 'HNRNPK', 'CHD2', 'SIN3A', 'E2F1', 'MXD1', 'HCFC1', 'ZBTB7A', 'HNRNPL', 'ZNF768', 'H3K4me2', 'ZNF605', 'POGK', 'ZBTB40', 'YEATS2', 'MNT', 'AGO2', 'THAP9', 'KDM5A', 'NFYC', 'NFYB', 'SP2', 'ZNF143', 'SUZ12', 'ADNP', 'RFX3', 'ATF7', 'REPIN1', 'ERF', 'KLF16', 'ZBTB7B', 'ZNF121', 'ZNF629', 'ZMYM2', 'HIC2', 'CBX5', 'CBFB', 'MEIS2', 'ARNTL', 'RAD51', 'NFYA', 'GMEB1', 'H3K9ac', 'DMAP1', 'TBP', 'DRAP1', 'KMT2A', 'KMT2B', 'ZNF350', 'ZNF788', 'SP4', 'PHF20', 'TFDP2', 'KDM5B', 'SMAD7', 'ZNF777', 'USF2', 'TFDP1', 'MAZ', 'YEATS4', 'SPEN', 'MYPOP', 'E2F8', 'SETDB1', 'RORA', 'CHD4', 'MEF2A', 'ATM', 'ZBTB37', 'ZKSCAN8', 'SMAD3', 'GATAD2B', 'NACC2', 'ZBTB34', 'ZNF12', 'ZNF614', 'MIER3', 'MEF2D', 'ZBTB20', 'ZBTB43', 'REST', 'MIER2', 'ZNF3', 'ZNF7', 'HNRNPUL1', 'FOXQ1', 'ZNF878', 'GATA4', 'FOXK2', 'ZNF710', 'MLX', 'ETS1', 'ZNF503', 'DLX6', 'SIX4', 'TBL1XR1', 'HMG20B', 'ONECUT2', 'EHMT2', 'SMARCC2', 'GPBP1L1', 'SCMH1', 'IKZF1', 'PHB2', 'SMARCE1', 'MBD4', 'TCF7', 'HMG20A', 'PROX1', 'H3K4me1', 'GATA2', 'NFIC', 'AHDC1', 'FOXP1', 'DPF2', 'NR2F2', 'ZFHX3', 'TRIM24', 'ONECUT1', 'CEBPA', 'FOXA1', 'FOXA2', 'HLF', 'ZNF384', 'H4K20me1', 'ZFP36', 'MAFK', 'MAFF', 'NFE2L2', 'ZNF644', 'SOX6', 'FOXA3', 'NFIL3', 'SOX13', 'ARID5B', 'GATAD2A', 'CEBPG', 'TEAD1', 'CEBPB', 'ARID3A', 'HNF4A', 'EP300', 'NR5A1', 'BCL6', 'SOX5', 'NCOA2', 'FOXP4', 'RARA', 'PAXIP1', 'TCF7L2', 'NR2F1', 'NFIA', 'LCOR', 'NR2F6', 'TEAD4', 'ZNF652', 'SMAD4', 'FOSL2', 'ZMYM4', 'TFAP4', 'SAP130', 'MAX', 'SALL1', 'TEAD3', 'ZNF687', 'XRCC5', 'NFE2', 'ZEB1', 'IKZF5', 'RBPJ', 'ZNF281', 'ATF2', 'TFE3', 'TARDBP', 'SP5', 'ATF3', 'ZNF460', 'ZMYM3', 'THRA', 'PITX1', 'MYBL2', 'FOXK1', 'IRF2', 'ZNF331', 'PBX2', 'LCORL', 'ETV4', 'ASH2L', 'NR2C2', 'KAT8', 'HNF1B', 'BHLHE40', 'ZNF865', 'NRF1', 'E2F4', 'MXI1', 'GATAD1', 'POLR2AphosphoS5', 'TAF1', 'PHF8', 'USF1', 'UBTF', 'ZNF501', 'MNX1', 'MYC', 'YY1', 'RBM39', 'ZNF770', 'ZNF48', 'ZNF263', 'H2AFZ', 'HNRNPLL', 'H3K4me3', 'RBFOX2', 'POLR2G', 'NONO', 'POLR2A', 'EGR1', 'ARID4A', 'ZFY', 'ZFP91', 'ZNF574', 'ELF1', 'HOXA3', 'THAP11', 'ZBTB33', 'HMGXB4', 'KDM2A', 'GABPB1', 'ZFX', 'PATZ1', 'ARID4B', 'POU2F1', 'LIN54', 'GABPA', 'MGA', 'CREB1', 'RFX1', 'MAFG', 'ATF1', 'PHF21A', 'ETV5', 'RXRA', 'TBX2', 'ZGPAT', 'HDAC1', 'H3K27ac', 'H3K79me2', 'CREM', 'MED1', 'SP1', 'ELF3', 'HDAC2', 'HNF4G', 'KDM1A', 'RXRB', 'ZNF219', 'POGZ', 'ZNF646', 'ZNF292', 'PPARG', 'ZNF217', 'FOXO1', 'THRB', 'CEBPD', 'RREB1', 'MIXL1', 'ZNF362', 'JUN', 'TBX3', 'ISL2', 'ZNF609', 'GFI1', 'HNF1A', 'ATF4', 'JUND', 'SKI', 'SKIL', 'HHEX', 'KDM6A', 'NCOR1', 'FOSL1', 'JUNB', 'RCOR2', 'FOXJ3', 'SIX1', 'RERE', 'HOMEZ', 'KAT2B', 'YBX1', 'H3K27me3', 'H3K36me3', 'EZH2', 'SYNCRIP', 'PLRG1', 'CTCF', 'RAD21', 'STAG1', 'SMC3', 'SRSF9', 'H3K9me3']
    ordered_index = [prot for prot in order if prot in enrichment_df.index]
    print(name)
    print(ordered_index)

    # ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
    # print(ordered_protein_names)

    ordered_log10_enrichment_df = log10_enrichment_df.loc[ordered_index, ordered_index]

    fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), gridspec_kw={'width_ratios': [1, 15]})
    ax_heatmap.set_aspect('equal')
    dendro = dendrogram(linkage_matrix_olo, labels=ordered_log10_enrichment_df.index, orientation='left', ax=ax_dendro)

    ax_dendro.invert_yaxis()  # Reverses the y-axis to make the dendrogram go from top to bottom
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
    plt.savefig(''+name, dpi=300, bbox_inches='tight')
