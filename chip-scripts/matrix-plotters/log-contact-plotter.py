"""
This plots a heatmap/contact matrix in log10 scale with ward's clustering.

Different clustering algs may easily be switched out 

The highlighting feature may be commented out, it has been useful to track previously
identified protein groups in new arrangements
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering
import os

enrichment_file_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/PrePPI/preppi-matrix.tsv'
enrichment_df = pd.read_csv(enrichment_file_path, sep='\t', index_col=0)

# OPTIONAL RESTRICT SET
subset = ['ZNF547', 'MXD4', 'ZNF563', 'FOXC1', 'ZNF407', 'NRL', 'NFKB2', 'ZNF430', 'GLYR1', 'TAF15', 'BORCS8', 'KDM5B', 'SALL2', 'ARNT2', 'POGK', 'ZNF605', 'ZNF350', 'SFPQ', 'ZNF580', 'ZBTB25', 'ZSCAN31', 'ZNF772', 'ZNF221', 'BAZ2A', 'ZNF556', 'HOXA5', 'HOXA10', 'ZNF792', 'ERF', 'KLF16', 'MEF2D', 'SMAD3', 'MIER3', 'ZNF614', 'IKZF5', 'ZNF331', 'RBPJ', 'FOXK1', 'ZGPAT', 'MYPOP', 'HDAC1', 'ASH2L', 'TARDBP', 'IRF2', 'REPIN1', 'KDM3A', 'E2F8', 'ZNF788', 'ZNF48', 'MNX1', 'CCDC6', 'ATF7', 'DRAP1', 'SMAD7', 'ZNF501', 'NONO', 'HOXA3', 'POU2F1', 'TAF1', 'POLR2A', 'ZFP91', 'TFDP1', 'MYC', 'KMT2B', 'KMT2A', 'AGO2', 'MAZ', 'RBM39', 'TFDP2', 'SP4', 'MXI1', 'YEATS2', 'THAP9', 'SIN3A', 'DMAP1', 'YEATS4', 'SPEN', 'SAP130', 'EGR1', 'ZFX', 'RBFOX2', 'NR2C2', 'CREB1', 'LIN54', 'ELF1', 'THAP11', 'GATAD1', 'PHF20', 'ZNF574', 'KDM2A', 'ARID4B', 'ZFY', 'POLR2G', 'PHF8', 'KAT8', 'E2F4', 'YY1', 'ARID4A', 'GMEB1', 'GABPA', 'GABPB1', 'NRF1', 'PATZ1', 'MAX', 'MGA', 'UBTF', 'USF1', 'ZNF687', 'HNRNPLL', 'HMGXB4', 'XRCC5','ZNF710', 'MIXL1', 'ZNF217', 'KDM1A', 'LCOR', 'ZNF219', 'HDAC2', 'PITX1', 'MED1', 'HOMEZ', 'RXRA', 'TFAP4', 'ZMYM4', 'POGZ', 'ZBTB7B', 'LCORL', 'PHF21A', 'HNF1B', 'ETV4', 'ATF2', 'TBX2', 'ELF3', 'MYBL2', 'SALL1', 'ETV5', 'SP1', 'CREM', 'SKI', 'HNF4G', 'RXRB', 'SMAD4', 'PAXIP1', 'NFIL3', 'TEAD1', 'TEAD4', 'CEBPG', 'NR2F6', 'HNF4A', 'CEBPB', 'JUND', 'TEAD3', 'ISL2', 'FOXJ3', 'FOXP4', 'ZNF609', 'GATAD2A', 'NCOR1', 'PROX1', 'FOXO1', 'HNF1A', 'ARID3A', 'HMG20A', 'AHDC1', 'NCOA2', 'EP300', 'FOXP1', 'ARID5B', 'TCF7L2', 'SOX5', 'SOX13', 'FOXA3', 'SOX6', 'RARA', 'BCL6', 'FOXA2', 'FOXA1', 'CEBPA', 'FOSL2','SMC3', 'STAG1', 'RAD21', 'CTCF']
enrichment_df = enrichment_df.loc[
    enrichment_df.index.intersection(subset),
    enrichment_df.columns.intersection(subset)
]

max = enrichment_df.max().max()
#log, replace 0 with a 1x-10
log10_enrichment_df = np.log10(enrichment_df.replace(0, 10 ** (-np.log10(max))))

#calc vmin vmax values
vmin = log10_enrichment_df.min().min()
vmax = log10_enrichment_df.max().max()
#vmin = -vmax

print(log10_enrichment_df.index.tolist())
print(log10_enrichment_df.columns.tolist())
print(len(log10_enrichment_df.index.tolist()))
print(len(log10_enrichment_df.columns.tolist()))

print(f"vmin: {vmin}, vmax: {vmax}")

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
name = os.path.basename(enrichment_file_path)
plt.savefig('subset-'+name[:-3] + 'png', dpi=300, bbox_inches='tight')