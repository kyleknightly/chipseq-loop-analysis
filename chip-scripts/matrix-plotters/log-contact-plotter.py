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

enrichment_file_path = '/mnt/altnas/work/Kyle.Knightly/looppi/no-empty-out/trans_enrichment_matrix.tsv'
enrichment_df = pd.read_csv(enrichment_file_path, sep='\t', index_col=0)

# OPTIONAL RESTRICT SET
# subset = ['STAG1', 'SMC3', 'RAD21', 'CTCF', 'ZNF710', 'DLX6', 'RCOR2', 'MED1', 'ZNF503', 'HOMEZ', 'RREB1', 'HNF4G', 'SKI', 'THRA', 'ZNF217', 'KDM1A', 'LCOR', 'ZNF219', 'RXRB', 'JUN', 'MIXL1', 'HDAC2', 'PITX1', 'ZNF292', 'RXRA', 'POGZ', 'ZMYM4', 'TFAP4', 'TEAD3', 'JUND', 'CEBPB', 'HNF4A', 'NR2F6', 'TEAD4', 'CEBPG', 'NFIL3', 'TEAD1', 'PAXIP1', 'SMAD4', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'RARA', 'SOX6', 'FOXA3', 'TCF7L2', 'ARID5B', 'SOX5', 'SOX13', 'FOXP1', 'BCL6', 'PROX1', 'FOXP4', 'ZNF609', 'PPARG', 'FOXO1', 'HNF1A', 'GATAD2A', 'NCOR1', 'GFI1', 'ARID3A', 'FOXJ3', 'ISL2', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC', 'HLF', 'H2AFZ', 'XRCC5', 'GABPA', 'GABPB1', 'NRF1', 'HMGXB4', 'HNRNPLL', 'ZNF687', 'USF1', 'MAX', 'MGA', 'UBTF', 'SAP130', 'EGR1', 'RBFOX2', 'ZFX', 'NR2C2', 'CREB1', 'PATZ1', 'ZNF574', 'THAP11', 'GATAD1', 'PHF20', 'ELF1', 'LIN54', 'NFYA', 'NFYB', 'DRAP1', 'SMAD7', 'ZNF501', 'HOXA3', 'POU2F1', 'NONO', 'TAF1', 'KMT2B', 'KMT2A', 'MYPOP', 'TFDP1', 'ZFP91', 'POLR2A', 'ARID4B', 'ZFY', 'POLR2G', 'PHF8', 'E2F4', 'KAT8', 'GMEB1', 'ARID4A', 'YY1', 'KDM2A', 'AGO2', 'SPEN', 'YEATS4', 'DMAP1', 'SIN3A', 'THAP9', 'YEATS2', 'MAZ', 'RBM39', 'TFDP2', 'SP4', 'MXI1', 'MYC', 'KDM3A', 'REPIN1', 'ATF7', 'MNX1', 'CCDC6', 'ZNF48', 'ZNF788', 'E2F8', 'KLF16', 'ERF', 'MBD1', 'CBFB', 'IKZF5', 'ZNF614', 'MIER3', 'SMAD3', 'MEF2D', 'GATAD2B', 'MEIS2', 'ZNF331', 'RBPJ', 'FOXK1', 'ZGPAT', 'SP5', 'CBX5', 'TARDBP', 'IRF2', 'HDAC1', 'ASH2L', 'ZBTB7B', 'LCORL', 'PHF21A', 'HNF1B', 'ELF3', 'TBX2', 'ETV4', 'ATF2', 'TFE3', 'SALL1', 'ETV5', 'SP1', 'MYBL2', 'CREM', 'TBP', 'MXD1', 'ZHX2', 'RFXAP', 'ZSCAN31', 'ZNF772', 'ZNF221', 'ZNF607', 'SALL2', 'ARNT2', 'KDM5B', 'BORCS8', 'TAF15', 'GLYR1', 'POGK', 'ZNF605', 'ZNF350', 'SFPQ', 'ZBTB25', 'ZNF580', 'ZNF691', 'ZBTB38', 'BAZ2A', 'ZNF556', 'ZNF414', 'HOXA5', 'ZNF511', 'ZNF792', 'MEIS1', 'HOXA10', 'MTA1', 'ZSCAN9', 'ZNF280D', 'BRD4', 'KLF12', 'MXD4', 'NRL', 'KAT7', 'ZNF747', 'ZMAT3', 'ZSCAN21', 'ZNF563', 'FOXC1', 'ZNF407', 'ZNF891', 'NFKB2', 'ZNF709', 'ZNF430', 'ZNF598', 'ZNF230', 'ZNF547', 'ZNF274', 'IRF9', 'ZNF543', 'CSRNP1', 'KLF6', 'ZFP90', 'BCL3', 'ZNF483', 'ZNF883']
# enrichment_df = enrichment_df.loc[
#     enrichment_df.index.intersection(subset),
#     enrichment_df.columns.intersection(subset)
# ]

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

linkage_matrix = linkage(enrichment_df, method='centroid')

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
# true_ticks = [10 ** tick for tick in log_ticks]

# # Format the tick labels
# # For values >= 1, display as integers; for small values, use scientific notation
# formatted_ticks = [
#     f"{t:.2e}" if t < 1 else f"{int(t):,}" for t in true_ticks
# ]

# Update colorbar with formatted labels
colorbar.set_ticks(log_ticks)
# colorbar.set_ticklabels(formatted_ticks)# Format as integers with commas

plt.subplots_adjust(wspace=0.07)
# plt.suptitle('Enrichment Heatmap of TF Overlaps with Clustering (Log10 Scale)')
# plt.xlabel('Transcription Factors')
# plt.ylabel('Transcription Factors')
name = os.path.basename(enrichment_file_path)
plt.savefig('no-empty-' + name[:-3] + 'png', dpi=300, bbox_inches='tight')