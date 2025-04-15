"""
This uses  new column-organized labeling

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

enrichment_file_path = '/mnt/altnas/work/Kyle.Knightly/enrichments-pseudo-anchor-contacts.tsv'
enrichment_df = pd.read_csv(enrichment_file_path, sep='\t', index_col=0)
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

# ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
# print(ordered_protein_names)

ordered_log10_enrichment_df = log10_enrichment_df.iloc[ordered_index, ordered_index]

fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(80, 50), gridspec_kw={'width_ratios': [1, 15]})

dendro = dendrogram(linkage_matrix, labels=ordered_log10_enrichment_df.index, orientation='left', ax=ax_dendro)
ax_dendro.invert_yaxis()  # Reverse the y-axis to make the dendrogram go from top to bottom
ax_dendro.set_xticks([])
ax_dendro.set_yticks([])

heatmap = sns.heatmap(
    ordered_log10_enrichment_df, 
    cmap='RdBu_r', # 'rdbu_r' "Spectral_r",    
    annot=False,      
    linewidths=0.05,   
    cbar_kws={'label': 'Log10 Enrichment Score'}, 
    center=0,  
    ax=ax_heatmap,
    xticklabels=True,
    yticklabels=False,
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

# ... previous code ...

# Set y-ticks for every row without labels
ytick_positions = list(range(len(ordered_log10_enrichment_df.index)))
adjusted_ytick_positions = [y + 0.5 for y in ytick_positions]


ax_heatmap.set_yticks(adjusted_ytick_positions)

ax_heatmap.set_yticklabels([''] * len(adjusted_ytick_positions))  # Remove default labels

# Adjust y-tick marks to match the manually placed labels
for i, tick in enumerate(ax_heatmap.yaxis.get_major_ticks()):
    if i % 3 == 0:
        tick.set_pad(30)  # Move tick mark to the right a lot
    elif i % 3 == 1:
        tick.set_pad(15)  # Move tick mark to the right a bit
    # Third tick stays at the default position (no need to set pad)

# Collect the labels in their original order
labels = ordered_log10_enrichment_df.index.tolist()

# Add custom labels for every third tick, combining three labels into one
for i in range(0, len(adjusted_ytick_positions), 3):
    if i + 2 < len(adjusted_ytick_positions):
        ax_heatmap.text(-0.5, i+1.75, labels[i], va='center', ha='right', fontsize=14)
        ax_heatmap.text(-20.5, i+1.75, labels[i + 1], va='center', ha='right', fontsize=14)
        ax_heatmap.text(-40.5, i+1.75, labels[i + 2], va='center', ha='right', fontsize=14)

# Access the colorbar
colorbar = heatmap.collections[0].colorbar
colorbar.ax.yaxis.label.set_size(30)  # Set the desired font size
colorbar.ax.tick_params(labelsize=30)

plt.subplots_adjust(wspace=0.5)
name = os.path.basename(enrichment_file_path)
plt.savefig('test-log10' + name[:-3] + 'png', dpi=300, bbox_inches='tight')