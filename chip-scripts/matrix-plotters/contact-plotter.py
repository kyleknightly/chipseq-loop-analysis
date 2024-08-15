"""
This plots a heatmap/contact matrix without log scale, usually for use in pure contacts, and with Ward's clustering
The highlighting feature may be commented out, it has been useful to track previously
identified protein groups in new arrangements
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import os

file_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/loop-ends/trans-probability-matrix.tsv'
df = pd.read_csv(file_path, sep='\t', index_col=0)
max = df.max().max()

#calc vmin vmax values
vmin = df.min().min()
vmax = df.max().max()
#vmin = -vmax
print(f"vmin: {vmin}, vmax: {vmax}")

linkage_matrix = linkage(df, method='ward')

ordered_index = leaves_list(linkage_matrix)

# ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
# print(ordered_protein_names)

ordered_df = df.iloc[ordered_index, ordered_index]

fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), gridspec_kw={'width_ratios': [1, 15]})

dendro = dendrogram(linkage_matrix, labels=ordered_df.index, orientation='left', ax=ax_dendro)
ax_dendro.invert_yaxis()  # Reverse the y-axis to make the dendrogram go from top to bottom
ax_dendro.set_xticks([])
ax_dendro.set_yticks([])

heatmap = sns.heatmap(
    ordered_df, 
    cmap='RdBu_r', # 'rdbu_r' "Spectral_r",    
    annot=False,      
    linewidths=0.05,   
    cbar_kws={'label': 'Log10 Enrichment Score'}, 
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
yticks = np.arange(len(ordered_df.index)) + 0.5
yticklabels = ordered_df.index
for i, label in enumerate(ax_heatmap.get_yticklabels()):
    if i % 2 == 0:
        label.set_x(-0.0005)  # Shift to the left
    else:
        label.set_x(-0.015)  # Shift further to the left

# Adjust tick lengths for staggered labels
ax_heatmap.yaxis.set_tick_params(which='both', length=0)
for i, tick in enumerate(ax_heatmap.yaxis.get_major_ticks()):
    if i % 2 == 0:
        tick.tick1line.set_markersize(2.75)  # Set longer tick length
    else:
        tick.tick1line.set_markersize(45)  # Set shorter tick length

#ax_heatmap.yaxis.tick_right()
#ax_heatmap.yaxis.set_label_position('right')
heatmap.set_xticks(np.arange(len(ordered_df.columns)) + 0.5)
heatmap.set_yticks(np.arange(len(ordered_df.index)) + 0.5)
heatmap.set_xticklabels(ordered_df.columns, rotation=90, fontsize=6)
heatmap.set_yticklabels(ordered_df.index, rotation=0, fontsize=9)

# Access the colorbar
colorbar = heatmap.collections[0].colorbar
# Set the font size for the colorbar label
colorbar.ax.yaxis.label.set_size(30)  # Set the desired font size
# Optional: set the font size for the colorbar ticks
colorbar.ax.tick_params(labelsize=30)

plt.subplots_adjust(wspace=0.05)
# plt.suptitle('Enrichment Heatmap of TF Overlaps with Clustering (Log10 Scale)')
# plt.xlabel('Transcription Factors')
# plt.ylabel('Transcription Factors')
name = os.path.basename(file_path)
plt.savefig('ward-' + name[:-3] + 'png', dpi=300, bbox_inches='tight')
