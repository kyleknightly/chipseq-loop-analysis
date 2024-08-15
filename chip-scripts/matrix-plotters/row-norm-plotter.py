"""
This normalizes each row of a matrix so that they sum to the same value,
and then plots it.

Fairly old, needs update
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import os

enrichment_file_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/same-anchor-interactions/anchor-contact-matrix.tsv'
df = pd.read_csv(enrichment_file_path, sep='\t', index_col=0)
#log, replace 0 with a 1x-10
# Define the target sum for each row
target_sum = 1  # You can change this value to your desired target sum

# Calculate the sum of each row
row_sums = df.sum(axis=1)

# Normalize each row to sum to 1
normalized_df = df.div(row_sums, axis=0)

# Scale each row to sum to the target value
normalized_df = normalized_df * target_sum

max = normalized_df.max().max()
min = normalized_df.min().min()

linkage_matrix = linkage(normalized_df, method='average')

ordered_index = leaves_list(linkage_matrix)

ordered_normalized_df = normalized_df.iloc[ordered_index, ordered_index]

fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 45), gridspec_kw={'width_ratios': [1, 15]})

dendro = dendrogram(linkage_matrix, labels=ordered_normalized_df.index, orientation='left', ax=ax_dendro)
ax_dendro.set_xticks([])
ax_dendro.set_yticks([])

heatmap = sns.heatmap(
    ordered_normalized_df, 
    cmap='RdBu_r', # "Spectral_r",    
    annot=False,      
    linewidths=0.05,   
    cbar_kws={'label': 'Log10 Enrichment Score'}, 
    center=0,  
    ax=ax_heatmap,
    xticklabels=True,
    yticklabels=True,
    vmin=min,  
    vmax=max   
)

#ax_heatmap.yaxis.tick_right()
#ax_heatmap.yaxis.set_label_position('right')
heatmap.set_xticks(np.arange(len(ordered_normalized_df.columns)) + 0.5)
heatmap.set_yticks(np.arange(len(ordered_normalized_df.index)) + 0.5)
heatmap.set_xticklabels(ordered_normalized_df.columns, rotation=90, fontsize=6)
heatmap.set_yticklabels(ordered_normalized_df.index, rotation=0, fontsize=6)

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
name = os.path.basename(enrichment_file_path)
plt.savefig('rownormed-' + name[:-3] + 'png', dpi=300, bbox_inches='tight')
