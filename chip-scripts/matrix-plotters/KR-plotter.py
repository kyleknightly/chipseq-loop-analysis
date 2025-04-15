"""
This normalizes a graph with KR normalization and plots it
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering
import os

file_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/loop-end-enrichments-pseudo-trans-contacts.tsv'
df = pd.read_csv(file_path, sep='\t', index_col=0)
max = df.max().max()

def kr_normalize(matrix, tol=1e-6, max_iters=1000):
    # Initialize the matrix as a numpy array
    A = matrix.values.astype(float)
    
    # Get the number of rows and columns
    n = A.shape[0]
    
    # Initialize row and column scalings
    r = np.ones(n)
    c = np.ones(n)
    
    # Perform KR normalization iteratively
    for iteration in range(max_iters):
        # Normalize rows
        r_new = 1 / np.dot(A, c)
        
        # Normalize columns
        c_new = 1 / np.dot(A.T, r_new)
        
        # Check for convergence
        if np.max(np.abs(r_new - r)) < tol and np.max(np.abs(c_new - c)) < tol:
            print(f"Converged in {iteration + 1} iterations.")
            break
        
        # Update row and column scalings
        r = r_new
        c = c_new
    else:
        print("KR normalization did not converge within the maximum number of iterations.")
    
    # Apply the scaling factors to the matrix
    D_r = np.diag(r)
    D_c = np.diag(c)
    
    normalized_matrix = np.dot(np.dot(D_r, A), D_c)
    
    # Return the result as a pandas DataFrame with the same index and columns as the original
    return pd.DataFrame(normalized_matrix, index=matrix.index, columns=matrix.columns)

df_normalized = kr_normalize(df)

#calc vmin vmax values
vmin = df_normalized.min().min()
vmax = df_normalized.max().max()

linkage_matrix = linkage(df_normalized, method='ward')

# Apply Optimal Leaf Ordering to the linkage matrix
linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, df_normalized)

# Get the ordered indices after optimal leaf ordering
ordered_index = leaves_list(linkage_matrix_olo)

# ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
# print(ordered_protein_names)

ordered_df_normalized = df_normalized.iloc[ordered_index, ordered_index]

fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), gridspec_kw={'width_ratios': [1, 15]})
ax_heatmap.set_aspect('equal')
dendro = dendrogram(linkage_matrix, labels=ordered_df_normalized.index, orientation='left', ax=ax_dendro)
ax_dendro.invert_yaxis()  # Reverse the y-axis to make the dendrogram go from top to bottom
ax_dendro.set_xticks([])
ax_dendro.set_yticks([])

heatmap = sns.heatmap(
    ordered_df_normalized, 
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
yticks = np.arange(len(ordered_df_normalized.index)) + 0.5
yticklabels = ordered_df_normalized.index
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

heatmap.set_xticks(np.arange(len(ordered_df_normalized.columns)) + 0.5)
heatmap.set_yticks(np.arange(len(ordered_df_normalized.index)) + 0.5)
heatmap.set_xticklabels(ordered_df_normalized.columns, rotation=90, fontsize=6)
heatmap.set_yticklabels(ordered_df_normalized.index, rotation=0, fontsize=9)

# Access the colorbar
colorbar = heatmap.collections[0].colorbar
# Set the font size for the colorbar label
colorbar.ax.yaxis.label.set_size(30)  # Set the desired font size
# Optional: set the font size for the colorbar ticks
colorbar.ax.tick_params(labelsize=30)

plt.subplots_adjust(wspace=0.07)
# plt.suptitle('Enrichment Heatmap of TF Overlaps with Clustering (Log10 Scale)')
# plt.xlabel('Transcription Factors')
# plt.ylabel('Transcription Factors')
name = os.path.basename(file_path)
plt.savefig('KR-' + name[:-3] + 'png', dpi=300, bbox_inches='tight')