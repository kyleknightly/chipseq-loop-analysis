"""
This plots a heatmap/contact matrix without log scale, usually for use in pure contacts, and with Ward's clustering
The highlighting feature may be commented out, it has been useful to track previously
identified protein groups in new arrangements
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering
import os

file_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/PrePPI/preppi-matrix.tsv'
df = pd.read_csv(file_path, sep='\t', index_col=0)
max = df.max().max()
# print("NaNs:\n", df[df.isna().any(axis=1)])

#log, replace 0 with a 1x-10
df = np.log10(df.replace(0, 10 ** (-np.log10(max)))) # OPTIONAL LOG VALUES
df.fillna(0, inplace=True)



#calc vmin vmax values
vmin = df.min().min()
vmax = df.max().max()
#vmin = -vmax
print(f"vmin: {vmin}, vmax: {vmax}")

linkage_matrix = linkage(df, method='ward')
# Apply Optimal Leaf Ordering to the linkage matrix
linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, df)

# Get the ordered indices after optimal leaf ordering
ordered_index = leaves_list(linkage_matrix_olo)
# ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
# print(ordered_protein_names)

desired_order_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/pipeline-running/geq50k_trans_enrichment_matrix'
desired_order_df = pd.read_csv(desired_order_file, sep='\t', index_col=0)
desired_order_linkage_matrix = linkage(desired_order_df, method='ward')

# Apply Optimal Leaf Ordering to the linkage matrix
desired_order_linkage_matrix_olo = optimal_leaf_ordering(desired_order_linkage_matrix, desired_order_df)

# Get the ordered indices after optimal leaf ordering
desired_order = desired_order_df.index[leaves_list(desired_order_linkage_matrix_olo)].tolist()
desired_order = [p for p in desired_order if p in df.index]
print(desired_order)
print('=======')
# desired_order=['GABPA', 'THAP11', 'ELF1', 'GABPB1', 'GMEB1', 'YEATS2', 'KAT8', 'PHF20', 'GATAD1', 'NFYC', 'NFYB', 'NFYA', 'SP1', 'SP4', 'DRAP1', 'TAF1', 'ZBTB25', 'MXD3', 'KLF6', 'YEATS4', 'KLF12', 'SPEN', 'MXD1', 'DMAP1', 'MYC', 'MXD4', 'HOXA3', 'THAP9', 'ZNF574', 'KDM2A', 'ARID4B', 'SAP130', 'ARID4A', 'PHF8', 'ZFY', 'ZFX', 'E2F4', 'TFDP2', 'TFDP1', 'LIN54', 'NR2C2', 'HMGXB3', 'KMT2B', 'KMT2A', 'ZBTB38', 'MTA1', 'ZNF598', 'ZFP91', 'ZNF350', 'NFKB2', 'KDM5B', 'MAZ', 'GLYR1', 'IRF9', 'ZNF274', 'ZNF883', 'ZNF547', 'ZNF407', 'HBP1', 'RFXAP', 'SALL2', 'ZFP90', 'HIVEP1', 'HMGXB4', 'ZNF691', 'ZNF230', 'ZNF501', 'SMAD7', 'MYPOP', 'UBTF', 'ZNF124', 'KDM3A', 'CBX5', 'KDM4B', 'ZNF792', 'ZNF511', 'KAT7', 'BCL3', 'ARNT2', 'ATF7', 'ZSCAN31', 'ZSCAN21', 'ZNF563', 'ZNF607', 'BORCS8', 'CSRNP1', 'BAZ2A', 'ZNF784', 'ZNF605', 'FOXC1', 'ZNF556', 'MNX1', 'MGA', 'NONO', 'POGK', 'TGIF2', 'POU2F1', 'ZNF483', 'ZNF580', 'ZNF772', 'ZNF747', 'ZNF430', 'ZMAT3', 'HOXA5', 'IRF2', 'ERF', 'USF1', 'ZNF788', 'EGR1', 'ZNF687', 'MAX', 'ZNF48', 'ZBED4', 'PATZ1', 'ZBTB7B', 'RXRA', 'ARNTL', 'CBFB', 'ZSCAN9', 'ZNF331', 'ZNF614', 'ZNF205', 'MEIS2', 'MEIS1', 'HDAC1', 'TFAP4', 'ZNF221', 'ZNF280D', 'ARID2', 'ZNF670', 'REPIN1', 'HOXA10', 'FBXL19', 'E2F8', 'CCDC6', 'ATF2', 'TFE3', 'HOXD1', 'MIER3', 'SMAD3', 'PHF21A', 'HNF1B', 'ZSCAN29', 'FOXK1', 'ETV5', 'ETV4', 'IKZF5', 'RBPJ', 'SP5', 'ZGPAT', 'MYBL2', 'CREB1', 'AHDC1', 'TCF7', 'FOXJ3', 'HMG20A', 'RCOR2', 'EP300', 'NCOA2', 'ISL2', 'GATA4', 'FOXP1', 'SOX5', 'SOX6', 'ARID5B', 'FOXA3', 'GATAD2A', 'BCL6', 'DPF2', 'FOXP4', 'TCF7L2', 'ZNF609', 'PROX1', 'ARID3A', 'RERE', 'GATA2', 'GFI1', 'ZNF644', 'HNF4A', 'TEAD1', 'MIXL1', 'FOXA1', 'FOXA2', 'CEBPA', 'NFIL3', 'PPARG', 'FOXO1', 'RARA', 'HMG20B', 'HNF1A', 'HOMEZ', 'GATAD2B', 'NACC2', 'ZNF710', 'PAXIP1', 'DLX6', 'ZNF503', 'MED1', 'NR5A1', 'CEBPG', 'LCOR', 'HDAC2', 'RXRB', 'KDM6A', 'SMAD4', 'THRA', 'LCORL', 'ZNF217', 'FOSL2', 'CEBPD', 'NFKBIZ', 'ZKSCAN8', 'ZNF219', 'PITX1', 'ZNF414', 'ZBTB20', 'SALL1', 'TARDBP', 'TBX2', 'SFPQ', 'MEF2D', 'ELF3', 'ZSCAN5A', 'POGZ', 'ZMYM4', 'ZNF292', 'RREB1', 'NFIA', 'HLF', 'ZNF362', 'TEAD4', 'SMC3', 'STAG1', 'RAD21', 'CTCF']
remaining = [p for p in df.index if p not in desired_order]
print(remaining)


ordered_df = df.loc[desired_order, desired_order]
# print(ordered_df.index[ordered_index].tolist())

fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), gridspec_kw={'width_ratios': [1, 15]})

# subset_df = desired_order_df.loc[ordered_df.index]  # or use df.loc[common]
# subset_linkage = linkage(subset_df, method='ward')
# subset_linkage_olo = optimal_leaf_ordering(subset_linkage, subset_df)

# dendro = dendrogram(subset_linkage_olo, labels=ordered_df.index, orientation='left', ax=ax_dendro)

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
# print(ordered_df.index.tolist())
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
plt.savefig('trans_geq50k_locked_' + name[:-3] + 'png', dpi=300, bbox_inches='tight')
