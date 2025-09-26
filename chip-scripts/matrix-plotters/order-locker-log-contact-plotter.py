"""
This allows you to plot a matrix 
in the order of another for the sake of comparison.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering
import os

enrichment_file_path = '/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/hepg2/trans_enrichment_matrix.tsv'
enrichment_df = pd.read_csv(enrichment_file_path, sep='\t', index_col=0)

desired_order=['ZSCAN9', 'HOXA5', 'ZNF511', 'BRD4', 'ZNF280D', 'ZNF598', 'BCL3', 'KLF6', 'ZFP90', 'ZBTB38', 'MBD1', 'ZNF543', 'ZNF547', 'ZNF274', 'ZNF883', 'ZNF230', 'ZNF483', 'KLF12', 'MXD4', 'ZNF430', 'ZNF563', 'FOXC1', 'ZNF891', 'ZNF407', 'NFKB2', 'ZNF772', 'ZSCAN21', 'IRF9', 'CSRNP1', 'ZMAT3', 'ZNF747', 'NRL', 'ZNF691', 'MEIS1', 'ZNF414', 'HOXA10', 'ZNF792', 'KAT7', 'BAZ2A', 'ZNF556', 'ZNF580', 'MTA1', 'ERF', 'IKZF5', 'ZNF614', 'SMAD3', 'MIER3', 'MEF2D', 'GATAD2B', 'MEIS2', 'CBFB', 'ZNF331', 'CBX5', 'KDM3A', 'TAF15', 'ZSCAN31', 'KDM5B', 'SALL2', 'BORCS8', 'GLYR1', 'RFXAP', 'ZNF607', 'ZBTB25', 'SFPQ', 'ZNF709', 'ZNF221', 'TBP', 'ARNT2', 'POGK', 'ZNF605', 'ZNF350', 'KLF16', 'ZNF48', 'E2F8', 'ZNF788', 'ZHX2', 'MXD1', 'AGO2', 'CTCF', 'RAD21', 'STAG1', 'SMC3', 'H2AFZ', 'GABPA', 'GABPB1', 'NRF1', 'HNRNPLL', 'ELF1', 'LIN54', 'CREB1', 'USF1', 'ZNF687', 'PHF20', 'GATAD1', 'THAP11', 'NFYA', 'NFYB', 'HMGXB4', 'ZNF574', 'KAT8', 'ARID4B', 'YY1', 'ARID4A', 'GMEB1', 'NR2C2', 'ZFX', 'ZFY', 'PHF8', 'RBFOX2', 'MAX', 'MGA', 'EGR1', 'SAP130', 'E2F4', 'POLR2G', 'POLR2A', 'KMT2A', 'KMT2B', 'UBTF', 'TFAP4', 'ZBTB7B', 'LCORL', 'FOXK1', 'ZGPAT', 'IRF2', 'RBPJ', 'TARDBP', 'SP5', 'MNX1', 'REPIN1', 'SMAD7', 'ATF7', 'CCDC6', 'SPEN', 'YEATS4', 'MAZ', 'RBM39', 'SIN3A', 'DMAP1', 'XRCC5', 'THAP9', 'MXI1', 'MYC', 'TFDP2', 'ZNF501', 'KDM2A', 'ZFP91', 'TFDP1', 'YEATS2', 'SP4', 'TAF1', 'HDAC1', 'MYPOP', 'NONO', 'HOXA3', 'POU2F1', 'PATZ1', 'SP1', 'ETV5', 'RXRA', 'SALL1', 'ZNF292', 'HDAC2', 'ZNF217', 'THRA', 'MYBL2', 'RREB1', 'ELF3', 'ATF2', 'ASH2L', 'TFE3', 'ETV4', 'TBX2', 'HNF1B', 'PHF21A', 'DRAP1', 'ZNF710', 'DLX6', 'MIXL1', 'PITX1', 'ZNF503', 'MED1', 'RCOR2', 'HOMEZ', 'KDM1A', 'LCOR', 'ZNF219', 'PPARG', 'FOXO1', 'RXRB', 'JUN', 'SKI', 'HNF4G', 'NCOR1', 'ARID3A', 'GFI1', 'GATAD2A', 'HNF1A', 'ZNF609', 'FOXJ3', 'ISL2', 'HLF', 'HNF4A', 'JUND', 'CEBPB', 'PAXIP1', 'NR2F6', 'PROX1', 'CEBPG', 'NFIL3', 'TEAD1', 'TEAD4', 'ZMYM4', 'POGZ', 'TEAD3', 'SMAD4', 'CREM', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'BCL6', 'TCF7L2', 'ARID5B', 'FOXP1', 'FOXP4', 'FOXA3', 'SOX13', 'SOX5', 'SOX6', 'RARA', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC']
enrichment_df = enrichment_df.loc[
    enrichment_df.index.intersection(desired_order),
    enrichment_df.columns.intersection(desired_order)
]
max = enrichment_df.max().max()
#log, replace 0 with a 1x-10
log10_enrichment_df = np.log10(enrichment_df.replace(0, 10 ** (-np.log10(max))))

#calc vmin vmax values
vmin = log10_enrichment_df.min().min()
vmax = log10_enrichment_df.max().max()
#vmin = -vmax

# print(log10_enrichment_df.index.tolist())
# print(log10_enrichment_df.columns.tolist())
# print(len(log10_enrichment_df.index.tolist()))
# print(len(log10_enrichment_df.columns.tolist()))

# print(f"vmin: {vmin}, vmax: {vmax}")

linkage_matrix = linkage(enrichment_df, method='ward')

# Apply Optimal Leaf Ordering to the linkage matrix
linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, enrichment_df)

# Get the ordered indices after optimal leaf ordering
ordered_index = leaves_list(linkage_matrix_olo)

ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
print(ordered_protein_names)

#geq10k cis enr order
# ['ETV4', 'SP5', 'ETV5', 'PATZ1', 'ELF3', 'TFAP4', 'HDAC1', 'HNF1B', 'ZMYM4', 'POGZ', 'ZBTB7B', 'RXRA', 'TFE3', 'ZNF414', 'SFPQ', 'ZNF280D', 'MEF2D', 'ZBTB20', 'ZNF670', 'THRA', 'RXRB', 'RARA', 'SMAD4', 'ZNF217', 'LCORL', 'PHF21A', 'TARDBP', 'TBX2', 'SALL1', 'ZNF219', 'MED1', 'LCOR', 'HDAC2', 'PITX1', 'ZNF710', 'GATAD2B', 'NACC2', 'SMC3', 'RAD21', 'CTCF', 'STAG1', 'GATA2', 'GFI1', 'HOMEZ', 'AHDC1', 'TEAD1', 'SOX5', 'SOX6', 'ARID3A', 'PAXIP1', 'GATAD2A', 'ARID5B', 'FOXA3', 'FOXA2', 'BCL6', 'FOXA1', 'PROX1', 'PPARG', 'HNF1A', 'TEAD4', 'RREB1', 'ZNF292', 'ZNF362', 'NR5A1', 'FOSL2', 'NFIL3', 'HNF4A', 'CEBPA', 'CEBPG', 'HLF', 'NFIA', 'ZNF644', 'CEBPD', 'KDM6A', 'DPF2', 'EP300', 'RCOR2', 'FOXO1', 'MIXL1', 'FOXJ3', 'ISL2', 'NCOA2', 'FOXP4', 'TCF7L2', 'FOXP1', 'HMG20A', 'ZNF609', 'TCF7', 'ZNF503', 'DLX6', 'GATA4', 'RERE', 'HMG20B', 'ZKSCAN8', 'NFKBIZ', 'KDM4B', 'ZNF124', 'ZNF511', 'ZSCAN9', 'SMAD3', 'MIER3', 'ZBTB25', 'ZNF792', 'ZNF580', 'KDM3A', 'CBX5', 'ZNF614', 'ZNF331', 'RBPJ', 'ERF', 'ZGPAT', 'IKZF5', 'IRF2', 'FOXK1', 'ATF2', 'CBFB', 'ZNF205', 'HMGXB3', 'ZSCAN5A', 'FBXL19', 'CSRNP1', 'BORCS8', 'ARID2', 'ZSCAN29', 'ARNTL', 'MEIS1', 'MTA1', 'ZNF691', 'ZNF556', 'ZNF483', 'ZNF230', 'ZNF747', 'ZNF772', 'ZNF788', 'MAZ', 'ZSCAN31', 'CCDC6', 'ZNF350', 'NONO', 'POU2F1', 'MNX1', 'E2F8', 'REPIN1', 'HOXA10', 'ZNF221', 'ZNF605', 'ZNF430', 'ZBED4', 'ZNF784', 'BCL3', 'ZNF547', 'ZNF407', 'ZNF883', 'ZNF274', 'ZNF598', 'HIVEP1', 'TAF1', 'DMAP1', 'YEATS4', 'MYC', 'ZNF48', 'KDM5B', 'MXD1', 'ATF7', 'ARNT2', 'POGK', 'FOXC1', 'ZMAT3', 'ZNF607', 'ZSCAN21', 'ZNF563', 'IRF9', 'SPEN', 'GLYR1', 'BAZ2A', 'ZBTB38', 'MXD4', 'MXD3', 'KLF6', 'KLF12', 'RFXAP', 'MEIS2', 'HOXA5', 'NFKB2', 'SALL2', 'ZFP90', 'HOXD1', 'THAP9', 'YEATS2', 'GATAD1', 'KAT8', 'ELF1', 'PHF20', 'KAT7', 'KMT2A', 'KMT2B', 'TGIF2', 'HBP1', 'SP1', 'TFDP1', 'TFDP2', 'E2F4', 'LIN54', 'MYBL2', 'SP4', 'DRAP1', 'USF1', 'THAP11', 'GABPB1', 'ARID4B', 'SAP130', 'ARID4A', 'KDM2A', 'HOXA3', 'ZNF574', 'EGR1', 'MGA', 'HMGXB4', 'ZFX', 'ZNF687', 'MAX', 'ZFP91', 'ZNF501', 'SMAD7', 'PHF8', 'ZFY', 'UBTF', 'MYPOP', 'GABPA', 'NR2C2', 'CREB1', 'GMEB1', 'NFYA', 'NFYC', 'NFYB']

ordered_log10_enrichment_df = log10_enrichment_df.loc[desired_order, desired_order]

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

plt.subplots_adjust(wspace=0.07)
# plt.suptitle('Enrichment Heatmap of TF Overlaps with Clustering (Log10 Scale)')
# plt.xlabel('Transcription Factors')
# plt.ylabel('Transcription Factors')
name = os.path.basename(enrichment_file_path)
plt.savefig('type-enr-order-' + name[:-3] + 'png', dpi=300, bbox_inches='tight')