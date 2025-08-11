import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data
file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/protein-pearson-scores.tsv'
df = pd.read_csv(file, sep='\t')

# Define your protein order here (replace with your actual list)
protein_order = ['STAG1', 'SMC3', 'RAD21', 'CTCF', 'ZNF710', 'DLX6', 'RCOR2', 'MED1', 'ZNF503', 'HOMEZ', 'RREB1', 'HNF4G', 'SKI', 'THRA', 'ZNF217', 'KDM1A', 'LCOR', 'ZNF219', 'RXRB', 'JUN', 'MIXL1', 'HDAC2', 'PITX1', 'ZNF292', 'RXRA', 'POGZ', 'ZMYM4', 'TFAP4', 'TEAD3', 'JUND', 'CEBPB', 'HNF4A', 'NR2F6', 'TEAD4', 'CEBPG', 'NFIL3', 'TEAD1', 'PAXIP1', 'SMAD4', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'RARA', 'SOX6', 'FOXA3', 'TCF7L2', 'ARID5B', 'SOX5', 'SOX13', 'FOXP1', 'BCL6', 'PROX1', 'FOXP4', 'ZNF609', 'PPARG', 'FOXO1', 'HNF1A', 'GATAD2A', 'NCOR1', 'GFI1', 'ARID3A', 'FOXJ3', 'ISL2', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC', 'HLF', 'H2AFZ', 'XRCC5', 'GABPA', 'GABPB1', 'NRF1', 'HMGXB4', 'HNRNPLL', 'ZNF687', 'USF1', 'MAX', 'MGA', 'UBTF', 'SAP130', 'EGR1', 'RBFOX2', 'ZFX', 'NR2C2', 'CREB1', 'PATZ1', 'ZNF574', 'THAP11', 'GATAD1', 'PHF20', 'ELF1', 'LIN54', 'NFYA', 'NFYB', 'DRAP1', 'SMAD7', 'ZNF501', 'HOXA3', 'POU2F1', 'NONO', 'TAF1', 'KMT2B', 'KMT2A', 'MYPOP', 'TFDP1', 'ZFP91', 'POLR2A', 'ARID4B', 'ZFY', 'POLR2G', 'PHF8', 'E2F4', 'KAT8', 'GMEB1', 'ARID4A', 'YY1', 'KDM2A', 'AGO2', 'SPEN', 'YEATS4', 'DMAP1', 'SIN3A', 'THAP9', 'YEATS2', 'MAZ', 'RBM39', 'TFDP2', 'SP4', 'MXI1', 'MYC', 'KDM3A', 'REPIN1', 'ATF7', 'MNX1', 'CCDC6', 'ZNF48', 'ZNF788', 'E2F8', 'KLF16', 'ERF', 'MBD1', 'CBFB', 'IKZF5', 'ZNF614', 'MIER3', 'SMAD3', 'MEF2D', 'GATAD2B', 'MEIS2', 'ZNF331', 'RBPJ', 'FOXK1', 'ZGPAT', 'SP5', 'CBX5', 'TARDBP', 'IRF2', 'HDAC1', 'ASH2L', 'ZBTB7B', 'LCORL', 'PHF21A', 'HNF1B', 'ELF3', 'TBX2', 'ETV4', 'ATF2', 'TFE3', 'SALL1', 'ETV5', 'SP1', 'MYBL2', 'CREM', 'TBP', 'MXD1', 'ZHX2', 'RFXAP', 'ZSCAN31', 'ZNF772', 'ZNF221', 'ZNF607', 'SALL2', 'ARNT2', 'KDM5B', 'BORCS8', 'TAF15', 'GLYR1', 'POGK', 'ZNF605', 'ZNF350', 'SFPQ', 'ZBTB25', 'ZNF580', 'ZNF691', 'ZBTB38', 'BAZ2A', 'ZNF556', 'ZNF414', 'HOXA5', 'ZNF511', 'ZNF792', 'MEIS1', 'HOXA10', 'MTA1', 'ZSCAN9', 'ZNF280D', 'BRD4', 'KLF12', 'MXD4', 'NRL', 'KAT7', 'ZNF747', 'ZMAT3', 'ZSCAN21', 'ZNF563', 'FOXC1', 'ZNF407', 'ZNF891', 'NFKB2', 'ZNF709', 'ZNF430', 'ZNF598', 'ZNF230', 'ZNF547', 'ZNF274', 'IRF9', 'ZNF543', 'CSRNP1', 'KLF6', 'ZFP90', 'BCL3', 'ZNF483', 'ZNF883']
# If protein_order is provided, filter and order the data
if protein_order:
    # Filter results to only include proteins in your order list
    plot_df = df[df['protein'].isin(protein_order)].copy()
    
    # Set protein as categorical with your specified order
    plot_df['protein'] = pd.Categorical(plot_df['protein'], categories=protein_order, ordered=True)
    plot_df = plot_df.sort_values('protein')
    
    print(f"Plotting {len(plot_df)} proteins from your specified order")
    print(f"Missing from order: {set(protein_order) - set(plot_df['protein'])}")
else:
    # If no order specified, use top 100 by correlation
    plot_df = df.dropna(subset=['pearson_correlation']).head(100).copy()
    print("No protein order specified, plotting top 100 by correlation")

# Create the bar chart
plt.figure(figsize=(40, 4))

# Create bars with color coding
colors = []
for corr in plot_df['pearson_correlation']:
    if pd.isna(corr):
        colors.append('gray')
    elif corr < 0:
        colors.append('red')
    elif corr > 0.5:
        colors.append('darkblue')
    else:
        colors.append('lightblue')

bars = plt.bar(range(len(plot_df)), plot_df['pearson_correlation'], color=colors)

# Formatting
plt.xlabel('Proteins', fontsize=14, fontweight='bold')
plt.ylabel('Pearson Correlation', fontsize=14, fontweight='bold')
plt.title('Protein-wise Pearson Correlations: Enrichment vs Control', fontsize=16, fontweight='bold')

# Set x-axis labels
plt.xticks(range(len(plot_df)), plot_df['protein'], rotation=90, fontsize=5)

# Add horizontal reference lines
plt.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=0.8)
plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, linewidth=0.6)
plt.axhline(y=-0.5, color='gray', linestyle='--', alpha=0.7, linewidth=0.6)

# Adjust layout and grid
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()

# Print some summary stats
print(f"\nSummary Statistics:")
print(f"Mean correlation: {plot_df['pearson_correlation'].mean():.4f}")
print(f"Median correlation: {plot_df['pearson_correlation'].median():.4f}")
print(f"Min correlation: {plot_df['pearson_correlation'].min():.4f}")
print(f"Max correlation: {plot_df['pearson_correlation'].max():.4f}")

# Save the plot
plot_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/prePPI-pearson-correlations.png'
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
plt.savefig(plot_file.replace('.png', '.pdf'), bbox_inches='tight')
print(f"\nPlot saved to: {plot_file}")

plt.show()