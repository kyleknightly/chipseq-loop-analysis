import pandas as pd
import ast
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.colors import LinearSegmentedColormap
import math

# # Define the colors for the colormap
# colors = [(0, 0, 0), (1, 0, 0)]  # white to red

# # Create the colormap
# white_red = LinearSegmentedColormap.from_list("white_red", colors)

# Load data
file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names=['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])

# Convert string representation of lists to actual tuples
df['prots1'] = df['prots1'].apply(ast.literal_eval).apply(tuple)
df['prots2'] = df['prots2'].apply(ast.literal_eval).apply(tuple)

#prots = list(set(df['prots1'].explode().tolist() + df['prots2'].explode().tolist()))

# List of proteins to plot
prots = ['MYC', 'ATAD3A', 'NR2C2', 'ZBTB34', 'ZNF225', 'DNMT1', 'GABPB1', 'ZNF317', 'ZNF7', 'ZNF646', 'ZNF550', 'MNX1', 'ZNF580', 'MZF1', 'ZBTB33', 'ZSCAN5A', 'MED13', 'MATR3', 'ZNF547', 'ZNF34', 'NRL', 'ZNF674', 'ZBTB2', 'THAP8', 'ATF3', 'MAFG', 'ZBTB39', 'ZNF865', 'ZNF827', 'ELK1', 'ATF2', 'RFX3', 'SIX1', 'ZNF431', 'ZNF160', 'IRX3', 'ETV4', 'MXD1', 'E2F1', 'ZFYVE20', 'NR2F1', 'TSC22D2', 'KDM6A', 'ZNF501', 'SSRP1', 'THAP4', 'HNF1A', 'ISL2', 'NCOA1', 'ZNF264', 'ZNF598', 'MIER2', 'GPN1', 'BRCA1', 'NPAS2', 'HNF4A', 'NFYC', 'POGK', 'FOXA1', 'ZFP41', 'ZSCAN29', 'GFI1', 'ZBTB7B', 'SNAPC5', 'DMAP1', 'ZHX2', 'ZNF383', 'RBAK', 'RFXAP', 'SMAD9', 'TEF', 'ZNF740', 'ZNF720', 'ZNF414', 'ZNF569', 'CENPBD1', 'THRA', 'HHEX', 'MGA', 'TOPORS', 'FOXJ3', 'MTERF4', 'RXRA', 'ZNF644', 'ZFAT', 'KLF9', 'SOX18', 'ZFP91', 'ZNF485', 'SOX6', 'ZNF773', 'YEATS2', 'ZBTB14', 'ERF', 'ZNF639', 'ZBTB24', 'FOXO4', 'MIXL1', 'GZF1', 'SMAD7', 'FOXK1', 'SNAI1', 'ETS1', 'ZNF75D', 'SPEN', 'HES4', 'ZMAT3', 'TARDBP', 'EEA1', 'NFAT5', 'SATB2', 'ZNF180', 'HMGXB4', 'CSRNP2', 'FOSL2', 'DRAP1', 'CREB3', 'GMEB1', 'STAG1', 'TCF12', 'HIC2', 'PROX1', 'TAF1', 'ZNF704', 'USF2', 'SP4', 'ZSCAN31', 'ZNF460', 'ZNF572', 'ARID4B', 'ZNF511', 'UBTF', 'PRMT3', 'TGIF2', 'ZNF335', 'ZNF33A', 'LCOR', 'SMAD1', 'SMAD3', 'ZBTB44', 'ZHX1', 'ZNF607', 'ZFP62', 'ZNF576', 'SP2', 'ZNF557', 'ZNF800', 'ZNF703', 'PPARG', 'HOMEZ', 'PATZ1', 'ZNF780A', 'RARG', 'ARNT2', 'SFPQ', 'ZNF12', 'ZMYM2', 'JUND', 'TBP', 'SAFB2', 'TFDP2', 'STAT6', 'ZNF616', 'ZNF784', 'SP1', 'ZNF512B', 'ZNF571', 'ZNF124', 'ZSCAN20', 'IRF9', 'PAF1', 'BRF2', 'ARNTL', 'ZNF691', 'E2F5', 'GATAD1', 'SMYD3', 'TRAFD1', 'NACC2', 'HLF', 'PRDM4', 'NFATC3', 'ZNF512', 'REXO4', 'E2F2', 'KDM1A', 'HIVEP1', 'ZNF329', 'GATAD2B', 'MYBL2', 'ZNF530', 'ADNP', 'BHLHA15', 'JDP2', 'TBX2', 'ARID2', 'PBX2', 'SALL2', 'AKAP8', 'ZNF556', 'CAMTA2', 'MAF1', 'KDM3A', 'ZNF670', 'ZSCAN9', 'JRK', 'TP53', 'NFYB', 'PRDM10', 'ZSCAN12', 'NKX3', 'FOXP1', 'ZNF574', 'EED', 'RFXANK', 'MYNN', 'KDM5B', 'TCF7', 'KMT2B', 'JARID2', 'TFE3', 'BCL3', 'TEAD1', 'LBX2', 'ZNF234', 'TFDP1', 'NFE2L1', 'ZMYM4', 'ZNF3', 'ATF6', 'ZNF136', 'SNAPC4', 'FOXC1', 'ELF3', 'LIN54', 'ZNF232', 'AHDC1', 'TIGD6', 'MAX', 'POU2F1', 'ZNF407', 'ARID4A', 'ZNF878', 'ZFY', 'ATF1', 'ETV5', 'ZNF33B', 'ZNF597', 'MXD3', 'ZBTB37', 'ZBTB49', 'RCOR2', 'CBX5', 'ZNF619', 'FUBP3', 'ZNF25', 'ZNF217', 'ZNF256', 'ZFP37', 'ZNF558', 'ZIK1', 'ZNF503', 'NFKBIZ', 'HNF1B', 'ZNF761', 'ZNF567', 'HINFP', 'CERS6', 'ZNF362', 'THRB', 'CREM', 'ZNF672', 'TIGD3', 'ZNF879', 'ONECUT2', 'ZSCAN30', 'CHCHD3', 'THAP11', 'ZNF274', 'ZXDC', 'ZZZ3', 'ZNF562', 'ZNF280B', 'NCOA5', 'ZNF746', 'HOXA9', 'TCF3', 'JUN', 'PITX1', 'ZBTB4', 'ZNF44', 'ZNF609', 'CBX2', 'MTERF2', 'ZNF224', 'ZNF263', 'ZBTB7A', 'ZNF713', 'KLF16', 'PHF21A', 'ZNF432', 'KAT8', 'CBFB', 'HDAC1', 'MTF1', 'ZNF841', 'ZNF564', 'ZNF326', 'BAZ2A', 'NFIA', 'TERF1', 'ZNF281', 'ZFP90', 'HOXA7', 'ZNF10', 'ZBTB1', 'NAIF1', 'PBX3', 'ZNF367', 'ZNF687', 'MYPOP', 'NFYA', 'FOXO1', 'ZNF710', 'GABPA', 'ZNF30', 'ZNF101', 'KDM4B', 'ZNF430', 'ZC3H8', 'ZNF697', 'ZBTB26', 'ZBTB8A', 'ZNF18', 'ZNF778', 'ZNF451', 'ZNF280D', 'DR1', 'NFIL3', 'ISX', 'ZFHX3', 'ZNF350', 'RBPJ', 'SP140L', 'FOXP4', 'ZBTB21', 'YEATS4', 'ZNF615', 'ZNF331', 'HMGA1', 'TFAP4', 'NR1H2', 'HSF2', 'EP300', 'ZNF563', 'ZNF527', 'IKZF5', 'DZIP1', 'ZNF121', 'ZNF770', 'KLF6', 'ZNF251', 'SALL1', 'ZNF296', 'FOXA3', 'KIAA2018', 'ZNF768', 'KLF15', 'SRY', 'AHR', 'HMG20B', 'MLXIP', 'DDIT3', 'FOSL1', 'ZNF318', 'REL', 'HDAC2', 'ZNF138', 'ZNF577', 'MED8', 'ZMYM3', 'ZNF775', 'MTF2', 'RARA', 'ZNF548', 'TCF7L2', 'ZBTB10', 'KAT7', 'RAD21', 'KMT2A', 'ZFP64', 'ZNF143', 'AKNA', 'ELF4', 'DPF2', 'ZBTB46', 'IRF1', 'ZNF219', 'MLX', 'ZNF441', 'CREB1', 'ZNF146', 'KDM2A', 'ELF1', 'ZBTB42', 'ZMAT5', 'NR0B2', 'ZBTB40', 'TEAD2', 'RELA', 'MEIS2', 'ZC3H13', 'CEBPD', 'PAX8', 'SMAD4', 'ZNF253', 'DLX6', 'KLF13', 'LRRFIP1', 'ESRRA', 'PREB', 'ZNF382', 'ZNF724P', 'IRF5', 'CENPT', 'GMEB2', 'ZNF337', 'IRF2', 'CENPX', 'HOXA10', 'ZNF614', 'ZBTB25', 'NONO', 'ETV6', 'ARID5B', 'ATRX', 'SKIL', 'TSC22D1', 'GPBP1L1', 'ARHGAP35', 'ZNF772', 'GLMP', 'BCL6', 'ZNF749', 'PHF20', 'THAP7', 'ZNF790', 'NR5A1', 'HNF4G', 'ZNF883', 'ZBTB20', 'PAXIP1', 'CTCF', 'KLF11', 'PHF8', 'TMF1', 'ZNF142', 'ZFP82', 'CSRNP1', 'ATF7', 'ZKSCAN8', 'ZNF490', 'ZFX', 'CC2D1A', 'SP110', 'ZNF608', 'DBP', 'PAWR', 'ZNF839', 'PIN1', 'SRF', 'NFKB2', 'FUBP1', 'GTF3A', 'FBXL19', 'ZNF513', 'RORA', 'ZNF605', 'ZBTB43', 'ZNF275', 'MXD4', 'ZNF552', 'ZNF526', 'ZNF776', 'TRIM24', 'ZNF260', 'ZNF221', 'ZKSCAN5', 'ZBTB3', 'GATA4', 'ZSCAN25', 'BATF2', 'ATF5', 'THYN1', 'RERE', 'ONECUT1', 'TOE1', 'ZFP14', 'ZNF891', 'ZNF48', 'ZNF292', 'HOXA5', 'SAP130', 'ZNF230', 'ZNF276', 'ZNF792', 'ZNF468', 'ZCCHC11', 'PLSCR1', 'ATF4', 'WIZ', 'MED1', 'ZNF83', 'GATAD2A', 'ZSCAN22', 'CEBPA', 'ZNF446', 'USF1', 'RXRB', 'ZUFSP', 'ZNF737', 'ZNF546', 'HMGXB3', 'ZNF747', 'ZNF629', 'ZC3H4', 'JUNB', 'MIER3', 'TEAD4', 'ZNF788', 'ZNF181', 'NRF1', 'ZHX3', 'ZNF709', 'ZFP36L1', 'ZNF484', 'ZSCAN21', 'SNAPC2', 'ZNF850', 'ZBED4', 'REPIN1', 'HOXA3', 'GATA2', 'ZNF816', 'MLLT10', 'ZNF781', 'ZNF343', 'ZNF235', 'NCOA2', 'ZNF766', 'RREB1', 'ZNF782', 'MEF2A', 'ZBED5', 'SOX5', 'NR3C1', 'RNF219', 'ZNF589', 'MYRF', 'CEBPG', 'LCORL', 'STAT5B', 'SIX4', 'ZNF205', 'GRHL1', 'ZNF660', 'BORCS8', 'ZNF707', 'POGZ', 'FOXQ1', 'NFIB', 'HBP1', 'ZNF510', 'ZNF496', 'HMG20A', 'ZNF354B', 'ARID3A', 'KLF12', 'AFF4', 'ZFP1', 'GLI4', 'ZNF483', 'MAZ', 'ZNF549', 'ZGPAT', 'DNMT3B', 'HOXD1', 'ZNF20', 'MTA3', 'FOXA2', 'BRD4', 'XBP1', 'FLYWCH1', 'DMTF1', 'SCMH1', 'E2F8', 'SP5', 'E2F4', 'NFE2', 'SETDB1', 'CREBL2', 'GLYR1', 'ZNF333', 'ZNF786', 'MTA1', 'ZNF678', 'AKAP8L', 'CCDC6', 'ELK4', 'ZFP36L2', 'MEF2D', 'IKZF4', 'ZNF17', 'THAP9', 'ZNF777', 'PRDM15', 'MEIS1', 'ZNF652', 'ZBTB38', 'SMC3', 'EGR1']  # Add more proteins as needed

# Determine the number of rows and columns for subplots (2 columns here)
n_cols = 24
n_rows = math.ceil(len(prots) / n_cols)  # Calculate rows based on number of proteins

# Create subplots with a grid layout
fig, axes = plt.subplots(n_rows, n_cols, figsize=(8*n_cols, 9 * n_rows))

# Flatten axes for easier indexing if it's a grid
axes = axes.flatten()

# Loop over each protein
for i, prot in enumerate(prots):
    print(i)
    # Initialize dictionaries to count upstream and downstream occurrences
    upstream_counts = defaultdict(int)
    downstream_counts = defaultdict(int)

    # Count occurrences for upstream and downstream, filtering for the current protein
    for index, row in df.iterrows():
        anchor_up = (row['chr1'], row['start1'], row['end1'])
        anchor_down = (row['chr2'], row['start2'], row['end2'])
        
        # Count only if the protein is in the respective protein list
        if prot in row['prots1']:
            upstream_counts[anchor_up] += 1
        if prot in row['prots2']:
            downstream_counts[anchor_down] += 1

    # Combine the counts into a single DataFrame for plotting
    combined_counts = pd.DataFrame(
        [(key, upstream_counts[key], downstream_counts[key]) for key in set(upstream_counts) | set(downstream_counts)],
        columns=['anchor', 'upstream_count', 'downstream_count']
    )

    # Count the number of occurrences for each (upstream_count, downstream_count) pair
    combined_counts['coordinate'] = combined_counts.apply(lambda row: (row['upstream_count'], row['downstream_count']), axis=1)
    coordinate_counts = combined_counts['coordinate'].value_counts().reset_index()
    coordinate_counts.columns = ['coordinate', 'count']

    # Separate the coordinate counts into x and y for plotting
    coordinate_counts['x'] = coordinate_counts['coordinate'].apply(lambda x: x[0])
    coordinate_counts['y'] = coordinate_counts['coordinate'].apply(lambda x: x[1])

    # Plot the scatter plot with a colormap based on the count
    scatter = axes[i].scatter(coordinate_counts['x'], coordinate_counts['y'], 
                              c=coordinate_counts['count'], cmap='plasma_r', s=50, edgecolor='none')

    # Set titles and labels
    axes[i].set_aspect('equal', adjustable='box')
    axes[i].set_box_aspect(1)
    axes[i].set_title(f'Upstream vs Downstream Anchors for {prot}')
    axes[i].set_xlabel('Upstream Count')
    axes[i].set_ylabel('Downstream Count')
    axes[i].grid(True)

    # Add a colorbar for each subplot, positioned next to the plot
    cbar = plt.colorbar(scatter, ax=axes[i], fraction=0.046, pad=0.04)
    cbar.set_label('Number of Anchors')

# Remove any unused subplots (if any)
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Adjust layout to prevent overlapping
plt.tight_layout()

# Save the entire figure
plt.savefig('new-mega-up-down-stream-anchors.png')
plt.show()
