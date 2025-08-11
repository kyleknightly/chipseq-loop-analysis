import pandas as pd
import seaborn as sns

looptypes = '/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/hepg2-ccre-loop-types.bedpe'
looptypes_df = pd.read_csv(looptypes, sep='\t', names = ['c1', 's1', 'e1', 't1', 'c2', 's2', 'e2', 't2'])

# Create long-format anchor regions with type
left = looptypes_df[["c1", "s1", "e1", "t1"]].copy()
left.columns = ["chrom", "start", "end", "type"]
right = looptypes_df[["c2", "s2", "e2", "t2"]].copy()
right.columns = ["chrom", "start", "end", "type"]

# Combine and deduplicate
anchortypes_df = pd.concat([left, right], ignore_index=True)
anchortypes_df.drop_duplicates(inplace=True)

# Optional: sort by genomic position for convenience
anchortypes_df = anchortypes_df.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)

# Build coord_key to join with anchors.map
anchortypes_df["coord_key"] = anchortypes_df["chrom"] + "_" + anchortypes_df["start"].astype(str) + "_" + anchortypes_df["end"].astype(str)

# Load anchors.map and join to get anchor_id
anchor_map_df = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchors.map", sep="\t", names=["coord_key", "anchor_id"])
anchor_type_df = anchortypes_df.merge(anchor_map_df, on="coord_key", how="inner")

# Load TF–anchor mapping
anchor_tf_hits_df = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/02_anchor_tf_hits.unique.tsv", sep="\t", names=["anchor_id", "TF"])

# Merge TFs with anchor types
merged = anchor_type_df.merge(anchor_tf_hits_df, on="anchor_id", how="inner")
# Collapse and clean types
type_map = {
    'CA': 'CA', 'CA-CTCF': 'CA', 'CA-H3K4me3': 'CA', 'CA-TF': 'CA',
    'PLS': 'PLS', 'dELS': 'dELS', 'pELS': 'pELS'
}

# Filter out Low-DNase
merged = merged[merged['type'] != 'Low-DNase'].copy()

# Map types to collapsed categories
merged['type'] = merged['type'].map(type_map)

# Drop rows with unmapped types (if any sneaked in)
merged = merged.dropna(subset=['type'])

# Debug: Check what we have before groupby
print("Merged dataframe columns:", merged.columns.tolist())
print("Merged dataframe shape:", merged.shape)
print("Sample merged data:")
print(merged.head())

# Count number of anchor types per TF
print("\nCreating pivot table...")
tf_anchor_counts_grouped = merged.groupby(["TF", "type"]).size()
print("Grouped data sample:")
print(tf_anchor_counts_grouped.head())

print("\nApplying unstack...")
tf_anchor_counts_unstacked = tf_anchor_counts_grouped.unstack(fill_value=0)
print("Unstacked data info:")
print("Type:", type(tf_anchor_counts_unstacked))
print("Index:", tf_anchor_counts_unstacked.index.name)
print("Columns:", tf_anchor_counts_unstacked.columns.tolist())
print("Shape:", tf_anchor_counts_unstacked.shape)

# Fix the naming conflict by renaming the index before reset_index
tf_anchor_counts_unstacked.index.name = 'Protein'
tf_anchor_counts = tf_anchor_counts_unstacked.reset_index()

print("\nTF Anchor Counts:")
print(tf_anchor_counts)

# Show summary statistics
print(f"\nSummary:")
print(f"Number of TFs: {len(tf_anchor_counts)}")
ccre_cols = [col for col in tf_anchor_counts.columns if col != 'Protein']
print(f"cCRE types found: {ccre_cols}")

# === INDIVIDUAL PROTEIN PLOTTING SECTION ===
import matplotlib.pyplot as plt
import numpy as np

# Define your protein order list
protein_order = ['STAG1', 'SMC3', 'RAD21', 'CTCF', 'ZNF710', 'DLX6', 'RCOR2', 'MED1', 'ZNF503', 'HOMEZ', 'RREB1', 'HNF4G', 'SKI', 'THRA', 'ZNF217', 'KDM1A', 'LCOR', 'ZNF219', 'RXRB', 'JUN', 'MIXL1', 'HDAC2', 'PITX1', 'ZNF292', 'RXRA', 'POGZ', 'ZMYM4', 'TFAP4', 'TEAD3', 'JUND', 'CEBPB', 'HNF4A', 'NR2F6', 'TEAD4', 'CEBPG', 'NFIL3', 'TEAD1', 'PAXIP1', 'SMAD4', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'RARA', 'SOX6', 'FOXA3', 'TCF7L2', 'ARID5B', 'SOX5', 'SOX13', 'FOXP1', 'BCL6', 'PROX1', 'FOXP4', 'ZNF609', 'PPARG', 'FOXO1', 'HNF1A', 'GATAD2A', 'NCOR1', 'GFI1', 'ARID3A', 'FOXJ3', 'ISL2', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC', 'HLF', 'H2AFZ', 'XRCC5', 'GABPA', 'GABPB1', 'NRF1', 'HMGXB4', 'HNRNPLL', 'ZNF687', 'USF1', 'MAX', 'MGA', 'UBTF', 'SAP130', 'EGR1', 'RBFOX2', 'ZFX', 'NR2C2', 'CREB1', 'PATZ1', 'ZNF574', 'THAP11', 'GATAD1', 'PHF20', 'ELF1', 'LIN54', 'NFYA', 'NFYB', 'DRAP1', 'SMAD7', 'ZNF501', 'HOXA3', 'POU2F1', 'NONO', 'TAF1', 'KMT2B', 'KMT2A', 'MYPOP', 'TFDP1', 'ZFP91', 'POLR2A', 'ARID4B', 'ZFY', 'POLR2G', 'PHF8', 'E2F4', 'KAT8', 'GMEB1', 'ARID4A', 'YY1', 'KDM2A', 'AGO2', 'SPEN', 'YEATS4', 'DMAP1', 'SIN3A', 'THAP9', 'YEATS2', 'MAZ', 'RBM39', 'TFDP2', 'SP4', 'MXI1', 'MYC', 'KDM3A', 'REPIN1', 'ATF7', 'MNX1', 'CCDC6', 'ZNF48', 'ZNF788', 'E2F8', 'KLF16', 'ERF', 'MBD1', 'CBFB', 'IKZF5', 'ZNF614', 'MIER3', 'SMAD3', 'MEF2D', 'GATAD2B', 'MEIS2', 'ZNF331', 'RBPJ', 'FOXK1', 'ZGPAT', 'SP5', 'CBX5', 'TARDBP', 'IRF2', 'HDAC1', 'ASH2L', 'ZBTB7B', 'LCORL', 'PHF21A', 'HNF1B', 'ELF3', 'TBX2', 'ETV4', 'ATF2', 'TFE3', 'SALL1', 'ETV5', 'SP1', 'MYBL2', 'CREM', 'TBP', 'MXD1', 'ZHX2', 'RFXAP', 'ZSCAN31', 'ZNF772', 'ZNF221', 'ZNF607', 'SALL2', 'ARNT2', 'KDM5B', 'BORCS8', 'TAF15', 'GLYR1', 'POGK', 'ZNF605', 'ZNF350', 'SFPQ', 'ZBTB25', 'ZNF580', 'ZNF691', 'ZBTB38', 'BAZ2A', 'ZNF556', 'ZNF414', 'HOXA5', 'ZNF511', 'ZNF792', 'MEIS1', 'HOXA10', 'MTA1', 'ZSCAN9', 'ZNF280D', 'BRD4', 'KLF12', 'MXD4', 'NRL', 'KAT7', 'ZNF747', 'ZMAT3', 'ZSCAN21', 'ZNF563', 'FOXC1', 'ZNF407', 'ZNF891', 'NFKB2', 'ZNF709', 'ZNF430', 'ZNF598', 'ZNF230', 'ZNF547', 'ZNF274', 'IRF9', 'ZNF543', 'CSRNP1', 'KLF6', 'ZFP90', 'BCL3', 'ZNF483', 'ZNF883']
protein_order.reverse()
# Filter and order the data
available_proteins = [p for p in protein_order if p in tf_anchor_counts['Protein'].values]
print(f"Plotting {len(available_proteins)} proteins: {available_proteins}")

if not available_proteins:
    print("No proteins from your list found in the data!")
    print(f"Available proteins include: {tf_anchor_counts['Protein'].head(10).tolist()}")
else:
    # Get data for the specified proteins in order
    protein_data = tf_anchor_counts[tf_anchor_counts['Protein'].isin(available_proteins)]
    protein_data = protein_data.set_index('Protein').loc[available_proteins].reset_index()

    # Get cCRE columns
    ccre_cols = [col for col in tf_anchor_counts.columns if col != 'Protein']
    print(f"cCRE types: {ccre_cols}")

    # Calculate proportions for each protein
    proportions_data = protein_data[['Protein'] + ccre_cols].copy()
    for ccre_col in ccre_cols:
        proportions_data[ccre_col] = proportions_data[ccre_col] / proportions_data[ccre_cols].sum(axis=1)

    # Replace any NaN values (from proteins with 0 total anchors) with 0
    proportions_data = proportions_data.fillna(0)

    print("\nProportions data:")
    print(proportions_data)

    # Create the plot - make it much wider for many proteins
    fig, ax = plt.subplots(figsize=(len(available_proteins) * 0.3, 4))  # Thinner bars: 0.3 instead of 0.8

    # Set up the bars
    x_spacing = 10  # or 2.5 — must be >= bar_width to avoid overlap
    x_pos = np.arange(len(available_proteins)) * x_spacing

    paired = plt.get_cmap("Paired")

    # Sample evenly spaced colors from it
    colors = [paired(i) for i in [0, 2, 4, 6]]

    # Create stacked bars with thinner width
    bottom = np.zeros(len(available_proteins))
    bars = []
    bar_width = 10  # Make bars thinner

    for i, ccre_type in enumerate(ccre_cols):
        values = proportions_data[ccre_type].values
        bar = ax.bar(x_pos, values, bottom=bottom, label=ccre_type,
                     color=colors[i], alpha=0.9, edgecolor='white', linewidth=0.2, width=bar_width)
        bars.append(bar)
        bottom += values

    # Customize the plot
    ax.set_xlabel('Proteins', fontsize=12, fontweight='bold')
    ax.set_ylabel('Proportion of Loop Anchors', fontsize=12, fontweight='bold')
    ax.set_title(f'cCRE Type Composition by Individual Protein ({len(available_proteins)} proteins)\n(Stacked Proportions)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(available_proteins, rotation=90, ha='center', fontsize=8)  # Vertical rotation, smaller font
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=20, 
          markerscale=2.0, handlelength=2.0, handletextpad=0.8)
    ax.grid(False)
    ax.set_ylim(0, 1.0)

    # Only add percentage labels for larger segments (>5%) and fewer proteins to avoid clutter
    if len(available_proteins) <= 20:  # Only add labels if not too many proteins
        bottom = np.zeros(len(available_proteins))
        for i, ccre_type in enumerate(ccre_cols):
            values = proportions_data[ccre_type].values
            for j, value in enumerate(values):
                if value > 0.05:  # Only label if > 5%
                    ax.text(j, bottom[j] + value/2, f'{value:.0%}', 
                           ha='center', va='center', fontweight='bold', fontsize=6)
            bottom += values

    plt.tight_layout()
    plt.savefig('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/protein_ccre_dists.png', 
                dpi=300, bbox_inches='tight')
    plt.show()

    # Print the numerical data
    print("\n" + "="*80)
    print("INDIVIDUAL PROTEIN cCRE TYPE PROPORTIONS")
    print("="*80)
    print(proportions_data.round(3))

    # Also show absolute counts
    print("\n" + "="*80)
    print("INDIVIDUAL PROTEIN cCRE TYPE ABSOLUTE COUNTS")
    print("="*80)
    abs_counts = protein_data[['Protein'] + ccre_cols]
    print(abs_counts)