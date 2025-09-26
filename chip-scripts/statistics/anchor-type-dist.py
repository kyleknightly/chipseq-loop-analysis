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
anchor_map_df = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/hepg2/anchors.map", sep="\t", names=["coord_key", "anchor_id"])
anchor_type_df = anchortypes_df.merge(anchor_map_df, on="coord_key", how="inner")

# Load TF–anchor mapping
anchor_tf_hits_df = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/hepg2/02_anchor_tf_hits.unique.tsv", sep="\t", names=["anchor_id", "TF"])

# Merge TFs with anchor types
merged = anchor_type_df.merge(anchor_tf_hits_df, on="anchor_id", how="inner")

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
protein_order = ['ZSCAN9', 'HOXA5', 'ZNF511', 'BRD4', 'ZNF280D', 'ZNF598', 'BCL3', 'KLF6', 'ZFP90', 'ZBTB38', 'MBD1', 'ZNF543', 'ZNF547', 'ZNF274', 'ZNF883', 'ZNF230', 'ZNF483', 'KLF12', 'MXD4', 'ZNF430', 'ZNF563', 'FOXC1', 'ZNF891', 'ZNF407', 'NFKB2', 'ZNF772', 'ZSCAN21', 'IRF9', 'CSRNP1', 'ZMAT3', 'ZNF747', 'NRL', 'ZNF691', 'MEIS1', 'ZNF414', 'HOXA10', 'ZNF792', 'KAT7', 'BAZ2A', 'ZNF556', 'ZNF580', 'MTA1', 'ERF', 'IKZF5', 'ZNF614', 'SMAD3', 'MIER3', 'MEF2D', 'GATAD2B', 'MEIS2', 'CBFB', 'ZNF331', 'CBX5', 'KDM3A', 'TAF15', 'ZSCAN31', 'KDM5B', 'SALL2', 'BORCS8', 'GLYR1', 'RFXAP', 'ZNF607', 'ZBTB25', 'SFPQ', 'ZNF709', 'ZNF221', 'TBP', 'ARNT2', 'POGK', 'ZNF605', 'ZNF350', 'KLF16', 'ZNF48', 'E2F8', 'ZNF788', 'ZHX2', 'MXD1', 'AGO2', 'CTCF', 'RAD21', 'STAG1', 'SMC3', 'H2AFZ', 'GABPA', 'GABPB1', 'NRF1', 'HNRNPLL', 'ELF1', 'LIN54', 'CREB1', 'USF1', 'ZNF687', 'PHF20', 'GATAD1', 'THAP11', 'NFYA', 'NFYB', 'HMGXB4', 'ZNF574', 'KAT8', 'ARID4B', 'YY1', 'ARID4A', 'GMEB1', 'NR2C2', 'ZFX', 'ZFY', 'PHF8', 'RBFOX2', 'MAX', 'MGA', 'EGR1', 'SAP130', 'E2F4', 'POLR2G', 'POLR2A', 'KMT2A', 'KMT2B', 'UBTF', 'TFAP4', 'ZBTB7B', 'LCORL', 'FOXK1', 'ZGPAT', 'IRF2', 'RBPJ', 'TARDBP', 'SP5', 'MNX1', 'REPIN1', 'SMAD7', 'ATF7', 'CCDC6', 'SPEN', 'YEATS4', 'MAZ', 'RBM39', 'SIN3A', 'DMAP1', 'XRCC5', 'THAP9', 'MXI1', 'MYC', 'TFDP2', 'ZNF501', 'KDM2A', 'ZFP91', 'TFDP1', 'YEATS2', 'SP4', 'TAF1', 'HDAC1', 'MYPOP', 'NONO', 'HOXA3', 'POU2F1', 'PATZ1', 'SP1', 'ETV5', 'RXRA', 'SALL1', 'ZNF292', 'HDAC2', 'ZNF217', 'THRA', 'MYBL2', 'RREB1', 'ELF3', 'ATF2', 'ASH2L', 'TFE3', 'ETV4', 'TBX2', 'HNF1B', 'PHF21A', 'DRAP1', 'ZNF710', 'DLX6', 'MIXL1', 'PITX1', 'ZNF503', 'MED1', 'RCOR2', 'HOMEZ', 'KDM1A', 'LCOR', 'ZNF219', 'PPARG', 'FOXO1', 'RXRB', 'JUN', 'SKI', 'HNF4G', 'NCOR1', 'ARID3A', 'GFI1', 'GATAD2A', 'HNF1A', 'ZNF609', 'FOXJ3', 'ISL2', 'HLF', 'HNF4A', 'JUND', 'CEBPB', 'PAXIP1', 'NR2F6', 'PROX1', 'CEBPG', 'NFIL3', 'TEAD1', 'TEAD4', 'ZMYM4', 'POGZ', 'TEAD3', 'SMAD4', 'CREM', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'BCL6', 'TCF7L2', 'ARID5B', 'FOXP1', 'FOXP4', 'FOXA3', 'SOX13', 'SOX5', 'SOX6', 'RARA', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC']
# protein_order.reverse()
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

    colors = sns.color_palette("Paired", n_colors=len(ccre_cols))

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
    plt.savefig('type-enr-ordered-protein_ccre_dists.png', 
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