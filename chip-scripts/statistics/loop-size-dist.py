import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# Load maps
print("Loading anchor and TF maps...")
anchor_map = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchors.map", sep="\t", header=None, names=["coord", "anchor_id"])
coord_to_anchor = dict(zip(anchor_map["coord"], anchor_map["anchor_id"]))

tf_map = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchor_tf_list.tsv", sep="\t", header=None, names=["anchor_id", "tf_list"])
anchor_to_tfs = {
    row["anchor_id"]: row["tf_list"].split("|")
    for _, row in tf_map.iterrows()
}

print(f"Loaded {len(coord_to_anchor)} anchor coordinates")
print(f"Loaded TF data for {len(anchor_to_tfs)} anchors")

# Function to get anchor ID from loop end
def coord_to_id(chrom, start, end):
    return coord_to_anchor.get(f"{chrom}_{start}_{end}", None)

# Track loop sizes per TF (keep all sizes, not averages)
tf_to_loop_sizes = defaultdict(list)

print("Processing BEDPE file...")
# Parse BEDPE
loop_count = 0
skipped_interchromosomal = 0
with open("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 6:
            continue
            
        chrom1, start1, end1, chrom2, start2, end2 = parts[:6]
        
        # Skip inter-chromosomal loops
        if chrom1 != chrom2:
            skipped_interchromosomal += 1
            continue
            
        a1 = coord_to_id(chrom1, start1, end1)
        a2 = coord_to_id(chrom2, start2, end2)
        
        if not a1 or not a2:
            continue
            
        # Get all TFs at both anchors
        tfs = set(anchor_to_tfs.get(a1, []) + anchor_to_tfs.get(a2, []))
        
        # Calculate loop span (distance between midpoints)
        midpoint1 = (int(start1) + int(end1)) // 2
        midpoint2 = (int(start2) + int(end2)) // 2
        span = abs(midpoint2 - midpoint1)
        
        # Add this loop size to all TFs present at either anchor
        for tf in tfs:
            tf_to_loop_sizes[tf].append(span)
            
        loop_count += 1
        if loop_count % 10000 == 0:
            print(f"Processed {loop_count} loops...")

print(f"Processed {loop_count} intra-chromosomal loops")
print(f"Skipped {skipped_interchromosomal} inter-chromosomal loops")
print(f"Found loop size data for {len(tf_to_loop_sizes)} TFs")

    # Convert to long-format DataFrame for easier plotting
print("Creating DataFrame...")
tf_loop_data = []
for tf, sizes in tf_to_loop_sizes.items():
    for size in sizes:
        # Filter out invalid sizes (0, negative, or very small values)
        if size > 0 and np.isfinite(size):
            tf_loop_data.append({"protein": tf, "size": size, "log10_size": np.log10(size)})

data_df = pd.DataFrame(tf_loop_data)

# Additional validation - remove any remaining NaN or inf values
print(f"Initial DataFrame shape: {data_df.shape}")
data_df = data_df.dropna()
data_df = data_df[np.isfinite(data_df['log10_size'])]
print(f"After removing invalid values: {data_df.shape}")
print(f"Created DataFrame with {len(data_df)} TF-loop size pairs")
print(f"Size range: {data_df['size'].min():,} - {data_df['size'].max():,} bp")

# Define your protein order list (same as original)
protein_order = ['STAG1', 'SMC3', 'RAD21', 'CTCF', 'ZNF710', 'DLX6', 'RCOR2', 'MED1', 'ZNF503', 'HOMEZ', 'RREB1', 'HNF4G', 'SKI', 'THRA', 'ZNF217', 'KDM1A', 'LCOR', 'ZNF219', 'RXRB', 'JUN', 'MIXL1', 'HDAC2', 'PITX1', 'ZNF292', 'RXRA', 'POGZ', 'ZMYM4', 'TFAP4', 'TEAD3', 'JUND', 'CEBPB', 'HNF4A', 'NR2F6', 'TEAD4', 'CEBPG', 'NFIL3', 'TEAD1', 'PAXIP1', 'SMAD4', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'RARA', 'SOX6', 'FOXA3', 'TCF7L2', 'ARID5B', 'SOX5', 'SOX13', 'FOXP1', 'BCL6', 'PROX1', 'FOXP4', 'ZNF609', 'PPARG', 'FOXO1', 'HNF1A', 'GATAD2A', 'NCOR1', 'GFI1', 'ARID3A', 'FOXJ3', 'ISL2', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC', 'HLF', 'H2AFZ', 'XRCC5', 'GABPA', 'GABPB1', 'NRF1', 'HMGXB4', 'HNRNPLL', 'ZNF687', 'USF1', 'MAX', 'MGA', 'UBTF', 'SAP130', 'EGR1', 'RBFOX2', 'ZFX', 'NR2C2', 'CREB1', 'PATZ1', 'ZNF574', 'THAP11', 'GATAD1', 'PHF20', 'ELF1', 'LIN54', 'NFYA', 'NFYB', 'DRAP1', 'SMAD7', 'ZNF501', 'HOXA3', 'POU2F1', 'NONO', 'TAF1', 'KMT2B', 'KMT2A', 'MYPOP', 'TFDP1', 'ZFP91', 'POLR2A', 'ARID4B', 'ZFY', 'POLR2G', 'PHF8', 'E2F4', 'KAT8', 'GMEB1', 'ARID4A', 'YY1', 'KDM2A', 'AGO2', 'SPEN', 'YEATS4', 'DMAP1', 'SIN3A', 'THAP9', 'YEATS2', 'MAZ', 'RBM39', 'TFDP2', 'SP4', 'MXI1', 'MYC', 'KDM3A', 'REPIN1', 'ATF7', 'MNX1', 'CCDC6', 'ZNF48', 'ZNF788', 'E2F8', 'KLF16', 'ERF', 'MBD1', 'CBFB', 'IKZF5', 'ZNF614', 'MIER3', 'SMAD3', 'MEF2D', 'GATAD2B', 'MEIS2', 'ZNF331', 'RBPJ', 'FOXK1', 'ZGPAT', 'SP5', 'CBX5', 'TARDBP', 'IRF2', 'HDAC1', 'ASH2L', 'ZBTB7B', 'LCORL', 'PHF21A', 'HNF1B', 'ELF3', 'TBX2', 'ETV4', 'ATF2', 'TFE3', 'SALL1', 'ETV5', 'SP1', 'MYBL2', 'CREM', 'TBP', 'MXD1', 'ZHX2', 'RFXAP', 'ZSCAN31', 'ZNF772', 'ZNF221', 'ZNF607', 'SALL2', 'ARNT2', 'KDM5B', 'BORCS8', 'TAF15', 'GLYR1', 'POGK', 'ZNF605', 'ZNF350', 'SFPQ', 'ZBTB25', 'ZNF580', 'ZNF691', 'ZBTB38', 'BAZ2A', 'ZNF556', 'ZNF414', 'HOXA5', 'ZNF511', 'ZNF792', 'MEIS1', 'HOXA10', 'MTA1', 'ZSCAN9', 'ZNF280D', 'BRD4', 'KLF12', 'MXD4', 'NRL', 'KAT7', 'ZNF747', 'ZMAT3', 'ZSCAN21', 'ZNF563', 'FOXC1', 'ZNF407', 'ZNF891', 'NFKB2', 'ZNF709', 'ZNF430', 'ZNF598', 'ZNF230', 'ZNF547', 'ZNF274', 'IRF9', 'ZNF543', 'CSRNP1', 'KLF6', 'ZFP90', 'BCL3', 'ZNF483', 'ZNF883']
protein_order.reverse()

# Filter to proteins available in the data
available_proteins = [p for p in protein_order if p in data_df['protein'].values]
print(f"Plotting {len(available_proteins)} proteins from the ordered list")

if not available_proteins:
    print("No proteins from your list found in the data!")
    print(f"Available proteins include: {data_df['protein'].unique()[:10]}")
else:
    # Filter data to only include proteins in our order list
    filtered_data = data_df[data_df['protein'].isin(available_proteins)].copy()
    
    # Create a categorical column to maintain protein order for plotting
    filtered_data['protein_cat'] = pd.Categorical(filtered_data['protein'], 
                                                  categories=available_proteins, 
                                                  ordered=True)
    
    # Sort by the categorical order
    filtered_data = filtered_data.sort_values('protein_cat')
    
    print(f"\nFiltered data summary:")
    print(f"Shape: {filtered_data.shape}")
    print(f"Size range: {filtered_data['size'].min():,} - {filtered_data['size'].max():,} bp")
    print(f"Log10 size range: {filtered_data['log10_size'].min():.2f} - {filtered_data['log10_size'].max():.2f}")
    
    # Print summary statistics per protein (both regular and log10)
    summary_stats = filtered_data.groupby('protein')['size'].agg(['count', 'mean', 'median', 'std']).round(0)
    summary_stats_log = filtered_data.groupby('protein')['log10_size'].agg(['mean', 'median', 'std']).round(3)
    summary_stats_log.columns = ['log10_mean', 'log10_median', 'log10_std']
    summary_combined = pd.concat([summary_stats, summary_stats_log], axis=1)
    print(f"\nTop 10 TFs by loop count:")
    print(summary_stats.sort_values('count', ascending=False).head(10))
    
    
    # Create the violin plot - maintain the long figure format
    fig, ax = plt.subplots(figsize=(len(available_proteins) * 0.3, 6))  # Keep the long format
    # Create violin plot using regular size data
    violin_parts = ax.violinplot([filtered_data[filtered_data['protein'] == protein]['size'].values 
                                  for protein in available_proteins],
                                 positions=range(len(available_proteins)),
                                 widths=0.8,
                                 showmeans=True,
                                 showmedians=True,
                                 showextrema=True)
    # Customize violin plot colors
    colors = sns.color_palette("viridis", n_colors=len(available_proteins))
    for i, pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(0.5)
    
    # Customize other elements
    violin_parts['cmeans'].set_color('red')
    violin_parts['cmeans'].set_linewidth(2)
    violin_parts['cmedians'].set_color('blue')
    violin_parts['cmedians'].set_linewidth(2)
    violin_parts['cbars'].set_color('black')
    violin_parts['cmaxes'].set_color('black')
    violin_parts['cmins'].set_color('black')
    
    # Customize the plot
    ax.set_xlabel('Transcription Factors', fontsize=12, fontweight='bold')
    ax.set_ylabel('Loop Size (bp)', fontsize=12, fontweight='bold')
    ax.set_title(f'Intra-chromosomal Loop Size Distribution by Transcription Factor ({len(available_proteins)} TFs)\n(Violin Plots)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(available_proteins)))
    ax.set_xticklabels(available_proteins, rotation=90, ha='center', fontsize=8)
    
    # Add a custom legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='red', lw=2, label='Mean'),
                       Line2D([0], [0], color='blue', lw=2, label='Median')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Format y-axis to show values in kb/Mb for readability
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}k' if x < 1000000 else f'{x/1000000:.1f}M'))
    ax.grid(True, alpha=0.3, axis='y')
    
    # Optional: Add count annotations above each violin
    for i, protein in enumerate(available_proteins):
        count = len(filtered_data[filtered_data['protein'] == protein])
        ax.text(i, ax.get_ylim()[1] * 0.95, f'n={count}', 
               ha='center', va='top', fontsize=6, rotation=90)
    
    plt.tight_layout()
    plt.savefig('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/tf_loop_size_distributions.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print numerical summary
    print("\n" + "="*80)
    print("TF INTRA-CHROMOSOMAL LOOP SIZE DISTRIBUTION SUMMARY")
    print("="*80)
    summary_detailed = filtered_data.groupby('protein')['size'].agg([
        'count', 'mean', 'median', 'std', 'min', 'max',
        lambda x: x.quantile(0.25),  # Q1
        lambda x: x.quantile(0.75)   # Q3
    ]).round(0)
    summary_detailed.columns = ['Count', 'Mean', 'Median', 'Std', 'Min', 'Max', 'Q1', 'Q3']
    
    # Format large numbers for readability
    for col in ['Mean', 'Median', 'Std', 'Min', 'Max', 'Q1', 'Q3']:
        summary_detailed[col] = summary_detailed[col].apply(lambda x: f'{int(x):,}')
    
    print(summary_detailed)
    
    # Save the full TF-loop size data to CSV for future use
    output_path = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/tf_loop_sizes_full.csv"
    data_df.to_csv(output_path, index=False)
    print(f"\nSaved full TF-loop size data to: {output_path}")
    print(f"This file contains {len(data_df)} TF-loop size pairs for {data_df['protein'].nunique()} TFs")