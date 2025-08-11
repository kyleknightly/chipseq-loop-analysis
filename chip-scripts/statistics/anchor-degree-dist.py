import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter

# Load anchor map
anchor_map = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchors.map", sep="\t", header=None, names=["coord", "anchor_id"])
coord_to_anchor = dict(zip(anchor_map["coord"], anchor_map["anchor_id"]))

# Load TF assignments
tf_map = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchor_tf_list.tsv", sep="\t", header=None, names=["anchor_id", "tf_list"])
anchor_to_tfs = {
    row["anchor_id"]: row["tf_list"].split("|")
    for _, row in tf_map.iterrows()
}

print(f"Loaded {len(coord_to_anchor)} anchors")
print(f"Loaded TFs for {len(anchor_to_tfs)} anchors")

# Function to get anchor ID from chrom/start/end
def coord_to_id(chrom, start, end):
    return coord_to_anchor.get(f"{chrom}_{start}_{end}", None)

# Step 1: Count degrees for each anchor
anchor_degrees = Counter()

with open("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        chrom1, start1, end1, chrom2, start2, end2 = parts[:6]

        # Skip inter-chromosomal loops
        if chrom1 != chrom2:
            continue

        a1 = coord_to_id(chrom1, start1, end1)
        a2 = coord_to_id(chrom2, start2, end2)

        if a1: anchor_degrees[a1] += 1
        if a2: anchor_degrees[a2] += 1

print(f"Counted degrees for {len(anchor_degrees)} anchors")

# Step 2: Assign anchor degrees to TFs
tf_to_anchor_degrees = defaultdict(list)

for anchor, degree in anchor_degrees.items():
    tfs = anchor_to_tfs.get(anchor, [])
    for tf in tfs:
        tf_to_anchor_degrees[tf].append(degree)

# Step 3: Create dataframe for plotting
anchor_df = []

for tf, degrees in tf_to_anchor_degrees.items():
    for d in degrees:
        anchor_df.append({
            'protein': tf,
            'degree': d,
            'log10_degree': np.log10(d) if d > 0 else 0
        })

data_df = pd.DataFrame(anchor_df)
print(f"Created anchor degree DataFrame: {data_df.shape}")

# Drop NaN/infinite
data_df = data_df.dropna()
data_df = data_df[np.isfinite(data_df['log10_degree'])]

# Step 4: Plot as violin plot
protein_order = ['STAG1', 'SMC3', 'RAD21', 'CTCF', 'ZNF710', 'DLX6', 'RCOR2', 'MED1', 'ZNF503', 'HOMEZ', 'RREB1', 'HNF4G', 'SKI', 'THRA', 'ZNF217', 'KDM1A', 'LCOR', 'ZNF219', 'RXRB', 'JUN', 'MIXL1', 'HDAC2', 'PITX1', 'ZNF292', 'RXRA', 'POGZ', 'ZMYM4', 'TFAP4', 'TEAD3', 'JUND', 'CEBPB', 'HNF4A', 'NR2F6', 'TEAD4', 'CEBPG', 'NFIL3', 'TEAD1', 'PAXIP1', 'SMAD4', 'FOSL2', 'CEBPA', 'FOXA1', 'FOXA2', 'RARA', 'SOX6', 'FOXA3', 'TCF7L2', 'ARID5B', 'SOX5', 'SOX13', 'FOXP1', 'BCL6', 'PROX1', 'FOXP4', 'ZNF609', 'PPARG', 'FOXO1', 'HNF1A', 'GATAD2A', 'NCOR1', 'GFI1', 'ARID3A', 'FOXJ3', 'ISL2', 'EP300', 'NCOA2', 'AHDC1', 'HMG20A', 'NFIC', 'HLF', 'H2AFZ', 'XRCC5', 'GABPA', 'GABPB1', 'NRF1', 'HMGXB4', 'HNRNPLL', 'ZNF687', 'USF1', 'MAX', 'MGA', 'UBTF', 'SAP130', 'EGR1', 'RBFOX2', 'ZFX', 'NR2C2', 'CREB1', 'PATZ1', 'ZNF574', 'THAP11', 'GATAD1', 'PHF20', 'ELF1', 'LIN54', 'NFYA', 'NFYB', 'DRAP1', 'SMAD7', 'ZNF501', 'HOXA3', 'POU2F1', 'NONO', 'TAF1', 'KMT2B', 'KMT2A', 'MYPOP', 'TFDP1', 'ZFP91', 'POLR2A', 'ARID4B', 'ZFY', 'POLR2G', 'PHF8', 'E2F4', 'KAT8', 'GMEB1', 'ARID4A', 'YY1', 'KDM2A', 'AGO2', 'SPEN', 'YEATS4', 'DMAP1', 'SIN3A', 'THAP9', 'YEATS2', 'MAZ', 'RBM39', 'TFDP2', 'SP4', 'MXI1', 'MYC', 'KDM3A', 'REPIN1', 'ATF7', 'MNX1', 'CCDC6', 'ZNF48', 'ZNF788', 'E2F8', 'KLF16', 'ERF', 'MBD1', 'CBFB', 'IKZF5', 'ZNF614', 'MIER3', 'SMAD3', 'MEF2D', 'GATAD2B', 'MEIS2', 'ZNF331', 'RBPJ', 'FOXK1', 'ZGPAT', 'SP5', 'CBX5', 'TARDBP', 'IRF2', 'HDAC1', 'ASH2L', 'ZBTB7B', 'LCORL', 'PHF21A', 'HNF1B', 'ELF3', 'TBX2', 'ETV4', 'ATF2', 'TFE3', 'SALL1', 'ETV5', 'SP1', 'MYBL2', 'CREM', 'TBP', 'MXD1', 'ZHX2', 'RFXAP', 'ZSCAN31', 'ZNF772', 'ZNF221', 'ZNF607', 'SALL2', 'ARNT2', 'KDM5B', 'BORCS8', 'TAF15', 'GLYR1', 'POGK', 'ZNF605', 'ZNF350', 'SFPQ', 'ZBTB25', 'ZNF580', 'ZNF691', 'ZBTB38', 'BAZ2A', 'ZNF556', 'ZNF414', 'HOXA5', 'ZNF511', 'ZNF792', 'MEIS1', 'HOXA10', 'MTA1', 'ZSCAN9', 'ZNF280D', 'BRD4', 'KLF12', 'MXD4', 'NRL', 'KAT7', 'ZNF747', 'ZMAT3', 'ZSCAN21', 'ZNF563', 'FOXC1', 'ZNF407', 'ZNF891', 'NFKB2', 'ZNF709', 'ZNF430', 'ZNF598', 'ZNF230', 'ZNF547', 'ZNF274', 'IRF9', 'ZNF543', 'CSRNP1', 'KLF6', 'ZFP90', 'BCL3', 'ZNF483', 'ZNF883']  # use your reversed TF list from earlier
available_proteins = [p for p in protein_order if p in data_df['protein'].values]

filtered_data = data_df[data_df['protein'].isin(available_proteins)].copy()
filtered_data['protein_cat'] = pd.Categorical(filtered_data['protein'], categories=available_proteins, ordered=True)
filtered_data = filtered_data.sort_values('protein_cat')

# Plot
fig, ax = plt.subplots(figsize=(len(available_proteins) * 0.3, 6))
violin_parts = ax.violinplot([filtered_data[filtered_data['protein'] == protein]['degree'].values 
                              for protein in available_proteins],
                             positions=range(len(available_proteins)),
                             widths=0.8,
                             showmeans=True,
                             showmedians=True,
                             showextrema=True)

# colors = sns.color_palette("mako", n_colors=len(available_proteins))
for i, pc in enumerate(violin_parts['bodies']):
    pc.set_facecolor('lightcyan')
    pc.set_alpha(0.7)
    pc.set_edgecolor('black')
    pc.set_linewidth(0.5)

violin_parts['cmeans'].set_color('red')
violin_parts['cmeans'].set_linewidth(5)
violin_parts['cmedians'].set_color('blue')
for part in ['cbars', 'cmaxes', 'cmins']:
    violin_parts[part].set_color('black')

ax.set_xlabel('Transcription Factors', fontsize=12, fontweight='bold')
ax.set_ylabel('Loop Anchor Degree (Loops per Anchor)', fontsize=12, fontweight='bold')
ax.set_title(f'Anchor Degree Distribution by Transcription Factor ({len(available_proteins)} TFs)', 
             fontsize=14, fontweight='bold')
ax.set_xticks(range(len(available_proteins)))
ax.set_xticklabels(available_proteins, rotation=90, ha='center', fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

for i, protein in enumerate(available_proteins):
    count = len(filtered_data[filtered_data['protein'] == protein])
    ax.text(i, ax.get_ylim()[1] * 0.95, f'n={count}', 
           ha='center', va='top', fontsize=6, rotation=90)

plt.tight_layout()
plt.savefig('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/tf_anchor_degree_distributions.png',
            dpi=300, bbox_inches='tight')
plt.show()
