import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === File paths ===
enrichment_matrix_file = '/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/trans_enrichment_matrix.tsv'
control_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/PrePPI/fulllist_totalLR_all.csv'

# === Load data ===
enrichment_matrix = pd.read_csv(enrichment_matrix_file, sep='\t', index_col=0)
control_df = pd.read_csv(control_file)

# === Standardize IDs ===
enrichment_matrix.index = enrichment_matrix.index.astype(str)
enrichment_matrix.columns = enrichment_matrix.columns.astype(str)
control_df['geneID1'] = control_df['geneID1'].astype(str)
control_df['geneID2'] = control_df['geneID2'].astype(str)

# === Debug: preview enrichment matrix and control file ===
print("🔍 Enrichment matrix shape:", enrichment_matrix.shape)
print("🔍 Enrichment sample:\n", enrichment_matrix.iloc[:5, :5])
print("🔍 Control CSV sample:\n", control_df.head())

# === Shared proteins ===
proteins = sorted(set(enrichment_matrix.index).intersection(enrichment_matrix.columns))
print("✅ Number of proteins in enrichment matrix:", len(proteins))

# === Build symmetric control matrix ===
control_matrix = pd.DataFrame(index=proteins, columns=proteins, dtype=float)

filled_count = 0
for _, row in control_df.iterrows():
    g1, g2, score = row['geneID1'], row['geneID2'], row['Score']
    if g1 in proteins and g2 in proteins:
        control_matrix.loc[g1, g2] = score
        control_matrix.loc[g2, g1] = score
        filled_count += 1

print(f"✅ Filled {filled_count} control matrix entries with CSV scores")

# === Spot check a few entries ===
spot_genes = proteins[:5]
for g1 in spot_genes:
    for g2 in spot_genes:
        if g1 != g2:
            print(f"Spot check ({g1}, {g2}): control={control_matrix.loc[g1, g2]}, enrich={enrichment_matrix.loc[g1, g2]}")

# === Extract matching upper triangle values ===
x_enrichment = []
y_control = []
pairs_used = 0

for i, p1 in enumerate(proteins):
    for j, p2 in enumerate(proteins):
        if i < j:
            enrich = enrichment_matrix.loc[p1, p2] if p1 in enrichment_matrix.index and p2 in enrichment_matrix.columns else np.nan
            control = control_matrix.loc[p1, p2] if p1 in control_matrix.index and p2 in control_matrix.columns else np.nan
            if pd.notna(enrich) and pd.notna(control):
                x_enrichment.append(enrich)
                y_control.append(control)
                pairs_used += 1

print(f"✅ Matched {pairs_used} valid protein pairs for scatter plot")
print("🔹 Enrichment value example:", x_enrichment[:5])
print("🔹 Control value example:", y_control[:5])

# === Convert to DataFrame for easier filtering
scatter_df = pd.DataFrame({
    'enrichment': x_enrichment,
    'control': y_control
})

# === Remove the single point with the highest control score
max_idx = scatter_df['control'].idxmax()
# === Clip control and enrichment values at the 99th percentile
control_clip = scatter_df['control'].quantile(0.99)
enrich_clip = scatter_df['enrichment'].quantile(0.99)

print(f"🚧 Clipping control values above {control_clip:.4f}")
print(f"🚧 Clipping enrichment values above {enrich_clip:.4f}")

scatter_df['control_clipped'] = scatter_df['control'].clip(upper=control_clip)
scatter_df['enrichment_clipped'] = scatter_df['enrichment'].clip(upper=enrich_clip)


plt.figure(figsize=(8, 6))
sns.scatterplot(data=scatter_df, x='enrichment', y='control', alpha=0.5, edgecolor=None)
plt.title("Enrichment vs. Control Score (Outlier Removed)")
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Co-enrichment Score")
plt.ylabel("CSV Pairwise Score")
plt.grid(True)
plt.tight_layout()
plt.savefig('enrichment-control-scatter.png')
