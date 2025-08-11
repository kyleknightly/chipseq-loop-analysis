import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from tqdm import tqdm

enrichment_matrix_file = '/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/trans_enrichment_matrix.tsv'
control_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/HINT/hint_binary_control_full.csv'

# Load matrices
enrichment_matrix = pd.read_csv(enrichment_matrix_file, sep='\t', index_col=0)
control_df = pd.read_csv(control_file)

# Ensure all identifiers are strings
enrichment_matrix.index = enrichment_matrix.index.astype(str)
enrichment_matrix.columns = enrichment_matrix.columns.astype(str)
control_df['geneID1'] = control_df['geneID1'].astype(str)
control_df['geneID2'] = control_df['geneID2'].astype(str)

# Construct symmetric control matrix
control_matrix = pd.DataFrame(index=enrichment_matrix.index, columns=enrichment_matrix.columns, dtype=float)
for _, row in control_df.iterrows():
    g1, g2, score = row['geneID1'], row['geneID2'], row['Score']
    if g1 in control_matrix.index and g2 in control_matrix.columns:
        control_matrix.loc[g1, g2] = score
        control_matrix.loc[g2, g1] = score

# Get common proteins
common = enrichment_matrix.index.intersection(control_matrix.index)
enrichment_matrix = enrichment_matrix.loc[common, common]
control_matrix = control_matrix.loc[common, common]

# Calculate Pearson correlation for each protein
protein_scores = []

for protein in tqdm(common, desc="Computing protein correlations"):
    enrichment_row = enrichment_matrix.loc[protein].values
    control_row = control_matrix.loc[protein].values
    
    # Mask NaNs
    mask = ~np.isnan(enrichment_row) & ~np.isnan(control_row)
    
    if np.sum(mask) >= 2:  # Need at least 2 points for correlation
        corr, p_value = pearsonr(enrichment_row[mask], control_row[mask])
        protein_scores.append({
            'protein': protein,
            'pearson_correlation': corr,
            'p_value': p_value,
            'n_valid_interactions': np.sum(mask)
        })
    else:
        protein_scores.append({
            'protein': protein,
            'pearson_correlation': np.nan,
            'p_value': np.nan,
            'n_valid_interactions': np.sum(mask)
        })

# Convert to DataFrame
results_df = pd.DataFrame(protein_scores)

# Sort by correlation (highest first, NaNs last)
results_df = results_df.sort_values('pearson_correlation', ascending=False, na_position='last')

# Display summary statistics
print(f"Total proteins analyzed: {len(results_df)}")
print(f"Proteins with valid correlations: {results_df['pearson_correlation'].notna().sum()}")
print(f"Mean correlation: {results_df['pearson_correlation'].mean():.4f}")
print(f"Median correlation: {results_df['pearson_correlation'].median():.4f}")
print(f"Std correlation: {results_df['pearson_correlation'].std():.4f}")

print("\nTop 10 correlations:")
print(results_df.head(10))

print("\nBottom 10 correlations:")
print(results_df.tail(10))

# Save results
output_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/HINT-pearson-scores.tsv'
results_df.to_csv(output_file, sep='\t', index=False)
print(f"\nResults saved to: {output_file}")