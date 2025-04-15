"""
This takes the file produced by community-colorer.py and finds all trans-community loops
designated by a color value of 0,0,0.

Then it calculates enrichments with anchor proportions of each protein and
creates a tsv
"""

import ast
import pandas as pd
from collections import Counter
import pickle

bedpe = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/colored-infomap.bedpe'
df = pd.read_csv(bedpe, sep='\t', names = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', '.1', '.2', '.3', '.4', 'color'])
print(df)
TFs_bedpe = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/paired-anchor-TFs.bed'
tf_df = pd.read_csv(TFs_bedpe, sep='\t', names = ['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'] )
tf_df['prots1'] = tf_df['prots1'].apply(ast.literal_eval)
tf_df['prots2'] = tf_df['prots2'].apply(ast.literal_eval)
tf_df['prots'] = tf_df.apply(lambda row: row['prots1'] + row['prots2'], axis=1)
print(tf_df)

# Merge the DataFrames on the specified columns
merged_df = pd.merge(df, tf_df, on=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])

# Filter the merged DataFrame to include only rows with '0,0,0' in the 'color' column
filtered_df = merged_df[merged_df['color'] == '0,0,0']

# Select the relevant columns for the new DataFrame
result_df = filtered_df[['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2', 'prots']]

print(result_df)

# Count the occurrences of each protein in the 'prots' column
protein_counts = Counter([protein for sublist in result_df['prots'] for protein in sublist])

# Convert the counts to a DataFrame
protein_counts_df = pd.DataFrame(protein_counts.items(), columns=['protein', 'count'])

with open('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/loop_proportions.pkl', 'rb') as pickle_file:
    proportions = pickle.load(pickle_file)

protein_counts_df['enrichments'] = protein_counts_df.apply(
    lambda row: row['count'] / (proportions[row['protein']] * len(filtered_df)), axis=1
)

protein_counts_df = protein_counts_df.sort_values(by='enrichments', ascending=False)
# Save the counts to a TSV file
output_tsv = 'trans-community-protein-counts.tsv'
protein_counts_df.to_csv(output_tsv, sep='\t', index=False)