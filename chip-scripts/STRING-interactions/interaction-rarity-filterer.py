"""
Removes proteins from an interaction list tsv that show up in too many lists
"""


import pandas as pd
import ast
from collections import Counter

# Load the DataFrame
infile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/controls/interaction-lists.tsv'
df = pd.read_csv(infile, sep='\t', header=0)
df['associated_proteins'] = df['associated_proteins'].apply(ast.literal_eval)

# Flatten the list of associated proteins and count occurrences
all_associated_proteins = [protein for sublist in df['associated_proteins'] for protein in sublist]
protein_counts = Counter(all_associated_proteins)

# Define the threshold n
n = 25  # Change this to your desired threshold

# Identify proteins to remove
proteins_to_remove = {protein for protein, count in protein_counts.items() if count >= n}

# Remove rows where the main protein is in proteins_to_remove
df_filtered = df[~df['protein'].isin(proteins_to_remove)]

# Remove proteins_to_remove from associated_proteins lists
df_filtered['associated_proteins'] = df_filtered['associated_proteins'].apply(lambda x: [protein for protein in x if protein not in proteins_to_remove])

# Reset index for the filtered DataFrame
df_filtered = df_filtered.reset_index(drop=True)

# Display the resulting DataFrame
print(df_filtered)

# Save the resulting DataFrame to a file
output_file = 'leq25-interaction_lists.tsv'
df_filtered.to_csv(output_file, sep='\t', index=False)
