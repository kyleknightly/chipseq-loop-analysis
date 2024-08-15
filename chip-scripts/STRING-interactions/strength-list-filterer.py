"""
This takes an interactions tsv (assumes col is already named combined_score)
keeps the 75th percentile
"""

import pandas as pd

# Load the data from the file
file_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/controls/interactions.tsv'  # Replace with your actual file path
df = pd.read_csv(file_path, sep='\t')

# Calculate the threshold for the top of combined scores
threshold = df['combined_score'].quantile(0.75)

# Filter the rows that are in the top of combined scores
filtered_df = df[df['combined_score'] >= threshold]

# Print the filtered DataFrame
filtered_df.to_csv('75percentile-interactions.tsv', sep='\t', index=False)
