"""
Takes a paired-anchor TFs bedpe and removes all loops >200mb in size
"""

import pandas as pd
import ast

file = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/flag-tag-catchseq/beds/paired-anchor-TFs.bed"
df = pd.read_csv(file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)

df['distance'] = df.apply(
    lambda row: abs(row['start1'] - row['start2']) if row['ch1'] == row['ch2'] else None,
    axis=1
)
print(len(df))
# Filter out loops larger than 200 Mb or spanning different chromosomes
condition = (df['distance'] <= 200000) | (df['distance'].isna())
filtered_df = df[condition]

# Identify and print rows removed
removed_rows = df[~condition]
print(f"Rows removed: {len(removed_rows)}")
print(removed_rows)
# filtered_df = df[(df['distance'] <= 200000) | (df['distance'].isna())]
# print(len(filtered_df))
# # Drop the distance column if it's not needed
# filtered_df = filtered_df.drop(columns=['distance'])

# Save the filtered file if needed
# filtered_df.to_csv("leq200mb-paired-anchor-TFs", sep="\t", index=False, header=False)