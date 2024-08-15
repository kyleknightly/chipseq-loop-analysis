"""
Converts a paired-anchor-TFs file to an anchor-TFs file
"""

import pandas as pd
import ast

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/paired-anchor-TFs.bed'

df = pd.read_csv(file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
#print(df)
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)

# Create DataFrames for the two sets of columns
df1 = df[['ch1', 'start1', 'end1', 'prots1']]
df2 = df[['ch2', 'start2', 'end2', 'prots2']]

# Rename the columns to match
df1.columns = ['ch', 'start', 'end', 'prots']
df2.columns = ['ch', 'start', 'end', 'prots']

# Concatenate the two DataFrames
combined_df = pd.concat([df1, df2])

combined_df['prots'] = combined_df['prots'].apply(tuple)

# Remove duplicates
combined_df = combined_df.drop_duplicates()

combined_df['prots'] = combined_df['prots'].apply(list)

# Print the resulting DataFrame
print(combined_df)

combined_df.to_csv('anchor-TFs.bed', sep='\t', header=False)