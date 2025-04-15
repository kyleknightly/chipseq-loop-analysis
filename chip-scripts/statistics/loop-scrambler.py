"""
This takes a paired anchor TF file and scrambles around the coanchors
This can produce control contacts for a random graph based on frequencies
"""

import pandas as pd
import ast


file = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names = ['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)

# Separate the DataFrame into two parts
first_four_columns = df[['chr1', 'start1', 'end1', 'prots1']]
last_four_columns = df[['chr2', 'start2', 'end2', 'prots2']]

# Shuffle the rows of the last four columns
shuffled_last_four_columns = last_four_columns.sample(frac=1).reset_index(drop=True)

# Concatenate the first four columns with the shuffled last four columns
result_df = pd.concat([first_four_columns, shuffled_last_four_columns], axis=1)

result_df.to_csv('2scrambled-paired-anchor-TFs.bedpe', sep='\t', index=False, header=False)