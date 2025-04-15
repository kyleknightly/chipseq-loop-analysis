"""
Removes any anchors from an anchor-TFs file that have CTCF
"""

import pandas as pd
import ast

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/same-anchor-interactions/anchor-TFs.bed'

df = pd.read_csv(file, sep='\t', names = ['ch', 'start', 'end', 'prots'])

df['prots'] = df['prots'].apply(ast.literal_eval)

print(df)
filtered_df = df[~df['prots'].apply(lambda x: 'CTCF' in x)]

print(filtered_df)
filtered_df.to_csv('non-ctcf-anchor-TFs.bed', sep='\t')