"""
Divides one matridx by another matching up indeces
Useful for enrichment with scrambled matrices
"""

import pandas as pd
import os

infile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/loop-ends/pseudo-trans-contacts-paired-anchor-TFs.bed'
df = pd.read_csv(infile, sep='\t', index_col = 0)
df = df.apply(pd.to_numeric, errors='coerce')
print(df)
divfile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/loop-scrambling/avg-scrambled-pseudo-trans-contacts-paired-anchor-TFs.tsv'
divdf = pd.read_csv(divfile, sep='\t', index_col = 0)
divdf = divdf.apply(pd.to_numeric, errors='coerce')
print(divdf)

# Initialize an empty DataFrame with the same shape
result_df = pd.DataFrame(index=df.index, columns=df.columns)

# Iterate through each element
for row in df.index:
    for col in df.columns:
        result_df.at[row, col] = df.at[row, col] / divdf.at[row, col]
name = os.path.basename(infile)

result_df.to_csv('avg-scrambled-enrichment-'+name, sep='\t', index=True)