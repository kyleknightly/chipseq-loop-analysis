"""
This takes an anchor-TFs file aand turns it into something usable by IGV
"""

import pandas as pd
import ast

infile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/anchor-TFs.bed'
df = pd.read_csv(infile, sep='\t', names=['chr', 'start', 'end', 'prots'])

# Uncomment and use ast.literal_eval to convert string representations of lists into actual lists
df['prots'] = df['prots'].apply(ast.literal_eval)

# Convert lists of proteins into comma-separated strings
df['prots'] = df['prots'].apply(lambda x: ','.join(x))

print(df)

df.to_csv('igv-anchor-TFs.bed', sep='\t', index=False, header=False)