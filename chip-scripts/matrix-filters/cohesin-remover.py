"""
This removes any short list of proteins from a matrix
"""

import os
import pandas as pd

remove = ['CTCF', 'RAD21', 'SMC3', 'STAG1']
matrix = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/matrices/geq10k-pseudo-trans-contacts.tsv'
df = pd.read_csv(matrix, sep='\t', header=0, index_col=0)
remove_rows = [item for item in remove if item in df.index]
remove_columns = [item for item in remove if item in df.columns]

# Remove the specified rows and columns
df = df.drop(index=remove_rows, columns=remove_columns)

name = os.path.basename(matrix)
df.to_csv('noCTCF-' + name, sep='\t', index=True)
