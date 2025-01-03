"""
This removes TFs from a matrix that don't show up at enough loop ends
It requiress a paired-anchor-TF file
n is that number
"""

import pandas as pd
import ast
import os

bed_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/paired-anchor-TFs.bed'

n=10000

df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
    #print(df)
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
#print(df)

proteins = set(df['prots1'].explode()).union(set(df['prots2'].explode()))
#proportion of anchors that each protein shows up at
protein_counts = {protein: 0 for protein in proteins}

for _, row in df.iterrows():
    proteins1 = set(row['prots1'])
    proteins2 = set(row['prots2'])
    for protein in proteins1:
        protein_counts[protein] += 1
    for protein in proteins2:
        protein_counts[protein] += 1

remove = []
for key, value in protein_counts.items():
    if value < n:
        remove.append(key)

print(remove)
print(len(remove))

matrix = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/pseudo-anchor-contacts.tsv'
df = pd.read_csv(matrix, sep='\t', header=0, index_col=0)
remove_rows = [item for item in remove if item in df.index]
remove_columns = [item for item in remove if item in df.columns]

# Remove the specified rows and columns
df = df.drop(index=remove_rows, columns=remove_columns)

name = os.path.basename(matrix)
df.to_csv('geq10k-' + name, sep='\t', index=True)