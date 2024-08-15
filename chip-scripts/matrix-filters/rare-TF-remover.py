"""
This removes TFs from a matrix that don't show up at enough anchors
n is that number
"""

import pandas as pd
import ast
import os

bed_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/paired-anchor-TFs.bed'

#how many friends makes you cool? (how many anchors does a protein need to be at to stay)
n=50

df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
    #print(df)
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
#print(df)

df1 = df[['ch1', 'start1', 'end1', 'prots1']]
df2 = df[['ch2', 'start2', 'end2', 'prots2']]
df2.columns = ['ch1', 'start1', 'end1', 'prots1']

comb_df = pd.concat([df1, df2])
comb_df['prots1'] = comb_df['prots1'].apply(tuple)
comb_df = comb_df.drop_duplicates()
comb_df['prots1'] = comb_df['prots1'].apply(list)
comb_df.reset_index(drop=True, inplace=True)

proteins = set(comb_df['prots1'].explode())
#proportion of anchors that each protein shows up at
protein_counts = {protein: 0 for protein in proteins}

for _, row in comb_df.iterrows():
    proteins1 = set(row['prots1'])
    for protein in proteins1:
        protein_counts[protein] += 1

remove = []
for key, value in protein_counts.items():
    if value <= n:
        remove.append(key)

print(remove)
print(len(remove))

matrix = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/same-anchor-interactions/non-ctcf-anchors/anchor-enrichments-pseudo-nonctcf-anchor-contact-matrix.tsv'
df = pd.read_csv(matrix, sep='\t', header=0, index_col=0)
remove_rows = [item for item in remove if item in df.index]
remove_columns = [item for item in remove if item in df.columns]

# Remove the specified rows and columns
df = df.drop(index=remove_rows, columns=remove_columns)

name = os.path.basename(matrix)
df.to_csv('geq50-' + name, sep='\t', index=True)