"""
This takes an anchor-TFs file and creates a contact matrix
(Cis anchor contacts)
"""

import pandas as pd
import ast

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/fake-anchors/subdivide/subdivisions-TFs.bed'

df = pd.read_csv(file, sep='\t', names = ['chr', 'start','end','prots'])

df['prots'] = df['prots'].apply(ast.literal_eval)

all_proteins = set()
for prots in df['prots']:
    all_proteins.update(prots)

all_proteins = sorted(all_proteins)

# Create an interaction matrix
matrix = pd.DataFrame(1, index=all_proteins, columns=all_proteins)

# Fill the interaction matrix
for prots in df['prots']:
    for i in range(len(prots)):
        for j in range(i, len(prots)):
            protein1 = prots[i]
            protein2 = prots[j]
            matrix.at[protein1, protein2] += 1
            matrix.at[protein2, protein1] += 1

print(matrix)  
matrix.to_csv('pseudo-subdivisions-contacts.tsv', sep='\t', index=True)