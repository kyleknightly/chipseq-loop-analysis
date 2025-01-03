"""
This takes an anchors-TF file and scrambles all of the proteins around,
preserving each anchor's original protein count

It does this n times, creating a contact matrix for each,
and then averaging them for use as a control
"""

import pandas as pd
import ast
import random

n=10000

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/anchor-TFs.bed'

df = pd.read_csv(file, sep='\t', names = ['ch', 'start', 'end', 'prots'])

df['prots'] = df['prots'].apply(ast.literal_eval)

df['count']=df['prots'].apply(lambda x: len(x))

print(df)

# print(len(all_prots))
# print(df['count'].sum())

matrices = []

for i in range(n):
    all_prots = df['prots'].explode().dropna().tolist()
    random.shuffle(all_prots)

    for index, row in df.iterrows():
        new_prots = [all_prots.pop(0) for i in range(row['count'])]
        df.at[index, 'prots'] = new_prots

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
    matrices.append(matrix)
    print(matrix)

avg_matrix = pd.DataFrame(0, index=matrices[0].index, columns=matrices[0].columns)

# Iterate over each DataFrame in matrices
for matrix in matrices:
    for row in avg_matrix.index:
        for col in avg_matrix.columns:
            avg_matrix.at[row, col] += matrix.at[row, col]

print(avg_matrix)
# Divide by the number of matrices to get the average
avg_matrix = avg_matrix / len(matrices)
print(avg_matrix)

avg_matrix.to_csv('meg-cis-shuf-contacts.tsv', sep='\t', header=True, index=True)