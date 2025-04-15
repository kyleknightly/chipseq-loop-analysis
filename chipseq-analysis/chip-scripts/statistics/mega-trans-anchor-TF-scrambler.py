"""
This takes a paired anchor TF file and scrambles around the proteins in each  n times
Then it produces a contact matrix for each, and takes the average of them
"""

import pandas as pd
import ast
import random

n=10000

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names = ['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
df['count1']=df['prots1'].apply(lambda x: len(x))
df['count2']=df['prots2'].apply(lambda x: len(x))

matrices = []
for i in range(n):
    print(i)
    all_prots1 = df['prots1'].explode().dropna().tolist()
    all_prots2 = df['prots2'].explode().dropna().tolist()
    random.shuffle(all_prots1)
    random.shuffle(all_prots2)
    #print(df)
    for index, row in df.iterrows():
        new_prots1 = [all_prots1.pop(0) for i in range(row['count1'])]
        new_prots2 = [all_prots2.pop(0) for i in range(row['count2'])]
        df.at[index, 'prots1'] = new_prots1
        df.at[index, 'prots2'] = new_prots2
    #print(df)
    all_proteins = set()
    for prots in df['prots1']:
        all_proteins.update(prots)
    for prots in df['prots2']:
        all_proteins.update(prots)

    all_proteins = sorted(all_proteins)
    
    #Initialize contact_matrix with a 1 in each (for analysis reasons)
    """
    EDIT STARTING VALUE (+1)
    """
    contact_matrix = pd.DataFrame(1, index=all_proteins, columns=all_proteins)

    # Step 4: Populate the contact map
    for _, row in df.iterrows():
        proteins1 = set(row['prots1'])
        proteins2 = set(row['prots2'])
    
        for protein in proteins1:
            for other_protein in proteins2:
                contact_matrix.at[protein, other_protein] += 1
        for protein in proteins2:
            for other_protein in proteins1:
                contact_matrix.at[protein, other_protein] += 1
    # Remove rows with no name
    contact_matrix = contact_matrix[contact_matrix.index.notna()]

    # Remove columns with no name
    contact_matrix = contact_matrix.loc[:, contact_matrix.columns.notna()]
    print(contact_matrix)
    matrices.append(contact_matrix)

avg_matrix = pd.DataFrame(0, index=matrices[0].index, columns=matrices[0].columns)

# Iterate over each DataFrame in matrices
for matrix in matrices:
    for row in avg_matrix.index:
        for col in avg_matrix.columns:
            avg_matrix.at[row, col] += matrix.at[row, col]

# Divide by the number of matrices to get the average
#print(avg_matrix)
avg_matrix = avg_matrix / len(matrices)
#print(avg_matrix)

avg_matrix.to_csv('meg-trans-shuf-contacts.tsv', sep='\t')
