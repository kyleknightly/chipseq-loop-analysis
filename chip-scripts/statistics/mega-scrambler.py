"""
This takes a paired anchor TF file and scrambles around the coanchors n times
Then it produces a contact matrix for each, and takes the average of them
"""

import pandas as pd
import ast

n=10

file = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names = ['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)

# Separate the DataFrame into two parts
first_four_columns = df[['chr1', 'start1', 'end1', 'prots1']]
last_four_columns = df[['chr2', 'start2', 'end2', 'prots2']]

matrices = []
for i in range(n):
    print(i)
    # Shuffle the rows of the last four columns
    shuffled_last_four_columns = last_four_columns.sample(frac=1).reset_index(drop=True)

    # Concatenate the first four columns with the shuffled last four columns
    result_df = pd.concat([first_four_columns, shuffled_last_four_columns], axis=1)

    proteins = set(result_df['prots1'].explode()).union(set(result_df['prots2'].explode()))
    
    #Initialize contact_matrix with a 1 in each (for analysis reasons)
    """
    EDIT STARTING VALUE (+1)
    """
    contact_matrix = pd.DataFrame(1, index=proteins, columns=proteins)

    # Step 4: Populate the contact map
    for _, row in result_df.iterrows():
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
print(avg_matrix)
avg_matrix = avg_matrix / len(matrices)
print(avg_matrix)

avg_matrix.to_csv('avg-scrambled-pseudo-trans-contacts-paired-anchor-TFs.tsv', sep='\t')
