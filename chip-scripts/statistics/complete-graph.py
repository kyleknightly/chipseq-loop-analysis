import pandas as pd
import numpy as np
import ast
from itertools import combinations

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names = ['ch', 'start', 'end', 'prots'])
df['prots'] = df['prots'].apply(ast.literal_eval)

# Get the unique proteins across all rows
all_proteins = set(protein for prots in df['prots'] for protein in prots)
protein_list = sorted(list(all_proteins))  # Sort to keep matrix consistent

# Initialize a contact matrix (size: number of unique proteins)
contact_matrix = pd.DataFrame(0, index=protein_list, columns=protein_list)

# Create contacts in the matrix
for i, row_i in df.iterrows():
    print(i)
    prots_i = row_i['prots']
    
    # Compare only with following rows to avoid double counting
    for j in range(i+1, len(df)):
        prots_j = df.loc[j, 'prots']
        
        # Record contacts for each combination of proteins between row i and row j
        for protein_i in prots_i:
            for protein_j in prots_j:
                contact_matrix.loc[protein_i, protein_j] += 1
                contact_matrix.loc[protein_j, protein_i] += 1  # Ensure symmetry

contact_matrix.to_csv('complete-graph-contacts.tsv', sep='\t', header=True, index=True)