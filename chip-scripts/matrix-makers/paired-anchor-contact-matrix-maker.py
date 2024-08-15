"""
This creates a contact matrix from a paired-anchor-TF file
This counts contacts as existing at EITHER anchor of a loop together
Cis AND trans
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import ast
import pickle
import os

def organize(contact_matrix):
    linkage_matrix = linkage(contact_matrix, method='average')
    dendro = leaves_list(linkage_matrix)
    
    # Reorder the DataFrame
    sorted_df = contact_matrix.iloc[dendro, :].iloc[:, dendro]
    return sorted_df

def contact_matrix(bed_file):
    
    # Pulls the anchors and proteins
    df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
    #print(df)
    df['prots1'] = df['prots1'].apply(ast.literal_eval)
    df['prots2'] = df['prots2'].apply(ast.literal_eval)
    print('DF MADE')
    
    proteins = set(df['prots1'].explode()).union(set(df['prots2'].explode()))
    
    #Initialize contact_matrix with a 1 in each (for analysis reasons)
    """
    EDIT STARTING VALUE (+1)
    """
    contact_matrix = pd.DataFrame(0, index=proteins, columns=proteins)

    # Step 4: Populate the contact map
    for _, row in df.iterrows():
        proteins1 = set(row['prots1'])
        proteins2 = set(row['prots2'])
    
        for protein in proteins1:
            for other_protein in proteins1.union(proteins2):
                contact_matrix.at[protein, other_protein] += 1
        for protein in proteins2:
            for other_protein in proteins1.union(proteins2):
                contact_matrix.at[protein, other_protein] += 1
    # Remove rows with no name
    contact_matrix = contact_matrix[contact_matrix.index.notna()]

    # Remove columns with no name
    contact_matrix = contact_matrix.loc[:, contact_matrix.columns.notna()]
    print('CONTACT MATRIX MADE')
    
    return contact_matrix

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/paired-anchor-TFs.bed'
contacts = contact_matrix(file)

organized_contacts=organize(contacts)
#print(organized_enrichments)
name = os.path.basename(file)
organized_contacts.to_csv('contacts-' + name, sep='\t', index=True)