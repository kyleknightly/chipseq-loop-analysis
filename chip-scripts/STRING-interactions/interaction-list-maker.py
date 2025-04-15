"""
Turns an interactions tsv into a tsv with each protein and a list of the proteins it interacts with
"""

import pandas as pd

infile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/STRING/groups/strong-interactions.tsv'

df = pd.read_csv(infile, sep='\t', header = 0)

# Create a dictionary to store the relationships
protein_dict = {}

# Populate the dictionary
for _, row in df.iterrows():
    prot1 = row['protein1']
    prot2 = row['protein2']
    
    if prot1 not in protein_dict:
        protein_dict[prot1] = set()
    if prot2 not in protein_dict:
        protein_dict[prot2] = set()
    
    protein_dict[prot1].add(prot2)
    protein_dict[prot2].add(prot1)

# Convert the dictionary to a DataFrame
protein_list = [(key, list(values)) for key, values in protein_dict.items()]
df_protein = pd.DataFrame(protein_list, columns=['protein', 'associated_proteins'])

df_protein.to_csv('strong-interaction-lists.tsv',sep='\t', index=False)
