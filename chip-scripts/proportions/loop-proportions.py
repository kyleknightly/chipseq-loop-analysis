"""
With a paired-anchor TF file it finds the proportion of loops for each protein
and creates a gpickle dict for use in enrichment
"""

import pandas as pd
import ast
import pickle
from collections import defaultdict

bed_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/paired-anchor-TFs.bed'

df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
    #print(df)
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
#print(df)

protein_counts = defaultdict(int)

# Iterate over the rows and count the occurrences of each protein
for _, row in df.iterrows():
    proteins = set(row['prots1']).union(set(row['prots2']))
    for protein in proteins:
        protein_counts[protein] += 1



# Step 2: Calculate the total number of paired anchors
total_loops = len(df)
print(len(df))

# Step 3: Compute the proportion for each protein
proportions = {protein: float(count) / total_loops for protein, count in protein_counts.items()}
print(proportions)

# with open('loop_proportions.pkl', 'wb') as pickle_file:
#     pickle.dump(proportions, pickle_file)