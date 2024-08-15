"""
This takes an anchor-TFs file and creates a dictionary of the proportion of anchors
with each protein, saved as a gpickle.
Used for enrichment
"""

import pandas as pd
import ast
import pickle

bed_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/fake-anchors/subdivide/subdivisions-TFs.bed'

df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1'])
    #print(df)
df['prots1'] = df['prots1'].apply(ast.literal_eval)
#print(df)


proteins = set(df['prots1'].explode())
#proportion of anchors that each protein shows up at
protein_counts = {protein: 0 for protein in proteins}

for _, row in df.iterrows():
    proteins1 = set(row['prots1'])
    for protein in proteins1:
        protein_counts[protein] += 1

# Step 2: Calculate the total number of paired anchors
total_anchors = len(df)
print(len(df))

# Step 3: Compute the proportion for each protein
proportions = {protein: float(count) / total_anchors for protein, count in protein_counts.items()}
print(proportions)

with open('subdivisions-proportions.pkl', 'wb') as pickle_file:
    pickle.dump(proportions, pickle_file)