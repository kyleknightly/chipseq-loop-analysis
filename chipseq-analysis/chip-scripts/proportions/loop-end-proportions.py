"""
This takes a paired-anchor TFs file and returns a gpickle dictionary with
the proportions of loop-ends that have each TF.
i.e. if an a protein is at one anchor, but that anchor is involved in 2 loops,
it has 2 loop ends.
"""

import pandas as pd
import ast
import pickle
from collections import defaultdict

bed_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/paired-anchor-all-merged-filtered-TFs.bed'

df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
    #print(df)
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
#print(df)

counts = defaultdict(int)

# Iterate over the rows and count the occurrences of each protein
for _, row in df.iterrows():
    for prot in row['prots1']:
        counts[prot] += 1
    for prot in row['prots2']:
        counts[prot] += 1



# Step 2: Calculate the total number of paired anchors
loop_ends = len(df)*2
print(len(df))
spot_check=['NFYA', 'NFYB','NFYC','FOSL1','CTCF','ZNF143']
for spot in spot_check:
    print(spot)
    print(counts[spot])
# Step 3: Compute the proportion for each protein
proportions = {protein: float(count) / loop_ends for protein, count in counts.items()}
# print(proportions)


with open('all-merged-filtered-loop_end_proportions.pkl', 'wb') as pickle_file:
    pickle.dump(proportions, pickle_file)