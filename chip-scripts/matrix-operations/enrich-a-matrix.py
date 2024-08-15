"""
This function takes a contact matrix and calculates enrichment for each cell

IMPORTANT:
Take care to look at the line that performs the division, the numbers may need to be adjusted
based on your data set
Generally, for anchor enrichments and anchor contacts you can use
number of anchors * proportionA * proportion B + 1 (pseudocount)
"""


import pandas as pd
import os
import pickle
#/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/anchor_proportions.pkl
#
with open('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/fake-anchors/subdivide/subdivisions-proportions.pkl', 'rb') as pickle_file:
    proportions = pickle.load(pickle_file)

contact_matrix = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/fake-anchors/subdivide/pseudo-subdivisions-contacts.tsv"
df = pd.read_csv(contact_matrix, sep='\t', header=0, index_col=0).astype(float)
print(df)
for row in df.index:
    for col in df.columns:
        if row in proportions and col in proportions:
            df.loc[row, col] /= (37825 * float(proportions[row]) * float(proportions[col])+1) #37836
            #anchors * (loops/anchor + 1(self)) * probability
print(df)

name = os.path.basename(contact_matrix)
df.to_csv('big-span-enrichments-' + name, sep='\t', index=True)
