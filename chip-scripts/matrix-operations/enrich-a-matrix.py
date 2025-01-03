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
with open('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/all-merged-filtered-loop_end_proportions.pkl', 'rb') as pickle_file:
    proportions = pickle.load(pickle_file)
spot_check=['NFYA', 'NFYB', 'NFYC']
# ['NFYA', 'NFYB','NFYC','FOSL1','CTCF','ZNF143']
# ['ZMYM4', 'ZNF219', 'SALL1', 'GATAD2A', 'ARID5B', 'ARID3A']
contact_matrix = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/geq10k-pseudo-all-merged-filtered-anchor-contacts.tsv"
df = pd.read_csv(contact_matrix, sep='\t', header=0, index_col=0).astype(float)
print(df)
for row in df.index:
    for col in df.columns:
        if row in proportions and col in proportions:
            # for spot in spot_check:
                # if row==col==spot:
                #     print(spot)
                #     print(df.loc[row,col])
                #     print(float(proportions[row]))
                #     print(252613 * 2 * float(proportions[row])* float(proportions[col])+1)
                #     print(df.loc[row,col]/(252613 * 2 * float(proportions[row])* float(proportions[col])+1))
            df.loc[row, col] /= (107619 * float(proportions[row]) * float(proportions[col])+1) 
            # df.loc[row, col] /= (252613 * 2 * float(proportions[row]) * float(proportions[col])+1) #LOOP ENDS
            #df.loc[row, col] /= (107619 * 4.636 * float(proportions[row]) * float(proportions[col])+1) #37836
            #anchors * (loops/anchor + 1(self)) * probability
# print(df)

name = os.path.basename(contact_matrix)
df.to_csv('enrichments-' + name, sep='\t', index=True)
