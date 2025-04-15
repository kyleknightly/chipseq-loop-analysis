"""
This takes an anchors-TF file and scrambles all of the proteins around,
preserving each anchor's original protein count
"""

import pandas as pd
import ast
import random

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/anchor-TFs.bed'

df = pd.read_csv(file, sep='\t', names = ['ch', 'start', 'end', 'prots'])

df['prots'] = df['prots'].apply(ast.literal_eval)

df['count']=df['prots'].apply(lambda x: len(x))

print(df)

all_prots = df['prots'].explode().dropna().tolist()

# print(len(all_prots))
# print(df['count'].sum())

random.shuffle(all_prots)
print(all_prots[:7])
print(all_prots[32])

for index, row in df.iterrows():
    new_prots = [all_prots.pop(0) for i in range(row['count'])]
    df.at[index, 'prots'] = new_prots

df = df.drop(columns=['count'])

df.to_csv('shuffled-anchor-TFs.bed', sep='\t', index=False, header=False)