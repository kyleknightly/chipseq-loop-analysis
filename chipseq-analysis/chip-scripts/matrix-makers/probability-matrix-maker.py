"""
This makes a matrix of P(A|B) for each protein using a paired-anchor-TFs file
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import ast
file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names = ['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
print('DF MADE')

proteins = set(df['prots1'].explode()).union(set(df['prots2'].explode()))
proteins = {p for p in proteins if pd.notna(p)}
#print(proteins)

matrix = pd.DataFrame(index=proteins, columns=proteins)


"""For each protein we count the loop ends it's at, and make a dict of counts of all the trans loop anchors"""
for prot in proteins:
    count = 0
    other_counts = {other_prot: 0 for other_prot in proteins}
    for _, row in df.iterrows():
        proteins1 = set(row['prots1'])
        proteins2 = set(row['prots2'])
        if prot in proteins1:
            count+=1
            for other_prot in proteins2:
                other_counts[other_prot]+=1
        if prot in proteins2:
            count+=1
            for other_prot in proteins1:
                other_counts[other_prot]+=1
    for other_prot, other_count in other_counts.items():
        matrix.at[prot, other_prot] = other_count/count
        #print(other_count/count)
    print(prot)

linkage_matrix = linkage(matrix, method='ward')
dendro = leaves_list(linkage_matrix)
sorted_matrix = matrix.iloc[dendro, :].iloc[:, dendro]

sorted_matrix.to_csv('trans-probability-matrix.tsv', sep='\t')