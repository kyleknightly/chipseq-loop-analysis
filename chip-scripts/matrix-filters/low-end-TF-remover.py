import pandas as pd
from collections import defaultdict
import ast
import os

"""
This finds TFs that appear at <n loop ends and removes them from a matrix.
"""

paired_anchor_tfs_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/paired-anchor-all-merged-filtered-TFs.bed'

paired_anchor_tfs = pd.read_csv(paired_anchor_tfs_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
paired_anchor_tfs['prots1'] = paired_anchor_tfs['prots1'].apply(ast.literal_eval)
paired_anchor_tfs['prots2'] = paired_anchor_tfs['prots2'].apply(ast.literal_eval)
print('paired-anchor-TFs created')

loop_end_counts = defaultdict(int)

for _, row in paired_anchor_tfs.iterrows():
    for prot in row['prots1']:
        loop_end_counts[prot] += 1
    for prot in row['prots2']:
        loop_end_counts[prot] += 1

matrix_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/pseudo-all-merged-filtered-anchor-contacts.tsv'
matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)

geq10k_proteins = [protein for protein, count in loop_end_counts.items() if count >= 10000]
matrix = matrix.loc[geq10k_proteins, geq10k_proteins]

name = os.path.basename(matrix_file)
matrix.to_csv('geq10k-' + name, sep='\t', index=True)