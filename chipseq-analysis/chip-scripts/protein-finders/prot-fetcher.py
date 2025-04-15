"""
If you have an anchor subset you can very quickly fetch the TFs from a complete anchor-TFs file
"""

import pandas as pd
import ast

anchor_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/same-anchor-interactions/non-ctcf-anchors/non-ctcf-anchors.bed'
tfs_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/same-anchor-interactions/anchor-TFs.bed'

anchor_df = pd.read_csv(anchor_file, sep = '\t', names = ['chr', 'start', 'end'])
tfs_df = pd.read_csv(tfs_file, sep = '\t', names = ['chr', 'start', 'end', 'prots'])
tfs_df['prots'] = tfs_df['prots'].apply(ast.literal_eval)

print(anchor_df)
print(tfs_df)

merged_df = pd.merge(anchor_df, tfs_df, on=['chr', 'start', 'end'], how='left')

print(merged_df)

merged_df.to_csv('non-ctcf-anchor-TFs.bed', sep='\t')