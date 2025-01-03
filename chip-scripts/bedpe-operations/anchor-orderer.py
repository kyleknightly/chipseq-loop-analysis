import pandas as pd
import ast
from collections import defaultdict

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)

anchors = pd.concat([
    df[['ch1', 'start1', 'end1']].rename(columns={'ch1': 'ch', 'start1': 'start', 'end1': 'end'}),
    df[['ch2', 'start2', 'end2']].rename(columns={'ch2': 'ch', 'start2': 'start', 'end2': 'end'})
]).drop_duplicates()

def ch_sort_key(ch):
    if ch[3:].isdigit():
        return int(ch[3:])  # Sort numeric chromosomes normally
    elif ch[3:] == 'X':
        return 1000  # Arbitrarily large value for 'X'
    elif ch[3:] == 'Y':
        return 1001  # 'Y' comes after 'X'

# Apply the custom sort
anchors['chr_sort'] = anchors['ch'].apply(ch_sort_key)

anchors = anchors.sort_values(['chr_sort', 'start', 'end']).drop('chr_sort', axis=1).reset_index(drop=True)

anchors['anchor_index'] = anchors.index
anchors.to_csv('anchors_sorted.tsv', sep='\t', index=False, header=True)