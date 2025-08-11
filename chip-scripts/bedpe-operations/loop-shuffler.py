import pandas as pd
import numpy as np

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe'
df = pd.read_csv(file, sep='\t', names = ['c1', 's1', 'e1', 'c2', 's2', 'e2'])

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe'
df = pd.read_csv(file, sep='\t', names=['c1', 's1', 'e1', 'c2', 's2', 'e2'])

# Extract second anchors and shuffle them
second_anchors = df[['c2', 's2', 'e2']].copy()
shuffled = second_anchors.sample(frac=1, random_state=42).reset_index(drop=True)

# Replace the original second anchor columns with the shuffled ones
df[['c2', 's2', 'e2']] = shuffled

df.to_csv('/mnt/altnas/work/Kyle.Knightly/looppi/random-hepg2/random-hepg2-loops.bedpe', sep='\t', header=False, index=False)