"""
This adds or subtracts a pseudo count from a contact matrix
"""
import pandas as pd

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/loop-ends/trans-contacts-paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', index_col = 0)
print(df)
df = df.apply(pd.to_numeric, errors='coerce')
df = df+1
print(df)
df.to_csv('pseudo-trans-contacts-paired-anchor-TFs.bed', sep='\t', index = True)