import pandas as pd

file_1 = '/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/rand-3/out/trans_enrichment_matrix.tsv'
file_2 = '/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/rand-3/out/trans_enrichment_matrix.tsv'
file_3 = '/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/rand-3/out/trans_enrichment_matrix.tsv'

df_1 = pd.read_csv(file_1, sep='\t', index_col = 0)
df_2 = pd.read_csv(file_2, sep='\t', index_col = 0)
df_3 = pd.read_csv(file_3, sep='\t', index_col = 0)

df_avg = (df_1 + df_2 + df_3) / 3

df_avg.to_csv('/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/new-avg-rand_enrichment_matrix.tsv', sep='\t', index=True)