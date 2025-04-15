"""
This will take a paired anchor TF file and turn it into an anchor-TF file with repeated anchors
"""
import pandas as pd
import ast
file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names = ['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)
print('DF MADE')
df_combined = pd.DataFrame({
    'chr': df['chr1'].append(df['chr2']).reset_index(drop=True),
    'start': df['start1'].append(df['start2']).reset_index(drop=True),
    'end': df['end1'].append(df['end2']).reset_index(drop=True),
    'prots': df['prots1'].append(df['prots2']).reset_index(drop=True)
})

print(df)
print(df_combined)
df_combined.to_csv('loop-end-TFs.bed', sep='\t', index = False, header = False)