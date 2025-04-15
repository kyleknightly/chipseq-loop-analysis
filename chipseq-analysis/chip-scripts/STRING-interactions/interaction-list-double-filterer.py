"""
This is taking an interactions lists file and an interactions file and intersecting them
"""

import pandas as pd
import ast

interactions = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/controls/75percentile-interactions.tsv'
interactions_lists = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/controls/leq25-interaction_lists.tsv'

interactions_df = df = pd.read_csv(interactions, sep='\t')
list_df = df = pd.read_csv(interactions_lists, sep='\t')
list_df['associated_proteins'] = list_df['associated_proteins'].apply(ast.literal_eval)

print(interactions_df)
print(list_df)

tuple_list = list(interactions_df[['protein1', 'protein2']].itertuples(index=False, name=None))
#print(tuple_list)


for index, row in list_df.iterrows():
    filtered_proteins = [protein2 for protein2 in row['associated_proteins'] if (row['protein'], protein2) in tuple_list]
    list_df.at[index, 'associated_proteins'] = filtered_proteins

print(list_df)
list_df = list_df[list_df['associated_proteins'].apply(lambda x: len(x) > 0)]

list_df.to_csv('leq25-75percentile-interaction_lists.tsv', sep = '\t', index = False)