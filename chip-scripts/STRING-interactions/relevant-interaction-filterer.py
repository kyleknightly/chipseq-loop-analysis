"""
Get paxdbs from https://www.uniprot.org/id-mapping
This takes a file of relevant proteins and intersects them with the physical links
from STRING
"""

import pandas as pd

infile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/controls/9606.protein.physical.links.v12.0.txt'
filter_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/controls/our-prot-paxdbs.tsv'

unfiltered = pd.read_csv(infile, sep=' ', header = 0)
filter = pd.read_csv(filter_file, sep='\t', header = 0)
print(filter)
paxdb_ids = set()
for ids in filter['PaxDb']:
    cleaned_ids = [id_.replace('-', '.') for id_ in str(ids).split(';')]
    paxdb_ids.update(cleaned_ids)
print(paxdb_ids)

df_filtered = unfiltered[
    unfiltered['protein1'].isin(paxdb_ids) &
    unfiltered['protein2'].isin(paxdb_ids)
]
# Iterate through each row and modify columns
for index, row in df_filtered.iterrows():
    df_filtered.at[index, 'protein1'] = str(row['protein1']).replace('.', '-')
    df_filtered.at[index, 'protein2'] = str(row['protein2']).replace('.', '-')


print(df_filtered)

# Remove ';' from the PaxDb column and create a mapping
filter['PaxDb'] = filter['PaxDb'].str.replace(';', '')
id_to_name_mapping = dict(zip(filter['PaxDb'], filter['From']))

# Replace protein IDs in the main DataFrame with protein names
df_filtered['protein1'] = df_filtered['protein1'].map(id_to_name_mapping).fillna(df_filtered['protein1'])
df_filtered['protein2'] = df_filtered['protein2'].map(id_to_name_mapping).fillna(df_filtered['protein2'])

print(df_filtered)
# Save the result to a new file
output_file = 'test.tsv'
df_filtered.to_csv(output_file, sep='\t', index=False)