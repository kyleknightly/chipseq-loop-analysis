import pandas as pd
import ast

file = '/mnt/altnas/work/Kyle.Knightly/anchor-star-ness.bed'
df = pd.read_csv(file, sep='\t', names=['chr','start','end', 'prots', 'starness'])
df['prots'] = df['prots'].apply(ast.literal_eval)

data = {}
data['none'] = []
data ['any'] = []

for index, row in df.iterrows():
    starness = row['starness']
    if row['prots']==[]:
        data['none'].append(starness)
    else:
        data['any'].append(starness)
    for protein in row['prots']:
        if protein not in data:
            data[protein] = []
        data[protein].append(starness)

tsv_data = []
for prot, vals in data.items():
    avg = sum(vals) / len(vals)
    tsv_data.append([prot, avg])
newdf = pd.DataFrame(tsv_data, columns=['Protein', 'Avg_Starness'])

newdf= newdf.sort_values(by='Avg_Starness')
newdf.to_csv('prot-avg-starness.tsv', sep='\t', header=False, index=False)