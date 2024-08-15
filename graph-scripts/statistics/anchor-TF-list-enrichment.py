"""
This is a widely applicable script

This takes a list of anchors and their TFs and returns each TF and how much it is enriched
in this list in comparison to it's prevelance in all OTHER anchors

Be sure to check the cols of your input file, this may need adjustment
"""

import pandas as pd
import ast
import pickle

file = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/semi_star_centers.bed'
df = pd.read_csv(file, sep='\t', names = ['chr', 'start', 'end', 'degree', 'prots'])
df['prots']=df['prots'].apply(ast.literal_eval)

# print(df)

proteins = set(df['prots'].explode())
protein_counts = {protein: 0 for protein in proteins if pd.notna(protein)}
for _, row in df.iterrows():
    proteins1 = set(row['prots'])
    for protein in proteins1:
        protein_counts[protein] += 1


allfile = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2-anchor-TFs.bed'
alldf = pd.read_csv(allfile, sep='\t', names = ['chr', 'start', 'end', 'prots'])
alldf['prots']=alldf['prots'].apply(ast.literal_eval)
# print(alldf)

# Extract the first three columns for comparison
alldf_subset = alldf[['chr', 'start', 'end']]
df_subset = df[['chr', 'start', 'end']]

# Perform a merge to find rows in alldf that are not in df
merged_df = alldf_subset.merge(df_subset, on=['chr', 'start', 'end'], how='left', indicator=True)

# Select rows that are only in alldf
otherdf = alldf[merged_df['_merge'] == 'left_only']

"""otherdf is anchors that aren't in our set"""

# print(otherdf)

other_protein_counts = {protein: 0 for protein in proteins if pd.notna(protein)}

for _, row in otherdf.iterrows():
    proteins = set(row['prots'])
    for protein in proteins1:
        other_protein_counts[protein] += 1

# Step 2: Calculate the total number of paired anchors
total_anchors = len(otherdf)

# Step 3: Compute the proportion for each protein
proportions = {protein: float(count) / total_anchors for protein, count in other_protein_counts.items()}


for protein, count in protein_counts.items():
    print(protein)
    print(count)
    print(proportions[protein])
    print(proportions[protein]*len(df))
    print((count+1)/(1+proportions[protein]*len(df)))
    print("===")
enrichment = {protein: ((count+1)/(1+proportions[protein]*len(df))) for protein, count in protein_counts.items()}
# print(enrichment)

# Sort the dictionary by values (descending order)
sorted_enrichment = sorted(enrichment.items(), key=lambda item: item[1], reverse=True)

# Write the sorted dictionary to a TSV file
with open('semi-star-enrichment.tsv', 'w') as f:
    for protein, value in sorted_enrichment:
        f.write(f"{protein}\t{value}\n")
