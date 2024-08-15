"""
This creates a tsv of every protein and the average degree of the anchors they inhabit
also runs 1 sample t-tests
"""

import pandas as pd
import ast
import networkx as nx
import scipy.stats as stats

# Load the TSV file into a DataFrame
degree_df = pd.read_csv('/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-degrees.bed', sep='\t', names = ['chr', 'start', 'end', 'degree'])
df = pd.read_csv('/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2-anchor-TFs.bed', sep='\t', names=['chr', 'start','end', 'prots'])
graph = nx.read_gpickle('/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle')

# Convert the `prots` column from strings to lists using ast.literal_eval
df['prots'] = df['prots'].apply(ast.literal_eval)

merged_df = pd.merge(df, degree_df[['chr', 'start', 'end', 'degree']], on=['chr', 'start', 'end'], how='left')
print(merged_df)

# Calculate the degree of each node
degrees = dict(graph.degree())
all_degrees= [degrees[node] for node in graph.nodes()]

# Create a dictionary to store the degrees for each protein
protein_degrees = {}

# Iterate through each row in the DataFrame
for index, row in merged_df.iterrows():
    degree = row['degree']
    for protein in row['prots']:
        if protein not in protein_degrees:
            protein_degrees[protein] = []
        protein_degrees[protein].append(degree)

# Initialize lists to store results
proteins = []
t_stats = []
p_values = []
avg_degrees = []

# Iterate through each protein in protein_degrees
for protein, degrees in protein_degrees.items():
    proteins.append(protein)
    avg_degrees.append(sum(degrees)/len(degrees))  # Store the average degree for the protein
    
    if len(degrees) > 1:
        t_stat, p_value = stats.ttest_ind(all_degrees, degrees)
        if t_stat<0:
            p_value = p_value/2
        else:
            p_value = 1 - (p_value / 2)
    else:
        t_stat, p_value = stats.ttest_1samp(all_degrees, degrees[0])
        if t_stat<0:
            p_value = p_value/2
        else:
            p_value = 1 - (p_value / 2)

    t_stats.append(t_stat)
    p_values.append(p_value)

# Create a DataFrame with the results
final_df = pd.DataFrame({
    'protein': proteins,
    'avg_degree': avg_degrees,
    't_stat': t_stats,
    'p_value': p_values
})
final_df= final_df.sort_values(by='p_value')
# Display the DataFrame
print(final_df)
final_df.to_csv('protein-degree-tests.tsv', sep='\t', index=False)