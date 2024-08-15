"""
This takes a networkx graph and creates a bed file of each anchor and its degree
"""

import networkx as nx
import pandas as pd

# Load the graph from a pickle file
graph = nx.read_gpickle('/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle')

# Calculate the degree of each node
degrees = dict(graph.degree())

data = []
for node in graph.nodes():
    chr, start, end = node
    degree = degrees[node]
    data.append([chr, start, end, degree])

# Convert to a DataFrame
df = pd.DataFrame(data)

# Save to TSV file
df.to_csv('anchor-degrees.bed', sep='\t', index=False, header=False)
