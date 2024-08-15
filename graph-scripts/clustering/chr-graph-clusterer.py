"""
This takes a networkx graph as a gpickle file and runs a 
clustering algorithm on it after first separating it by chromosome

It returns a bed file of each anchor and its cluster
"""

import networkx as nx
from community import community_louvain
import os
import pandas as pd

# Load the graph from gpickle file
graph = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle'
G = nx.read_gpickle(graph)

# Function to split the graph into subgraphs based on the 'chr' attribute
def split_graph_by_attribute(graph, attribute):
    subgraphs = {}
    for node, data in graph.nodes(data=True):
        attr_value = data.get(attribute)
        if attr_value is not None:
            if attr_value not in subgraphs:
                subgraphs[attr_value] = graph.subgraph([])
            subgraphs[attr_value] = nx.compose(subgraphs[attr_value], graph.subgraph([node] + list(graph.neighbors(node))))
    return subgraphs

# Split the graph based on the 'chr' attribute
subgraphs = split_graph_by_attribute(G, 'chr')

# Function to apply clustering to each subgraph
def cluster_subgraphs(subgraphs):
    clustered_subgraphs = {}
    for attr_value, subgraph in subgraphs.items():
        # Applying Louvain method for community detection
        partition = community_louvain.best_partition(subgraph)
        clustered_subgraphs[attr_value] = partition
    return clustered_subgraphs

# Cluster each subgraph
clustered_subgraphs = cluster_subgraphs(subgraphs)

# Prepare the BED data
bed_data = []

for chr_value, partition in clustered_subgraphs.items():
    for node, community in partition.items():
        chr_name, start, end = node
        bed_data.append([chr_name, start, end, community])

# Convert the data to a DataFrame
bed_df = pd.DataFrame(bed_data, columns=['chr', 'start', 'end', 'community'])

# Save to a BED file
bed_df.to_csv('communities.bed', sep='\t', header=False, index=False)