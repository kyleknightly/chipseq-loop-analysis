"""
This takes a networkx graph as a gpickle file and runs a 
clustering algorithm on it

It returns a bed file of each anchor and its cluster
"""
import networkx as nx
from community import community_louvain
import os
import pandas as pd
import leidenalg as la
import igraph as ig

# Load the graph from gpickle file
graph = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle'
G_nx = nx.read_gpickle(graph)
G_ig = ig.Graph.from_networkx(G_nx)

def find_leiden_communities(G):
    temp_graph = ig.Graph.from_networkx(G) if isinstance(G, nx.Graph) else G
    leiden_communities = la.find_partition(temp_graph,
                                            partition_type=la.ModularityVertexPartition, 
                                            #weights='weight',
                                            seed=0)
    return leiden_communities
    
def find_infomap_communities(G, weight_attribute=None):
    temp_graph = ig.Graph.from_networkx(G) if isinstance(G, nx.Graph) else G
    infomap_communities = temp_graph.community_infomap(edge_weights=weight_attribute)
    return infomap_communities

leiden_communities = find_leiden_communities(G_ig)
infomap_communities = find_infomap_communities(G_ig)

# Create a BED file from the community assignments
def create_bed_file(communities, nodes, output_file):
    with open(output_file, 'w') as f:
        for node, community in zip(nodes, communities.membership):
            chr, start, end = node
            f.write(f"{chr}\t{start}\t{end}\t{community}\n")


nodes = list(G_nx.nodes)

create_bed_file(infomap_communities, nodes, 'infomap_communities.bed')