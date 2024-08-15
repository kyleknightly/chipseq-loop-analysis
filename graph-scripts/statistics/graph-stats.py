"""
This prints summary statistics of a networkx graph,
clustering coefficient, degree distribution, connectivity, etc.
"""

import networkx as nx
import os

# Function to calculate summary statistics of a graph
def graph_summary_stats(G):
    stats = {}
    
    # Number of nodes and edges
    stats['number_of_nodes'] = G.number_of_nodes()
    stats['number_of_edges'] = G.number_of_edges()
    
    # Average degree
    degrees = dict(G.degree()).values()
    stats['average_degree'] = sum(degrees) / len(degrees)
    
    # Connectivity measures
    #stats['average_node_connectivity'] = nx.average_node_connectivity(G)
    
    # Average clustering coefficient
    stats['average_clustering_coefficient'] = nx.average_clustering(G)
    
    # Degree distribution (optional)
    stats['degree_distribution'] = dict(G.degree())
    
    return stats

# Path to the gpickle file
file_path = "/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle"

# Load the graph from the gpickle file
G = nx.read_gpickle(file_path)

# Calculate summary statistics
summary_stats = graph_summary_stats(G)

# Print summary statistics
for stat, value in summary_stats.items():
    print(f"{stat}: {value}")


