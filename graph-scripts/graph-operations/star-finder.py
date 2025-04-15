"""
This will rank nodes by their star-ness
Specifically, it will take the set of co-anchors, and return the fraction
of possible edges that are present
"""

import networkx as nx

G = nx.read_gpickle('/mnt/altnas/work/Kyle.Knightly/leq200mb-hepg2-anchor-graph.gpickle')

# Dictionary to store the connectedness fraction for each node
connectedness = {}

# Iterate through each node in the graph
for node in G.nodes:
    # Get the neighbors of the node
    neighbors = list(G.neighbors(node))
    
    # Calculate the number of neighbors (N)
    n_neighbors = len(neighbors)
    
    # If the node has less than 2 neighbors, connectedness is 0 by definition (no possible edges)
    if n_neighbors < 2:
        connectedness[node] = (0.0, G.nodes[node].get('prots', None))
        continue
    
    # Subgraph induced by the neighbors of the node
    neighbor_subgraph = G.subgraph(neighbors)
    
    # Calculate the number of actual edges between the neighbors (A)
    actual_edges = neighbor_subgraph.number_of_edges()
    
    # Calculate the number of possible edges in a complete subgraph (N)
    possible_edges = n_neighbors * (n_neighbors - 1) / 2  # Combination of n_neighbors choose 2
    
    # Compute the connectedness fraction: A / N
    connectedness[node] = (actual_edges / possible_edges, G.nodes[node].get('prots', None))

# print(connectedness)
with open('leq200mb-anchor-star-ness.bed', 'w') as f:
    for key, value in connectedness.items():
        # print(key)
        # print(value)
        chrom, start, end = key
        
        connectedness, prots = value
        f.write(f"{chrom}\t{start}\t{end}\t{prots}\t{connectedness}\n")