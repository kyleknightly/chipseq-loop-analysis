import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict

# Path to the gpickle file
file_path = "/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2/hepg2-anchor-graph.gpickle"

# Load the graph from the gpickle file
G = nx.read_gpickle(file_path)

kb = 50  # loop size threshold

# Identify edges to remove based on the size attribute
edges_to_remove = [
    (u, v) for u, v, attr in G.edges(data=True) if attr.get('size', 0) > (kb*1000)
]
print(len(edges_to_remove))
# Remove the identified edges from the graph
G.remove_edges_from(edges_to_remove)

nodes_to_remove = []
for node in G.nodes:
    neighbors = list(G.neighbors(node))
    n_neighbors = len(neighbors)
    if n_neighbors==0:
        nodes_to_remove.append(node)

G.remove_nodes_from(nodes_to_remove)


for node in G.nodes:
    # Get the 'prots' attribute, default to an empty list if not present
    prots = G.nodes[node].get('prots', [])
    
    # Ensure 'prots' is a list (handle cases where it might not be)
    if not isinstance(prots, list):
        prots = []
    
    # Append "none" or "any" based on the current state of the list
    if not prots:  # If the list is empty
        prots.append("none")
    else:  # If the list is non-empty
        prots.append("any")
    
    # Update the 'prots' attribute in the node
    G.nodes[node]['prots'] = prots



# Dictionary to store the connectedness fraction for each node
coeffs = {}
degs = defaultdict(list)
# Iterate through each node in the graph
for node in G.nodes:
    # Get the neighbors of the node
    neighbors = list(G.neighbors(node))
    # Calculate the number of neighbors (N)
    n_neighbors = len(neighbors)
    # If the node has less than 2 neighbors, connectedness is 0 by definition (no possible edges)
    if n_neighbors > 1:
        for protein in G.nodes[node].get('prots', None):
            degs[protein].append(n_neighbors)

         # Subgraph induced by the neighbors of the node
        neighbor_subgraph = G.subgraph(neighbors)
        
        # Calculate the number of actual edges between the neighbors (A)
        actual_edges = neighbor_subgraph.number_of_edges()
        
        # Calculate the number of possible edges in a complete subgraph (N)
        possible_edges = n_neighbors * (n_neighbors - 1) / 2  # Combination of n_neighbors choose 2
        
        # Compute the connectedness fraction: A / N
        coeffs[node] = (actual_edges / possible_edges, G.nodes[node].get('prots', None))

protein_totals = defaultdict(float)  # Sum of coefficients
protein_counts = defaultdict(int)    # Count of appearances

avg_degs = {protein: sum(degs[protein])/len(degs[protein]) for protein in degs.keys()}
# print(avg_degs)


# Iterate through the coeffs dictionary
for node, (coefficient, prots) in coeffs.items():
    if prots:  # Ensure the protein list is not empty or None
        for protein in prots:
            protein_totals[protein] += coefficient
            protein_counts[protein] += 1

# Calculate the average coefficient for each protein
avg_coeffs = {
    protein: protein_totals[protein] / protein_counts[protein]
    for protein in protein_totals
}

plt.figure(figsize=(30, 20))

common_proteins = set(avg_degs.keys()) & set(avg_coeffs.keys())

# Define the proteins you want to highlight with red dots
highlight_proteins = ['CTCF', 'STAG1', 'RAD21', 'SMC3']  # Replace with your proteins

# Prepare data for scatter plot with color coding
colors = ['red' if protein in highlight_proteins else 'blue' for protein in common_proteins]

# Extract degree and coefficient values for common proteins
x_values = [avg_degs[protein] for protein in common_proteins]
y_values = [avg_coeffs[protein] for protein in common_proteins]

plt.scatter(x_values, y_values, c=colors, alpha=0.7)

for protein, x, y in zip(common_proteins, x_values, y_values,):
    plt.text(x, y, protein, fontsize=8, ha='right', va='bottom')

plt.xlim(0, 9)  # Fixed x-axis range
plt.ylim(0, 0.41)  # Fixed y-axis range

# Customize the plot
plt.xlabel('Average Degree')
plt.ylabel('Average Clustering Coefficient')
plt.title(' ')
plt.grid(True)


plt.savefig('leq'+str(kb)+'kb-stariness-degree.png')