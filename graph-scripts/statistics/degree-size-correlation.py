"""
This creates a scatterplot/box plot to show correlation between the degree of anchors
and the size of their loops
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

file = "/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle"
G = nx.read_gpickle(file)

# Step 1: Calculate node degrees
node_degrees = dict(G.degree())

# Step 2: Collect edge sizes associated with each node
node_edge_sizes = {node: [] for node in G.nodes()}

for u, v, data in G.edges(data=True):
    size = data.get('size', 1)  # Default to 1 if 'size' is not specified
    node_edge_sizes[u].append(size)
    node_edge_sizes[v].append(size)

# Step 3: Calculate the average edge size for each node
node_avg_edge_sizes = {node: np.mean(sizes) if sizes else 0 for node, sizes in node_edge_sizes.items()}

# Group average edge sizes by node degree
degree_to_sizes = {}

for node, degree in node_degrees.items():
    if degree not in degree_to_sizes:
        degree_to_sizes[degree] = []
    degree_to_sizes[degree].append(node_avg_edge_sizes[node])

# Prepare data for box plot
degrees = sorted(degree_to_sizes.keys())
boxplot_data = [degree_to_sizes[degree] for degree in degrees]


# Customize the appearance of the outliers
flierprops = dict(marker='.', color='blue', alpha=0.5, markersize=5, linestyle='none')

# Create the box plot
plt.figure(figsize=(10, 5))
plt.boxplot(boxplot_data, labels=degrees)
plt.title("Box Plot of Average Edge Sizes by Node Degree")
plt.xlabel("Node Degree")
plt.ylabel("Average Edge Size")
plt.grid(True)

# Save the figure
plt.savefig('degree-size-boxplot.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
