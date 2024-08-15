"""
This plots a networkx graph of the genome, it is not very good
"""

import networkx as nx
import matplotlib.pyplot as plt
import os


print("started")


file_path = "/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle"


def plot_graph(G):
    plt.figure(figsize=(80, 80))
    pos = nx.spring_layout(G)  # Use spring layout for better visualization
    nx.draw(G, pos, with_labels=False, node_size=50, alpha = 0.7, width=0.5, node_color="skyblue", edge_color="gray")
    name = os.path.basename(file_path)
    plt.savefig(name[:-7]+'png')
"""
def plot_colored_graph(G):
    node_colors = [data['color'] for _, data in G.nodes(data=True)]
    plt.figure(figsize=(80, 80))
    pos = nx.spring_layout(G, k=0.01)  # Use spring layout for better visualization
    nx.draw(G, pos, with_labels=False, node_size=30, alpha = 0.5, width=0.5, node_color=node_colors, edge_color="gray")
    name = os.path.basename(file_path)
    plt.savefig(name[:-7]+'png')
"""

# Load the graph
G = nx.read_gpickle(file_path)
print("PICKLE READ")
plot_graph(G)
print("NORMAL GRAPH MADE")
