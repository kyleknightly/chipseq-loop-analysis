

import pandas as pd
import ast
import networkx as nx
import scipy.stats as stats
import matplotlib.pyplot as plt

# Path to the gpickle file
file_path = "/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2/hepg2-anchor-graph.gpickle"

# Load the graph from the gpickle file
G = nx.read_gpickle(file_path)
# G = nx.erdos_renyi_graph(50000,0.0001)
print('graph made')

degrees = [deg for _, deg in G.degree()]

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(degrees, bins=range(min(degrees), max(degrees) + 2), edgecolor='black')
ax.set_xlabel('Degree')
ax.set_ylabel('Nodes')
ax.set_yscale('log')
plt.savefig("log10_np_degree_distribution.png", dpi=300, bbox_inches='tight')