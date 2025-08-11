import networkx as nx
from collections import defaultdict

# Load your graph
G = nx.read_gpickle('/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/directed/directed-dEP-graph.gpickle')

print(type(G))
filename = 'directed-dEP-loops.txt'
unique_nodes = set()
for u, v in G.edges():
    unique_nodes.update([u, v])

node_to_id = {node: i+1 for i, node in enumerate(sorted(unique_nodes))}


with open(filename, "w") as f:
    for u, v in G.edges():
        if u != v:  # avoid self-loops
            f.write(f"{node_to_id[u]} {node_to_id[v]} 1\n")

# Optional: Write node mapping
with open(f"directed_dEP_node_map.txt", "w") as f:
        for node, node_id in node_to_id.items():
            ch, start, end = node
            f.write(f"{node_id}\t{ch}:{start}-{end}\n")