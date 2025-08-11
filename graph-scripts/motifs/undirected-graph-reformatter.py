import networkx as nx
from collections import defaultdict

# Load your graph
G = nx.read_gpickle('/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/hepg2-anchor-graph.gpickle')

# Group edges by chromosome (still include both directions)
chromosome_edges = defaultdict(list)
for u, v in G.edges():
    chr_u, chr_v = u[0], v[0]
    chromosome_edges[chr_u].append((u, v))
    if chr_u != chr_v:
        chromosome_edges[chr_v].append((u, v))

# Now for each chromosome, assign compact node IDs and write files
for chrom_id, edges in chromosome_edges.items():
    print(f"Processing {chrom_id} with {len(edges)} edges")

    # Build per-chromosome node ID map
    unique_nodes = set()
    for u, v in edges:
        unique_nodes.update([u, v])
    node_to_id = {node: i+1 for i, node in enumerate(sorted(unique_nodes))}

    # Write .txt for mfinder
    filename = f"{chrom_id}_loops.txt"
    with open(filename, "w") as f:
        for u, v in edges:
            if u != v:  # avoid self-loops
                f.write(f"{node_to_id[u]} {node_to_id[v]} 1\n")

    # Optional: Write node mapping
    with open(f"{chrom_id}_node_map.txt", "w") as f:
        for node, node_id in node_to_id.items():
            ch, start, end = node
            f.write(f"{node_id}\t{ch}:{start}-{end}\n")