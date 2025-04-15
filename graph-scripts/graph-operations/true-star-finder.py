"""
This finds star graph centers in a networkx graph, high degree nodes whose neighbors aren't connected directly

It then creates a bed file with each of these anchors and their proteins
"""

import networkx as nx

G = nx.read_gpickle('/mnt/altnas/work/Kyle.Knightly/anchor-graph/LCL-loops/LCL-anchor-graph.gpickle')

def find_stars(graph):
    star_centers = []
    for node in graph.nodes:
        neighbors = list(graph.neighbors(node))
        degree = graph.degree[node]
        
        # Check if the node's degree equals the size of its connected component minus one
        component_size = len(nx.node_connected_component(graph, node))
        if degree == component_size - 1:
            # Ensure all neighbors have a degree of 1
            if all(graph.degree[neighbor] == 1 for neighbor in neighbors):
                star_centers.append((node, degree))
    
    # Sort the star centers by degree in descending order
    star_centers.sort(key=lambda x: x[1], reverse=True)
    
    return star_centers

def find_semi_stars(graph):
    star_centers = []
    other_centers = []
    for node in graph.nodes:
        neighbors = list(graph.neighbors(node))
        degree = graph.degree[node]
        prots = graph.nodes[node].get('prots', None)
        
        # Check if any of the neighbors are directly connected to each other
        is_star = True
        con_count = 0
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                if graph.has_edge(neighbors[i], neighbors[j]):
                    is_star = False
                    con_count+=1
        if is_star:
            star_centers.append((node, degree, prots))
        else:
            other_centers.append((node, degree, con_count, prots))
    
    # Sort the star centers by degree in descending order
    star_centers.sort(key=lambda x: x[1], reverse=True)
    other_centers.sort(key=lambda x: x[2], reverse=True)
    
    return star_centers, other_centers


semi_stars, non_stars = find_semi_stars(G)

with open('LCL_semi_star_centers.bed', 'w') as f:
    for node, degree, prots in semi_stars:
        chrom, start, end = node
        f.write(f"{chrom}\t{start}\t{end}\t{degree}\t{prots}\n")

with open('LCL_non_star_centers.bed', 'w') as f:
    for node, degree, con_count, prots in non_stars:
        chrom, start, end = node
        f.write(f"{chrom}\t{start}\t{end}\t{degree}\t{con_count}\t{prots}\n")




