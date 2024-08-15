"""
This takes a loop list and creates a networkx graph
You can edit the names parameter of the df definition to include proteins
In that case, uncomment the ast.literal_eval lines and the prot assignment in nodes and edges

Each node is a tuple of (chr, start, end), aand can hold other properties like span, proteins
Each edge can also hold properties like size, proteins at either anchor
"""

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import ast

# Define the function to read the BED file and create the graph

bed_file = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/LCL-loops/loopList_CTCF_noCTCF.bedpe'
df = pd.read_csv(bed_file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'ch2', 'start2', 'end2', 'ctcf'])
# df['prots1'] = df['prots1'].apply(ast.literal_eval)
# df['prots2'] = df['prots2'].apply(ast.literal_eval)
print('DF MADE')


G = nx.Graph()
# Add nodes and edges to the graph
for _, row in df.iterrows():
    node1 = (row['ch1'], row['start1'], row['end1'])
    node2 = (row['ch2'], row['start2'], row['end2'])
    
    G.add_node(node1, 
            #    prots = row['prots1'], 
               chr = row['ch1'])
    G.add_node(node2, 
            #    prots = row['prots2'], 
               chr = row['ch2'])
    G.add_edge(node1, node2,
                # prots = list(set(row['prots1']) | set(row['prots2'])), 
                size = abs(int(row['start2'])-int(row['end1']))
                )
    #print(abs(int(row['start2'])-int(row['end1'])))
# Display basic information about the graph
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")

nx.write_gpickle(G, 'LCL-anchor-graph.gpickle')

