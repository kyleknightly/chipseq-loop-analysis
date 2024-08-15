"""
This plots a networkx graph of the genome by first splitting it into chromosomes and arranging
those in a circle. Again, this is not very good
"""

import networkx as nx
import matplotlib.pyplot as plt
import os
import numpy as np

print("started")

file_path = "/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-graph.gpickle"
G = nx.read_gpickle(file_path)

# Extract the chr attribute from each node
chr_groups = {}
for node, data in G.nodes(data=True):
    chr_value = data['chr']
    if chr_value not in chr_groups:
        chr_groups[chr_value] = []
    chr_groups[chr_value].append(node)
print(chr_groups.keys())
# Create a spring layout for each chr group
group_pos = {}
group_sizes = {}
for chr_value, nodes in chr_groups.items():
    subgraph = G.subgraph(nodes)
    sub_pos = nx.spring_layout(subgraph)
    group_pos[chr_value] = sub_pos
    group_sizes[chr_value] = len(nodes)

chromosome_order = [f'chr{num}' for num in range(1, 23)]
chromosome_order.append('chrX')
print(chromosome_order)

# Calculate total size to determine angle distribution
total_size = sum(group_sizes.values())

# Arrange the groups in a circle with angles proportional to their size
cumulative_angle = np.pi/2
group_centers = {}
radius = 10  # Base radius of the circle
num_groups = len(chr_groups)
print(list(enumerate(chromosome_order)))
for i, chr_value in enumerate(chromosome_order):
    angle = (np.pi/2)- ((2 * np.pi * i) / num_groups) # Adjusted to start from 12 o'clock
    print (chr_value, angle)
    center_x = radius * np.cos(angle)
    center_y = radius * np.sin(angle)
    group_centers[chr_value] = (center_x, center_y)
print(group_centers)

# for chr_value, size in group_sizes.items():
#     # Calculate the angle step for this group
#     angle_step = 2 * np.pi / num_groups
#     cumulative_angle += angle_step  # Place the center of the group in the middle of its allocated angle
    
#     center_x = radius * np.cos(cumulative_angle)
#     center_y = radius * np.sin(cumulative_angle)
#     group_centers[chr_value] = (center_x, center_y)
    
#     cumulative_angle += angle_step  # Move to the next group's angle

# Calculate the scaling factor for each subgraph based on its size
max_size = max(group_sizes.values())
scale_factors = {chr_value: np.sqrt(size / max_size) for chr_value, size in group_sizes.items()}

# Offset each group's positions based on their circular arrangement and scale them
pos = {}
for chr_value, nodes in chr_groups.items():
    center_x, center_y = group_centers[chr_value]
    sub_pos = group_pos[chr_value]
    scale = scale_factors[chr_value]
    for node, (x, y) in sub_pos.items():
        pos[node] = (center_x + scale * x, center_y + scale * y)

for node, data in G.nodes(data=True):
    data['color'] = 'red' if 'CTCF' in data['prots'] else 'blue'

node_colors = [data['color'] for _, data in G.nodes(data=True)]

# Plotting
plt.figure(figsize=(80, 80))
nx.draw(G, pos, with_labels=False, node_size=10, alpha=0.7, node_color=node_colors, edge_color="gray", width=0.5)
print
for chr_value in group_centers.keys():
    print(chr_value, group_centers[chr_value])
    center_x, center_y = group_centers[chr_value]
    plt.text(center_x, center_y+1.5, chr_value, fontsize=12, ha='center', va='center', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
# Save the plot
name = os.path.basename(file_path)
plt.savefig('1grouped-' + name[:-7] + 'png')
plt.show()
