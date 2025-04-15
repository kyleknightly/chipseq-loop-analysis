import pandas as pd
import ast
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.colors import LinearSegmentedColormap

# Define the colors for the colormap
colors = [(0, 0, 0), (1, 0, 0)]  # white to red

# Create the colormap
white_red = LinearSegmentedColormap.from_list("white_red", colors)

# Load data
file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep='\t', names=['chr1', 'start1', 'end1', 'prots1', 'chr2', 'start2', 'end2', 'prots2'])

# Convert string representation of lists to actual tuples
df['prots1'] = df['prots1'].apply(ast.literal_eval).apply(tuple)
df['prots2'] = df['prots2'].apply(ast.literal_eval).apply(tuple)

# Initialize dictionaries to count upstream and downstream occurrences
upstream_counts = defaultdict(int)
downstream_counts = defaultdict(int)

prot = "CTCF"

# Count occurrences for upstream and downstream, filtering for 'CTCF'
for index, row in df.iterrows():
    anchor_up = (row['chr1'], row['start1'], row['end1'])
    anchor_down = (row['chr2'], row['start2'], row['end2'])
    
    # Count only if 'CTCF' is in the respective protein list
    if prot in row['prots1']:
        upstream_counts[anchor_up] += 1
    if prot in row['prots2']:
        downstream_counts[anchor_down] += 1

# Combine the counts into a single DataFrame for plotting
combined_counts = pd.DataFrame(
    [(key, upstream_counts[key], downstream_counts[key]) for key in set(upstream_counts) | set(downstream_counts)],
    columns=['anchor', 'upstream_count', 'downstream_count']
)

# Count the number of occurrences for each (upstream_count, downstream_count) pair
combined_counts['coordinate'] = combined_counts.apply(lambda row: (row['upstream_count'], row['downstream_count']), axis=1)
coordinate_counts = combined_counts['coordinate'].value_counts().reset_index()
coordinate_counts.columns = ['coordinate', 'count']

# Separate the coordinate counts into x and y for plotting
coordinate_counts['x'] = coordinate_counts['coordinate'].apply(lambda x: x[0])
coordinate_counts['y'] = coordinate_counts['coordinate'].apply(lambda x: x[1])

# Plot the scatter plot with a colormap based on the count
plt.figure(figsize=(8, 8))
scatter = plt.scatter(coordinate_counts['x'], coordinate_counts['y'], 
                      c=coordinate_counts['count'], cmap='plasma_r', s=50, edgecolor='none')
plt.colorbar(scatter, label='Number of Anchors')
plt.xlabel('Upstream Count')
plt.ylabel('Downstream Count')
plt.title('Scatter Plot of Upstream vs Downstream Anchors')
plt.grid(True)
plt.savefig(prot+'-up-down-stream-anchors-plasma.png')
plt.show()
