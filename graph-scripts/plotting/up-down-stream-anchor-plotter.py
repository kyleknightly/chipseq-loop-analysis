"""
This plots a scatter plot where each point is an anchor,
The point on the plot is given by how many times it is up or downstream in a loop
It also will color CTCF anchors red
"""

import pandas as pd
import ast
import matplotlib.pyplot as plt

file = '/mnt/altnas/work/Kyle.Knightly/hepg2-loopList_CTCF_noCTCF.bedpe'
df = pd.read_csv(file, sep='\t', names = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'class'])
# df['prots1'] = df['prots1'].apply(ast.literal_eval)
# df['prots2'] = df['prots2'].apply(ast.literal_eval)

# Convert lists to sets in the 'prots1' and 'prots2' columns
# df['prots1'] = df['prots1'].apply(tuple)
# df['prots2'] = df['prots2'].apply(tuple)

anchors = {}
for idx, row in df.iterrows():
    if (row['chr1'], row['start1'], row['end1']) not in anchors.keys():
        if row['class']==1 or row['class']==2:
            anchors[(row['chr1'], row['start1'], row['end1'])]=1
    if (row['chr2'], row['start2'], row['end2']) not in anchors.keys():
        if row['class']==1 or row['class']==3:
            anchors[(row['chr2'], row['start2'], row['end2'])]=1
        

# Step 1: Count occurrences for each set
count1 = df.groupby(['chr1', 'start1', 'end1']).size().reset_index(name='count_as_1')
count2 = df.groupby(['chr2', 'start2', 'end2',]).size().reset_index(name='count_as_2')

# Renaming columns for consistency
count1.rename(columns={'chr1': 'chr', 'start1': 'start', 'end1': 'end'}, inplace=True)
count2.rename(columns={'chr2': 'chr', 'start2': 'start', 'end2': 'end'}, inplace=True)

# Step 2: Merge the counts
merged_counts = pd.merge(count1, count2, on=['chr', 'start', 'end'], how='outer').fillna(0)

# Generate the key column for lookup
merged_counts['key'] = list(zip(merged_counts['chr'], merged_counts['start'], merged_counts['end']))

def get_color(key):
    value = anchors.get(key, None)
    if value == 1:
        return 'red'
    elif value == 2:
        return 'yellow'
    elif value == 3:
        return 'orange'
    elif value == 4:
        return 'green'
    else:
        return 'blue'  # Default color if the key is not found in the anchors

merged_counts['color'] = merged_counts['key'].apply(get_color)

# Step 4: Plot the scatterplot with colors
plt.figure(figsize=(6, 6))
plt.scatter(merged_counts['count_as_1'], merged_counts['count_as_2'], c=merged_counts['color'], alpha=0.3)
plt.xlabel('Upstream')
plt.ylabel('Downstream')
plt.title('Scatterplot of Up vs Down-stream anchors, CTCF in red')
plt.grid(True)
plt.savefig('up-down-stream-anchors.png')

