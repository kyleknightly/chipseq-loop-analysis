"""
The goal of this is to take a paired-anchor-TFs file, and get an average of the
skips (how many anchors are passed up) each time a protein is at an anchor.
"""

import pandas as pd
import ast
from collections import defaultdict
import matplotlib.pyplot as plt 

file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/paired-anchor-TFs.bed'
df = pd.read_csv(file, sep="\t", header=None, names = ['ch1', 'start1', 'end1', 'prots1', 'ch2', 'start2', 'end2', 'prots2'])
df['prots1'] = df['prots1'].apply(ast.literal_eval)
df['prots2'] = df['prots2'].apply(ast.literal_eval)

# anchorfile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/anchors_sorted.tsv'
# anchordf = pd.read_csv(anchorfile, sep='\t', header =0)
# # Initialize looped_anchors column as an empty list
# anchordf['looped_anchors'] = [[] for _ in range(len(anchordf))]
# print(anchordf)

# def find_anchor_index(ch, start, end, anchor_df):
#     result = anchor_df[(anchor_df['ch'] == ch) & (anchor_df['start'] == start) & (anchor_df['end'] == end)]
#     if not result.empty:
#         return result.index[0]  # Return the index of the found anchor
#     else:
#         return None


# for index, row in df.iterrows():
#     print(index/len(df))
#     # Find the index of the first anchor in the anchordf DataFrame
#     index_anchor1 = find_anchor_index(row['ch1'], row['start1'], row['end1'], anchordf)
    
#     # Find the index of the second anchor in the anchordf DataFrame
#     index_anchor2 = find_anchor_index(row['ch2'], row['start2'], row['end2'], anchordf)
    
#     if index_anchor1 is not None and index_anchor2 is not None:
#         # Add the second anchor's index to the first anchor's looped_anchors column
#         anchordf.at[index_anchor1, 'looped_anchors'].append(index_anchor2)
        
#         # Add the first anchor's index to the second anchor's looped_anchors column
#         anchordf.at[index_anchor2, 'looped_anchors'].append(index_anchor1)

# anchordf['looped_anchors'] = anchordf['looped_anchors'].apply(lambda x: ','.join(map(str, x)) if x else '')

# anchordf.to_csv('looped-anchor-ids.tsv',sep='\t', index=False, header=True)
# print(anchordf)

anchordf = pd.read_csv('/mnt/altnas/work/Kyle.Knightly/looped-anchor-ids.tsv', sep='\t', header=0)
anchordf['looped_anchors'] = anchordf['looped_anchors'].apply(lambda x: list(map(int, x.split(','))))
def calculate_difference(row):
    max_anchor = max(row['looped_anchors'])  # Maximum value in looped_anchors
    anchor_index = row['anchor_index']       # anchor_index
    looped_len = len(row['looped_anchors'])  # Length of looped_anchors
    return abs(abs(max_anchor - anchor_index) - looped_len) 

anchordf['skips'] = anchordf.apply(calculate_difference, axis=1)
print(anchordf)

"""
plt.hist(anchordf['skips'], bins=70, edgecolor='black')

# Add titles and labels
plt.title('Skips per Anchor')
plt.xlabel('Skips')
plt.ylabel('Frequency')
plt.savefig('skips_per_anchor.png')
"""
# Step 1: Unpair the columns into a new dataframe
df1 = df[['ch1', 'start1', 'end1', 'prots1']].rename(columns={'ch1': 'ch', 'start1': 'start', 'end1': 'end', 'prots1': 'prots'})
df2 = df[['ch2', 'start2', 'end2', 'prots2']].rename(columns={'ch2': 'ch', 'start2': 'start', 'end2': 'end', 'prots2': 'prots'})

# Concatenate df1 and df2 to form a new unpaired dataframe
new_df = pd.concat([df1, df2], ignore_index=True)

# Step 2: Merge with the anchor dataframe on 'ch', 'start', and 'end' columns to get the 'skips' value
merged_df = pd.merge(new_df, anchordf[['ch', 'start', 'end', 'skips']], on=['ch', 'start', 'end'], how='left')

# Step 1: Explode the 'prots' column so that each protein gets its own row
df_exploded = merged_df.explode('prots')

# Step 2: Group by 'prots' and calculate the average 'skips' value for each protein
protein_avg_skips = df_exploded.groupby('prots')['skips'].mean().reset_index()

# Rename the columns for clarity
protein_avg_skips.columns = ['protein', 'average_skips']

protein_avg_skips_sorted = protein_avg_skips.sort_values(by='average_skips', ascending=False)

# Step 2: Plot a histogram of the average skips
plt.figure(figsize=(75, 10))
plt.bar(protein_avg_skips_sorted['protein'], protein_avg_skips_sorted['average_skips'])

# Add titles and labels
plt.title('Average Skips for Each Protein (Sorted)')
plt.xlabel('Proteins')
plt.ylabel('Average Skips')
plt.xticks(rotation=90, fontsize=6)  # Rotate x-axis labels for better readability

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label cutoff
plt.savefig('prot-avg-skips.png')
