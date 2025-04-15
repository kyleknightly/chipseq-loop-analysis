"""
This takes a bed file of anchors and their degrees and the loop list, and
creates a bedpe for use in Juicebox where the degrees are colored with a cmap
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


file = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2-loops.bed'
bedpe = pd.read_csv(file, sep=' ', names = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
groupfile = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/anchor-degrees.bed'
anchors = pd.read_csv(groupfile, sep='\t', names = ['chr', 'start', 'end', 'degree'])
# print(bedpe)
# print(anchors)
# Merge bedpe with anchors on the first set of coordinates
bedpe_merged_1 = pd.merge(bedpe, anchors, how='left', left_on=['chr1', 'start1', 'end1'], right_on=['chr', 'start', 'end'])
bedpe_merged_1 = bedpe_merged_1.rename(columns={'degree': 'degree1'}).drop(columns=['chr', 'start', 'end'])

# Merge the result with anchors on the second set of coordinates
bedpe_merged_2 = pd.merge(bedpe_merged_1, anchors, how='left', left_on=['chr2', 'start2', 'end2'], right_on=['chr', 'start', 'end'])
bedpe_merged_2 = bedpe_merged_2.rename(columns={'degree': 'degree2'}).drop(columns=['chr', 'start', 'end'])

#Average degree
bedpe_merged_2['avg_degree'] = bedpe_merged_2[['degree1', 'degree2']].mean(axis=1)
print(bedpe_merged_2)

# Normalize the avg_degree column to range between 0 and 1
norm = matplotlib.colors.Normalize(vmin=bedpe_merged_2['avg_degree'].min(), vmax=bedpe_merged_2['avg_degree'].max())

# Map the normalized values to colormap
cmap = plt.cm.rainbow
bedpe_merged_2['color'] = bedpe_merged_2['avg_degree'].apply(lambda x: cmap(norm(x)))

# Convert RGBA to RGB and format as "R,G,B"
bedpe_merged_2['color'] = bedpe_merged_2['color'].apply(lambda x: f"{int(x[0]*255)},{int(x[1]*255)},{int(x[2]*255)}")




new_df = bedpe_merged_2[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]
# Create a DataFrame with four columns containing periods
periods_df = pd.DataFrame('.', index=new_df.index, columns=['.', '.', '.', '.'])
# Add the mod_group column
mod_group_df = bedpe_merged_2[['color']]

final = pd.concat([new_df, periods_df, mod_group_df], axis=1)

# Display the final DataFrame
print(final)

# Optionally, save the resulting DataFrame to a file
final.to_csv('colored-degrees.bedpe', sep='\t', header=False, index=False)