"""
This takes a bed file with each anchor's group and the genome loop file, and assigns a color to each loop
If a loop is between two different groups it is colored black

The output is a bedpe for use in Juicebox
"""

import pandas as pd

colors = {
    'cyan': "0,255,255",
    'deep_blue': "0,0,255",
    'dark_green': "0,128,0",
    'purple': "128,0,128",
    'orange': "255,165,0",
    'magenta': "255,0,255",
    'yellow': "255,255,0",
    'gray': "128,128,128",
    'lime': "50,205,50",
    'sky_blue': "135,206,235",
    'olive': "128,128,0",
    'teal': "0,128,128",
    'navy': "0,0,128",
    'maroon': "128,0,0",
    'peach': "255,218,185",
    'coral': "255,127,80",
    'mint': "189,252,201",
    'lavender': "230,230,250",
    'rose': "255,228,225",
    'beige': "245,245,220",
    'mustard': "255,219,88",
    'turquoise': "64,224,208",
    'viridian': "64,130,109",
    'cerulean': "42,82,190",
    'auburn': "165,42,42",
    'amber': "255,191,0",
    'chartreuse': "127,255,0",
    'azure': "240,255,255",
    'violet': "238,130,238",
    'indigo': "75,0,130"
}

file = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2-loops.bed'
bedpe = pd.read_csv(file, sep=' ', names = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
groupfile = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/infomap_communities.bed'
anchors = pd.read_csv(groupfile, sep='\t', names = ['chr', 'start', 'end', 'group'])
#print(bedpe)
#print(anchors)
# Merge bedpe with anchors on the first set of coordinates
bedpe_merged_1 = pd.merge(bedpe, anchors, how='left', left_on=['chr1', 'start1', 'end1'], right_on=['chr', 'start', 'end'])
bedpe_merged_1 = bedpe_merged_1.rename(columns={'group': 'group1'}).drop(columns=['chr', 'start', 'end'])

# Merge the result with anchors on the second set of coordinates
bedpe_merged_2 = pd.merge(bedpe_merged_1, anchors, how='left', left_on=['chr2', 'start2', 'end2'], right_on=['chr', 'start', 'end'])
bedpe_merged_2 = bedpe_merged_2.rename(columns={'group': 'group2'}).drop(columns=['chr', 'start', 'end'])

#print(bedpe_merged_2)
# Determine if groups are the same or different
bedpe_merged_2['groups'] = bedpe_merged_2.apply(lambda row: row['group1'] if row['group1'] == row['group2'] else 'DIFFERENT', axis=1)

# Drop intermediate group columns
bedpe_final = bedpe_merged_2.drop(columns=['group1', 'group2'])

# Display the final DataFrame
#print(bedpe_final)

# Function to map group number to dictionary value using modulo
def map_group_to_value(number, dict):
    modulo_value = number % len(dict)
    #print(modulo_value)
    #print(list(dict.keys()))
    return dict[list(dict.keys())[modulo_value]]

# Add the new column based on the groups modulo operation
bedpe_final['mod_group'] = bedpe_final['groups'].apply(lambda x: "0,0,0" if x == "DIFFERENT" else map_group_to_value(int(x), colors))

print(bedpe_final)

new_df = bedpe_final[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]
# Create a DataFrame with four columns containing periods
periods_df = pd.DataFrame('.', index=new_df.index, columns=['.', '.', '.', '.'])
# Add the mod_group column
mod_group_df = bedpe_final[['mod_group']]

final = pd.concat([new_df, periods_df, mod_group_df], axis=1)

# Display the final DataFrame
print(final)

# Optionally, save the resulting DataFrame to a file
final.to_csv('colored-infomap.bedpe', sep='\t', header=False, index=False)