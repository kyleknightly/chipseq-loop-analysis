import pandas as pd
import matplotlib.pyplot as plt

degfile = '/mnt/altnas/work/Kyle.Knightly/leq200mb-protein-degree-tests.tsv'
degdf= pd.read_csv(degfile, sep='\t', header=0)
# print(degdf)

starfile = '/mnt/altnas/work/Kyle.Knightly/leq200mb-prot-avg-starness.tsv'
stardf = pd.read_csv(starfile, sep='\t', names = ['protein', 'coef'])
# print(skipdf)

merged_df = pd.merge(degdf, stardf, on='protein')

# Create the scatter plot
plt.figure(figsize=(10, 6))
plt.rcParams['font.family'] = 'Times New Roman'

######
# Define three lists of proteins (you can replace these with your actual lists)
list_red = ['STAG1', 'CTCF', 'RAD21', 'SMC3']  # Proteins to color red
list_green = ['ZNF146', 'SCMH1', 'CENPX', 'ZMAT5', 'MLLT10', 'ZC3H4']  # Proteins to color green
list_blue = ['none']  # Proteins to color blue

# Create the scatter plot
plt.figure(figsize=(30, 20))

# Plot each group with its specific color
# plt.scatter(merged_df[merged_df['protein'].isin(list_red)]['avg_degree'], 
#             merged_df[merged_df['protein'].isin(list_red)]['coef'], 
#             color='red', label='Red Group', alpha=1)

# plt.scatter(merged_df[merged_df['protein'].isin(list_green)]['avg_degree'], 
#             merged_df[merged_df['protein'].isin(list_green)]['coef'], 
#             color='green', label='Green Group', alpha=0.7)

# plt.scatter(merged_df[merged_df['protein'].isin(list_blue)]['avg_degree'], 
#             merged_df[merged_df['protein'].isin(list_blue)]['coef'], 
#             color='blue', label='Blue Group', alpha=0.7)

# # Plot all other proteins in black
# other_proteins = merged_df[~merged_df['protein'].isin(list_red + list_green + list_blue)]
# plt.scatter(other_proteins['avg_degree'], other_proteins['coef'], 
#             color='black', label='Other Proteins', alpha=0.7)


######
plt.scatter(merged_df['avg_degree'], merged_df['coef'], alpha=0.7)
plt.scatter(merged_df[merged_df['protein'].isin(list_red)]['avg_degree'], 
            merged_df[merged_df['protein'].isin(list_red)]['coef'], 
            color='red', label='Red Group', alpha=1)

# Add labels for each point
for i, row in merged_df.iterrows():
    plt.text(row['avg_degree'], row['coef'], row['protein'], fontsize=9, ha='right')

# Customize the plot
plt.xlabel('Average Degree')
plt.ylabel('Average Clustering Coefficient')
plt.title(' ')
plt.grid(True)


plt.savefig('leq200mb-stariness-degree.png')