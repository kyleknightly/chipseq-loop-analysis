import pandas as pd
import matplotlib.pyplot as plt

degfile = '/mnt/altnas/work/Kyle.Knightly/anchor-graph/hepg2/protein-degree-tests.tsv'
degdf= pd.read_csv(degfile, sep='\t', header=0)
# print(degdf)

skipfile = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/stats/protein-avg-skips.tsv'
skipdf = pd.read_csv(skipfile, sep='\t', header=0)
# print(skipdf)

merged_df = pd.merge(degdf, skipdf, on='protein')

# Create the scatter plot
plt.figure(figsize=(30, 20))
plt.scatter(merged_df['avg_degree'], merged_df['skips'], alpha=0.7)

# Add labels for each point
for i, row in merged_df.iterrows():
    plt.text(row['avg_degree'], row['skips'], row['protein'], fontsize=9, ha='right')

# Customize the plot
plt.xlabel('Average Degree')
plt.ylabel('Skips')
plt.title('Skippiness and Average Degree of Each Protein')
plt.grid(True)


plt.savefig('skips-degree.png')