"""
This code aims to create a figure for print on 8.5x11" paper in portrait
It will plot half of the matrix, rotated 45 degrees with the diagonal running up and down
It will stagger labels in and out, and plot a dendogram
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering
import os
from io import BytesIO
from PIL import Image

file = '/mnt/altnas/work/Kyle.Knightly/enrichments-pseudo-trans-contacts.tsv'
in_df = pd.read_csv(file, sep='\t', index_col=0)
max = in_df.max().max()
#log, replace 0 with a 1x-10
log_df = np.log10(in_df.replace(0, 10 ** (-np.log10(max))))

#calc vmin vmax values
vmin = log_df.min().min()
vmax = log_df.max().max()
#vmin = -vmax

linkage_matrix = linkage(in_df, method='ward')

# Apply Optimal Leaf Ordering to the linkage matrix
linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, in_df)

# Get the ordered indices after optimal leaf ordering
ordered_index = leaves_list(linkage_matrix_olo)

df = log_df.iloc[ordered_index, ordered_index]

# Mask for the upper triangular part (you can switch this to `np.triu` for the upper part)
mask = np.tril(np.ones(df.shape), k=-1)  # k=1 keeps the diagonal, set k=0 to remove it

# Create the heatmap
fig, ax = plt.subplots(figsize=(8.5, 11))

sns.heatmap(df, 
            mask=mask, 
            square=True, 
            cmap="RdBu_r", 
            annot=False, 
            cbar=False, 
            xticklabels=False, 
            yticklabels=False)

plt.axis('off')  # Turn off the axis for the plot to remove extra space

# Save the plot to a buffer as an image
buf = BytesIO()
plt.savefig(buf, format='png', bbox_inches='tight')
buf.seek(0)

# Open the image from the buffer and rotate it 45 degrees
im = Image.open(buf)
rotated_im = im.rotate(-45, expand=True)

# Create a new figure for displaying the rotated image
plt.figure(figsize=(8.5, 11))
plt.imshow(rotated_im)
plt.axis('off')  # Turn off axis for the rotated image
# Save the figure
plt.savefig('enrichments-pseudo-trans-contacts.png', dpi=800, bbox_inches='tight')

