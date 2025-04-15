"""
Creates a spearman matrix out of any other matrix
"""

import pandas as pd
import os

matrix = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/removing-rare-TFs/geq250-geq50-anchor-enrichments-pseudo-contacts-paired-anchor-TFs.bed'

# Load the contact matrix from a TSV file
contact_matrix = pd.read_csv(matrix, sep='\t', index_col=0)

# Calculate the Spearman correlation matrix
spearman_matrix = contact_matrix.corr(method='spearman')
name = os.path.basename(matrix)
# Save the Spearman correlation matrix to a new TSV file
spearman_matrix.to_csv('spearman-' + name, sep='\t')

