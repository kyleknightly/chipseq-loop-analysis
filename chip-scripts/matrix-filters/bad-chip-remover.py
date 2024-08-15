"""
Takes a matrix (heatmap of protein interactions of some sort) and removes proteins
with <peak_filter peaks in it's file.

This iterates through a chipdir, which holds the peaks for each protein
It expects the names to be "PROTEIN-whatever.bed"
"""

import os
import pandas as pd

peak_filter=250
chipdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/chipbin10beds_merged'
peak_counts = {}
for bed_filename in os.listdir(chipdir):
    if bed_filename.endswith('.bed'): #verify filetype
        bed_filepath = os.path.join(chipdir, bed_filename) #get path
        with open(bed_filepath, 'r') as bed_file: #open
            num = sum(1 for line in bed_file)
            peak_counts[bed_filename.split('-')[0]] = num

remove = []
for key, value in peak_counts.items():
    if value <= peak_filter:
        remove.append(key)

print(remove)
print(len(remove))

matrix = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/loop-ends/loop-end-enrichments-pseudo-trans-contacts-paired-anchor-TFs.bed'
df = pd.read_csv(matrix, sep='\t', header=0, index_col=0)
remove_rows = [item for item in remove if item in df.index]
remove_columns = [item for item in remove if item in df.columns]

# Remove the specified rows and columns
df = df.drop(index=remove_rows, columns=remove_columns)

name = os.path.basename(matrix)
df.to_csv('geq250-' + name, sep='\t', index=True)
