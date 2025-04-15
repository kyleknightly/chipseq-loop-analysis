"""
To use as a control and a substitute to using anchors, this subdivides the genome into sections of a certain length
"""

import pandas as pd

lengths = {'chr1': 248956422, 'chr10': 133520353, 'chr11': 135074598, 'chr12': 133181459, 'chr13': 114226953, 'chr14': 105352939, 'chr15': 101724648, 'chr16': 90019594, 'chr17': 82917884, 'chr18': 80247771, 'chr19': 58575455, 'chr2': 241878248, 'chr20': 64120656, 'chr21': 46286522, 'chr22': 50691732, 'chr3': 197943179, 'chr4': 189904718, 'chr5': 181261319, 'chr6': 170554363, 'chr7': 158977292, 'chr8': 144997812, 'chr9': 137618973, 'chrX': 155149106}
total_length = sum(lengths.values())
span_length = 100000

bed_entries = []
# Generate spans for each chromosome
for chr_name, chr_length in lengths.items():
    for start in range(1, chr_length + 1, span_length):
        end = min(start + span_length - 1, chr_length)
        bed_entries.append([chr_name, start, end])

# Convert to DataFrame
bed_df = pd.DataFrame(bed_entries, columns=['chromosome', 'start', 'end'])

# Save to BED file
bed_df.to_csv('subdivisions.bed', sep='\t', header=False, index=False)