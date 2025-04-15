"""
To use as a control, this creates fake anchors evenly distributed around the genome
"""

import pandas as pd

lengths = {'chr1': 248956422, 'chr10': 133520353, 'chr11': 135074598, 'chr12': 133181459, 'chr13': 114226953, 'chr14': 105352939, 'chr15': 101724648, 'chr16': 90019594, 'chr17': 82917884, 'chr18': 80247771, 'chr19': 58575455, 'chr2': 241878248, 'chr20': 64120656, 'chr21': 46286522, 'chr22': 50691732, 'chr3': 197943179, 'chr4': 189904718, 'chr5': 181261319, 'chr6': 170554363, 'chr7': 158977292, 'chr8': 144997812, 'chr9': 137618973, 'chrX': 155149106}
total_length = sum(lengths.values())
total_anchors = 10000
chr_proportions = {chr: lengths[chr]/total_length for chr in lengths.keys()}
anchors_per_chr = {chr: int(chr_proportions[chr]*total_anchors) for chr in chr_proportions.keys()}
print(anchors_per_chr)
span_length = 1000

# Dictionary to store the spans for each chromosome
spans_dict = {}

for chr_name, chr_length in lengths.items():
    num_spans = anchors_per_chr[chr_name]
    
    # Calculate the separation distance between the start of each span
    separation_distance = (chr_length - span_length * num_spans) // (num_spans + 1)
    
    spans = []
    start = separation_distance  # Start the first span after the initial separation

    for i in range(num_spans):
        end = start + span_length
        spans.append((start, end))
        start = end + separation_distance  # Move to the next start position
    
    spans_dict[chr_name] = spans

bed_data = []
for chr_name, spans in spans_dict.items():
    for start, end in spans:
        bed_data.append([chr_name, start, end])

bed_df = pd.DataFrame(bed_data, columns=['chr', 'start', 'end'])

bed_df.to_csv('big-spans.bed', sep='\t', header=False, index=False)
