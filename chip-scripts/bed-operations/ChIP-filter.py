"""
This goes through ChIP seq experiments and uses their scores to filter for top n%
"""

import os
import numpy as np

chipdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/chip-tracks'
outdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/new-filtered-tracks'

# Ensure the output directory exists (create if it doesn't)
os.makedirs(outdir, exist_ok=True)

# Iterate through all files in chipdir
for filename in os.listdir(chipdir):
    # Filter to files ending with narrowPeak or bed (adjust as needed)
    if not (filename.endswith('.narrowPeak') or filename.endswith('.bed')):
        continue
    
    input_path = os.path.join(chipdir, filename)
    lines = []
    scores = []

    # Read the file, store lines and corresponding score (5th column)
    with open(input_path, 'r') as f:
        for line in f:
            # Skip empty/comment lines
            if not line.strip() or line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            # The "score" in narrowPeak is column index 4 (5th column)
            print(parts[6])
            try:
                score = float(parts[6])
                print(score)
            except ValueError:
                print('failed')
                # If we fail to parse the score, just skip this line
                continue

            lines.append(parts)
            scores.append(score)

    # If file was empty or no valid scores, skip
    if not scores:
        continue

    # Calculate the 90th percentile (top 10%)
    cutoff = np.percentile(scores, 90)

    # Filter lines that meet or exceed the cutoff
    filtered_lines = [l for l in lines if float(l[4]) >= cutoff]

    # Write filtered lines to output directory
    output_path = os.path.join(outdir, filename)
    with open(output_path, 'w') as out_f:
        for parts in filtered_lines:
            out_f.write('\t'.join(parts) + '\n')