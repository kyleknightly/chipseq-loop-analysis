import os
import re
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === CONFIGURE THIS ===
directory = "/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/sorted-loops"

# === REGEX to match motif result lines ===
motif_pattern = re.compile(r'^(\d+)\s+[\d.]+\s+[\d.+-]+\s+[\d.+-]+\s+[\d.]+\s+\d+\s+[\d.]+')

# === Parse motif IDs per chromosome ===
motif_sets = defaultdict(set)

for filename in os.listdir(directory):
    if filename.startswith("chr") and filename.endswith("_OUT.txt"):
        chrom = filename.split("_")[0]  # e.g., "chr1"
        filepath = os.path.join(directory, filename)
        with open(filepath, 'r') as file:
            for line in file:
                match = motif_pattern.match(line.strip())
                if match:
                    motif_id = int(match.group(1))
                    motif_sets[chrom].add(motif_id)

# === Summary metrics ===
all_chroms = sorted(motif_sets)
all_motifs = set.union(*motif_sets.values())
common_motifs = set.intersection(*motif_sets.values())
unique_counts = {c: len(motif_sets[c]) for c in all_chroms}

print(f"Total motif IDs found: {len(all_motifs)}")
print(f"Motif IDs common to all chromosomes: {len(common_motifs)}")
print(f"Unique motif counts per chromosome:\n{unique_counts}")

# === Make presence/absence matrix ===
presence_matrix = pd.DataFrame(index=sorted(all_motifs), columns=all_chroms, data=0)

for chrom in all_chroms:
    presence_matrix.loc[presence_matrix.index.isin(motif_sets[chrom]), chrom] = 1

# === Save matrix to file (optional) ===
presence_matrix.to_csv("motif_presence_matrix.csv")

# === Visualize as heatmap ===
plt.figure(figsize=(14, 10))
sns.heatmap(presence_matrix, cmap="Greens", cbar_kws={"label": "Motif Present"})
plt.title("Motif ID Presence Across Chromosomes")
plt.xlabel("Chromosome")
plt.ylabel("Motif ID")
plt.tight_layout()
plt.savefig("motif_presence_heatmap.png", dpi=300)
