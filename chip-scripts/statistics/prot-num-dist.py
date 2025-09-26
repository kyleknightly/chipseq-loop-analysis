# --- configuration you can edit ---
FILE = "/mnt/altnas/work/Kyle.Knightly/looppi/out/annotated_anchors.bed"                     # path to your TSV
UNIQUE_PROTEINS_PER_ROW = False        # True = de-duplicate proteins within a row
SAVE_PNG = False                        # set False to skip plotting
OUT_TSV = None                         # None => auto name like input_protein_list_sizes.tsv
OUT_PNG = 'loop-anchor-protein-list-sizes.png'                         # None => auto name like input_protein_list_sizes.png
# -----------------------------------

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Read TSV: dtype=str + keep_default_na=False so "NA"/"." stay as strings, not NaN
df = pd.read_csv(FILE, sep="\t", dtype=str, keep_default_na=False)

# The proteins list is assumed to be the LAST column
proteins_col = df.columns[-1]
s = df[proteins_col].astype(str).str.strip()

MISSING_MARKERS = {"", "NA", "N/A", ".", "-", "NONE", "NULL"}

def count_list(cell: str) -> int:
    up = cell.strip().upper()
    if up in MISSING_MARKERS:
        return 0
    toks = [t.strip() for t in cell.split("|") if t.strip() != ""]
    if UNIQUE_PROTEINS_PER_ROW:
        toks = list(set(toks))
    return len(toks)

sizes = s.apply(count_list)

# Build a tidy distribution table: ensure keys 0..max appear even if missing
if sizes.empty:
    dist = pd.Series([0], index=[0], name="count")
else:
    dist = sizes.value_counts().sort_index()
    dist = dist.reindex(range(0, sizes.max() + 1), fill_value=0)
dist.index.name = "list_size"
dist.name = "count"

# Print to console
total = 0
print("list_size\tcount")
for k, v in dist.items():
    total+=v
    print(f"{k}\t{v}")
print(total)

# Save counts TSV
in_path = Path(FILE)
if OUT_TSV is None:
    OUT_TSV = str(in_path.with_suffix("")) + "_protein_list_sizes.tsv"
dist.to_csv(OUT_TSV, sep="\t")

# Plot histogram PNG
if SAVE_PNG:
    if OUT_PNG is None:
        OUT_PNG = str(in_path.with_suffix("")) + "_protein_list_sizes.png"
    plt.figure(figsize = (40,5))
    plt.bar(dist.index, dist.values)
    plt.yscale("log")  
    plt.xlabel("# proteins in final-column list")
    plt.ylabel("# rows")
    plt.title("Distribution of protein-list sizes")
    plt.xticks(dist.index, fontsize = 5, rotation=90)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    plt.close()

    print(f"\nSaved counts to: {OUT_TSV}\nSaved histogram to: {OUT_PNG}")
else:
    print(f"\nSaved counts to: {OUT_TSV}")
