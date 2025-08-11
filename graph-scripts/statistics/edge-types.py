import pandas as pd
from collections import Counter
from itertools import product
import pprint

# Load BEDPE file
bedpe_path = '/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/hepg2-ccre-loop-types.bedpe'
cols = ["chr1", "start1", "end1", "type1", "chr2", "start2", "end2", "type2"]
df = pd.read_csv(bedpe_path, sep="\t", names=cols)

# # Normalize type strings (strip any prefix like "CA-")
# df["type1"] = df["type1"].str.replace(r'^.*?-', '', regex=True)
# df["type2"] = df["type2"].str.replace(r'^.*?-', '', regex=True)

# Total number of loops
num_loops = len(df)

# Get marginal frequencies
left_counts = df["type1"].value_counts()
right_counts = df["type2"].value_counts()

P_left = (left_counts / num_loops).to_dict()
P_right = (right_counts / num_loops).to_dict()

# Compute expected counts under independence
all_types = set(P_left) | set(P_right)
expected_counts = {}

for a, b in product(all_types, repeat=2):
    if a > b:
        continue  # only consider unordered pairs once
    p = P_left.get(a, 0) * P_right.get(b, 0)
    if a != b:
        p += P_left.get(b, 0) * P_right.get(a, 0)
    expected_counts[(a, b)] = num_loops * p

# Compute observed counts
def unordered_pair(a, b):
    return tuple(sorted([a, b]))

df["edge_type"] = df.apply(lambda row: unordered_pair(row["type1"], row["type2"]), axis=1)
observed_counts = df["edge_type"].value_counts().to_dict()

# Compute enrichment
enrichment = {}
for pair in expected_counts:
    obs = observed_counts.get(pair, 0)
    exp = expected_counts[pair]
    enrichment[pair] = obs / exp if exp > 0 else float('nan')

# Display summary
print(f"\nTotal loops: {num_loops}")
print("\nTop 20 Enrichments (Observed / Expected):")
for pair, enr in sorted(enrichment.items(), key=lambda x: -x[1])[:20]:
    obs = observed_counts.get(pair, 0)
    exp = expected_counts[pair]
    print(f"{pair[0]}–{pair[1]}: Enrichment = {enr:.2f}  (Obs = {obs}, Exp = {exp:.1f})")

# Optional: Save full enrichment table
out_df = pd.DataFrame([
    {"type1": a, "type2": b, "observed": observed_counts.get((a, b), 0),
     "expected": expected_counts.get((a, b), 0), "enrichment": enrichment.get((a, b), float('nan'))}
    for (a, b) in expected_counts
])
out_df = out_df.sort_values("enrichment", ascending=False)
out_df.to_csv("anchor_type_enrichment.tsv", sep="\t", index=False)
print("\nSaved enrichment table to anchor_type_enrichment.tsv")
