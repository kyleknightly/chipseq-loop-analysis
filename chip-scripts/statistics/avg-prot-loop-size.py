import pandas as pd
from collections import defaultdict

# Load maps
anchor_map = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchors.map", sep="\t", header=None, names=["coord", "anchor_id"])
coord_to_anchor = dict(zip(anchor_map["coord"], anchor_map["anchor_id"]))

tf_map = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/looppi/hepg2/anchor_tf_list.tsv", sep="\t", header=None, names=["anchor_id", "tf_list"])
anchor_to_tfs = {
    row["anchor_id"]: row["tf_list"].split("|")
    for _, row in tf_map.iterrows()
}

# Function to get anchor ID from loop end
def coord_to_id(chrom, start, end):
    return coord_to_anchor.get(f"{chrom}_{start}_{end}", None)

# Track loop sizes per TF
tf_to_loop_sizes = defaultdict(list)

# Parse BEDPE
with open("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe") as f:
    for line in f:
        chrom1, start1, end1, chrom2, start2, end2 = line.strip().split()[:6]
        a1 = coord_to_id(chrom1, start1, end1)
        a2 = coord_to_id(chrom2, start2, end2)
        if not a1 or not a2:
            continue
        tfs = set(anchor_to_tfs.get(a1, []) + anchor_to_tfs.get(a2, []))
        midpoint1 = (int(start1) + int(end1)) // 2
        midpoint2 = (int(start2) + int(end2)) // 2
        span = abs(midpoint2 - midpoint1)
        for tf in tfs:
            tf_to_loop_sizes[tf].append(span)

# Compute averages
tf_avg_sizes = {
    tf: sum(sizes)/len(sizes)
    for tf, sizes in tf_to_loop_sizes.items()
}

# Output
df = pd.DataFrame(list(tf_avg_sizes.items()), columns=["TF", "average_loop_size"])
df.to_csv("tf_average_loop_sizes.tsv", sep="\t", index=False)