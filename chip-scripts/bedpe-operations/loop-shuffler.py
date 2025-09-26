#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

# --------- inputs/outputs ----------
in_file  = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe'
out_file = '/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/rand-3/rand-3.bedpe'
seed = None  # set None for nondeterministic

# --------- load ---------
df = pd.read_csv(in_file, sep='\t', names=['c1','s1','e1','c2','s2','e2'], index_col=False)

# Represent each anchor by a stable ID string; also keep a mapping to coordinates
def anchor_id(c, s, e):
    return f'{c}:{int(s)}-{int(e)}'

left_ids  = df.apply(lambda r: anchor_id(r.c1, r.s1, r.e1), axis=1).to_numpy()
right_ids = df.apply(lambda r: anchor_id(r.c2, r.s2, r.e2), axis=1).to_numpy()

# Universe of anchors (nodes)
all_ids = np.concatenate([left_ids, right_ids])
# Degree of each anchor = number of times it appears across both ends
unique_ids, counts = np.unique(all_ids, return_counts=True)
id2deg = dict(zip(unique_ids, counts))
ids = unique_ids

# Sanity: total stubs must be even and equal to 2 * |E|
num_edges = len(df)
total_stubs = counts.sum()
assert total_stubs == 2 * num_edges, "Total stubs must equal 2 * number of edges."

# Build the stub list: one entry per stub containing the anchor id
stubs = np.repeat(ids, counts)  # length = total_stubs

rng = np.random.default_rng(seed)

def random_matching_no_self(stubs, max_restarts=50):
    """
    Random perfect matching of stubs, avoiding self-loops.
    Start from a random permutation and fix same-anchor pairs by swapping.
    """
    n = len(stubs)
    assert n % 2 == 0
    for _ in range(max_restarts):
        perm = stubs.copy()
        rng.shuffle(perm)

        # Pair as (0,1), (2,3), ...
        i = 0
        ok = True
        while i < n:
            a = perm[i]
            b = perm[i+1]
            if a == b:
                # Find a later stub to swap with perm[i+1] that avoids creating self-pairs in BOTH pairs
                swapped = False
                for j in range(i+2, n):
                    # Try swap perm[i+1] <-> perm[j]
                    c = perm[j]
                    if (c != a) and ( (j % 2 == 0 and perm[j+1] != b) or (j % 2 == 1 and perm[j-1] != b) ):
                        # Do the swap
                        perm[i+1], perm[j] = perm[j], perm[i+1]
                        swapped = True
                        break
                if not swapped:
                    ok = False
                    break
            i += 2

        if ok:
            return perm  # a valid pairing with no self-loops

    raise RuntimeError("Failed to construct a self-loop-free matching after several restarts. "
                       "Graph may be too constrained; consider allowing self-loops or increasing max_restarts.")

# Get a random perfect matching of the stubs with no self-loops
perm = random_matching_no_self(stubs)

# Build the new edge list from the paired permutation
pairs = [(perm[i], perm[i+1]) for i in range(0, len(perm), 2)]

# Optional: keep output symmetric order (e.g., as given). We just write as-is.
# Convert anchor IDs back to BEDPE rows by splitting the "chr:start-end" string
def id_to_fields(aid):
    chrom, rest = aid.split(':', 1)
    s, e = rest.split('-', 1)
    return chrom, int(s), int(e)

new_rows = []
for a_id, b_id in pairs:
    c1, s1, e1 = id_to_fields(a_id)
    c2, s2, e2 = id_to_fields(b_id)
    new_rows.append((c1, s1, e1, c2, s2, e2))

out = pd.DataFrame(new_rows, columns=['c1','s1','e1','c2','s2','e2'])
out.to_csv(out_file, sep='\t', header=False, index=False)
print(f"Wrote configuration-model rewired loops (no self-loops) to:\n  {out_file}")