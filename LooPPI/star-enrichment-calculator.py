#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from tqdm.auto import tqdm

tqdm.pandas(desc="Computing (vectorized)")

# -----------------------------
# Config / input paths
# -----------------------------
indir = '/mnt/altnas/work/Kyle.Knightly/looppi/no-empty-out'
annot_loop_file = indir + 'annotated_loops.bedpe'
annot_anchor_file = indir + 'annotated_anchors.bed'
trans_count_file = indir + 'trans-contacts.tsv'
cis_count_file   = indir + 'cis_contacts.tsv'
prot_end_props_file    = indir + 'loop_end_proportions.tsv'
prot_anchor_props_file = indir + 'anchor_proportions.tsv'

out_cis   = indir + 'cis_enrichment_matrix.tsv'
out_trans = indir + 'trans_enrichment_matrix.tsv'

# -----------------------------
# Pseudocounts / eps
# -----------------------------
pseudo  = 1.0     # change if you want a different prior
epsilon = 1e-6    # floor for logs/products; keeps values strictly >0

def safe_div(numer, denom, p=pseudo):
    """Compute (numer+p)/(denom+p) elementwise with nan/inf safety."""
    numer = np.nan_to_num(numer, nan=0.0, posinf=0.0, neginf=0.0)
    denom = np.nan_to_num(denom, nan=0.0, posinf=0.0, neginf=0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (numer + p) / (denom + p)
    return out

# -----------------------------
# Load data
# -----------------------------
cis_count_long   = pd.read_csv(cis_count_file,   sep='\t', names=['p1','p2','count'], index_col=False)
trans_count_long = pd.read_csv(trans_count_file, sep='\t', names=['p1','p2','count'], index_col=False)

prot_anchor_props = pd.read_csv(prot_anchor_props_file, sep='\t', index_col=0, names=['count','prop'])
prot_end_props    = pd.read_csv(prot_end_props_file,    sep='\t', index_col=0, names=['count','prop'])

annot_loops = pd.read_csv(
    annot_loop_file, sep='\t',
    names=['c1','s1','e1','id1','t1','p1','c2','s2','e2','id2','t2','p2'], index_col=False
)

# -----------------------------
# Parse pipe-separated protein lists
# -----------------------------
def parse_list_col(x):
    if pd.isna(x) or x == "":
        return []
    return [t.strip() for t in str(x).split("|") if t.strip()]

annot_loops['p1'] = annot_loops['p1'].apply(parse_list_col)
annot_loops['p2'] = annot_loops['p2'].apply(parse_list_col)

# -----------------------------
# Proteins
# -----------------------------
prots = prot_end_props.index.tolist()
P = len(prots)
prot2idx = {p:i for i,p in enumerate(prots)}

# Align and sanitize anchor proportions (P_anchor), fill missing with 0
anchor_props = prot_anchor_props.reindex(prots)['prop'].fillna(0.0).to_numpy()

# Number of anchors (|V|)
anchors = len(pd.read_csv(annot_anchor_file, sep='\t', header=None))

# -----------------------------
# Helper functions
# -----------------------------
def attach_indices(df):
    df = df.copy()
    df['i1'] = df['p1'].map(prot2idx)
    df['i2'] = df['p2'].map(prot2idx)
    return df.dropna(subset=['i1','i2']).astype({'i1':int,'i2':int})

def keep_upper(df):
    return df.loc[df['i1'] <= df['i2']].copy()

cis_count_long   = attach_indices(cis_count_long)
trans_count_long = attach_indices(trans_count_long)
cis_count_long_u   = keep_upper(cis_count_long)
trans_count_long_u = keep_upper(trans_count_long)

# -----------------------------
# CIS enrichment (with pseudocount)
# -----------------------------
# E_V expected from anchor proportions
exp_cis_full = anchors * (anchor_props[:,None] * anchor_props[None,:])  # P x P

cis_counts = np.zeros((P,P), dtype=float)
if not cis_count_long_u.empty:
    g = cis_count_long_u.groupby(['i1','i2'], as_index=False)['count'].sum()
    cis_counts[g['i1'].to_numpy(), g['i2'].to_numpy()] = g['count'].to_numpy()

# Enrichment with pseudocounts
cis_enrich_upper = safe_div(cis_counts, exp_cis_full, p=pseudo)

cis_matrix = cis_enrich_upper + cis_enrich_upper.T
np.fill_diagonal(cis_matrix, np.diag(cis_enrich_upper))

# -----------------------------
# TRANS expected using CIS enrichment (epsilon-guarded logs)
# -----------------------------
# Replace non-finite with epsilon and floor to epsilon so log is defined
cis_arr_eps = np.clip(np.nan_to_num(cis_matrix, nan=epsilon, posinf=epsilon, neginf=epsilon), epsilon, None)

# Convert loop p1/p2 protein lists into index lists once
loop_A, loop_B = [], []
for _, row in annot_loops.iterrows():
    A = [prot2idx[p] for p in row['p1'] if p in prot2idx]
    B = [prot2idx[p] for p in row['p2'] if p in prot2idx]
    if len(A) == 0 and len(B) == 0:
        continue
    loop_A.append(np.array(A, dtype=int))
    loop_B.append(np.array(B, dtype=int))

trans_exp_upper = np.zeros((P,P), dtype=float)

for A_idx, B_idx in tqdm(zip(loop_A, loop_B), total=len(loop_A), desc="Computing trans expected (epsilon)"):
    # Compute log-space products to avoid collapse; zero out self-terms
    if A_idx.size > 0:
        subA = np.log(cis_arr_eps[:, A_idx])
        for pos, j in enumerate(A_idx):
            subA[j, pos] = 0.0
        vA = np.exp(subA.sum(axis=1))
    else:
        vA = np.ones(P)

    if B_idx.size > 0:
        subB = np.log(cis_arr_eps[:, B_idx])
        for pos, j in enumerate(B_idx):
            subB[j, pos] = 0.0
        vB = np.exp(subB.sum(axis=1))
    else:
        vB = np.ones(P)

    M = (vA[:,None] * vB[None,:]) + (vB[:,None] * vA[None,:])
    diag_fix = (vA * vB)
    np.fill_diagonal(M, np.diag(M) - diag_fix)

    iu = np.triu_indices(P)
    trans_exp_upper[iu] += M[iu]

# -----------------------------
# TRANS enrichment (with pseudocount)
# -----------------------------
trans_counts = np.zeros((P,P), dtype=float)
if not trans_count_long_u.empty:
    g = trans_count_long_u.groupby(['i1','i2'], as_index=False)['count'].sum()
    trans_counts[g['i1'].to_numpy(), g['i2'].to_numpy()] = g['count'].to_numpy()

trans_enrich_upper = safe_div(trans_counts, trans_exp_upper, p=pseudo)

trans_matrix = trans_enrich_upper + trans_enrich_upper.T
np.fill_diagonal(trans_matrix, np.diag(trans_enrich_upper))

# -----------------------------
# Save outputs
# -----------------------------
cis_df = pd.DataFrame(cis_matrix, index=prots, columns=prots)
trans_df = pd.DataFrame(trans_matrix, index=prots, columns=prots)

cis_df.to_csv(out_cis,   sep='\t', float_format="%.10g")
trans_df.to_csv(out_trans, sep='\t', float_format="%.10g")

print(f"Pseudocount used: {pseudo}")
print(f"Saved:\n  {out_cis}\n  {out_trans}")
