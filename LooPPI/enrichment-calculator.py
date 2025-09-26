#!/usr/bin/env python3
"""
Simple CIS/TRANS Enrichment from /out tables
===========================================

Inputs (produced by your pipeline, no headers unless noted):
  - annotated_loops.bedpe              (one row per loop)
  - annotated_anchors.bed              (one row per anchor)
  - cis_contacts.tsv                   (p1 \t p2 \t count)
  - trans-contacts.tsv                 (p1 \t p2 \t count)
  - loop_end_proportions.tsv           (protein \t count \t prop)
  - anchor_proportions.tsv             (protein \t count \t prop)

Expectations:
  TRANS:  E[p,q] = N_loops * (π_end[p] * π_end[q] * (2 if p≠q else 1))
  CIS:    E[p,q] = N_anchors * (π_anchor[p] * π_anchor[q])        for p≠q
          E[p,p] = N_anchors * (π_anchor[p])                      for p=p

Enrichments use a pseudocount of 1 added to both numerator and denominator:
  enrich = (obs + 1) / (exp + 1)

Outputs (TSV):
  - cis_expected.tsv,   trans_expected.tsv
  - cis_enrichment.tsv, trans_enrichment.tsv

Edit CONFIG and Run.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd

# =========================
# CONFIG
# =========================
CONFIG = {
    # Base directory with the input tables
    "indir": "/mnt/altnas/work/Kyle.Knightly/looppi/out/",  # change if needed

    # Filenames (relative to indir)
    "annot_loops": "annotated_loops.bedpe",
    "annot_anchors": "annotated_anchors.bed",
    "cis_counts": "cis_contacts.tsv",
    "trans_counts": "trans-contacts.tsv",
    "loop_end_props": "loop_end_proportions.tsv",
    "anchor_props": "anchor_proportions.tsv",

    # Outputs (relative to indir)
    "out_cis_exp": "cis_expected.tsv",
    "out_trans_exp": "trans_expected.tsv",
    "out_cis_enr": "cis_enrichment.tsv",
    "out_trans_enr": "trans_enrichment.tsv",

    # Pseudocount to add to obs and exp before dividing
    "pseudocount": 1.0,
}

# =========================
# Helpers
# =========================

def _read_props_tsv(path: Path) -> pd.DataFrame:
    """Read a (protein, count, prop) table with or without header.
    Returns df indexed by protein with columns ['count','prop'] (floats).
    """
    # Try headerless first
    df = pd.read_csv(path, sep="\t", header=None, names=["protein","count","prop"]) \
           .set_index("protein")
    # If the first row looks like a header, re-read with header
    if not pd.to_numeric(df["prop"], errors="coerce").notna().all():
        df = pd.read_csv(path, sep="\t").set_index(df.columns[0])
        # best-effort normalize column names
        cols = {c.lower(): c for c in df.columns}
        count_col = cols.get("count") or list(df.columns)[0]
        prop_col  = cols.get("prop") or list(df.columns)[-1]
        df = df[[count_col, prop_col]].rename(columns={count_col:"count", prop_col:"prop"})
    # ensure numeric
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0.0)
    df["prop"]  = pd.to_numeric(df["prop"],  errors="coerce").fillna(0.0)
    return df


def _read_counts_long(path: Path) -> pd.DataFrame:
    """Read (p1, p2, count) with or without header; return normalized with p1<=p2."""
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=["p1","p2","count"])
        # if header accidentally parsed as row
        if not pd.to_numeric(df["count"], errors="coerce").notna().all():
            raise ValueError
    except Exception:
        df = pd.read_csv(path, sep="\t")
        # guess column names
        cols = {c.lower(): c for c in df.columns}
        p1c = cols.get("p1") or list(df.columns)[0]
        p2c = cols.get("p2") or list(df.columns)[1]
        cc  = cols.get("count") or list(df.columns)[2]
        df = df[[p1c,p2c,cc]].rename(columns={p1c:"p1", p2c:"p2", cc:"count"})
    # normalize order to upper triangle
    df = df.dropna(subset=["p1","p2"]).copy()
    df["p1"] = df["p1"].astype(str)
    df["p2"] = df["p2"].astype(str)
    a = df["p1"].values
    b = df["p2"].values
    swap = a > b
    a2 = np.where(swap, b, a)
    b2 = np.where(swap, a, b)
    df["p1"], df["p2"] = a2, b2
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0.0)
    # aggregate duplicates
    df = df.groupby(["p1","p2"], as_index=False)["count"].sum()
    return df


def _sym_from_upper(df_long: pd.DataFrame, prots: List[str]) -> np.ndarray:
    idx = {p:i for i,p in enumerate(prots)}
    P = len(prots)
    M = np.zeros((P,P), dtype=float)
    if df_long.empty:
        return M
    i = df_long["p1"].map(idx)
    j = df_long["p2"].map(idx)
    mask = i.notna() & j.notna()
    if not mask.any():
        return M
    I = i[mask].astype(int).to_numpy()
    J = j[mask].astype(int).to_numpy()
    V = df_long.loc[mask, "count"].to_numpy(dtype=float)
    M[I,J] = V
    # reflect to full symmetry; keep diagonal as-is
    M = M + M.T
    np.fill_diagonal(M, np.diag(M)/2.0)  # undo doubling on diag
    return M


# =========================
# Main
# =========================

def main():
    indir = Path(CONFIG["indir"]).expanduser()
    # Load counts
    cis_long   = _read_counts_long(indir / CONFIG["cis_counts"])  # p1,p2,count with p1<=p2
    trans_long = _read_counts_long(indir / CONFIG["trans_counts"])  # p1,p2,count with p1<=p2

    # Load proportions
    end_props_df    = _read_props_tsv(indir / CONFIG["loop_end_props"])   # π_end
    anchor_props_df = _read_props_tsv(indir / CONFIG["anchor_props"])     # π_anchor

    # Universe of proteins
    prots = sorted(set(end_props_df.index) | set(anchor_props_df.index) |
                   set(cis_long["p1"]) | set(cis_long["p2"]) |
                   set(trans_long["p1"]) | set(trans_long["p2"]))
    P = len(prots)

    # Observed matrices (symmetric)
    O_cis   = _sym_from_upper(cis_long, prots)
    O_trans = _sym_from_upper(trans_long, prots)

    # Loop & anchor counts
    try:
        N_loops = sum(1 for _ in open(indir / CONFIG["annot_loops"], "r"))
    except FileNotFoundError:
        # Fallback to sum of trans counts divided by average pair multiplicity is messy; require file ideally
        raise FileNotFoundError("annotated_loops.bedpe not found; needed to compute N_loops.")
    try:
        N_anchors = sum(1 for _ in open(indir / CONFIG["annot_anchors"], "r"))
    except FileNotFoundError:
        raise FileNotFoundError("annotated_anchors.bed not found; needed to compute N_anchors.")

    # Align proportions
    pi_end    = end_props_df.reindex(prots)["prop"].fillna(0.0).to_numpy(float)
    pi_anchor = anchor_props_df.reindex(prots)["prop"].fillna(0.0).to_numpy(float)

    # ===============
    # Expectations
    # ===============
    # TRANS: off-diag 2*N*pi_p*pi_q ; diag N*pi_p^2
    E_trans = (2.0 * N_loops) * (pi_end[:,None] * pi_end[None,:])
    np.fill_diagonal(E_trans, N_loops * (pi_end**2))

    # CIS: off-diag N*pi_p*pi_q ; diag N*pi_p
    E_cis = (N_anchors) * (pi_anchor[:,None] * pi_anchor[None,:])
    np.fill_diagonal(E_cis, N_anchors * pi_anchor)

    # ===============
    # Enrichments with pseudocount
    # ===============
    p = float(CONFIG["pseudocount"]) if CONFIG.get("pseudocount") is not None else 1.0
    CIS_enr   = (O_cis   + p) / (E_cis   + p)
    TRANS_enr = (O_trans + p) / (E_trans + p)

    # ===============
    # Save
    # ===============
    index = pd.Index(prots, name="protein")
    def _to_df(M: np.ndarray) -> pd.DataFrame:
        return pd.DataFrame(M, index=index, columns=prots)

    out_cis_exp   = indir / CONFIG["out_cis_exp"]
    out_trans_exp = indir / CONFIG["out_trans_exp"]
    out_cis_enr   = indir / CONFIG["out_cis_enr"]
    out_trans_enr = indir / CONFIG["out_trans_enr"]

    _to_df(E_cis  ).to_csv(out_cis_exp,   sep='\t', float_format='%.10g')
    _to_df(E_trans).to_csv(out_trans_exp, sep='\t', float_format='%.10g')
    _to_df(CIS_enr).to_csv(out_cis_enr,   sep='\t', float_format='%.10g')
    _to_df(TRANS_enr).to_csv(out_trans_enr, sep='\t', float_format='%.10g')

    # Console summary
    print("=== Simple CIS/TRANS Enrichment ===")
    print(f"Proteins: {P}")
    print(f"Loops (N): {N_loops} ; Anchors (|V|): {N_anchors}")
    print(f"Sum obs CIS:   {O_cis.sum():.0f} ; Sum exp CIS:   {E_cis.sum():.2f}")
    print(f"Sum obs TRANS: {O_trans.sum():.0f} ; Sum exp TRANS: {E_trans.sum():.2f}")
    print(f"Wrote:\n  {out_cis_exp}\n  {out_trans_exp}\n  {out_cis_enr}\n  {out_trans_enr}")


if __name__ == "__main__":
    main()
