#!/usr/bin/env python3
"""
Type–Type Contact Enrichment for Multi‑Labeled Networks
======================================================

Given a NetworkX graph whose nodes have a list[str] attribute `labels`, this script:
  • Counts the observed number of edges connecting each pair of label types (a,b).
  • Computes the expected number under a configuration‑style null using *typed stubs*:
        let k_a = total number of stubs carrying type a (sum of degrees of nodes
                that include label a). Let m = number of edges.
      Then
        E[E_ab] = k_a k_b / (2m - 1)                for a ≠ b
        E[E_aa] = k_a (k_a − 1) / (4m − 2)          for a = b
     (Matches the derivation in your screenshot.)
  • Returns/saves observed, expected, and enrichment matrices (obs/exp), symmetric
    with rows/columns = types.
  • (Optional) Verifies the expectation formulas via Monte‑Carlo random stub pairings.
  • (Optional) Plots a clustered heatmap of log10(enrichment), similar to your example.

Notes on multi‑labels:
  - If a node has multiple labels, its degree contributes to *each* label’s stub count.
    Likewise, an edge between two multi‑labeled nodes contributes to all label‑pair
    combinations across its endpoints. This mirrors the typed‑stub counting in the
    expectation so obs/exp are comparable (though
    Σ_{a≤b} E[E_ab] can exceed m when labels overlap).

Edit the CONFIG block and hit Run—no CLI args needed.
"""

from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional
import csv
import numpy as np
import pandas as pd
import networkx as nx

# Optional plotting & clustering deps
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering

# =========================
# CONFIGURATION (edit me!)
# =========================
CONFIG = {
    # Path to a graph written by your generator (NetworkX gpickle with node attr `labels`).
    # Example: "./out/graph.gpickle" produced by make_multilabeled_er_config.py
    "in_graph_path": "./out/graph.gpickle",

    # If True, include unlabeled nodes as a distinct type name (e.g., "UNLABELED").
    # If False, edges involving unlabeled nodes do not contribute to any type pair.
    "include_unlabeled": False,
    "unlabeled_name": "UNLABELED",

    # Monte‑Carlo verification of expectation (random stub pairings)
    "mc_verify": True,
    "mc_reps": 500,        # number of random pairings to estimate E[E_ab]
    "mc_seed": 12345,

    # Output options
    "save_outputs": True,
    "out_dir": "./out_type_type",
    # CSV names
    "obs_csv": "observed_type_pairs.csv",
    "exp_csv": "expected_type_pairs.csv",
    "enr_csv": "enrichment_type_pairs.csv",
    "mc_csv": "mc_expected_type_pairs.csv",  # written only if mc_verify=True
    "diff_csv": "mc_minus_expected.csv",           # MC - analytic expected
    "relerr_csv": "relative_error_mc_vs_expected.csv",  # (MC-Exp)/Exp
    "print_mc_vs_expected_stats": True,

    # Plotting options (log10 heatmap with clustering)
    "plot_heatmap": True,
    "subset": None,        # list of type names to restrict to (or None)
    "cmap": "RdBu_r",
    "figsize": (30, 24),
    "save_png": "enrichment_heatmap.png",
    # Dendrogram/ordering
    "linkage_method": "ward",   # 'ward', 'single', 'complete', 'average', 'centroid', ...
    "metric": "euclidean",
}

# -----------------------
# Helpers
# -----------------------

def get_types(G: nx.Graph, include_unlabeled: bool, unlabeled_name: str) -> List[str]:
    types = set()
    for u in G.nodes():
        labs = G.nodes[u].get("labels", []) or []
        if labs:
            for a in labs:
                types.add(str(a))
        elif include_unlabeled:
            types.add(unlabeled_name)
    return sorted(types)


def node_labels_resolved(G: nx.Graph, u, include_unlabeled: bool, unlabeled_name: str) -> List[str]:
    labs = G.nodes[u].get("labels", []) or []
    if labs:
        return [str(x) for x in labs]
    return [unlabeled_name] if include_unlabeled else []


def compute_stub_counts(
    G: nx.Graph,
    types: Sequence[str],
    include_unlabeled: bool,
    unlabeled_name: str,
) -> Dict[str, int]:
    tset = set(types)
    k: Dict[str, int] = {t: 0 for t in types}
    for u, deg in G.degree():
        labs = node_labels_resolved(G, u, include_unlabeled, unlabeled_name)
        for a in labs:
            if a in tset:
                k[a] += int(deg)
    return k


def count_observed_pairs(
    G: nx.Graph,
    types: Sequence[str],
    include_unlabeled: bool,
    unlabeled_name: str,
) -> Dict[Tuple[str, str], int]:
    tset = set(types)
    obs: Dict[Tuple[str, str], int] = {(a, b): 0 for a in types for b in types if a <= b}

    for u, v in G.edges():
        Au = node_labels_resolved(G, u, include_unlabeled, unlabeled_name)
        Bv = node_labels_resolved(G, v, include_unlabeled, unlabeled_name)
        if not Au or not Bv:
            continue
        for a in Au:
            if a not in tset:
                continue
            for b in Bv:
                if b not in tset:
                    continue
                x, y = (a, b) if a <= b else (b, a)
                obs[(x, y)] += 1
    return obs


def expected_pairs_from_stubs(k: Dict[str, int], m: int) -> Dict[Tuple[str, str], float]:
    exp: Dict[Tuple[str, str], float] = {}
    denom_ab = max(1, 2 * m - 1)  # avoid zero for tiny graphs
    denom_aa = max(1, 4 * m - 2)
    types = list(k.keys())
    for i, a in enumerate(types):
        for j, b in enumerate(types):
            if a > b:
                continue  # fill only a≤b
            ka, kb = k[a], k[b]
            if a == b:
                exp[(a, a)] = (ka * (ka - 1)) / denom_aa
            else:
                exp[(a, b)] = (ka * kb) / denom_ab
    return exp


def dict_to_matrix(d: Dict[Tuple[str, str], float], types: Sequence[str]) -> np.ndarray:
    n = len(types)
    M = np.zeros((n, n), dtype=float)
    idx = {t: i for i, t in enumerate(types)}
    for (a, b), val in d.items():
        i, j = idx[a], idx[b]
        M[i, j] = val
        M[j, i] = val
    return M


def as_dataframe(M: np.ndarray, types: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame(M, index=types, columns=types)


# -----------------------
# Monte‑Carlo verification
# -----------------------

def mc_expected_from_graph(
    G: nx.Graph,
    types: Sequence[str],
    include_unlabeled: bool,
    unlabeled_name: str,
    reps: int = 1000,
    seed: int = 0,
) -> Dict[Tuple[str, str], float]:
    """Monte‑Carlo estimate of E[E_ab] using *true* stubs (length = 2m).

    Unlike the earlier approach that duplicated stubs per label (wrong for
    multi‑labeled nodes), this treats each physical stub once and attaches the
    *set of labels* for its node. For each random perfect matching of stubs,
    we increment all label‑pair combinations across paired stubs.
    """
    rng = np.random.default_rng(seed)
    types = list(types)
    tset = set(types)

    # Build 2m physical stubs, each carrying the node's label list (filtered to known types)
    stub_labels: List[List[str]] = []
    for u, deg in G.degree():
        labs = node_labels_resolved(G, u, include_unlabeled, unlabeled_name)
        labs = [a for a in labs if a in tset]
        for _ in range(int(deg)):
            stub_labels.append(labs)

    L = len(stub_labels)
    m = G.number_of_edges()
    if L != 2 * m:
        raise ValueError(
            f"Internal inconsistency: built {L} stubs but graph has 2m = {2*m}."
        )

    acc: Dict[Tuple[str, str], float] = {(a, b): 0.0 for a in types for b in types if a <= b}

    for _ in range(int(reps)):
        perm = rng.permutation(L)
        for i in range(0, L, 2):
            A = stub_labels[perm[i]]
            B = stub_labels[perm[i + 1]]
            if not A or not B:
                continue
            for a in A:
                for b in B:
                    x, y = (a, b) if a <= b else (b, a)
                    acc[(x, y)] += 1.0

    for key in acc:
        acc[key] /= float(reps)
    return acc

# -----------------------
# Plotting
# -----------------------

def plot_enrichment_heatmap(
    df_enr: pd.DataFrame,
    subset: Optional[List[str]] = None,
    cmap: str = "RdBu_r",
    figsize: Tuple[int, int] = (30, 24),
    linkage_method: str = "ward",
    metric: str = "euclidean",
    save_path: Optional[Path] = None,
) -> None:
    # Restrict to subset if provided
    if subset is not None:
        df_enr = df_enr.loc[
            df_enr.index.intersection(subset), df_enr.columns.intersection(subset)
        ]

    # Log10 transform around 1. Replace nonpositive or NaN with a small epsilon.
    X = df_enr.copy().astype(float)
    pos = X > 0
    if pos.values.any():
        min_pos = np.nanmin(X[pos].values)
        eps = min_pos / 10.0
    else:
        eps = 1e-6
    X = np.log10(np.where(pos, X, eps))
    X = pd.DataFrame(X, index=df_enr.index, columns=df_enr.columns)

    # Cluster order (rows = cols)
    # linkage operates on observations (rows). Ensure symmetry by applying to X
    # and reindex both axes with the same ordering.
    Z = linkage(X.fillna(0.0), method=linkage_method, metric=metric)
    Z = optimal_leaf_ordering(Z, X.fillna(0.0))
    order = leaves_list(Z)
    Xo = X.iloc[order, order]

    # Plot with dendrogram on the left and heatmap on the right
    fig, (ax_dendro, ax_heat) = plt.subplots(
        1, 2, figsize=figsize, gridspec_kw={"width_ratios": [1, 12]}
    )
    from scipy.cluster.hierarchy import dendrogram

    dendrogram(Z, labels=Xo.index, orientation="left", ax=ax_dendro)
    ax_dendro.invert_yaxis()
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])

    sns.heatmap(
        Xo,
        cmap=cmap,
        center=0.0,  # log10(1) = 0
        linewidths=0.05,
        square=True,
        ax=ax_heat,
        cbar_kws={"label": "log10(enrichment)"},
    )

    plt.subplots_adjust(wspace=0.05)
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    else:
        plt.show()


# -----------------------
# Main routine
# -----------------------

def main() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    in_path = Path(CONFIG["in_graph_path"])
    if not in_path.exists():
        raise FileNotFoundError(
            f"Graph not found: {in_path}. Set CONFIG['in_graph_path'] correctly."
        )
    G: nx.Graph = nx.read_gpickle(in_path)

    # Basic counts
    m = G.number_of_edges()
    if m == 0:
        raise ValueError("Graph has no edges—cannot compute contact enrichments.")

    types = get_types(G, CONFIG["include_unlabeled"], CONFIG["unlabeled_name"])
    if not types:
        raise ValueError(
            "No types found. Ensure nodes have a 'labels' list or enable include_unlabeled."
        )

    # Stub counts and observed
    k = compute_stub_counts(G, types, CONFIG["include_unlabeled"], CONFIG["unlabeled_name"])
    obs_dict = count_observed_pairs(
        G, types, CONFIG["include_unlabeled"], CONFIG["unlabeled_name"]
    )
    exp_dict = expected_pairs_from_stubs(k, m)

    # Optional Monte‑Carlo validation of expectation
    if CONFIG.get("mc_verify", False):
        mc = mc_expected_from_graph(
            G,
            types,
            CONFIG["include_unlabeled"],
            CONFIG["unlabeled_name"],
            reps=int(CONFIG.get("mc_reps", 500)),
            seed=int(CONFIG.get("mc_seed", 0))
        )
    else:
        mc = None

    # Matrices
    Obs = dict_to_matrix(obs_dict, types)
    Exp = dict_to_matrix(exp_dict, types)

    # Enrichment (obs/exp); guard zero expectations
    with np.errstate(divide="ignore", invalid="ignore"):
        Enr = Obs / Exp
        Enr[~np.isfinite(Enr)] = np.nan  # inf/NaN -> NaN

    df_obs = as_dataframe(Obs, types)
    df_exp = as_dataframe(Exp, types)
    df_enr = as_dataframe(Enr, types)

    # Optional save
    if CONFIG["save_outputs"]:
        outdir = Path(CONFIG["out_dir"])
        outdir.mkdir(parents=True, exist_ok=True)
        df_obs.to_csv(outdir / CONFIG["obs_csv"])
        df_exp.to_csv(outdir / CONFIG["exp_csv"])
        df_enr.to_csv(outdir / CONFIG["enr_csv"])
        if mc is not None:
            # write MC expected matrix
            mcM = dict_to_matrix(mc, types)
            df_mc = pd.DataFrame(mcM, index=types, columns=types)
            df_mc.to_csv(outdir / CONFIG["mc_csv"])

            # Differences: MC - Analytic
            diffM = mcM - Exp
            df_diff = pd.DataFrame(diffM, index=types, columns=types)
            df_diff.to_csv(outdir / CONFIG["diff_csv"])  # MC - Expected

            # Relative error: (MC - Exp)/Exp  (NaN where Exp==0)
            with np.errstate(divide="ignore", invalid="ignore"):
                relM = diffM / Exp
                relM[~np.isfinite(relM)] = np.nan
            df_rel = pd.DataFrame(relM, index=types, columns=types)
            df_rel.to_csv(outdir / CONFIG["relerr_csv"])
# Quick on-screen summary
    total_obs_pairs = df_obs.values[np.triu_indices(len(types))].sum()
    total_exp_pairs = df_exp.values[np.triu_indices(len(types))].sum()
    print("=== Type–Type Contact Enrichment ===")
    print(f"Graph: {in_path}")
    print(f"Nodes: {G.number_of_nodes()}, Edges: {m}")
    print(f"#Types: {len(types)} (include_unlabeled={CONFIG['include_unlabeled']})")
    print(f"Sum of observed counts over a≤b: {total_obs_pairs:.1f}")
    print(
        f"Sum of expected counts over a≤b: {total_exp_pairs:.3f}  (can exceed m with overlapping labels)"
    )

    # Light sanity check of the expected formulas (non‑overlapping labels only):
    if all(len(G.nodes[u].get("labels", []) or []) in (0, 1) for u in G.nodes()):
        ratio = total_exp_pairs / m if m > 0 else float("nan")
        print(
            f"(Sanity) Expected total / m ≈ {ratio:.3f}  [should be ~1.0 when labels are disjoint]"
        )
    else:
        print(
            "(Info) Nodes are multi‑labeled; total expected over a≤b can legitimately exceed m."
        )

    # MC vs Expected stats (upper‑triangle, finite entries)
    if mc is not None and CONFIG.get("print_mc_vs_expected_stats", True):
        mcM = dict_to_matrix(mc, types)
        diffM = mcM - Exp
        with np.errstate(divide="ignore", invalid="ignore"):
            relM = diffM / Exp
        n = len(types)
        tri = np.triu_indices(n)
        diff_vals = diffM[tri].astype(float)
        rel_vals = relM[tri].astype(float)
        # Mask to finite entries
        mask_diff = np.isfinite(diff_vals)
        mask_rel = np.isfinite(rel_vals)
        if mask_diff.any():
            mad = np.nanmean(np.abs(diff_vals[mask_diff]))
            mx = np.nanmax(np.abs(diff_vals[mask_diff]))
            print(f"MC−Exp | mean abs diff: {mad:.4g}, max abs diff: {mx:.4g}")
        if mask_rel.any():
            med_rel = np.nanmedian(np.abs(rel_vals[mask_rel]))
            p95_rel = np.nanpercentile(np.abs(rel_vals[mask_rel]), 95)
            print(f"MC−Exp | median abs rel err: {med_rel:.3%}, 95th: {p95_rel:.3%}")

    # Plot heatmap if requested
    if CONFIG.get("plot_heatmap", False):
        save_path = Path(CONFIG["out_dir"]) / CONFIG["save_png"] if CONFIG["save_outputs"] else None
        plot_enrichment_heatmap(
            df_enr,
            subset=CONFIG.get("subset"),
            cmap=CONFIG.get("cmap", "RdBu_r"),
            figsize=CONFIG.get("figsize", (30, 24)),
            linkage_method=CONFIG.get("linkage_method", "ward"),
            metric=CONFIG.get("metric", "euclidean"),
            save_path=save_path,
        )

    return df_obs, df_exp, df_enr


if __name__ == "__main__":
    main()
