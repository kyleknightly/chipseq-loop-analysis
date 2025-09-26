"""
Erdös Renyi null model for multi-label networks.
Input N, P, |L|, p(L), and p(empty)
Returns an NP graph with multiple labels depending on function p, and a fraction of empty nodes
"""

#!/usr/bin/env python3
"""
make_multilabeled_er_config.py
--------------------------------
Configuration-driven ("parameter block") script that you can **just run**
(no command-line flags).

It generates a null-model *multi-labeled* network:
  1) Builds an Erdős-Rényi G(n, p) graph (undirected).
  2) Creates L labels with assignment probabilities that decay by label index i,
     using either an exponential or power-law scheme.
  3) Forces a user-set fraction of nodes to have **no labels** at all.
  4) Attaches a list[str] node attribute `labels`.
  5) Optionally saves outputs (NetworkX gpickle + CSVs).

Outputs (if SAVE_OUTPUTS=True):
  - graph.gpickle  : NetworkX graph with node attribute "labels" (list[str])
  - node_labels.csv: node id, semicolon-separated label list
  - edges.csv      : edge list (u, v)
"""

import csv
import math
import random
from pathlib import Path
from typing import List, Sequence, Tuple

import networkx as nx

# =========================
# CONFIGURATION (edit me!)
# =========================
CONFIG = {
    # Graph size and density
    "n": 500,                 # number of nodes
    "p": 0.01,               # edge probability for G(n,p)

    # Labels
    "L": 20,                 # number of labels (L1..L20)
    "label_prefix": "L",    # prefix for label names

    # Decay scheme for label assignment probabilities by index i=1..L
    # Choose one: "exp" or "power"
    "scheme": "exp",
    # For scheme == "exp":    p_i = p0 * exp(-lam * (i-1))
    "p0": 0.6,
    "lam": 0.01,
    # For scheme == "power":  p_i = p0 / (i ** alpha)
    "alpha": 1.0,

    # Fraction of nodes forced to have NO labels at all (kept empty)
    "unlabeled_frac": 0.30,

    # Random seed for reproducibility (used for both labels + ER graph)
    "seed": 42,

    # Output options
    "save_outputs": True,
    "out_dir": "./out",
}


# ===============
# CORE FUNCTIONS
# ===============

def compute_label_probs(
    L: int,
    scheme: str = "exp",
    p0: float = 0.6,
    lam: float = 0.25,
    alpha: float = 1.0,
) -> List[float]:
    """Compute per‑label assignment probabilities for labels i=1..L.

    Two decay schemes:
      - "exp"   : p_i = p0 * exp(-lam * (i-1))
      - "power" : p_i = p0 / (i ** alpha)
    Values are clipped to [0, 1].
    """
    probs: List[float] = []
    for i in range(1, L + 1):
        if scheme == "exp":
            p = CONFIG["p0"] * math.exp(-CONFIG["lam"] * (i - 1))
        elif scheme == "power":
            p = CONFIG["p0"] / (i ** CONFIG["alpha"]) if i > 0 else 0.0
        else:
            raise ValueError("Unknown scheme. Use 'exp' or 'power'.")
        probs.append(max(0.0, min(1.0, p)))
    return probs


def assign_labels_to_nodes(
    n: int,
    label_probs: Sequence[float],
    unlabeled_frac: float,
    rng: random.Random,
    label_prefix: str = "L",
) -> List[List[str]]:
    """Assign *multi‑labels* to n nodes.

    Steps:
      1) Choose a set (fraction `unlabeled_frac`) of nodes to be **forced unlabeled**.
      2) For the remaining nodes, assign each label i independently with
         probability `label_probs[i-1]`.
      3) Nodes not forced‑unlabeled can still end up empty by chance.

    Returns: list of lists of label strings, length n.
    """
    if not (0.0 <= unlabeled_frac <= 1.0):
        raise ValueError("unlabeled_frac must be in [0, 1].")

    L = len(label_probs)
    labels_per_node: List[List[str]] = [[] for _ in range(n)]

    # Forced‑unlabeled set
    k_unlabeled = int(round(unlabeled_frac * n))
    forced_unlabeled = set(rng.sample(range(n), k_unlabeled)) if k_unlabeled > 0 else set()

    for u in range(n):
        if u in forced_unlabeled:
            continue  # stays empty
        node_labels: List[str] = []
        for i in range(1, L + 1):
            if rng.random() < label_probs[i - 1]:
                node_labels.append(f"{label_prefix}{i}")
        labels_per_node[u] = sorted(node_labels)

    return labels_per_node


def build_er_graph(n: int, p: float, seed: int) -> nx.Graph:
    """Build an undirected Erdős–Rényi G(n, p) graph."""
    return nx.gnp_random_graph(n=n, p=p, seed=seed, directed=False)


def attach_labels(G: nx.Graph, labels_per_node: Sequence[Sequence[str]]) -> None:
    """Attach labels to the graph as node attribute 'labels' (list[str])."""
    if G.number_of_nodes() != len(labels_per_node):
        raise ValueError("labels_per_node length must match number of nodes in G.")
    for u in G.nodes():
        G.nodes[u]["labels"] = list(labels_per_node[u])


def save_outputs(
    G: nx.Graph,
    out_dir: str,
    graph_name: str = "graph.gpickle",
    node_labels_name: str = "node_labels.csv",
    edges_name: str = "edges.csv",
) -> Tuple[Path, Path, Path]:
    """Save GPickle + CSVs to `out_dir`.

    Returns: (gpickle_path, node_labels_csv, edges_csv)
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    gp_path = out_path / graph_name
    nx.write_gpickle(G, gp_path)

    nl_path = out_path / node_labels_name
    with nl_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node", "labels"])  # header
        for u in G.nodes():
            labs = G.nodes[u].get("labels", [])
            writer.writerow([u, ";".join(labs)])

    e_path = out_path / edges_name
    with e_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["u", "v"])  # header
        for u, v in G.edges():
            writer.writerow([u, v])

    return gp_path, nl_path, e_path


# ==========
#   MAIN
# ==========

def main() -> nx.Graph:
    rng = random.Random(CONFIG["seed"])

    # 1) Label probabilities
    label_probs = compute_label_probs(
        L=CONFIG["L"],
        scheme=CONFIG["scheme"],
        p0=CONFIG["p0"],
        lam=CONFIG["lam"],
        alpha=CONFIG["alpha"],
    )

    # 2) Assign labels
    labels_per_node = assign_labels_to_nodes(
        n=CONFIG["n"],
        label_probs=label_probs,
        unlabeled_frac=CONFIG["unlabeled_frac"],
        rng=rng,
        label_prefix=CONFIG["label_prefix"],
    )

    # 3) Build ER graph
    G = build_er_graph(n=CONFIG["n"], p=CONFIG["p"], seed=CONFIG["seed"])

    # 4) Attach labels
    attach_labels(G, labels_per_node)

    # 5) Optionally save
    if CONFIG["save_outputs"]:
        gp, nl, ed = save_outputs(G, CONFIG["out_dir"])

    # 6) Summary
    total_labels = sum(len(Ls) for Ls in labels_per_node)
    num_labeled_nodes = sum(1 for Ls in labels_per_node if len(Ls) > 0)
    avg_labels_per_node = total_labels / CONFIG["n"] if CONFIG["n"] > 0 else 0.0
    avg_labels_per_labeled_node = (
        total_labels / num_labeled_nodes if num_labeled_nodes > 0 else 0.0
    )

    # Label counts (top 10)
    counts = {}
    for Ls in labels_per_node:
        for lab in Ls:
            counts[lab] = counts.get(lab, 0) + 1
    top10 = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:10]

    print("=== Multi‑labeled ER Graph Summary ===")
    print(f"Nodes: {CONFIG['n']}, Edges: {G.number_of_edges()}, p: {CONFIG['p']}")
    print(f"Labels: {CONFIG['L']} (scheme={CONFIG['scheme']}, p0={CONFIG['p0']}, lam={CONFIG['lam']}, alpha={CONFIG['alpha']})")
    print(f"Forced unlabeled fraction: {CONFIG['unlabeled_frac']}")
    print(f"Avg labels per node (incl. unlabeled): {avg_labels_per_node:.3f}")
    print(f"Avg labels per *labeled* node: {avg_labels_per_labeled_node:.3f}")
    print("Top‑10 label counts:", top10)
    if CONFIG["save_outputs"]:
        print("Saved files to:")
        print("  -", Path(CONFIG["out_dir"]) / "graph.gpickle")
        print("  -", Path(CONFIG["out_dir"]) / "node_labels.csv")
        print("  -", Path(CONFIG["out_dir"]) / "edges.csv")

    return G


if __name__ == "__main__":
    main()
