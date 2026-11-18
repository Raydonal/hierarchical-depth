"""Descriptive topology of the three real networks used in the article
(email-Eu-core, OpenFlights, C. elegans): node/edge counts, density, mean
degree, reciprocity, global clustering coefficient (transitivity), and
degree assortativity. These are standard network-science descriptive
statistics, independent of the cascade/mixed-Poisson estimator -- reported so
a reader can place each network in the usual structural-statistics landscape
before reading what the depth estimator adds on top.

Every network is reduced to a simple undirected graph (directions folded,
multi-edges merged, self-loops dropped). For C. elegans that graph is the
union of chemical synapses and gap junctions among neurons, as produced by
``scripts/data_loaders.load_celegans`` (R/Rp mirror rows and the NMJ
pseudo-node removed).

Writes ``results/network_topology.json`` (read by ``make_tables.py``) and
prints the same numbers.

Run from the repo root:

    PYTHONPATH=. python3 scripts/network_topology.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts import common  # noqa: E402
from scripts.data_loaders import (  # noqa: E402
    celegans_edges,
    openflights_edges,
)
from scripts.data_loaders import email_edges as email_eu_core_edges  # noqa: E402


def _undirected_adjacency(directed_edges):
    """``directed_edges``: iterable of (u, v). Returns (nodes, adj) with adj a
    dict[node] -> set(neighbors), self-loops dropped, both directions folded
    into one undirected simple graph."""
    adj = defaultdict(set)
    nodes = set()
    for u, v in directed_edges:
        if u == v:
            continue
        adj[u].add(v)
        adj[v].add(u)
        nodes.add(u)
        nodes.add(v)
    return nodes, adj


def _reciprocity(directed_edges):
    edge_set = set(directed_edges)
    if not edge_set:
        return float("nan")
    reciprocated = sum(1 for (u, v) in edge_set if (v, u) in edge_set)
    return reciprocated / len(edge_set)


def _clustering_and_assortativity(nodes, adj):
    """Global (transitivity) clustering coefficient and Pearson degree
    assortativity, both computed directly from the adjacency dict -- no
    external graph library. Transitivity = 3 * triangles / connected triples
    (= closed triples / connected triples);
    assortativity is the Pearson correlation of (deg(u), deg(v)) over the
    edge list, each edge counted once in each direction (the standard
    definition, Newman 2002)."""
    degree = {u: len(adj[u]) for u in nodes}
    triangles = 0
    triples = 0
    for u in nodes:
        neighbors = adj[u]
        d = len(neighbors)
        triples += d * (d - 1) // 2
        neighbors = list(neighbors)
        for i in range(len(neighbors)):
            ni = adj[neighbors[i]]
            for j in range(i + 1, len(neighbors)):
                if neighbors[j] in ni:
                    triangles += 1
    # ``triangles`` was accumulated over every vertex, so each triangle has
    # already been counted three times (once per corner): it *is* the number
    # of closed triples. An earlier version multiplied by 3 again and reported
    # transitivities three times too large.
    transitivity = (triangles / triples) if triples else float("nan")

    x, y = [], []
    seen = set()
    for u in nodes:
        for v in adj[u]:
            if (v, u) in seen:
                continue
            seen.add((u, v))
            x.append(degree[u])
            y.append(degree[v])
            x.append(degree[v])
            y.append(degree[u])
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)
    assortativity = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else float("nan")
    return transitivity, assortativity


def _summary(name, directed_edges, weighted=False):
    directed_edges = list(directed_edges)
    nodes, adj = _undirected_adjacency(directed_edges)
    n = len(nodes)
    m = sum(len(v) for v in adj.values()) // 2
    density = 2 * m / (n * (n - 1)) if n > 1 else float("nan")
    degrees = np.array([len(adj[u]) for u in nodes], dtype=float)
    transitivity, assortativity = _clustering_and_assortativity(nodes, adj)
    reciprocity = _reciprocity(directed_edges)
    out = {
        "n_nodes": n,
        "n_edges_undirected": m,
        "n_edges_directed": len(set(directed_edges)),
        "density": density,
        "mean_degree": float(degrees.mean()),
        "degree_std": float(degrees.std()),
        "reciprocity": reciprocity,
        "global_clustering": transitivity,
        "degree_assortativity": assortativity,
    }
    print(f"-- {name} --")
    for k, v in out.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")
    return out


def main():
    networks = {
        "email-Eu-core": _summary("email-Eu-core", email_eu_core_edges()),
        "OpenFlights": _summary("OpenFlights", openflights_edges()),
        "C. elegans": _summary("C. elegans", celegans_edges()),
    }
    common.write_json("network_topology", {
        "metadata": common.metadata(
            "scripts/network_topology.py", None,
            note="Deterministic (no random numbers). Simple undirected graphs; "
                 "C. elegans = union of chemical synapses and gap junctions "
                 "among neurons, NMJ pseudo-node and R/Rp mirror rows removed. "
                 "Reciprocity for C. elegans is not meaningful (gap junctions "
                 "enter as single pairs)."),
        "networks": networks,
    })


if __name__ == "__main__":
    main()
