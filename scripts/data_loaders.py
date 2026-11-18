"""Shared loaders for the four empirical data sets.

Every analysis, topology and figure script reads the raw files through this
module, so that a data-processing decision is made (and documented) exactly
once.

C. elegans (``data/celegans-NeuronConnect.csv``, WormAtlas NeuronConnect.xls)
-----------------------------------------------------------------------------
The file lists every contact from the point of view of *both* partners:

* ``S``/``Sp`` rows: Neuron 1 sends a chemical synapse to Neuron 2
  (``Sp`` = polyadic). ``R``/``Rp`` rows are the mirror images of those same
  synapses seen from the receiving neuron, so they are redundant and are
  dropped. A (pre, post) pair can appear once as ``S`` and once as ``Sp``;
  the two contact counts are added.
* ``EJ`` rows (gap junctions) are undirected and listed in both directions;
  each unordered pair is kept once. If the two directions disagree on ``Nbr``
  the maximum is used and the number of disagreements is reported.
* ``NMJ`` rows all point to a single pseudo-node called ``NMJ`` (neuromuscular
  junctions); it is not a neuron and those rows are dropped.
* Self-loops are dropped.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"

CELEGANS_PATH = DATA_DIR / "celegans-NeuronConnect.csv"
EMAIL_PATH = DATA_DIR / "email-Eu-core.txt"
EMAIL_LABELS_PATH = DATA_DIR / "email-Eu-core-department-labels.txt"
AIRPORTS_PATH = DATA_DIR / "openflights-airports.dat"
ROUTES_PATH = DATA_DIR / "openflights-routes.dat"
TURBULENCE_PATH = DATA_DIR / "turbulence-mendeley-cylinder50mm-16ms.txt"


# --------------------------------------------------------------------------
# C. elegans
# --------------------------------------------------------------------------


@dataclass
class CelegansData:
    neurons: list            # sorted names, NMJ pseudo-node excluded
    chemical: dict           # (pre, post) -> number of synaptic contacts
    gap: dict                # (u, v) with u < v -> number of gap junctions
    out_deg: np.ndarray      # weighted chemical out-degree, aligned with neurons
    in_deg: np.ndarray       # weighted chemical in-degree
    gap_deg: np.ndarray      # weighted gap-junction degree
    total_deg: np.ndarray    # out + in + gap
    info: dict = field(default_factory=dict)


def load_celegans(verbose: bool = False) -> CelegansData:
    chemical, gap_dir = {}, {}
    counts = {"rows": 0, "S_Sp_rows": 0, "R_Rp_rows_dropped": 0, "EJ_rows": 0,
              "NMJ_rows_dropped": 0, "self_loops_dropped": 0}
    names = set()
    with CELEGANS_PATH.open() as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            if len(row) < 4:
                continue
            n1, n2, typ = row[0].strip(), row[1].strip(), row[2].strip()
            try:
                nbr = int(row[3])
            except ValueError:
                continue
            counts["rows"] += 1
            if typ == "NMJ" or n1 == "NMJ" or n2 == "NMJ":
                counts["NMJ_rows_dropped"] += 1
                continue
            names.update((n1, n2))
            if typ in ("R", "Rp"):
                counts["R_Rp_rows_dropped"] += 1
                continue
            if n1 == n2:
                counts["self_loops_dropped"] += 1
                continue
            if typ in ("S", "Sp"):
                counts["S_Sp_rows"] += 1
                chemical[(n1, n2)] = chemical.get((n1, n2), 0) + nbr
            elif typ == "EJ":
                counts["EJ_rows"] += 1
                gap_dir[(n1, n2)] = gap_dir.get((n1, n2), 0) + nbr

    gap, disagreements, one_sided = {}, 0, 0
    for (u, v), w in gap_dir.items():
        key = (u, v) if u < v else (v, u)
        if key in gap:
            continue
        w_rev = gap_dir.get((v, u))
        if w_rev is None:
            one_sided += 1
            w_rev = w
        elif w_rev != w:
            disagreements += 1
        gap[key] = max(w, w_rev)

    neurons = sorted(names)
    idx = {name: i for i, name in enumerate(neurons)}
    n = len(neurons)
    out_deg = np.zeros(n, dtype=int)
    in_deg = np.zeros(n, dtype=int)
    gap_deg = np.zeros(n, dtype=int)
    for (u, v), w in chemical.items():
        out_deg[idx[u]] += w
        in_deg[idx[v]] += w
    for (u, v), w in gap.items():
        gap_deg[idx[u]] += w
        gap_deg[idx[v]] += w

    info = dict(counts)
    info.update({
        "n_neurons": n,
        "n_chemical_directed_pairs": len(chemical),
        "n_chemical_contacts": int(sum(chemical.values())),
        "n_gap_junction_pairs": len(gap),
        "n_gap_junction_contacts": int(sum(gap.values())),
        "gap_junction_direction_disagreements": disagreements,
        "gap_junction_one_sided_pairs": one_sided,
    })
    if verbose:
        print("C. elegans loader:")
        for k, v in info.items():
            print(f"  {k}: {v}")
        print(f"  gap-junction pairs whose two directions disagree on Nbr: "
              f"{disagreements} (max used)")
    return CelegansData(neurons, chemical, gap, out_deg, in_deg, gap_deg,
                        out_deg + in_deg + gap_deg, info)


def celegans_edges():
    """Connections among neurons (chemical, as directed pairs, followed by gap
    junctions as one pair each) -- the union used for the topology table."""
    data = load_celegans()
    return list(data.chemical) + list(data.gap)


# --------------------------------------------------------------------------
# email-Eu-core
# --------------------------------------------------------------------------


def email_edges():
    edges = np.loadtxt(EMAIL_PATH, dtype=int)
    return [(int(u), int(v)) for u, v in edges]


def email_departments():
    dept = {}
    with EMAIL_LABELS_PATH.open() as fh:
        for line in fh:
            node, d = line.split()
            dept[int(node)] = int(d)
    return dept


def email_degrees():
    edges = np.loadtxt(EMAIL_PATH, dtype=int)
    src, dst = edges[:, 0], edges[:, 1]
    n_nodes = int(max(src.max(), dst.max())) + 1
    out_deg = np.bincount(src, minlength=n_nodes)
    in_deg = np.bincount(dst, minlength=n_nodes)
    return {"out-degree": out_deg, "in-degree": in_deg, "total-degree": out_deg + in_deg}


# --------------------------------------------------------------------------
# OpenFlights
# --------------------------------------------------------------------------


def openflights_edges():
    """One (source id, destination id) pair per row of the route table (a pair
    served by several airlines appears several times)."""
    edges = []
    with ROUTES_PATH.open(encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if len(row) < 6:
                continue
            src_id, dst_id = row[3], row[5]
            if src_id == "\\N" or dst_id == "\\N":
                continue
            edges.append((src_id, dst_id))
    return edges


def openflights_airports():
    """``{airport id: (country, lon, lat)}``."""
    out = {}
    with AIRPORTS_PATH.open(encoding="utf-8") as fh:
        for row in csv.reader(fh):
            out[row[0]] = (row[3], float(row[7]), float(row[6]))
    return out


def openflights_degrees():
    out_deg, in_deg = {}, {}
    for src_id, dst_id in openflights_edges():
        out_deg[src_id] = out_deg.get(src_id, 0) + 1
        in_deg[dst_id] = in_deg.get(dst_id, 0) + 1
    ids = sorted(set(out_deg) | set(in_deg))
    out_arr = np.array([out_deg.get(i, 0) for i in ids])
    in_arr = np.array([in_deg.get(i, 0) for i in ids])
    return ids, {"out-degree": out_arr, "in-degree": in_arr,
                 "total-degree": out_arr + in_arr}


# --------------------------------------------------------------------------
# Turbulence
# --------------------------------------------------------------------------


def load_velocity():
    data = np.loadtxt(TURBULENCE_PATH)
    t, v = data[:, 0], data[:, 1]
    return t, v, float(t[1] - t[0])
