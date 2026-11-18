"""SUPERSEDED by ``make_data_figures.py`` (Figs. 5a-c); kept only so the old
``fig_real_networks`` can be rebuilt with ``--legacy``.

Draw the actual connection structure of the three real networks used in
the article -- not a schematic, the real graphs -- one panel each, plus a
descriptive-statistics table (``figures/network_topology.json``, written by
``network_topology.py``, which this script also regenerates if missing).

Panel choice is data-driven, not decorative: OpenFlights carries real
geographic coordinates, so it is drawn on its actual map; email-Eu-core
carries a real (flat) department partition, so its adjacency matrix is
sorted by department to show the block structure directly; C. elegans has
neither, so its adjacency matrix is sorted by degree, which is the
structure the power-law fit in Section 6 is about (a few hub neurons
carrying most of the connections).

Run from the repo root:

    PYTHONPATH=. .venv/bin/python3 scripts/network_topology.py
    PYTHONPATH=. .venv/bin/python3 scripts/make_network_figures.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.analyze_openflights import COUNTRY_TO_CONTINENT  # noqa: E402
from scripts.data_loaders import celegans_edges, openflights_edges  # noqa: E402
from scripts.data_loaders import email_edges as email_eu_core_edges  # noqa: E402

DATA_DIR = REPO_ROOT / "data"
FIGURES = REPO_ROOT / "figures"

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10.5,
    "legend.fontsize": 7.5,
})
CONTINENT_COLORS = {
    "Africa": "#DD8452", "Asia": "#C44E52", "Europe": "#4C72B0",
    "North America": "#55A868", "Oceania": "#8172B2", "South America": "#937860",
    "Antarctica": "0.5", "Unknown": "0.7",
}


def panel_openflights(ax):
    airport_xy, airport_country = {}, {}
    with (DATA_DIR / "openflights-airports.dat").open(encoding="utf-8") as fh:
        for row in csv.reader(fh):
            airport_id, country = row[0], row[3]
            lat, lon = float(row[6]), float(row[7])
            airport_xy[airport_id] = (lon, lat)
            airport_country[airport_id] = country

    edges = openflights_edges()
    rng = np.random.default_rng(0)
    sample = [edges[i] for i in rng.choice(len(edges), size=6000, replace=False)]
    for src, dst in sample:
        if src not in airport_xy or dst not in airport_xy:
            continue
        (x0, y0), (x1, y1) = airport_xy[src], airport_xy[dst]
        ax.plot([x0, x1], [y0, y1], color="0.6", lw=0.15, alpha=0.25, zorder=1)

    used = {u for e in edges for u in e if u in airport_xy}
    for continent in CONTINENT_COLORS:
        xs, ys = [], []
        for aid in used:
            country = airport_country.get(aid, "?")
            cont = COUNTRY_TO_CONTINENT.get(country, "Unknown")
            if cont == continent:
                x, y = airport_xy[aid]
                xs.append(x)
                ys.append(y)
        if xs:
            ax.scatter(xs, ys, s=2.5, color=CONTINENT_COLORS[continent],
                       label=continent, zorder=2, linewidths=0)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(f"OpenFlights: {len(used)} of 3330 routed airports have coordinates "
                 f"({len(sample)} routes shown)", fontsize=8.5)
    ax.legend(loc="lower left", markerscale=4, frameon=False, ncol=2, handletextpad=0.2,
              columnspacing=0.6, fontsize=6.5)


def panel_email_eu_core(ax):
    edges = email_eu_core_edges()
    dept = {}
    with (DATA_DIR / "email-Eu-core-department-labels.txt").open() as fh:
        for line in fh:
            node, d = line.split()
            dept[int(node)] = int(d)
    n_nodes = max(max(u, v) for u, v in edges) + 1
    order = sorted(range(n_nodes), key=lambda v: (dept.get(v, -1), v))
    rank = {v: i for i, v in enumerate(order)}

    A = np.zeros((n_nodes, n_nodes), dtype=bool)
    for u, v in edges:
        A[rank[u], rank[v]] = True
        A[rank[v], rank[u]] = True
    ax.imshow(A, cmap="Greys", vmax=1, interpolation="nearest", aspect="equal")

    dept_sorted = [dept.get(v, -1) for v in order]
    boundaries = [i for i in range(1, n_nodes) if dept_sorted[i] != dept_sorted[i - 1]]
    for b in boundaries:
        ax.axhline(b - 0.5, color=(0.85, 0.3, 0.3), lw=0.15, alpha=0.5)
        ax.axvline(b - 0.5, color=(0.85, 0.3, 0.3), lw=0.15, alpha=0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"email-Eu-core: {n_nodes} nodes, sorted by department (42 blocks)", fontsize=9.5)


def panel_celegans(ax):
    edges = celegans_edges()
    deg = {}
    for u, v in edges:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    nodes = sorted(deg, key=lambda v: -deg[v])
    rank = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    A = np.zeros((n, n), dtype=bool)
    for u, v in edges:
        A[rank[u], rank[v]] = True
        A[rank[v], rank[u]] = True
    ax.imshow(A, cmap="Greys", vmax=1, interpolation="nearest", aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"C. elegans: {n} neurons, sorted by degree (hubs first)", fontsize=9.5)


def main():
    if "--legacy" not in sys.argv:
        raise SystemExit(
            "This script is superseded by scripts/make_theory_figures.py and "
            "scripts/make_data_figures.py (`make figures`). Its old outputs in "
            "figures/ are kept untouched because the manuscript may still "
            "reference them; pass --legacy to regenerate them anyway.")
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0), constrained_layout=True)
    panel_openflights(axes[0])
    panel_email_eu_core(axes[1])
    panel_celegans(axes[2])
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"fig_real_networks.{ext}", dpi=220)
    plt.close(fig)
    print("wrote fig_real_networks (.pdf + .png) to", FIGURES)


if __name__ == "__main__":
    main()
