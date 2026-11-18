#!/usr/bin/env python3
"""Figures 5 and 6 of the manuscript (real networks, model fits, turbulence).

    fig5a_openflights.pdf     route map drawn from the route table alone
    fig5b_email.pdf           adjacency matrix sorted by department
    fig5c_celegans.pdf        adjacency matrix sorted by total degree
    fig6a_email_fit.pdf       degree CCDF vs fitted mixed-Poisson model
    fig6b_openflights_fit.pdf
    fig6c_celegans_fit.pdf
    fig6d_turbulence.pdf      depth-test p-value vs lag, full and thinned

Figure 6 reads ``results/analyze_*.json`` (run the analysis scripts first);
raw data come through ``scripts/data_loaders.py``. Style: ``scripts/figstyle.py``.

Model shown in Fig. 6(a-c): the homogeneous mixed-Poisson model of the
selected class at the *tested* null, i.e. depth ``N_round``, rate
``lambda_hat`` and the shape ``a`` refitted under ``H0: N = N_round`` (the
JSON field ``a_hat_H0``) -- the parameter values the depth test simulates
from. The band is the pointwise 95% envelope of the CCDF over 500 simulated
sequences of the same size as the data.

Run from the repo root::

    PYTHONPATH=. python3 scripts/make_data_figures.py [fig5a fig5b ... fig6d]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts import figstyle  # noqa: E402

figstyle.apply()

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from htheory import network  # noqa: E402
from scripts import common, data_loaders  # noqa: E402
from scripts.analyze_openflights import COUNTRY_TO_CONTINENT  # noqa: E402

FULL, HALF = figstyle.FULL_WIDTH, figstyle.HALF_WIDTH
OI = figstyle.OKABE_ITO
N_ENVELOPE = 500
ENVELOPE_SEED = 20260905

CONTINENTS = {            # colour, marker
    "Africa": (OI["orange"], "o"),
    "Asia": (OI["vermillion"], "s"),
    "Europe": (OI["blue"], "^"),
    "North America": (OI["green"], "D"),
    "Oceania": (OI["purple"], "v"),
    "South America": (OI["sky"], "P"),
}


# --------------------------------------------------------------------------
# Figure 5
# --------------------------------------------------------------------------


def fig5a_openflights():
    airports = data_loaders.openflights_airports()
    edges = data_loaders.openflights_edges()
    fig, ax = plt.subplots(figsize=(FULL, 2.6))
    rng = np.random.default_rng(0)
    sample = [edges[i] for i in rng.choice(len(edges), size=6000, replace=False)]
    segs = [[airports[u][1:], airports[v][1:]] for u, v in sample
            if u in airports and v in airports]
    ax.add_collection(LineCollection(segs, colors="0.55", linewidths=0.2, alpha=0.18,
                                     rasterized=True, zorder=1))
    used = sorted({u for e in edges for u in e if u in airports})
    cont = np.array([COUNTRY_TO_CONTINENT.get(airports[u][0], "Other") for u in used])
    xy = np.array([airports[u][1:] for u in used])
    other = ~np.isin(cont, list(CONTINENTS))
    if other.any():
        ax.scatter(xy[other, 0], xy[other, 1], s=3, color="0.5", marker="x",
                   linewidths=0.4, rasterized=True, zorder=2)
    handles = []
    for name, (color, marker) in CONTINENTS.items():
        m = cont == name
        ax.scatter(xy[m, 0], xy[m, 1], s=3.2, color=color, marker=marker,
                   linewidths=0, rasterized=True, zorder=2)
        handles.append(Line2D([], [], ls="none", marker=marker, color=color, ms=4,
                              label=name))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-58, 84)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.legend(handles=handles, loc="outside lower center", ncols=6, fontsize=8,
               handletextpad=0.0, columnspacing=0.9, borderpad=0.1)
    figstyle.save(fig, "fig5a_openflights")
    print(f"  {len(used)} airports with coordinates, {len(segs)} sampled routes drawn")


def _matrix_axes(ax, n, label):
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect("equal")
    ticks = [0, 500, 1000] if n > 600 else [0, 100, 200]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel(label)
    ax.set_ylabel("vertex")
    ax.spines[["top", "right"]].set_visible(True)


def fig5b_email():
    edges = data_loaders.email_edges()
    dept = data_loaders.email_departments()
    n = max(max(u, v) for u, v in edges) + 1
    order = sorted(range(n), key=lambda v: (dept.get(v, -1), v))
    rank = {v: i for i, v in enumerate(order)}
    img = np.ones((n, n, 3))
    # Within-department blocks shaded light blue (Okabe-Ito sky blue, very
    # light) instead of drawing 41 boundary lines over the data.
    dept_sorted = np.array([dept.get(v, -1) for v in order])
    same = dept_sorted[:, None] == dept_sorted[None, :]
    img[same] = (0.84, 0.92, 0.97)
    for u, v in edges:
        img[rank[u], rank[v]] = 0.0          # row = sender, column = recipient
    fig, ax = plt.subplots(figsize=(HALF, HALF))
    ax.imshow(img, interpolation="nearest", rasterized=True, zorder=1)
    _matrix_axes(ax, n, "vertex (sorted by department)")
    figstyle.save(fig, "fig5b_email")


def fig5c_celegans():
    data = data_loaders.load_celegans()
    keep = np.nonzero(data.total_deg > 0)[0]
    order = keep[np.argsort(-data.total_deg[keep], kind="stable")]
    rank = {data.neurons[i]: r for r, i in enumerate(order)}
    n = len(order)
    img = np.ones((n, n, 3))
    black = np.array([0.0, 0.0, 0.0])
    orange = np.array([int(OI["orange"][i:i + 2], 16) / 255 for i in (1, 3, 5)])
    for (u, v) in data.chemical:
        img[rank[u], rank[v]] = black        # row = presynaptic, column = postsynaptic
    for (u, v) in data.gap:                  # symmetric, drawn on top
        img[rank[u], rank[v]] = orange
        img[rank[v], rank[u]] = orange
    fig, ax = plt.subplots(figsize=(HALF, HALF))
    ax.imshow(img, interpolation="nearest", rasterized=True)
    _matrix_axes(ax, n, "neuron (sorted by total degree)")
    ax.set_ylabel("neuron")
    fig.legend(handles=[Patch(fc="black", ec="none", label="chemical synapse"),
                        Patch(fc=OI["orange"], ec="none", label="gap junction")],
               loc="outside upper center", ncols=2, fontsize=8, handlelength=1.0,
               handleheight=0.8, columnspacing=1.2, borderpad=0.1)
    figstyle.save(fig, "fig5c_celegans")


# --------------------------------------------------------------------------
# Figure 6
# --------------------------------------------------------------------------


def _ccdf(D):
    D = np.sort(np.asarray(D))
    ks, first = np.unique(D, return_index=True)
    keep = ks >= 1
    return ks[keep], (1.0 - first / D.size)[keep]


def _envelope(cls, s, n, k_grid, rng):
    """Pointwise 2.5/50/97.5% of P(D >= k) over simulated sequences of size n."""
    a_vec = [cls["a_hat_H0"]] * cls["N_round"]
    out = np.empty((N_ENVELOPE, k_grid.size))
    for b in range(N_ENVELOPE):
        D = np.sort(network.sample_degrees(a_vec, s, cls["lambda_hat"], rng, size=n))
        out[b] = 1.0 - np.searchsorted(D, k_grid, side="left") / n
    return np.percentile(out, [2.5, 50, 97.5], axis=0)


def _series_by_name(result, name):
    return next(r for r in result["series"] if r["series"] == name)


def _fit_panel(stem, result_name, raw, names, legend_loc="lower left"):
    """``raw``: {series name: degree array as analysed}."""
    result = common.read_json(result_name)
    rng = np.random.default_rng(ENVELOPE_SEED)
    fig, ax = plt.subplots(figsize=(HALF, 2.4))
    colors = [OI["black"], OI["blue"], OI["vermillion"]]
    handles = []
    kmax = max(int(raw[nm].max()) for nm in names)
    floor = 0.5 / max(raw[nm].size for nm in names)
    for i, nm in enumerate(names):
        rec = _series_by_name(result, nm)
        D = raw[nm]
        assert D.size == rec["n"], (nm, D.size, rec["n"])
        color, marker = colors[i], figstyle.MARKERS[i]
        key = rec["selected_class"]
        short = nm.replace("-degree", "")
        if key is None:
            label = f"{short}: no admissible fit"
        else:
            cls = rec["classes"][key]
            s = float(key)
            k_grid = np.unique(np.round(np.geomspace(1, 3 * kmax, 160)).astype(int))
            lo, med, hi = _envelope(cls, s, D.size, k_grid, rng)
            ax.fill_between(k_grid, np.maximum(lo, floor / 10), np.maximum(hi, floor / 10),
                            color=color, alpha=0.18, lw=0, zorder=1)
            ax.plot(k_grid, np.maximum(med, floor / 10), color=color, lw=0.9,
                    ls=figstyle.LINESTYLES[i], zorder=2)
            s_txt = "1/2" if key == "0.5" else "1"
            label = (f"{short}: $s={s_txt}$, $N={cls['N_round']}$" if len(names) > 1
                     else f"model: $s={s_txt}$, $N={cls['N_round']}$")
        ks, cc = _ccdf(D)
        ax.plot(ks, cc, ls="none", marker=marker, ms=2.6, mfc="none" if i else color,
                mec=color, mew=0.6, zorder=3)
        if len(names) > 1:
            handles.append(Line2D([], [], color=color, ls=figstyle.LINESTYLES[i]
                                  if key else "none", marker=marker, ms=3.5,
                                  mfc="none" if i else color, mec=color, lw=0.9,
                                  label=label))
        else:
            handles.append(Line2D([], [], color=color, ls="none", marker=marker,
                                  ms=3.5, label=f"data ({short}-degree)"))
            if key is not None:
                handles.append(Line2D([], [], color=color, lw=0.9, label=label))
                handles.append(Patch(fc=color, alpha=0.18, ec="none",
                                     label="95% envelope"))
            else:
                handles.append(Line2D([], [], ls="none", label="no admissible fit"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.8, 3 * kmax)
    ax.set_ylim(floor / (4 if len(names) > 1 else 1), 1.6)   # room for the legend
    ax.set_xlabel(r"degree $k$")
    ax.set_ylabel(r"$P(D\geq k)$")
    ax.legend(handles=handles, loc=legend_loc, fontsize=8, handlelength=2.0,
              handletextpad=0.5)
    figstyle.save(fig, stem)


def fig6a_email_fit():
    _fit_panel("fig6a_email_fit", "analyze_email_eu_core",
               data_loaders.email_degrees(), ["out-degree"])


def fig6b_openflights_fit():
    _, deg = data_loaders.openflights_degrees()
    _fit_panel("fig6b_openflights_fit", "analyze_openflights",
               {k: v[v > 0] for k, v in deg.items()}, ["out-degree"])


def fig6c_celegans_fit():
    data = data_loaders.load_celegans()
    raw = {"out-degree": data.out_deg, "in-degree": data.in_deg,
           "total-degree": data.total_deg}
    _fit_panel("fig6c_celegans_fit", "analyze_celegans",
               {k: v[v > 0] for k, v in raw.items()},
               ["out-degree", "in-degree", "total-degree"])


def fig6d_turbulence():
    result = common.read_json("analyze_turbulence")
    fig, ax = plt.subplots(figsize=(HALF, 2.4))
    styles = {"full": dict(color=OI["black"], marker="o", mfc=OI["black"], ls="-"),
              "thinned": dict(color=OI["vermillion"], marker="s", mfc="none", ls="--")}
    floor = None
    for variant, st in styles.items():
        rows = [r for r in result["series"] if r["variant"] == variant
                and r["selected_class"] is not None]
        rows.sort(key=lambda r: r["lag_ms"])
        lag = np.array([r["lag_ms"] for r in rows])
        cls = [r["classes"][r["selected_class"]] for r in rows]
        p = np.array([c["p_value"] for c in cls], dtype=float)
        floor = 1.0 / cls[0]["n_boot"]
        censored = p < floor
        p_plot = np.where(censored, floor, p)
        ax.plot(lag, p_plot, color=st["color"], ls=st["ls"], lw=0.9, zorder=2)
        ax.plot(lag[~censored], p_plot[~censored], ls="none", marker=st["marker"],
                color=st["color"], mfc=st["mfc"], mew=0.8, ms=4.5, zorder=3)
        # p = 0 out of n_boot draws: shown at the resolution limit, arrow down
        ax.plot(lag[censored], p_plot[censored], ls="none", marker=st["marker"],
                color=st["color"], mfc=st["mfc"], mew=0.8, ms=4.5, zorder=3)
        for x in lag[censored]:
            ax.annotate("", xy=(x, floor * 0.55), xytext=(x, floor * 0.93),
                        arrowprops=dict(arrowstyle="-|>", color=st["color"], lw=0.7,
                                        mutation_scale=5, shrinkA=0, shrinkB=0),
                        zorder=1)
    ax.axhline(0.05, color=figstyle.GREY, ls=(0, (4, 2)), lw=0.8, zorder=1)
    ax.annotate("0.05", xy=(1.0, 0.05), xycoords=("axes fraction", "data"),
                xytext=(-1, 2), textcoords="offset points", ha="right", va="bottom",
                fontsize=8, color="0.3")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.4, 150)
    ax.set_ylim(floor * 0.45, 1.6)
    ax.set_xlabel(r"lag $\tau$ (ms)")
    ax.set_ylabel(r"$p$-value of $H_0\!:N=N_{\mathrm{rd}}$")
    handles = [Line2D([], [], color=st["color"], ls=st["ls"], lw=0.9, marker=st["marker"],
                      mfc=st["mfc"], mew=0.8, ms=4.5,
                      label="full sample" if v == "full" else "thinned")
               for v, st in styles.items()]
    handles.append(Line2D([], [], ls="none", marker=r"$\downarrow$", color="0.2", ms=6,
                          label=rf"$p<{floor:g}$"))
    fig.legend(handles=handles, loc="outside upper center", ncols=3, fontsize=8,
               handlelength=1.8, handletextpad=0.4, columnspacing=0.8, borderpad=0.1)
    figstyle.save(fig, "fig6d_turbulence")


FIGS = {"fig5a": fig5a_openflights, "fig5b": fig5b_email, "fig5c": fig5c_celegans,
        "fig6a": fig6a_email_fit, "fig6b": fig6b_openflights_fit,
        "fig6c": fig6c_celegans_fit, "fig6d": fig6d_turbulence}


def main(argv=None):
    wanted = (argv if argv is not None else sys.argv[1:]) or list(FIGS)
    for key in wanted:
        FIGS[key]()


if __name__ == "__main__":
    main()
