#!/usr/bin/env python3
"""Figures 1-4 of the manuscript (theory, identifiability, network model).

    fig1_cascade.pdf          simulated cascade + B_N(E)
    fig2_entropy.pdf          parametric entropy curves, both classes
    fig3_identifiability.pdf  rho(a) = kappa3/kappa2 + size/power of the depth test
    fig4_network_model.pdf    the hierarchical mixed-Poisson network model

Figure 3(b) reads ``results/compute_article_numbers.json``; nothing else here
depends on the analysis scripts. Style: ``scripts/figstyle.py``.

Run from the repo root::

    PYTHONPATH=. python3 scripts/make_theory_figures.py [fig1 fig2 fig3 fig4]
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
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from scipy import special  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from htheory import network  # noqa: E402
from htheory.background import sample_conditional  # noqa: E402
from htheory.entropy import parametric_curve  # noqa: E402
from htheory.statfactor import B  # noqa: E402
from scripts import common  # noqa: E402

FULL = figstyle.FULL_WIDTH
OI = figstyle.OKABE_ITO


# --------------------------------------------------------------------------
# Figure 1
# --------------------------------------------------------------------------


def fig1_cascade():
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.3), width_ratios=[1.0, 1.0])

    # (a) Ancestral Monte Carlo simulation of the chain T_0 -> ... -> T_3
    # (s = 1/2, a_i = 2, T_0 = 1 fixed): the marginal of log T_i widens with i.
    ax = axes[0]
    rng = np.random.default_rng(0)
    n_paths, a_level, depths = 20_000, 2.0, [0, 1, 2, 3]
    T = np.full(n_paths, 1.0)
    logT = {}
    for level in depths[1:]:
        T = sample_conditional(T, a_level, 0.5, rng)
        logT[level] = np.log(T)
    lo = min(np.quantile(x, 0.001) for x in logT.values())
    hi = max(np.quantile(x, 0.999) for x in logT.values())
    grid = np.linspace(lo, hi, 500)
    colors = figstyle.depth_colors(len(depths))
    step = 1.0
    offsets = {d: (len(depths) - 1 - d) * step for d in depths}
    peak = max(gaussian_kde(logT[d])(grid).max() for d in depths[1:])
    for d in depths:
        base = offsets[d]
        ax.axhline(base, color="0.8", lw=0.5, zorder=1)
        if d == 0:
            # T_0 is a constant: a point mass at log T_0 = 0, drawn as a line.
            ax.plot([0, 0], [base, base + 0.85 * step], color=colors[d], lw=1.6,
                    solid_capstyle="butt", zorder=3)
        else:
            dens = gaussian_kde(logT[d])(grid) / peak * 0.85 * step
            ax.fill_between(grid, base, base + dens, color=colors[d], alpha=0.35,
                            lw=0, zorder=2)
            ax.plot(grid, base + dens, color=colors[d], lw=1.2, zorder=3)
        ax.annotate(f"$T_{d}$", xy=(0.0, base), xycoords=("axes fraction", "data"),
                    xytext=(2, 3), textcoords="offset points", ha="left", va="bottom")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.02, offsets[0] + 0.95 * step)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(r"$\ln T_i$")
    ax.set_ylabel("density (offset)")
    figstyle.panel_label(ax, "(a)")

    # (b) B_N(E), a = 1, s = 1/2: close in the bulk, separated in the tail.
    ax = axes[1]
    E = np.linspace(0.0, 10.0, 240)
    marks = np.linspace(20, len(E) - 40, 5).astype(int)
    for N in range(4):
        st = figstyle.depth_style(N, 4)
        ax.plot(E, B(E, [1.0] * N, 0.5), color=st["color"], ls=st["linestyle"],
                marker=st["marker"], markevery=[int(m) + 6 * N for m in marks], ms=3.5,
                label=f"$N={N}$")
    ax.set_yscale("log")
    ax.set_xlim(0, 10)
    ax.set_xlabel(r"$\beta_0 E$")
    ax.set_ylabel(r"$B_N(E)$")
    ax.legend(loc="lower left", handlelength=3.0)
    figstyle.panel_label(ax, "(b)")
    figstyle.save(fig, "fig1_cascade")


# --------------------------------------------------------------------------
# Figure 2
# --------------------------------------------------------------------------


def entropy_curves(s, a=1.0, n_max=5, n_points=320):
    """``{N: (B, s_tilde)}`` -- the logic of the original Figs. 1-2 script."""
    E_max = 400.0 if s == 1.0 else 90.0   # the power-law tail needs more room in E
    return {N: parametric_curve([a] * N, s, n=n_points, E_max=E_max)
            for N in range(n_max + 1)}


def fig2_entropy(n_points=320):
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.55), sharey=False)
    n_depth = 6
    for ax, s, lab, cls in zip(axes, (0.5, 1.0), ("(a)", "(b)"),
                               (r"$s=1/2$", r"$s=1$")):
        curves = entropy_curves(s, n_points=n_points)
        for N, (Bv, sv) in curves.items():
            st = figstyle.depth_style(N, n_depth)
            order = np.argsort(Bv)
            Bv, sv = np.asarray(Bv)[order], np.asarray(sv)[order]
            # one marker per curve, at its maximum, so depth is not colour-only
            ax.plot(Bv, sv, color=st["color"], ls=st["linestyle"])
            j = int(np.argmax(sv))
            ax.plot(Bv[j], sv[j], color=st["color"], marker=st["marker"], ls="none",
                    ms=4.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$B$")
        ax.set_ylabel(r"$\tilde{s}(B)$")
        figstyle.panel_label(ax, lab)
        ax.annotate(cls, xy=(1.0, 1.0), xycoords="axes fraction", xytext=(0, 3),
                    textcoords="offset points", ha="right", va="bottom")
    handles = [Line2D([], [], label=f"$N={N}$", ms=4.5,
                      **figstyle.depth_style(N, n_depth)) for N in range(n_depth)]
    fig.legend(handles=handles, loc="outside lower center", ncols=6,
               columnspacing=1.0, handlelength=2.6, handletextpad=0.5)
    figstyle.save(fig, "fig2_entropy")


# --------------------------------------------------------------------------
# Figure 3
# --------------------------------------------------------------------------


def fig3_identifiability():
    numbers = common.read_json("compute_article_numbers")
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.4))

    # (a) rho(a) = kappa_3/kappa_2 of one level; monotone in a for each class.
    ax = axes[0]
    a = np.geomspace(0.05, 50, 300)
    rho_half = special.polygamma(2, a) / special.polygamma(1, a)             # < 0
    rho_one = -special.polygamma(2, a + 1) / special.polygamma(1, a + 1)     # > 0
    ax.plot(a, -rho_half, color=OI["blue"], ls="-",
            label=r"$s=1/2$:  $-\rho(a)$")
    ax.plot(a, rho_one, color=OI["vermillion"], ls="--",
            label=r"$s=1$:  $+\rho(a)$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.05, 50)
    ax.set_xlabel(r"$a$")
    ax.set_ylabel(r"$|\rho(a)|=|\kappa_3/\kappa_2|$")
    ax.legend(loc="upper right")
    figstyle.panel_label(ax, "(a)")

    # (b) rejection rate of H0: N = 1 at level 0.05, by true depth.
    ax = axes[1]
    pbd = numbers["power_by_depth"]
    sizes = pbd["sample_sizes"]
    depth_order = sorted(pbd["by_N_true"], key=int)
    cols = [OI["black"], OI["blue"], OI["vermillion"], OI["green"]]
    handles = []
    for i, key in enumerate(depth_order):
        rate = np.array([pbd["by_N_true"][key][str(n)]["rate"] for n in sizes])
        se = np.array([pbd["by_N_true"][key][str(n)]["se"] for n in sizes])
        ax.errorbar(sizes, rate, yerr=1.96 * se, color=cols[i],
                    ls=figstyle.LINESTYLES[i], marker=figstyle.MARKERS[i],
                    capsize=2, elinewidth=0.8, capthick=0.8)
        handles.append(Line2D([], [], color=cols[i], ls=figstyle.LINESTYLES[i],
                              marker=figstyle.MARKERS[i],
                              label=rf"$N_{{\mathrm{{true}}}}={key}$"))
    ax.axhline(0.05, color=figstyle.GREY, ls=(0, (4, 2)), lw=0.8)
    handles.append(Line2D([], [], color=figstyle.GREY, ls=(0, (4, 2)), lw=0.8,
                          label="nominal level"))
    ax.set_xscale("log")
    ax.set_xlim(70, 7000)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel(r"sample size $n$")
    ax.set_ylabel(r"rejection rate of $H_0\!:N=1$")
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.45),
              handlelength=2.8)
    figstyle.panel_label(ax, "(b)")
    figstyle.save(fig, "fig3_identifiability")


# --------------------------------------------------------------------------
# Figure 4
# --------------------------------------------------------------------------

LEAF_COLORS = [OI["blue"], OI["sky"], OI["vermillion"], OI["orange"]]   # level-2
LEAF_MARKERS = ["o", "s", "^", "D"]
L1_COLORS = [OI["blue"], OI["vermillion"]]                              # level-1


def _box(ax, x, y, text, color, w=0.2, h=0.13, face=None, fontsize=9):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.0,rounding_size=0.03",
                                fc="white" if face is None else face, ec=color,
                                lw=1.0, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, zorder=4)


def _arrow(ax, xy0, xy1, color="0.35"):
    ax.annotate("", xy=xy1, xytext=xy0, zorder=2,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0.7,
                                shrinkA=0, shrinkB=0, mutation_scale=6))


def _panel_tree(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    y_root, y1, y2, y_leaf = 0.93, 0.70, 0.47, 0.27
    h = 0.13
    x2 = [0.275, 0.485, 0.695, 0.905]
    x1 = [0.5 * (x2[0] + x2[1]), 0.5 * (x2[2] + x2[3])]
    x0 = 0.5 * (x1[0] + x1[1])
    _box(ax, x0, y_root, r"$\lambda$", "0.2", w=0.14, h=0.11)
    for c, x in enumerate(x1):
        _box(ax, x, y1, rf"$X_{{1,{c + 1}}}$", L1_COLORS[c], w=0.18, h=h)
        _arrow(ax, (x0, y_root - 0.055), (x, y1 + h / 2))
    for c, x in enumerate(x2):
        _box(ax, x, y2, rf"$X_{{2,{c + 1}}}$", LEAF_COLORS[c], w=0.18, h=h)
        _arrow(ax, (x1[c // 2], y1 - h / 2), (x, y2 + h / 2))
        for dx in (-0.06, 0.0, 0.06):
            _arrow(ax, (x, y2 - h / 2), (x + dx, y_leaf + 0.04))
            ax.plot(x + dx, y_leaf, marker=LEAF_MARKERS[c], color=LEAF_COLORS[c],
                    ms=4.5, mec="0.2", mew=0.4, zorder=4)
    for y, lab in ((y_root, "rate"), (y1, "level 1"), (y2, "level 2"),
                   (y_leaf, "vertices")):
        ax.text(0.0, y, lab, ha="left", va="center", fontsize=8, color="0.3")
    ax.text(0.5, 0.115, r"$W_v=X_{1,c_1(v)}\,X_{2,c_2(v)}$", ha="center", va="center")
    ax.text(0.5, -0.02, r"$D_v\mid W_v\sim\mathrm{Poisson}(\lambda W_v)$",
            ha="center", va="center")


def _panel_graph(ax, seed=12, per_leaf=14, lam=6.0, a=2.0):
    """One realization of the N = 2 model: propensities W_v = X_1 X_2 shared
    within a community, edges Chung-Lu with p_uv proportional to W_u W_v,
    denser within a level-2 community than within a level-1 community than
    across level-1 communities. The seed was picked for legibility only (no
    community so sparse or so dense that it hides the others)."""
    rng = np.random.default_rng(seed)
    X1 = sample_conditional(np.ones(2), a, 0.5, rng)
    X2 = sample_conditional(np.ones(4), a, 0.5, rng)
    leaf = np.repeat(np.arange(4), per_leaf)
    W = X1[leaf // 2] * X2[leaf]
    n = leaf.size
    same_leaf = leaf[:, None] == leaf[None, :]
    same_l1 = (leaf[:, None] // 2) == (leaf[None, :] // 2)
    theta = np.where(same_leaf, 1.0, np.where(same_l1, 0.22, 0.04))
    P = theta * np.outer(W, W)
    np.fill_diagonal(P, 0.0)
    P *= lam * W.mean() * n / P.sum()          # mean degree = lam * mean(W)
    P = np.clip(P, 0, 1)
    U = rng.random((n, n))
    A = np.triu(U < P, k=1)
    A = A | A.T
    deg = A.sum(axis=1)

    centers = np.array([[-1.25, 0.62], [-1.25, -0.62], [1.25, 0.62], [1.25, -0.62]])
    r = 0.50 * np.sqrt(rng.random(n))
    phi = 2 * np.pi * rng.random(n)
    xy = centers[leaf] + np.c_[1.25 * r * np.cos(phi), r * np.sin(phi)]

    iu, ju = np.nonzero(np.triu(A, k=1))
    for i, j in zip(iu, ju):
        within = leaf[i] == leaf[j]
        ax.plot(xy[[i, j], 0], xy[[i, j], 1], color="0.35" if within else "0.6",
                lw=0.35 if within else 0.3, alpha=0.55 if within else 0.45, zorder=1)
    for c in range(4):
        m = leaf == c
        ax.scatter(xy[m, 0], xy[m, 1], s=4 + 1.6 * deg[m], marker=LEAF_MARKERS[c],
                   color=LEAF_COLORS[c], edgecolors="0.15", linewidths=0.35, zorder=3)
        cx, cy = centers[c]
        ax.text(cx, cy + (0.62 if cy > 0 else -0.62),
                rf"$X_{{1,{c // 2 + 1}}}X_{{2,{c + 1}}}={W[m][0]:.2f}$",
                ha="center", va="bottom" if cy > 0 else "top", fontsize=8)
    ax.set_xlim(-2.05, 2.05)
    ax.set_ylim(-1.55, 1.55)
    ax.axis("off")


def _panel_degree_distributions(ax, a=2.0, lam=20.0, n=100_000, seed=20260826):
    rng = np.random.default_rng(seed)
    for N in range(4):
        D = network.sample_degrees([a] * N, 0.5, lam, rng, size=n)
        vmr = D.var() / D.mean()
        k = np.sort(D)
        ks, first = np.unique(k, return_index=True)
        ccdf = 1.0 - first / n                      # P(D >= k)
        keep = ks >= 1
        ks, ccdf = ks[keep], ccdf[keep]
        st = figstyle.depth_style(N, 4)
        marks = np.unique(np.searchsorted(ks, np.geomspace(2, ks.max(), 9)))
        marks = marks[marks < ks.size]
        ax.plot(ks, ccdf, color=st["color"], ls=st["linestyle"], marker=st["marker"],
                markevery=[int(m) for m in marks], ms=3.5,
                label=rf"$N={N}$  (VMR$\,=\,${vmr:.1f})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 2000)
    ax.set_ylim(0.5 / n, 2)
    ax.set_xlabel(r"degree $k$")
    ax.set_ylabel(r"$P(D\geq k)$")
    ax.legend(loc="lower left", handlelength=3.0)


def fig4_network_model():
    fig = plt.figure(figsize=(FULL, 3.75))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1.0], width_ratios=[1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    _panel_tree(ax_a)
    _panel_graph(ax_b)
    _panel_degree_distributions(ax_c)
    figstyle.panel_label(ax_a, "(a)")
    figstyle.panel_label(ax_b, "(b)")
    figstyle.panel_label(ax_c, "(c)")
    figstyle.save(fig, "fig4_network_model")


FIGS = {"fig1": fig1_cascade, "fig2": fig2_entropy, "fig3": fig3_identifiability,
        "fig4": fig4_network_model}


def main(argv=None):
    wanted = (argv if argv is not None else sys.argv[1:]) or list(FIGS)
    for key in wanted:
        FIGS[key]()


if __name__ == "__main__":
    main()
