#!/usr/bin/env python3
"""Recompute the Section 4 numbers the article quotes (log-cumulant Monte
Carlo check and the depth-test power curve), writing them to
``results/compute_article_numbers.json`` for the figure and table scripts and
for the prose to cite directly rather than from memory. Every rejection rate
is stored with its Monte Carlo standard error, sqrt(p (1 - p) / R), R being
the number of replicates.

Run from the repo root::

    PYTHONPATH=. python3 scripts/compute_article_numbers.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from htheory import identifiability as ident  # noqa: E402
from htheory.background import sample_TN  # noqa: E402

from scripts import common  # noqa: E402

SEED = 20260826
ALPHA = 0.05


def _with_se(powers, n_reps):
    """``{n: rate}`` -> ``{str(n): {"rate", "se", "n_reps"}}``."""
    return {str(n): {"rate": float(p), "se": float(np.sqrt(p * (1 - p) / n_reps)),
                     "n_reps": n_reps} for n, p in powers.items()}


def main():
    rng = np.random.default_rng(SEED)

    # Monte Carlo cross-check of Theorem 3's closed-form log-cumulants,
    # heterogeneous three-level gamma cascade.
    a = [1.5, 3.0, 0.7]
    T = sample_TN(a, 0.5, rng, size=4_000_000)
    kstats = [float(v) for v in ident.empirical_log_cumulants(T, r_max=4)]
    theory = [float(v) for v in ident.log_cumulants_theory(a, 0.5, r_max=4)]

    # kappa_3/kappa_2 as a function of a, the injective map the homogeneous
    # estimator inverts (Theorem 3, Step 4).
    a_grid = np.geomspace(0.05, 50, 200)
    ratio, ratio_s1 = [], []
    for av in a_grid:
        k2, k3 = ident.log_cumulants_theory([av], 0.5, r_max=3)
        ratio.append(float(k3 / k2))
        k2, k3 = ident.log_cumulants_theory([av], 1.0, r_max=3)
        ratio_s1.append(float(k3 / k2))

    # Power curve: true depth N=3, a=2, tested against the shallow null N0=1,
    # gamma branch -- same parameters as the draft's existing table.
    powers = ident.power_curve(
        a_true=[2.0, 2.0, 2.0], s=0.5, N0=1, rng=rng,
        sample_sizes=(100, 300, 1000, 5000), n_reps=200, n_boot=300,
    )

    # Full factorial design: tested null fixed at N0=1 (the depth every real
    # application in Sec. 6 fits under the gamma branch), true depth crossed
    # against sample size. N_true=1 is the null itself, so its "power" is the
    # test's empirical size -- it should track the nominal 0.05 rather than
    # showing a spurious trend; N_true=2,3,5 show power rising both with n
    # and with the depth gap N_true - N0. a=2 at every level throughout, so
    # depth is the only thing varying across rows.
    depth_grid = [1, 2, 3, 5]
    power_by_depth = {}
    for n_true in depth_grid:
        p = ident.power_curve(
            a_true=[2.0] * n_true, s=0.5, N0=1, rng=rng,
            sample_sizes=(100, 300, 1000, 5000), n_reps=100, n_boot=200,
        )
        power_by_depth[str(n_true)] = p
        print(f"N_true={n_true}: {p}")

    out = {
        "metadata": common.metadata(
            "scripts/compute_article_numbers.py", SEED, alpha=ALPHA,
            note="Simulation study (not an empirical depth test): the "
                 "replicate/bootstrap counts below are those of the study design."),
        "cumulant_check": {"a": a, "s": 0.5, "n": 4_000_000,
                           "empirical": kstats, "theory": theory},
        "kappa_ratio_curve": {"a": a_grid.tolist(), "ratio": ratio,
                              "ratio_s1": ratio_s1},
        "power_curve": {"N_true": 3, "N0": 1, "a": 2.0, "s": 0.5, "n_boot": 300,
                        "by_n": _with_se(powers, 200)},
        "power_by_depth": {"N0": 1, "a": 2.0, "s": 0.5, "n_boot": 200,
                           "sample_sizes": [100, 300, 1000, 5000],
                           "by_N_true": {k: _with_se(v, 100)
                                         for k, v in power_by_depth.items()}},
    }
    common.write_json("compute_article_numbers", out)
    print("cumulant check:", out["cumulant_check"])


if __name__ == "__main__":
    main()
