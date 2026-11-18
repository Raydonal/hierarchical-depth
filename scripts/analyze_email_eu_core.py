"""Recover the depth estimate on the email-Eu-core network.

Analyzes the out-/in-/total-degree of the SNAP email-Eu-core
network (1005 nodes, 25571 directed edges, no self-loops filtered), fit under the
homogeneous mixed-Poisson model for both universality classes, plus a formal
bootstrap p-value (``network.test_depth_degrees``) for whichever branch's fit
succeeds.

Run from the repo root:

    PYTHONPATH=. .venv/bin/python3 scripts/analyze_email_eu_core.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts import common, data_loaders  # noqa: E402

SEED = 20260905
SYSTEM = "email-Eu-core"


def main():
    rng_test = np.random.default_rng(SEED)
    rng_ci = np.random.default_rng([SEED, 1])
    series = []
    for name, D in data_loaders.email_degrees().items():
        # All 1005 vertices are kept, including those with degree 0 in the
        # series at hand: a zero is a legitimate mixed-Poisson outcome.
        series.append(common.analyze_degree_series(
            SYSTEM, name, D, rng_test, rng_ci,
            n_zero=int(np.sum(D == 0)), zero_policy="kept"))
    common.write_json("analyze_email_eu_core", {
        "metadata": common.metadata(
            "scripts/analyze_email_eu_core.py", SEED,
            rng_streams={"depth_test": "default_rng(SEED)",
                         "bootstrap_ci": "default_rng([SEED, 1])"},
            n_boot_test=common.N_BOOT_TEST, n_boot_ci=common.N_BOOT_CI,
            note=common.SELECTION_NOTE),
        "series": series,
    })


if __name__ == "__main__":
    main()
