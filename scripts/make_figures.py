#!/usr/bin/env python3
"""Regenerate the parametric entropy curves ``(B(E), s(E))`` for ``a = 1`` and
``N = 0..5``, both universality classes (the original paper's Figs. 1 and 2).

Thin wrapper: the curves are now drawn as the two panels of
``figures/fig2_entropy.pdf`` by ``make_theory_figures.fig2_entropy`` in the
manuscript's figure style. The earlier one-file-per-class outputs
(``fig1_stretched_exponential``, ``fig2_power_law``) are no longer
regenerated; the existing files are left in place.

Run from the repo root::

    PYTHONPATH=. python3 scripts/make_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.make_theory_figures import fig2_entropy  # noqa: E402


def main():
    fig2_entropy()


if __name__ == "__main__":
    main()
