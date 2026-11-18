"""Numerical companion to *Hierarchical entropies* (Physica A **699** (2026) 131899)
and to the follow-on article on the identifiability of the depth (in preparation).

The modules follow the paper's logical chain and may only import backwards along
it::

    meijer -> background -> statfactor -> entropy -> maxent
                          -> identifiability
                          -> network

``meijer``          robust/cached Meijer G evaluation (every closed form is one).
``background``      the reservoir temperature law ``f_N(T_N)``, Eq. (33), and the
                    SDE hierarchy, Eq. (29), that produces it.
``statfactor``      the generalized statistical factor ``B_N^{(s)}(E; a)``, Eq. (34).
``entropy``         the entropy generator ``s~_N^{(s)}(E)``, Eq. (39), and ``S[p]``.
``maxent``          maximization of ``S[p]`` under the standard constraints.
``identifiability``  Theorem 3 (``docs/theorems.md``): log-cumulant identifiability
                    of the depth ``N``, point estimation, and a bootstrap depth test.
``network``         the same cascade as a mixed-Poisson network model, with a
                    depth estimator from the degree sequence alone.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [
    "meijer",
    "background",
    "statfactor",
    "entropy",
    "maxent",
    "identifiability",
    "network",
]
