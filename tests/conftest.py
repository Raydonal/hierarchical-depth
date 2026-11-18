"""Shared fixtures and parameter sets for the test suite.

The parameter sets are chosen to exercise the corners that actually break:
``N = 0`` (must fall out of Eq. (34), not be special-cased), integer ``a``
(coincident Meijer G parameters), ``a < 1`` (the ``s = 1`` shift of Eq. (38)
drives ``a - 1`` non-positive), and unequal ``a_i`` (catches anything that
silently assumes ``a_i = a``).
"""

from __future__ import annotations

import numpy as np
import pytest

CLASSES = [0.5, 1.0]

# (a, label). Kept short: each entry costs Meijer G evaluations in every test.
PARAM_SETS = [
    ([], "N=0"),
    ([1.0], "N=1 a=1"),
    ([0.5], "N=1 a=1/2"),
    ([2.0], "N=1 a=2"),
    ([1.0, 1.0], "N=2 a=1"),
    ([1.5, 2.0, 0.8], "N=3 unequal"),
]

A_VECTORS = [a for a, _ in PARAM_SETS]
A_IDS = [label for _, label in PARAM_SETS]


@pytest.fixture
def rng():
    return np.random.default_rng(20260826)


@pytest.fixture(params=A_VECTORS, ids=A_IDS)
def a_vec(request):
    return request.param


@pytest.fixture(params=CLASSES, ids=["s=1/2", "s=1"])
def s_class(request):
    return request.param
