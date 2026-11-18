"""The Meijer G wrapper: known values, the underflow contract, and caching."""

from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from htheory.meijer import (
    UNDERFLOW_Z,
    MeijerGFailure,
    _meijerg_cached,
    meijerg,
    vectorize,
    with_dps,
)


def test_g_1_0_0_1_is_the_exponential():
    """``G^{1,0}_{0,1}(z | - ; 0) = exp(-z)`` -- the identity that makes N = 0
    fall out of Eq. (34) without a special case."""
    for z in (0.1, 1.0, 7.5):
        assert float(meijerg([[], []], [[0], []], z)) == pytest.approx(np.exp(-z), rel=1e-14)


def test_g_2_0_0_2_is_the_bessel_k_form():
    """``G^{2,0}_{0,2}(x | - ; a, 0) = 2 x^{a/2} K_a(2 sqrt x)``, Eq. (46)."""
    a = 1.3
    for x in (0.2, 1.0, 5.0):
        want = 2 * mp.mpf(x) ** (a / 2) * mp.besselk(a, 2 * mp.sqrt(x))
        assert float(meijerg([[], []], [[a, 0], []], x)) == pytest.approx(float(want), rel=1e-12)


def test_g_1_1_1_1_is_the_power_law():
    """``G^{1,1}_{1,1}(x | 1-a ; 0) = Gamma(a) (1+x)^{-a}``, the identity behind
    Eq. (57)."""
    a = 2.4
    for x in (0.3, 2.0, 20.0):
        want = float(mp.gamma(a)) * (1 + x) ** (-a)
        assert float(meijerg([[1 - a], []], [[0], []], x)) == pytest.approx(want, rel=1e-12)


def test_nonconvergence_far_out_returns_zero_when_opted_in():
    """Past ``underflow_above`` a failure means the value has underflowed."""
    huge = meijerg([[], []], [[1.0, 1.0, 1.0, 0], []], 1e12, underflow_above=UNDERFLOW_Z)
    assert huge >= 0
    assert huge < 1e-30


def test_nonconvergence_raises_by_default():
    """Silence is not the default. A caller that has not said where its
    underflow regime begins gets an exception, not a plausible-looking zero."""
    with pytest.raises(MeijerGFailure):
        meijerg([[-1.0, -1.0], []], [[0.0], []], 1e-4)


def test_nonconvergence_below_the_threshold_raises():
    """The tripwire.

    ``G^{1,2}_{2,1}(1e-4 | -1, -1 ; 0)`` is worth about 0.9996 -- mpmath cannot
    evaluate it at any working precision, because the parameters are
    degenerate. Returning 0 for it (which this wrapper used to do
    unconditionally) is not a small error, it is the opposite answer. Below
    ``underflow_above`` the wrapper must refuse rather than guess.
    """
    with pytest.raises(MeijerGFailure, match="not negligible"):
        meijerg([[-1.0, -1.0], []], [[0.0], []], 1e-4, underflow_above=UNDERFLOW_Z)


def test_argument_must_be_positive():
    with pytest.raises(ValueError, match="must be positive"):
        meijerg([[], []], [[0], []], 0.0)
    with pytest.raises(ValueError, match="must be positive"):
        meijerg([[], []], [[0], []], -1.0)


def test_results_are_cached_on_the_exact_argument():
    _meijerg_cached.cache_clear()
    meijerg([[], []], [[0.7, 0], []], 1.25)
    first = _meijerg_cached.cache_info()
    meijerg([[], []], [[0.7, 0], []], 1.25)
    second = _meijerg_cached.cache_info()
    assert second.hits == first.hits + 1
    assert second.misses == first.misses


def test_vectorize_preserves_shape_and_scalarness():
    f = vectorize(lambda x: mp.mpf(x) ** 2)
    assert isinstance(f(3.0), float)
    out = f(np.array([[1.0, 2.0], [3.0, 4.0]]))
    assert out.shape == (2, 2)
    assert out[1, 1] == pytest.approx(16.0)


def test_with_dps_raises_working_precision():
    @with_dps(50)
    def digits():
        return mp.mp.dps

    assert digits() == 50
