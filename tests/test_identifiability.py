"""Theorem 3: the cascade depth ``N`` and shapes ``a`` are identifiable from
``T_N`` alone, via the log-cumulants.

:func:`test_depth_test_does_not_reject_a_true_null` and
:func:`test_depth_test_rejects_a_false_null` are the two directions of the
check on the original paper's conclusion
that the depth is invisible from the foreground distribution is checked here
to be false as an estimation claim, not just as a plot.
"""

from __future__ import annotations

import numpy as np
import pytest

from htheory import identifiability as ident
from htheory.background import sample_TN


def test_log_cumulants_are_zero_for_N_equals_zero(s_class):
    """``N = 0``: ``T_N = T_0`` deterministically, no fluctuation at any order."""
    out = ident.log_cumulants_theory([], s_class, r_max=4)
    assert out == pytest.approx(np.zeros(3))


def test_log_cumulants_are_additive_across_levels(s_class):
    """Cumulants of independent summands add -- the fact the whole module
    rests on (see the derivation in :mod:`htheory.identifiability`)."""
    single = ident.log_cumulants_theory([1.7], s_class, r_max=4)
    triple = ident.log_cumulants_theory([1.7, 1.7, 1.7], s_class, r_max=4)
    assert triple == pytest.approx(3 * single)

    mixed = ident.log_cumulants_theory([1.2, 0.6, 2.5], s_class, r_max=4)
    parts = sum(
        ident.log_cumulants_theory([x], s_class, r_max=4) for x in (1.2, 0.6, 2.5)
    )
    assert mixed == pytest.approx(parts)


def test_fit_homogeneous_recovers_exact_theoretical_cumulants(s_class):
    """No sampling noise: feed the closed-form ``(kappa2, kappa3)`` for a known
    ``(N, a)`` back into the estimator and recover both to solver precision.

    This isolates the estimator's algebra from Monte Carlo variance -- if this
    fails, the bug is in :func:`fit_homogeneous` or :func:`_ratio_theory`, not
    in sampling.
    """
    for N, a in [(1, 1.0), (2, 1.0), (3, 0.7), (5, 3.0)]:
        k2, k3 = ident.log_cumulants_theory([a] * N, s_class, r_max=3)
        fit = ident.fit_homogeneous(k2, k3, s_class)
        assert fit.a == pytest.approx(a, rel=1e-6)
        assert fit.N == pytest.approx(N, rel=1e-6)
        assert fit.N_round == N


@pytest.mark.slow
def test_log_cumulants_match_montecarlo(a_vec, s_class, rng):
    """Cross-validates the closed form against ancestral sampling, at the
    precision the sum-rule-style tests elsewhere in the suite use."""
    if len(a_vec) == 0:
        pytest.skip("N = 0 has no fluctuation to estimate")
    T = sample_TN(a_vec, s_class, rng, size=2_000_000)
    theory = ident.log_cumulants_theory(a_vec, s_class, r_max=3)
    emp = ident.empirical_log_cumulants(T, r_max=3)
    assert emp == pytest.approx(theory, rel=0.03, abs=0.02)


@pytest.mark.slow
def test_fit_homogeneous_on_a_monte_carlo_sample(s_class, rng):
    N, a = 3, 2.0
    T = sample_TN([a] * N, s_class, rng, size=500_000)
    k2, k3 = ident.empirical_log_cumulants(T, r_max=3)
    fit = ident.fit_homogeneous(k2, k3, s_class)
    assert fit.a == pytest.approx(a, rel=0.1)
    assert fit.N == pytest.approx(N, rel=0.1)


def test_depth_test_does_not_reject_a_true_null(s_class, rng):
    N0, a = 1, 2.0
    T = sample_TN([a] * N0, s_class, rng, size=3000)
    result = ident.test_depth(T, N0, s_class, rng, n_boot=400)
    assert result.p_value > 0.05


def test_depth_test_rejects_a_false_null(s_class, rng):
    """True depth 3 tested against the null ``N0 = 1``: a large third
    log-cumulant discrepancy the bootstrap null distribution cannot absorb."""
    T = sample_TN([2.0, 2.0, 2.0], s_class, rng, size=3000)
    result = ident.test_depth(T, 1, s_class, rng, n_boot=400)
    assert result.p_value < 0.05


@pytest.mark.slow
def test_power_curve_increases_with_sample_size(s_class, rng):
    """The "defeito 2" evidence block: power against a true ``N = 3`` cascade
    at the shallow null ``N0 = 1`` should climb toward 1 as ``n`` grows."""
    powers = ident.power_curve(
        a_true=[2.0, 2.0, 2.0],
        s=s_class,
        N0=1,
        rng=rng,
        sample_sizes=(100, 1000),
        n_reps=60,
        n_boot=200,
    )
    assert powers[1000] >= powers[100]
    assert powers[1000] > 0.7
