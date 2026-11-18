"""The network application, article draft §05: mixed-Poisson degrees generated
by the same cascade as the reservoir temperature, and the depth estimated
from the degree sequence alone -- without ever observing community labels.
"""

from __future__ import annotations

import numpy as np
import pytest

from htheory import network as net


def test_log_factorial_moment_theory_is_additive_across_levels(s_class):
    single = net.log_factorial_moment_theory([1.7], s_class, k=2)
    triple = net.log_factorial_moment_theory([1.7, 1.7, 1.7], s_class, k=2)
    assert triple == pytest.approx(3 * single)


def test_log_factorial_moment_theory_is_zero_for_N_equals_zero(s_class):
    """No levels: ``W = 1`` deterministically, so ``E[W^k] = 1`` for every ``k``."""
    assert net.log_factorial_moment_theory([], s_class, k=2) == pytest.approx(0.0)


def test_s1_branch_rejects_k_too_large_for_a():
    """``E[W^k]`` diverges in the ``s = 1`` branch once ``k >= a + 1`` -- a
    heavy-tail fact about the model, not a numerical bug."""
    with pytest.raises(ValueError, match="k = "):
        net.log_factorial_moment_theory([1.5], 1.0, k=3)


def test_mean_degree_equals_lambda_exactly(s_class, rng):
    """``E[W] = 1`` by the cascade's own mean-preservation, Eq. (30) of the
    original paper, so ``E[D] = lambda`` with no correction term."""
    lam = 12.0
    D = net.sample_degrees([1.5, 2.0], s_class, lam, rng, size=300_000)
    assert D.mean() == pytest.approx(lam, rel=0.02)


def test_empirical_log_factorial_moments_match_theory(s_class, rng):
    # a = 3.0, not 2.0: the s = 1 branch needs a > k - 1 = 2 for E[W^3] to
    # exist at all (see _level_log_moment), so a smaller a would make this
    # test's own theoretical reference undefined for that branch.
    a, N, lam = 3.0, 2, 20.0
    D = net.sample_degrees([a] * N, s_class, lam, rng, size=300_000)
    lam_hat, L = net.empirical_log_factorial_moments(D, k_max=3)
    theory = [net.log_factorial_moment_theory([a] * N, s_class, k) for k in (2, 3)]
    assert lam_hat == pytest.approx(lam, rel=0.02)
    assert L == pytest.approx(theory, rel=0.15, abs=0.05)


@pytest.mark.slow
def test_fit_homogeneous_network_recovers_known_depth(s_class, rng):
    """The identifiability chain closes: depth recovered from degrees alone,
    the claim the original paper's Eq. (33)/Fig. 1-2 said was impossible."""
    a, N, lam = 3.0, 3, 15.0
    D = net.sample_degrees([a] * N, s_class, lam, rng, size=2_000_000)
    fit = net.fit_homogeneous_network(D, s_class)
    assert fit.N_round == N
    assert fit.a == pytest.approx(a, rel=0.15)
    assert fit.lam == pytest.approx(lam, rel=0.05)


@pytest.mark.slow
def test_test_depth_degrees_does_not_reject_the_true_depth(rng):
    """``fit_homogeneous_network``'s point estimate had no formal p-value;
    under H0 (the simulated depth *is* N0) the bootstrap test should not
    reject at the nominal level, up to Monte-Carlo noise from a single draw.
    """
    a, N0, lam = 2.5, 2, 20.0
    D = net.sample_degrees([a] * N0, 0.5, lam, rng, size=200_000)
    result = net.test_depth_degrees(D, N0, 0.5, rng, n_boot=500)
    assert result.p_value > 0.05
    assert result.a_hat == pytest.approx(a, rel=0.2)


@pytest.mark.slow
def test_test_depth_degrees_rejects_the_wrong_depth(rng):
    """Testing ``H0: N = 1`` against data actually generated at depth 3 must
    reject most of the time -- otherwise the test has no power to do the one
    thing it exists for."""
    a, N_true, lam = 2.5, 3, 20.0
    D = net.sample_degrees([a] * N_true, 0.5, lam, rng, size=200_000)
    result = net.test_depth_degrees(D, 1, 0.5, rng, n_boot=500)
    assert result.p_value < 0.05


def test_fit_homogeneous_network_distinguishes_shallow_from_deep(rng):
    """Two degree sequences with matched mean but different depth must not be
    fit to the same N -- otherwise the estimator would be useless for exactly
    the claim it is built to support."""
    s = 0.5
    D_shallow = net.sample_degrees([3.0], s, 15.0, rng, size=800_000)
    D_deep = net.sample_degrees([3.0, 3.0, 3.0], s, 15.0, rng, size=800_000)
    fit_shallow = net.fit_homogeneous_network(D_shallow, s)
    fit_deep = net.fit_homogeneous_network(D_deep, s)
    assert fit_shallow.N_round == 1
    assert fit_deep.N_round == 3
