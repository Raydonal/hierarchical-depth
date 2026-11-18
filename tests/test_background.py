"""The fluctuating background: conditionals, the marginal, and the SDE.

The load-bearing test here is :func:`test_eq33_branch_assignment_is_corrected`,
which pins the erratum: the branch labels of Eq. (33) as printed are crossed.
"""

from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest
from scipy import stats

from htheory import background as bg
from htheory.meijer import meijerg


def test_conditional_is_gamma_for_s_half():
    """Eq. (31): Gamma(shape a_i, rate a_i/T_{i-1})."""
    T = np.array([0.3, 1.0, 2.5])
    got = bg.conditional_pdf(T, 1.7, 1.2, 0.5)
    want = stats.gamma.pdf(T, 1.2, scale=1.7 / 1.2)
    assert got == pytest.approx(want, rel=1e-12)


def test_conditional_is_inverse_gamma_for_s_one():
    """Eq. (32): inverse-gamma(shape a_i + 1, scale a_i T_{i-1})."""
    T = np.array([0.3, 1.0, 2.5])
    got = bg.conditional_pdf(T, 1.7, 1.2, 1.0)
    want = stats.invgamma.pdf(T, 2.2, scale=1.2 * 1.7)
    assert got == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_both_conditionals_have_mean_T_prev(s):
    """The linear drift of Eq. (29) puts every level at the mean temperature of
    the one above it, so <T_i> = T_0 for all i."""
    T_prev = 1.7
    grid = np.geomspace(1e-4, 400.0, 20000)
    pdf = bg.conditional_pdf(grid, T_prev, 1.2, s)
    assert np.trapezoid(grid * pdf, grid) == pytest.approx(T_prev, rel=2e-3)


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_marginal_matches_the_conditional_at_N_equals_1(s):
    """For N = 1 the marginal of Eq. (33) *is* the conditional, Eqs. (31)/(32),
    with T_{i-1} = T_0."""
    a = 1.7
    T = np.array([0.4, 1.0, 2.5])
    assert bg.f_N(T, [a], s) == pytest.approx(bg.conditional_pdf(T, 1.0, a, s), rel=1e-10)


def test_eq33_branch_assignment_is_corrected():
    """The erratum, pinned.

    Eq. (33) as printed labels the inverse-gamma (``G^{0,N}_{N,0}``) row
    ``s = 1/2`` and the gamma (``G^{N,0}_{0,N}``) row ``s = 1``.  That is
    backwards: the gamma chain is the ``s = 1/2`` conditional, Eq. (31).  This
    test fails if the assignment is ever flipped back.
    """
    a, T = 1.7, 0.9
    # The gamma row of Eq. (33), written out: w/(T0 Gamma(a)) G^{N,0}_{0,N}(w T | - ; a-1).
    gamma_row = a / mp.gamma(a) * meijerg([[], []], [[a - 1], []], a * T)
    assert bg.f_N(T, [a], 0.5) == pytest.approx(float(gamma_row), rel=1e-12)
    assert bg.f_N(T, [a], 0.5) == pytest.approx(stats.gamma.pdf(T, a, scale=1.0 / a), rel=1e-12)
    assert bg.f_N(T, [a], 1.0) == pytest.approx(
        stats.invgamma.pdf(T, a + 1, scale=a), rel=1e-12
    )
    # And the two rows are genuinely different functions, so the labels matter.
    assert bg.f_N(T, [a], 0.5) != pytest.approx(bg.f_N(T, [a], 1.0), rel=1e-3)


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_ancestral_sampling_reproduces_the_marginal(s, rng):
    """Eq. (30) by Monte Carlo against Eq. (33): the check that never touches a
    Meijer G on the sampling side."""
    a = [1.5, 2.0, 0.8]
    draws = bg.sample_TN(a, s, rng, size=400_000)
    assert draws.mean() == pytest.approx(1.0, abs=0.03)
    # Compare the empirical density with f_N on a few interior points.
    for T, width in ((0.5, 0.1), (1.0, 0.2), (2.0, 0.3)):
        emp = np.mean(np.abs(draws - T) < width / 2) / width
        assert emp == pytest.approx(float(bg.f_N(T, a, s)), rel=0.06)


def test_sample_conditional_draws_independently_per_parent(rng):
    """Each parent temperature must get its own variate.

    Letting numpy broadcast instead would hand a whole array a single gamma
    draw in the ``s = 1`` branch, where ``T_prev`` enters the scale rather than
    a numpy parameter -- a bug that silently collapses the cascade to one
    random factor and leaves the ``s = 1/2`` branch looking fine.
    """
    parents = np.full(5000, 1.0)
    for s in (0.5, 1.0):
        draws = bg.sample_conditional(parents, 1.3, s, rng)
        assert draws.shape == parents.shape
        assert np.unique(draws).size == draws.size


def test_kappa_from_a_inverts_the_definition():
    a, gamma = np.array([1.5, 2.0]), np.array([10.0, 100.0])
    kappa = bg.kappa_from_a(a, gamma)
    assert 2 * gamma / kappa**2 == pytest.approx(a)


def test_N_equals_zero_has_no_density():
    with pytest.raises(ValueError, match="delta"):
        bg.f_N(1.0, [], 0.5)


def test_a_must_be_positive():
    with pytest.raises(ValueError, match="positive"):
        bg.f_N(1.0, [-0.5], 0.5)


@pytest.mark.slow
@pytest.mark.parametrize("s", [0.5, 1.0])
def test_sde_stationary_law_is_the_conditional(s, rng):
    """Euler-Maruyama on Eq. (29) for a single level, against Eq. (31)/(32).

    One level only: the wide separation of time scales that Eq. (29) assumes
    (gamma_i >> gamma_{i-1}) makes a faithful multi-level simulation far more
    expensive than the marginal is worth here, and the chain is checked instead
    by ancestral sampling above.
    """
    a = 2.0
    paths = bg.simulate_sde([a], s, gamma=[50.0], rng=rng, dt=5e-4, n_steps=120_000, n_paths=16)
    draws = paths.ravel()
    assert draws.mean() == pytest.approx(1.0, rel=0.05)
    for T, width in ((0.6, 0.2), (1.0, 0.2)):
        emp = np.mean(np.abs(draws - T) < width / 2) / width
        assert emp == pytest.approx(float(bg.conditional_pdf(T, 1.0, a, s)), rel=0.15)


def test_sde_rejects_a_step_too_coarse_for_gamma():
    with pytest.raises(ValueError, match="too coarse"):
        bg.simulate_sde([2.0], 0.5, gamma=[100.0], rng=np.random.default_rng(0), dt=0.1)
