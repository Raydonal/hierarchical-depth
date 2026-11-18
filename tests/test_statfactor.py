"""The generalized statistical factor: invariants, closed forms, cross-checks.

Two invariants anchor the package; anything that breaks either is wrong:

1. the flat-energy sum rule, Eq. (27) -- ``int_0^inf B_N dE = <T_N> = 1/beta0``;
2. the boundary conditions, Eq. (7), with ``B'(E) < 0``, Eq. (43).
"""

from __future__ import annotations

import numpy as np
import pytest

from htheory import statfactor as sf
from htheory.entropy import s_tilde


def test_N_zero_is_the_boltzmann_factor():
    """Eq. (36), which must *fall out* of Eq. (34) with empty parameter
    vectors, not be special-cased."""
    E = np.array([0.0, 0.5, 2.0, 9.0])
    for beta0 in (1.0, 0.4, 3.0):
        assert sf.B(E, [], 0.5, beta0) == pytest.approx(np.exp(-beta0 * E), rel=1e-12)
        assert sf.B(E, [], 1.0, beta0) == pytest.approx(np.exp(-beta0 * E), rel=1e-12)


def test_boundary_conditions(a_vec, s_class):
    """Eq. (7): B(0) = 1 and B(inf) = 0."""
    assert float(sf.B(0.0, a_vec, s_class)) == pytest.approx(1.0, rel=1e-12)
    assert float(sf.B(1e4, a_vec, s_class)) < 1e-3
    assert float(sf.B(1e4, a_vec, s_class)) >= 0.0


def test_B_is_strictly_decreasing(a_vec, s_class):
    """Eq. (43). This is what makes s''(B) < 0, Eq. (44), with no closed form."""
    E = np.geomspace(1e-3, 200.0, 40)
    b = np.asarray(sf.B(E, a_vec, s_class), dtype=float)
    assert np.all(np.diff(b) < 0)
    assert np.all(np.asarray(sf.B_prime(E[:12], a_vec, s_class), dtype=float) < 0)


def test_sum_rule_via_B_tilde(a_vec, s_class):
    """Eq. (27) in its cheapest form.

    ``Btilde(0) = beta0 int_0^inf B dE``, Eq. (37), so the sum rule is exactly
    ``Btilde(0) = 1`` -- and by Eq. (38) that is ``B(0; a + 3 - 4s) = 1``, which
    is Eq. (7) for the shifted parameters. No quadrature needed.
    """
    for beta0 in (1.0, 0.7):
        assert float(sf.B_tilde(0.0, a_vec, s_class, beta0)) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.slow
def test_sum_rule_by_quadrature(a_vec, s_class):
    """Eq. (27) the honest way: integrate Eq. (34) numerically.

    The power-law class has an ``E^{-a-1}`` tail, so the truncated integral is
    compared with the exact remainder ``Btilde(E_max)/beta0`` rather than with
    ``1/beta0`` alone.
    """
    E_max = 400.0
    grid = np.concatenate([[0.0], np.geomspace(1e-6, E_max, 900)])
    b = np.asarray(sf.B(grid, a_vec, s_class), dtype=float)
    head = float(np.trapezoid(b, grid))
    tail = float(sf.B_tilde(E_max, a_vec, s_class))
    assert head + tail == pytest.approx(1.0, rel=2e-4)


def test_B_tilde_shifts_the_vector_but_not_omega(a_vec, s_class):
    """Eq. (38), and the one bookkeeping error the paper warns about.

    ``Btilde(E; a) = B(E; a + 3 - 4s)`` with ``omega = prod a_i`` held at its
    *unshifted* value. Recomputing omega from the shifted vector gives a
    different function; this pins that it is not done.
    """
    if not a_vec:
        pytest.skip("N = 0 has no vector to shift")
    a = np.asarray(a_vec, dtype=float)
    shifted = a + sf.shift(s_class)
    # Through the private entry point: for s = 1 the shift can drive entries to
    # zero or below, which the public B() rightly refuses as a reservoir shape
    # even though it is a fine Meijer G parameter.
    correct = sf._B(1.3, shifted, s_class, 1.0, sf.omega_of(a))
    assert float(sf.B_tilde(1.3, a, s_class)) == pytest.approx(float(correct), rel=1e-12)

    wrong_omega = float(np.prod(shifted))
    if wrong_omega <= 0:
        # For s = 1 the shift can send prod(a) to zero or negative, so the
        # wrong reading is not merely a different number -- it has no value at
        # all. Nothing left to compare.
        return
    wrong = sf._B(1.3, shifted, s_class, 1.0, wrong_omega)
    assert float(correct) != pytest.approx(float(wrong), rel=1e-6)


def test_shift_values():
    assert sf.shift(0.5) == 1
    assert sf.shift(1) == -1
    with pytest.raises(ValueError, match="s must be"):
        sf.shift(0.75)


def test_B_tilde_is_the_integral_of_B(a_vec, s_class):
    """Eq. (37) checked directly: Btilde(E1) - Btilde(E2) = beta0 int_E1^E2 B."""
    E1, E2 = 0.4, 6.0
    grid = np.linspace(E1, E2, 1200)
    integral = float(np.trapezoid(np.asarray(sf.B(grid, a_vec, s_class), dtype=float), grid))
    delta = float(sf.B_tilde(E1, a_vec, s_class)) - float(sf.B_tilde(E2, a_vec, s_class))
    assert delta == pytest.approx(integral, rel=1e-5)


# --- the N = 1 closed forms, Section 5.2 ------------------------------------


def test_K_distribution_closed_form():
    """Eq. (47): B_1^{(1/2)} is the K-distribution."""
    E = np.array([0.0, 0.3, 1.0, 4.0])
    for a in (0.7, 1.0, 2.5):
        assert sf.B(E, [a], 0.5) == pytest.approx(sf.B_K_distribution(E, a), rel=1e-10)


def test_stretched_exponential_at_a_one_half():
    """Eq. (52): at a = 1/2 the K-distribution collapses to exp(-sqrt(2 beta0 E))."""
    E = np.array([0.0, 0.3, 1.0, 4.0, 20.0])
    assert sf.B(E, [0.5], 0.5) == pytest.approx(sf.B_stretched_exact(E), rel=1e-10)
    assert sf.B_K_distribution(E, 0.5) == pytest.approx(sf.B_stretched_exact(E), rel=1e-10)


def test_q_exponential_closed_form():
    """Eq. (57): B_1^{(1)} = (1 + beta0 E/omega)^{-a-1}, with omega = a."""
    E = np.array([0.0, 0.3, 1.0, 4.0, 50.0])
    for a in (0.6, 1.0, 2.0):
        assert sf.B(E, [a], 1.0) == pytest.approx(sf.B_q_exponential(E, a), rel=1e-10)


def test_q_exponential_B_tilde_closed_form():
    """Eq. (58), and Eq. (59): Btilde = B^{a/(a+1)}."""
    E = np.array([0.0, 0.3, 1.0, 4.0, 50.0])
    for a in (0.6, 1.0, 2.0):
        want = sf.B_tilde_q_exponential(E, a)
        assert sf.B_tilde(E, [a], 1.0) == pytest.approx(want, rel=1e-10)
        assert want == pytest.approx(np.power(sf.B_q_exponential(E, a), a / (a + 1.0)), rel=1e-12)


def test_beta0_only_sets_the_energy_scale(a_vec, s_class):
    """B depends on E and beta0 only through beta0 E."""
    assert float(sf.B(2.0, a_vec, s_class, beta0=0.5)) == pytest.approx(
        float(sf.B(1.0, a_vec, s_class, beta0=1.0)), rel=1e-12
    )


def test_tail_exponent_theta():
    """theta = sum_i (a_i + 1/2 - N)/(N+1), Eq. (35)."""
    assert sf.tail_exponent_theta([1.0, 1.0]) == pytest.approx((1.5 - 2 + 1.5 - 2) / 3)


# --- fully independent cross-validation ------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("s", [0.5, 1.0])
def test_compounding_reproduces_eq34_by_quadrature(s):
    """Eq. (24) integrated against Eq. (33), versus the closed form Eq. (34).

    This is the check that the corrected branch assignment is the right one:
    the printed labels give the *other* class here.

    The s = 1 tail of f_N decays like T^(-min a_i - 1) = T^(-1.8) here; at the
    default cutoff T_max = 400 the omitted mass is 2.5e-5, i.e. 4e-4 relative
    to B(3.0) = 0.059, so the cutoff is raised to 4e4 (residual ~1e-7).
    """
    a = [1.5, 2.0, 0.8]
    for E in (0.3, 1.0, 3.0):
        q = sf.B_quadrature(E, a, s, T_max=4e4)
        assert q == pytest.approx(float(sf.B(E, a, s)), rel=2e-4)


@pytest.mark.slow
@pytest.mark.parametrize("s", [0.5, 1.0])
def test_montecarlo_reproduces_eq34(s, rng):
    """Eq. (24) by ancestral Monte Carlo -- no Meijer G on the sampling side."""
    a = [1.5, 2.0, 0.8]
    for E in (0.3, 1.0, 3.0):
        est, err = sf.B_montecarlo(E, a, s, rng, n=400_000)
        assert abs(est - float(sf.B(E, a, s))) < 5 * err


def test_no_spurious_zeros_at_small_energy(a_vec, s_class):
    """Regression: B must stay near 1 as E -> 0, for both classes.

    The power-law branch of Eq. (34) is a Meijer G that mpmath cannot evaluate
    at small argument when the a_i differ by integers, and the wrapper used to
    answer 0 for it -- so B went from 0.9976 to 0.0 with no warning, and the
    sum rule quietly lost 4e-4 of its mass. This walks the small-energy end
    that the other tests skip.
    """
    E = np.geomspace(1e-9, 1e-1, 25)
    b = np.asarray(sf.B(E, a_vec, s_class), dtype=float)
    assert np.all(b > 0)
    assert np.all(np.diff(b) < 0)
    assert b[0] > 0.999
    assert b[0] < 1.0


@pytest.mark.parametrize("a", [[1.0], [1.0, 1.0], [1.5, 2.0, 0.8]])
def test_powerlaw_series_and_meijer_g_agree_where_both_work(a):
    """The small-argument series is the *same function* as Eq. (34), not an
    approximation to it: where mpmath can also evaluate the G, they agree to
    full working precision."""
    omega = sf.omega_of(a)
    for E in (0.02, 0.05, 0.2):
        value, converged = sf._B_powerlaw_series(E, np.asarray(a), 1.0, omega)
        if not converged:
            continue
        assert float(value) == pytest.approx(float(sf.B(E, a, 1.0)), rel=1e-12)


def test_energy_must_be_nonnegative(a_vec, s_class):
    with pytest.raises(ValueError, match="nonnegative"):
        sf.B(-1.0, a_vec, s_class)


def test_s_tilde_vanishes_at_both_ends(a_vec, s_class):
    """Eq. (39) with Eq. (7): s~(0) = -1 + 1 = 0, the gauge s(1) = 0 of Eq. (14)."""
    assert float(s_tilde(0.0, a_vec, s_class)) == pytest.approx(0.0, abs=1e-12)
