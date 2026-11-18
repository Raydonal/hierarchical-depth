"""The entropy generator: Eq. (39), inversion, and the exactly solvable cases."""

from __future__ import annotations

import numpy as np
import pytest

from htheory import entropy as en
from htheory import statfactor as sf


def test_s_tilde_is_eq39(a_vec, s_class):
    """s~(E) = (beta0 E - 1) B(E; a) + B(E; a + 3 - 4s), Eq. (39)."""
    E = np.array([0.0, 0.4, 2.0, 9.0])
    want = (E - 1.0) * np.asarray(sf.B(E, a_vec, s_class), dtype=float) + np.asarray(
        sf.B_tilde(E, a_vec, s_class), dtype=float
    )
    assert en.s_tilde(E, a_vec, s_class) == pytest.approx(want, rel=1e-12, abs=1e-15)


def test_s_tilde_boundary_values(a_vec, s_class):
    """s~(0) = 0 by the gauge alpha = -1, Eqs. (14)-(16); s~(inf) = 0."""
    assert float(en.s_tilde(0.0, a_vec, s_class)) == pytest.approx(0.0, abs=1e-12)
    assert float(en.s_tilde(1e5, a_vec, s_class)) == pytest.approx(0.0, abs=1e-2)


def test_s_tilde_is_positive_in_between(a_vec, s_class):
    E = np.geomspace(1e-3, 50.0, 30)
    assert np.all(np.asarray(en.s_tilde(E, a_vec, s_class), dtype=float) > 0)


def test_E_of_B_inverts_B(a_vec, s_class):
    """B is strictly decreasing, Eq. (43), so the inverse is well posed."""
    E = np.array([0.05, 0.5, 3.0, 25.0])
    b = np.asarray(sf.B(E, a_vec, s_class), dtype=float)
    assert en.E_of_B(b, a_vec, s_class) == pytest.approx(E, rel=1e-8)
    assert float(en.E_of_B(1.0, a_vec, s_class)) == 0.0


def test_E_of_B_rejects_values_outside_the_domain(a_vec, s_class):
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="must lie in"):
            en.E_of_B(bad, a_vec, s_class)


def test_s_prime_of_B_is_the_inverse_function(a_vec, s_class):
    """s'(B) = alpha + beta0 E(B) = beta0 E(B) - 1, Eqs. (12) and (16).

    Checked against a finite difference of s(B) itself, so it is the *relation*
    being tested and not just the formula being restated.
    """
    b = 0.35
    h = 1e-5
    fd = (float(en.s_of_B(b + h, a_vec, s_class)) - float(en.s_of_B(b - h, a_vec, s_class))) / (2 * h)
    assert float(en.s_prime_of_B(b, a_vec, s_class)) == pytest.approx(fd, rel=1e-5)


def test_s_of_B_is_concave(a_vec, s_class):
    """Eq. (44). s''(B) = beta0/B'(E) < 0 for every admissible (N, a)."""
    b = np.linspace(0.02, 0.98, 40)
    s_vals = np.asarray(en.s_of_B(b, a_vec, s_class), dtype=float)
    second = np.diff(s_vals, 2)
    assert np.all(second < 0)


def test_parametric_curve_matches_pointwise_evaluation(a_vec, s_class):
    b, s_vals = en.parametric_curve(a_vec, s_class, n=25, E_max=30.0)
    assert np.all(np.diff(b) < 0)
    assert s_vals == pytest.approx(np.asarray(en.s_of_B(b, a_vec, s_class), dtype=float), rel=1e-6, abs=1e-9)


# --- the exactly solvable cases --------------------------------------------


def test_boltzmann_shannon_recovered():
    """Eq. (21): the N = 0 member is s(B) = -B ln B, exactly."""
    b = np.array([0.05, 0.3, 0.7, 1.0])
    assert en.s_of_B(b, [], 0.5) == pytest.approx(en.s_boltzmann_shannon(b), rel=1e-8, abs=1e-12)
    assert en.s_of_B(b, [], 1.0) == pytest.approx(en.s_boltzmann_shannon(b), rel=1e-8, abs=1e-12)


@pytest.mark.parametrize("a", [0.6, 1.0, 2.0, 3.5])
def test_tsallis_kernel_recovered_exactly(a):
    """Eq. (63): for N = 1, s = 1 the construction closes in elementary form."""
    b = np.array([0.05, 0.3, 0.7, 1.0])
    assert en.s_of_B(b, [a], 1.0) == pytest.approx(en.s_tsallis_kernel(b, a), rel=1e-8, abs=1e-12)


@pytest.mark.parametrize("a", [0.6, 1.0, 2.0, 3.5])
def test_the_two_tsallis_indices_are_distinct(a):
    """Eq. (65). q is the *entropy* index and is below 1; q_stat is the
    *distribution* index and is above 1. Conflating them is the standard error
    this test exists to prevent."""
    q, q_stat = en.q_indices(a)
    assert q == pytest.approx(a / (a + 1.0))
    assert q_stat == pytest.approx((a + 2.0) / (a + 1.0))
    assert q < 1.0 < q_stat
    assert q + q_stat == pytest.approx(2.0)


@pytest.mark.parametrize("a", [0.6, 2.0])
def test_simplex_entropy_is_tsallis_for_N1_powerlaw(a):
    """Eqs. (64) and (67): summed over the simplex, the N = 1 power-law kernel
    is the Tsallis entropy -- up to the positive factor and linear term of
    Eq. (68), both irrelevant to the maximizer."""
    E = np.linspace(0.1, 6.0, 9)
    b = np.asarray(sf.B(E, [a], 1.0), dtype=float)
    Z = float(b.sum())
    p = b / Z
    got = en.S_simplex(p, [a], 1.0, Z=Z)
    q, _ = en.q_indices(a)
    # Eq. (68): Z^{-1} s(Zp) = (a+1)(Z^{q-1} p^q - p).
    want = (a + 1.0) * (Z ** (q - 1.0) * float(np.sum(p**q)) - 1.0)
    assert got == pytest.approx(want, rel=1e-6)


def test_tsallis_entropy_formula():
    p = np.array([0.5, 0.3, 0.2])
    a = 2.0
    q, _ = en.q_indices(a)
    assert en.tsallis_entropy(p, a) == pytest.approx((1 - np.sum(p**q)) / (q - 1))
    assert en.tsallis_entropy(p, a) == pytest.approx((a + 1) * (np.sum(p**q) - 1))


def test_S_simplex_requires_Z():
    with pytest.raises(ValueError, match="Eq. \\(28\\)"):
        en.S_simplex(np.array([0.5, 0.5]), [1.0], 0.5)


def test_S_simplex_is_shannon_minus_log_Z_for_N0():
    """Eq. (22): S[p] = -sum p ln p - ln Z."""
    E = np.linspace(0.1, 6.0, 9)
    b = np.exp(-E)
    Z = float(b.sum())
    p = b / Z
    want = float(-np.sum(p * np.log(p)) - np.log(Z))
    assert en.S_simplex(p, [], 0.5, Z=Z) == pytest.approx(want, rel=1e-7)


# --- the tabulated interpolants --------------------------------------------


def test_interpolants_agree_with_the_exact_functions(a_vec, s_class):
    """Five digits *relative*, across ten decades of B, on a 300-point grid.

    Relative accuracy is the point. s~ vanishes at both ends, so a spline
    through it directly holds absolute error while relative error degrades to
    ~1e-2 in the tail -- fine for a plot, useless for a maximization that lives
    there. Tabulating the bounded ratio u = s~/(B(1 + beta0 E)) instead keeps
    relative error tracking E(B)'s own ~1e-6. This test is what pins that: it
    compares relative, with no absolute escape hatch.
    """
    interp = en.make_interpolants(a_vec, s_class, n=300, B_min=1e-10)
    b = np.array([0.9, 0.5, 0.1, 1e-3, 1e-6, 1e-9])
    assert interp.E_of_B(b) == pytest.approx(np.asarray(en.E_of_B(b, a_vec, s_class), dtype=float), rel=1e-4)
    assert interp.s_of_B(b) == pytest.approx(np.asarray(en.s_of_B(b, a_vec, s_class), dtype=float), rel=1e-4)


def test_interpolants_saturate_instead_of_extrapolating(a_vec, s_class):
    """Outside the tabulated range the splines must clamp.

    PCHIP extrapolates, and extrapolating log(1 + beta0 E) a few decades past
    the last knot overflows on the way back -- a NaN in the middle of an
    optimizer rather than a large number.
    """
    interp = en.make_interpolants(a_vec, s_class, n=120, B_min=1e-8)
    assert np.isfinite(interp.E_of_B(1e-40))
    assert interp.E_of_B(1e-40) == pytest.approx(interp.E_of_B(interp.B_min))
    assert np.isfinite(interp.s_of_B(1e-40))


def test_make_E_of_B_is_the_interpolant_alone(a_vec, s_class):
    f = en.make_E_of_B(a_vec, s_class, n=120, B_min=1e-8)
    assert f(0.5) == pytest.approx(en.make_interpolants(a_vec, s_class, n=120, B_min=1e-8).E_of_B(0.5))
