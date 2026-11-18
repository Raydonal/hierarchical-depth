"""The generalized statistical factor ``B_N^{(s)}(E; a)``, Eq. (34).

``B_N(E)`` is the Boltzmann factor averaged over the innermost reservoir
temperature, Eq. (24)::

    B_N(E) = int_0^inf dT_N exp(-E/T_N) f_N(T_N)

with ``f_N`` from :mod:`htheory.background`.  Compounding the two chains of
Eq. (33) gives the two closed forms of Eq. (34), both Meijer G.

Two invariants anchor the whole package and are checked in
``tests/test_statfactor.py``:

* the flat-energy sum rule, Eq. (27): ``int_0^inf B_N(E) dE = <T_N> = 1/beta0``
  for every ``N``, every ``a`` and both classes -- equivalently ``Btilde(0) = 1``.
  This is *not* the partition function: ``Z = sum_i B(E_i)``, Eq. (28), needs the
  density of states, and the two coincide only under the flat measure rho(E) = 1;
* the boundary conditions, Eq. (7): ``B(0) = 1``, ``B(inf) = 0``, together with
  ``B'(E) < 0``, Eq. (43), which is what makes ``s''(B) < 0``, Eq. (44).

``N = 0`` needs no special case: with empty parameter vectors ``omega = 1`` and
Eq. (34) collapses to ``G^{1,0}_{0,1}(beta0 E | - ; 0) = exp(-beta0 E)``, the
Boltzmann factor of Eq. (36).
"""

from __future__ import annotations

import numpy as np
import mpmath as mp

from .meijer import UNDERFLOW_Z, MeijerGFailure, meijerg, vectorize
from .background import _as_a, _omega, f_N, sample_TN

__all__ = [
    "B",
    "B_tilde",
    "B_prime",
    "shift",
    "omega_of",
    "B_quadrature",
    "B_montecarlo",
    "B_K_distribution",
    "B_stretched_exact",
    "B_q_exponential",
    "B_tilde_q_exponential",
    "tail_exponent_theta",
]


def omega_of(a) -> float:
    """``omega = prod a_i``, Eq. (33). The scalar that Eq. (38) leaves alone."""
    return _omega(_as_a(a))


def shift(s) -> int:
    """The parameter shift of Eq. (38): ``a -> a + 3 - 4s``.

    ``+1`` for the stretched-exponential class, ``-1`` for the power-law class.
    It moves the *vector* ``a`` -- and with it ``Gamma(a)`` -- but never the
    scalar ``omega``; see :func:`B_tilde`.
    """
    if s == 0.5:
        return 1
    if s == 1:
        return -1
    raise ValueError(f"s must be 1/2 or 1, got {s}")


def _B_powerlaw_series(E, a, beta0, omega, tol: float = 1e-18, kmax: int = 80):
    """``B_N^{(1)}`` from the inverse moments of ``T_N``, for small ``beta0 E``.

    For the inverse-gamma cascade ``T_N = omega T_0 / prod Y_i`` with
    ``Y_i ~ Gamma(a_i + 1, 1)``, so the *inverse* moments are elementary::

        E[T_N^{-k}] = prod_i Gamma(a_i + 1 + k) / Gamma(a_i + 1) / (omega T_0)^k

    and expanding ``exp(-E/T_N)`` term by term in Eq. (24) gives

        B = sum_k (-z)^k / k! * prod_i Gamma(a_i + 1 + k)/Gamma(a_i + 1),
        z = beta0 E / omega.

    The series is asymptotic, not convergent: the coefficients grow like
    ``(k!)^N``, so for ``N >= 2`` the terms shrink, reach a minimum, and then
    diverge. Summed to that minimum it is extremely accurate at small ``z`` --
    which is exactly the regime where mpmath cannot evaluate the Meijer G of
    Eq. (34) at all, at any working precision, because the parameters
    ``(-a, 0)`` are degenerate whenever the ``a_i`` differ by integers.

    Returns ``(value, converged)``. ``converged`` is False once the terms start
    growing, and the value must then be discarded in favour of the G.
    """
    z = mp.mpf(beta0) * E / omega
    term = mp.mpf(1)
    total = mp.mpf(1)
    smallest = None
    for k in range(1, kmax):
        term *= -z / k
        for x in a:
            # Gamma(a + 1 + k)/Gamma(a + k) = a + k, applied k times in sequence.
            term *= mp.mpf(x) + k
        total += term
        size = abs(term)
        if size < tol * abs(total):
            return total, True
        if smallest is not None and size > smallest:
            return total, False
        smallest = size
    return total, False


def _B_scalar(E, a, s, beta0, omega):
    if E < 0:
        raise ValueError(f"energy must be nonnegative, got {E}")
    if E == 0:
        # Eq. (7). The z -> 0 limit of both G functions of Eq. (34) is the
        # leading term z^0, whose coefficient is exactly the reciprocal of the
        # prefactor, so B(0) = 1 for every N and a.
        return mp.mpf(1)
    if s == 0.5:
        # G^{N+1,0}_{0,N+1}( beta0 omega E | - ; a, 0 ) / Gamma(a)
        pref = mp.mpf(1)
        for x in a:
            pref /= mp.gamma(mp.mpf(x))
        b_upper = [mp.mpf(x) for x in a] + [mp.mpf(0)]
        return pref * meijerg(
            [[], []], [b_upper, []], mp.mpf(beta0) * omega * E,
            underflow_above=UNDERFLOW_Z,
        )
    if s == 1:
        # Small argument first: mpmath cannot evaluate this G there at all when
        # the a_i differ by integers, and used to return 0 for it -- turning
        # B = 0.9976 into 0.0 without a word. The series is exact in that
        # regime and cheaper; past it the G takes over. The two agree to full
        # precision where both work, which tests/test_statfactor.py pins.
        value, converged = _B_powerlaw_series(E, a, beta0, omega)
        if converged:
            return value
        # G^{1,N}_{N,1}( beta0 E / omega | -a ; 0 ) / Gamma(a + 1)
        pref = mp.mpf(1)
        for x in a:
            pref /= mp.gamma(mp.mpf(x) + 1)
        a_upper = [-mp.mpf(x) for x in a]
        return pref * meijerg(
            [a_upper, []], [[mp.mpf(0)], []], mp.mpf(beta0) * E / omega,
            underflow_above=UNDERFLOW_Z,
        )
    raise ValueError(f"s must be 1/2 or 1, got {s}")


def B(E, a, s, beta0: float = 1.0, omega: float | None = None):
    """``B_N^{(s)}(E; a)``, Eq. (34).

    ``omega`` defaults to ``prod a_i`` but can be overridden.  That override is
    the whole point of the signature: Eq. (38) shifts the vector ``a`` while
    holding ``omega`` at its *unshifted* value, so :func:`B_tilde` must be able
    to pass the two independently.  Recomputing ``omega`` from the shifted
    vector is the single easiest mistake to make here.
    """
    return _B(E, _as_a(a), s, beta0, omega)


def _B(E, a: np.ndarray, s, beta0, omega):
    """:func:`B` on an already-coerced ``a``, with no positivity check."""
    if omega is None:
        omega = _omega(a)
    return vectorize(lambda e: _B_scalar(e, a, s, beta0, omega))(E)


def B_tilde(E, a, s, beta0: float = 1.0):
    """``Btilde_N^{(s)}(E) = beta0 int_E^inf B_N(E') dE'``, Eqs. (37)-(38).

    Eq. (38) says this is the same function with ``a -> a + 3 - 4s`` and
    ``omega`` untouched.  Structurally it is a size-biasing of the reservoir
    temperature: ``beta0 int_E^inf B = E[exp(-E/T*)]`` with ``T*`` the
    length-biased version of ``T_N``, and size-biasing moves a gamma's shape
    while leaving its rate alone -- which is exactly why ``omega`` does not move.
    """
    a = _as_a(a)
    # The shifted vector is not a reservoir shape any more -- for s = 1 and
    # a_i <= 1 it is zero or negative -- so it bypasses the positivity check.
    return _B(E, a + shift(s), s, beta0, _omega(a))


def B_prime(E, a, s, beta0: float = 1.0, omega: float | None = None):
    """``dB/dE``, Eq. (43): strictly negative, hence ``s''(B) < 0``, Eq. (44).

    Eq. (43) writes this as ``-int_0^inf (dT/T) exp(-E/T) f_N(T)``, manifestly
    negative; numerically it is cheaper to differentiate Eq. (34) directly, and
    accurate enough for the monotonicity checks that use it.
    """
    a = _as_a(a)
    if omega is None:
        omega = _omega(a)

    def d(e):
        h = max(1e-6, 1e-6 * e)
        return mp.diff(lambda x: _B_scalar(float(x), a, s, beta0, omega), mp.mpf(e), h=h)

    return vectorize(d)(E)


# --------------------------------------------------------------------------
# Independent evaluations, for cross-validation only
# --------------------------------------------------------------------------


def B_quadrature(E, a, s, beta0: float = 1.0, T_max: float = 400.0):
    """Eq. (24) by direct quadrature against :func:`~htheory.background.f_N`.

    Slow -- every integrand evaluation is an mpmath Meijer G -- and used only to
    check Eq. (34).  ``T_max`` (in units of ``T0``) truncates the upper tail,
    which is heavy for ``s = 1``: the missing mass is ``O(T_max^{-min a_i})``
    there, and exponentially small for ``s = 1/2``.
    """
    a = _as_a(a)
    T0 = 1.0 / beta0

    def integrand(T):
        return mp.e ** (-mp.mpf(E) / T) * mp.mpf(float(f_N(float(T), a, s, T0)))

    return float(mp.quad(integrand, [0, 0.25 * T0, T0, 4 * T0, T_max * T0]))


def B_montecarlo(E, a, s, rng, beta0: float = 1.0, n: int = 200_000):
    """Eq. (24) by ancestral Monte Carlo: ``mean(exp(-E/T_N))``.

    Returns ``(estimate, standard_error)``.  The cheapest fully independent
    check that the corrected branch assignment of Eq. (33) is the right one,
    since it never touches a Meijer G.
    """
    T = sample_TN(a, s, rng, T0=1.0 / beta0, size=n)
    vals = np.exp(-np.asarray(E, dtype=float)[..., None] / T)
    est = vals.mean(axis=-1)
    err = vals.std(axis=-1, ddof=1) / np.sqrt(n)
    if np.ndim(E) == 0:
        return float(est), float(err)
    return est, err


# --------------------------------------------------------------------------
# The N = 1 closed forms, Section 5.2 -- the only exactly solvable cases
# --------------------------------------------------------------------------


def B_K_distribution(E, a: float, beta0: float = 1.0):
    """``B_1^{(1/2)}(E; a) = 2 (w b0 E)^{a/2} K_a(2 sqrt(w b0 E)) / Gamma(a)``,
    Eq. (47), with ``omega = a``.  The canonical K-distribution.

    Neither ``E(B)`` nor ``Btilde(B)`` is available in closed form here -- only
    the factor itself is.
    """

    def scalar(e):
        x = mp.mpf(beta0) * a * e
        if x == 0:
            return mp.mpf(1)
        return 2 * x ** (mp.mpf(a) / 2) * mp.besselk(a, 2 * mp.sqrt(x)) / mp.gamma(a)

    return vectorize(scalar)(E)


def B_stretched_exact(E, beta0: float = 1.0):
    """``B_1^{(1/2)}(E; 1/2) = exp(-sqrt(2 beta0 E))``, Eq. (52).

    The one point where the K-distribution collapses to a pure stretched
    exponential, via ``K_{1/2}(z) = sqrt(pi/2z) exp(-z)``, Eq. (51).
    """
    return np.exp(-np.sqrt(2.0 * beta0 * np.asarray(E, dtype=float)))


def B_q_exponential(E, a: float, beta0: float = 1.0):
    """``B_1^{(1)}(E; a) = (1 + beta0 E / omega)^{-a-1}``, Eq. (57), ``omega = a``.

    A ``q_stat``-exponential with ``q_stat = (a + 2)/(a + 1) > 1``, Eq. (66) --
    *not* the entropy index ``q = a/(a+1) < 1``; see
    :func:`~htheory.entropy.q_indices`.
    """
    return np.power(1.0 + beta0 * np.asarray(E, dtype=float) / a, -a - 1.0)


def B_tilde_q_exponential(E, a: float, beta0: float = 1.0):
    """``Btilde_1^{(1)}(E; a) = (1 + beta0 E / omega)^{-a}``, Eq. (58).

    Equivalently ``B^{a/(a+1)}``, Eq. (59) -- the only nontrivial case where
    ``Btilde`` is an explicit function of ``B`` rather than of ``E``.
    """
    return np.power(1.0 + beta0 * np.asarray(E, dtype=float) / a, -a)


def tail_exponent_theta(a) -> float:
    """``theta = sum_i (a_i + 1/2 - N)/(N + 1)``: the algebraic prefactor
    exponent of the stretched-exponential tail, Eq. (35)."""
    a = _as_a(a)
    N = a.size
    return float(np.sum(a + 0.5 - N) / (N + 1))
