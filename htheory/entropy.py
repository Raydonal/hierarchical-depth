"""The hierarchical entropies: the generator ``s~_N^{(s)}(E)`` and ``S[p]``.

The inverse maximum-entropy construction of Section 3 runs::

    s'(B) = alpha + beta0 E(B)                          Eq. (12)
    s~(E) = (beta0 E - 1) B(E) + Btilde(E)              Eq. (18)
    S_p[p] = sum_i Z^{-1} s(Z p_i)   (+ irrelevant linear term)   Eq. (10)

with the gauge ``alpha = -1`` fixed by the flat-energy sum rule, Eqs. (15)-(16).
Specializing Eq. (18) to the two hierarchical classes gives Eq. (39), which is
what :func:`s_tilde` evaluates -- and Eq. (38) makes the second term the *same*
Meijer G as the first, with ``a`` shifted by ``3 - 4s``.

The map ``E -> (B(E), s~(E))`` is a parametric curve; that is how Figs. 1 and 2
are drawn, and it is all that is available in general, because ``E(B)`` has no
closed form except in the three exactly solvable cases (BG; ``s = 1/2, a = 1/2``;
the whole ``s = 1`` family).  :func:`E_of_B` inverts ``B`` numerically, which
:func:`s_of_B` and :mod:`htheory.maxent` use.

Concavity, Eq. (44), needs none of that: ``s''(B) = beta0 / B'(E) < 0`` follows
from ``B'(E) < 0``, Eq. (43), for every admissible ``(N, a)`` and both classes.
"""

from __future__ import annotations

import numpy as np
import mpmath as mp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .meijer import vectorize
from .background import _as_a, _omega
from .statfactor import B, B_tilde, B_prime, shift

__all__ = [
    "s_tilde",
    "E_of_B",
    "make_E_of_B",
    "make_interpolants",
    "s_of_B",
    "s_prime_of_B",
    "parametric_curve",
    "S_simplex",
    "s_boltzmann_shannon",
    "s_tsallis_kernel",
    "q_indices",
    "tsallis_entropy",
]


def s_tilde(E, a, s, beta0: float = 1.0):
    """``s~_N^{(s)}(E) = (beta0 E - 1) B_N(E; a) + B_N(E; a + 3 - 4s)``, Eq. (39).

    Written out per class this is Eq. (40) for ``s = 1/2`` and Eq. (41) for
    ``s = 1``; both are the same two Meijer G calls with shifted ``a``.

    Boundary values follow from Eq. (7): ``s~(0) = -1 + 1 = 0`` -- the gauge
    ``s(1) = 0`` of Eq. (14) -- and ``s~(inf) = 0``, since ``beta0 E B -> 0``
    (exponentially for ``s = 1/2``; as ``E^{-a}`` for ``s = 1``).
    """
    E_arr = np.asarray(E, dtype=float)
    return (beta0 * E_arr - 1.0) * B(E_arr, a, s, beta0) + B_tilde(E_arr, a, s, beta0)


def parametric_curve(a, s, beta0: float = 1.0, E=None, n: int = 400, E_max: float = 60.0):
    """The pair ``(B(E), s~(E))`` on an energy grid -- the curve of Figs. 1, 2.

    The parametric form is not a convenience: outside the three solvable cases
    ``E(B)`` simply has no closed form, so ``s(B)`` exists only as this curve or
    through :func:`E_of_B`.  The default grid is geometric in ``1 + beta0 E`` so
    that both the ``B ~ 1`` knee and the tail are resolved.
    """
    if E is None:
        E = (np.geomspace(1.0, 1.0 + beta0 * E_max, n) - 1.0) / beta0
    E = np.asarray(E, dtype=float)
    return B(E, a, s, beta0), s_tilde(E, a, s, beta0)


# --------------------------------------------------------------------------
# Numerical inversion E(B)
# --------------------------------------------------------------------------


def _E_of_B_scalar(b, a, s, beta0, xtol):
    if not 0.0 < b <= 1.0:
        raise ValueError(f"B must lie in (0, 1], got {b}")
    if b == 1.0:
        return 0.0
    # B is strictly decreasing, Eq. (43), so bracket by doubling.
    hi = 1.0 / beta0
    for _ in range(200):
        if float(B(hi, a, s, beta0)) < b:
            break
        hi *= 2.0
    else:  # pragma: no cover - would mean B does not vanish, contradicting Eq. (7)
        raise RuntimeError(f"could not bracket E(B) for B = {b}")
    return brentq(lambda e: float(B(e, a, s, beta0)) - b, 0.0, hi, xtol=xtol, rtol=1e-13)


def E_of_B(Bval, a, s, beta0: float = 1.0, xtol: float = 1e-12):
    """Invert Eq. (34) for ``E(B)`` by bisection.

    Well posed for every ``(N, a)`` and both classes because ``B`` is strictly
    decreasing from ``B(0) = 1`` to ``B(inf) = 0``, Eqs. (7) and (43).
    """
    a = _as_a(a)
    return vectorize(lambda b: _E_of_B_scalar(b, a, s, beta0, xtol))(Bval)


class _Interpolants:
    """Tabulated ``E(B)`` and ``s(B)`` sharing one energy grid.

    Both are exact formulas evaluated on the grid -- ``E`` by construction,
    ``s`` from Eq. (39) -- and only the *interpolation* between grid points is
    approximate.  Monotone (PCHIP) splines are used so that ``E(B)`` keeps the
    strict decrease of Eq. (43) and never manufactures a non-physical inversion.
    """

    def __init__(self, grid, b, s_vals, beta0):
        keep = b > 0
        x = np.log(b[keep])[::-1]
        rising = np.concatenate([[True], np.diff(x) > 0])
        self.x = x[rising]
        self.beta0 = beta0
        self._lo, self._hi = float(self.x[0]), float(self.x[-1])
        self._E = PchipInterpolator(self.x, np.log1p(beta0 * grid[keep])[::-1][rising])
        # s~ is *not* tabulated directly. It vanishes at both ends -- s~(1) = 0
        # by the gauge of Eq. (14), s~(0) = 0 because B does -- so a spline
        # through s~ itself keeps absolute accuracy and loses all relative
        # accuracy in the tail, which is exactly where the maximization needs
        # it. What is tabulated is the bounded ratio
        #
        #     u(B) = s~ / (B (1 + beta0 E)),
        #
        # which runs from 0 at B = 1 to 1 as B -> 0 (there s~ -> beta0 E B by
        # Eq. (39)), and s~ is rebuilt from it. Relative error then tracks the
        # ~1e-6 of E(B) instead of degrading with B.
        e_grid = (beta0 * grid[keep])[::-1][rising]
        s_grid = s_vals[keep][::-1][rising]
        self._u = PchipInterpolator(self.x, s_grid / (np.exp(self.x) * (1.0 + e_grid)))
        self.B_min = float(np.exp(self.x[0]))

    def _logB(self, Bval):
        """``log B`` clipped to the tabulated range.

        PCHIP extrapolates, and extrapolating ``log(1 + beta0 E)`` a few decades
        past the last knot overflows the exponential on the way back.  An
        optimizer *will* probe outside the range, so clamping there is the
        difference between a saturating value and a NaN.
        """
        return np.clip(np.log(np.asarray(Bval, dtype=float)), self._lo, self._hi)

    def E_of_B(self, Bval):
        """``E(B)``, the inverse of Eq. (34)."""
        out = np.maximum(np.expm1(self._E(self._logB(Bval))) / self.beta0, 0.0)
        return out if np.ndim(Bval) else float(out)

    def s_of_B(self, Bval):
        """``s(B) = s~(E(B))``, Eq. (39) read as a function of the factor."""
        logb = self._logB(Bval)
        out = self._u(logb) * np.exp(logb) * (1.0 + self.beta0 * self.E_of_B(np.exp(logb)))
        return out if np.ndim(Bval) else float(out)


def make_interpolants(a, s, beta0: float = 1.0, n: int = 300, B_min: float = 1e-12):
    """Tabulate ``E(B)`` and ``s(B)`` once, on a shared grid.

    :func:`E_of_B` costs a bracketed bisection -- tens of Meijer G evaluations
    -- per point, and :func:`s_of_B` costs that plus two more.  That is fine for
    a handful of values and hopeless inside an optimizer loop.  This pays for
    ``n`` grid evaluations once.

    Valid on ``[B_min, 1]``; below that the splines extrapolate and should not
    be trusted, so keep ``B_min`` under anything the caller will ask for.
    """
    a = _as_a(a)
    E_max = 1.0 / beta0
    for _ in range(400):
        if float(B(E_max, a, s, beta0)) < B_min:
            break
        E_max *= 2.0
    else:  # pragma: no cover - contradicts B(inf) = 0, Eq. (7)
        raise RuntimeError("B does not reach B_min; check Eq. (34)")

    # Log-spaced in E, not in 1 + beta0 E: for N >= 1 the factor leaves E = 0
    # non-analytically (B = 1 - O(E^{1/(N+1)}) in the stretched-exponential
    # class), so a grid that merely starts at zero resolves the B ~ 1 knee --
    # the high-probability end, the one the maximization cares about -- to only
    # a digit or two.
    grid = np.concatenate([[0.0], np.geomspace(1e-9 / beta0, E_max, n)])
    b = np.asarray(B(grid, a, s, beta0), dtype=float)
    s_vals = np.asarray(s_tilde(grid, a, s, beta0), dtype=float)
    return _Interpolants(grid, b, s_vals, beta0)


def make_E_of_B(a, s, beta0: float = 1.0, n: int = 300, B_min: float = 1e-12):
    """``E(B)`` alone, as a fast callable. See :func:`make_interpolants`."""
    return make_interpolants(a, s, beta0, n, B_min).E_of_B


def s_of_B(Bval, a, s, beta0: float = 1.0):
    """``s(B) = s~(E(B))``: the entropy density as a function of the factor."""
    return s_tilde(E_of_B(Bval, a, s, beta0), a, s, beta0)


def s_prime_of_B(Bval, a, s, beta0: float = 1.0):
    """``s'(B) = alpha + beta0 E(B) = beta0 E(B) - 1``, Eqs. (12) and (16).

    Exact, not a finite difference: the inverse function is the derivative here.
    This is what makes the maximization in :mod:`htheory.maxent` cheap -- the
    gradient of ``S[p]`` costs one inversion per point and no differentiation.
    """
    return beta0 * np.asarray(E_of_B(Bval, a, s, beta0), dtype=float) - 1.0


# --------------------------------------------------------------------------
# The entropy on the probability simplex, Eq. (10)
# --------------------------------------------------------------------------


def S_simplex(p, a, s, beta0: float = 1.0, Z: float | None = None, c: float = 0.0):
    """``S_p[p] = sum_i Z^{-1} s(Z p_i) + c sum_i p_i``, Eq. (10).

    The linear term is irrelevant under normalization -- it shifts ``alpha`` by
    ``c`` and leaves ``beta0`` alone, Eq. (11) -- and defaults to zero.

    ``Z`` is the partition function of Eq. (28), ``sum_i B(E_i)``.  It is *not*
    the ``1/beta0`` of the flat-energy sum rule, Eq. (27); the two coincide only
    under the flat measure ``rho(E) = 1``.  It must be supplied because the
    energy spectrum is not knowable from ``(a, s)`` alone.
    """
    p = np.asarray(p, dtype=float)
    if Z is None:
        raise ValueError("Z = sum_i B(E_i), Eq. (28), must be supplied")
    x = Z * p
    inside = x > 0
    out = np.zeros_like(x)
    out[inside] = np.asarray(s_of_B(x[inside], a, s, beta0), dtype=float) / Z
    return float(np.sum(out) + c * np.sum(p))


# --------------------------------------------------------------------------
# The exactly solvable cases
# --------------------------------------------------------------------------


def s_boltzmann_shannon(Bval):
    """``s(B) = -B ln B``, Eq. (21): the ``N = 0`` member, recovered exactly."""
    x = np.asarray(Bval, dtype=float)
    return np.where(x > 0, -x * np.log(np.where(x > 0, x, 1.0)), 0.0)


def s_tsallis_kernel(Bval, a: float):
    """``s(B) = -(a + 1)(B - B^{a/(1+a)})``, Eq. (63).

    The ``N = 1``, ``s = 1`` case, where the construction closes in elementary
    form: Eq. (59) gives ``Btilde = B^{a/(a+1)}`` and Eq. (60) gives ``E(B)``.
    Summed over the simplex this is exactly the Tsallis entropy, Eqs. (64), (67).
    """
    x = np.asarray(Bval, dtype=float)
    q, _ = q_indices(a)
    return -(a + 1.0) * (x - np.power(x, q))


def q_indices(a: float):
    """``(q, q_stat)`` from Eq. (65): ``q = a/(a+1) < 1``, ``q_stat = 2 - q > 1``.

    The two indices must not be conflated.  ``q`` is the *entropy* index -- the
    exponent in ``S_q``, obtained with the standard linear mean-energy
    constraint.  ``q_stat = (a+2)/(a+1)`` is the *distribution* index, the one
    labelling the heavy-tailed factor of Eq. (66).  They are additive duals.
    """
    q = a / (a + 1.0)
    q_stat = 2.0 - q
    return q, q_stat


def tsallis_entropy(p, a: float):
    """``S_q[p] = (1 - sum_i p_i^q)/(q - 1) = (a+1)(sum_i p_i^{a/(a+1)} - 1)``,
    Eqs. (64) and (67), with the *entropy* index ``q = a/(a+1)``."""
    p = np.asarray(p, dtype=float)
    q, _ = q_indices(a)
    return float((1.0 - np.sum(np.power(p, q))) / (q - 1.0))
