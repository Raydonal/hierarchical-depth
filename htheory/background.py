"""The fluctuating background: the reservoir temperature law and its dynamics.

Three descriptions of the same object, in increasing order of directness:

* :func:`f_N` -- the closed form, Eq. (33), a Meijer G;
* :func:`sample_TN` -- ancestral sampling down the chain of conditionals,
  Eqs. (31)/(32), which is what Eq. (33) is the marginal of;
* :func:`simulate_sde` -- Euler-Maruyama on the SDE hierarchy, Eq. (29), whose
  stationary Fokker-Planck solution (with ``T_{i-1}`` held fixed) *is* the
  conditional.

The three agree only if the branch labels of Eq. (33) are read as corrected --
see :data:`ERRATUM_EQ33`.  ``tests/test_background.py`` pins that.
"""

from __future__ import annotations

import numpy as np
import mpmath as mp
from scipy import stats

from .meijer import UNDERFLOW_Z, UNDERFLOW_Z_SMALL, meijerg, vectorize

__all__ = [
    "ERRATUM_EQ33",
    "conditional_pdf",
    "sample_conditional",
    "f_N",
    "sample_TN",
    "simulate_sde",
    "kappa_from_a",
]

ERRATUM_EQ33 = """\
Eq. (33) as printed has its two branch labels crossed: the row labelled s = 1/2
is the form that compounds to the power-law factor, and the row labelled s = 1
the one that compounds to the stretched exponential -- inconsistent with
Eqs. (31)/(32) and (34).  The corrected assignment, implemented here, is

    s = 1/2 (gamma chain)     f_N = w/(T0 G(a))     G^{N,0}_{0,N}( w T_N/T0 | - ; a-1  )
    s = 1   (inverse-gamma)   f_N = 1/(T0 w G(a+1)) G^{0,N}_{N,0}( T_N/(w T0) | -a-1 ; - )

verified to 8 digits against quadrature of Eq. (24) and against Monte-Carlo
ancestral sampling.\
"""


def _as_a(a, require_positive: bool = True) -> np.ndarray:
    """Coerce the shape parameter to a 1-D float array (``N = 0`` gives size 0).

    ``require_positive`` is on for the physical parameter ``a_i = 2 gamma_i /
    kappa_i**2``, which is positive by construction.  It is off only for the
    *shifted* vector ``a + 3 - 4s`` of Eq. (38), where the ``s = 1`` shift is
    ``-1`` and drives entries to zero or below whenever ``a_i <= 1`` -- still a
    perfectly admissible Meijer G parameter, just no longer a reservoir shape.
    """
    arr = np.atleast_1d(np.asarray(a, dtype=float))
    if arr.ndim != 1:
        raise ValueError(f"a must be a scalar or 1-D sequence, got shape {arr.shape}")
    if require_positive and arr.size and np.any(arr <= 0):
        raise ValueError(f"all a_i must be positive, got {arr}")
    return arr


def _omega(a: np.ndarray) -> float:
    """The scalar ``omega = prod a_i``; 1 for the empty hierarchy, Eq. (33)."""
    return float(np.prod(a)) if a.size else 1.0


def kappa_from_a(a, gamma):
    """Invert ``a_i = 2 gamma_i / kappa_i**2`` for the noise amplitudes.

    The pair ``(gamma_i, kappa_i)`` in Eq. (29) enters the stationary law only
    through this combination, so a hierarchy is specified either way.
    """
    a = _as_a(a)
    gamma = np.broadcast_to(np.asarray(gamma, dtype=float), a.shape)
    return np.sqrt(2.0 * gamma / a)


# --------------------------------------------------------------------------
# The conditionals, Eqs. (31) and (32)
# --------------------------------------------------------------------------


def conditional_pdf(T, T_prev, a_i, s):
    """``f(T_i | T_{i-1})``, Eq. (31) for ``s = 1/2`` and Eq. (32) for ``s = 1``.

    Both have mean ``T_{i-1}``: the linear drift of Eq. (29) drives every level
    to the temperature of the one above it, so ``<T_i> = T_0`` for all ``i``.
    """
    T = np.asarray(T, dtype=float)
    if s == 0.5:
        # Gamma(shape a_i, rate a_i/T_prev).
        return stats.gamma.pdf(T, a_i, scale=T_prev / a_i)
    if s == 1:
        # Inverse-gamma(shape a_i + 1, scale a_i * T_prev).
        return stats.invgamma.pdf(T, a_i + 1.0, scale=a_i * T_prev)
    raise ValueError(f"s must be 1/2 or 1, got {s}")


def sample_conditional(T_prev, a_i, s, rng, size=None):
    """Draw ``T_i ~ f(T_i | T_{i-1})``, Eqs. (31)/(32).

    ``size`` defaults to the shape of ``T_prev``, so that an array of parent
    temperatures gets an independent draw each: leaving it to broadcast would
    silently hand the whole array a *single* gamma variate in the ``s = 1``
    branch, where ``T_prev`` enters the scale rather than a numpy parameter.
    """
    T_prev = np.asarray(T_prev, dtype=float)
    if size is None and T_prev.ndim:
        size = T_prev.shape
    if s == 0.5:
        # Gamma(a_i, rate a_i/T_prev).
        return rng.gamma(a_i, 1.0, size=size) * T_prev / a_i
    if s == 1:
        # Inverse-gamma(a_i + 1, scale a_i T_prev) = scale / Gamma(a_i + 1, 1).
        return a_i * T_prev / rng.gamma(a_i + 1.0, 1.0, size=size)
    raise ValueError(f"s must be 1/2 or 1, got {s}")


# --------------------------------------------------------------------------
# The marginal, Eq. (33)  -- corrected branch assignment, see ERRATUM_EQ33
# --------------------------------------------------------------------------


def _f_N_scalar(T, a, s, T0):
    a = np.asarray(a, dtype=float)
    if T <= 0:
        return mp.mpf(0)
    if a.size == 0:
        raise ValueError("f_N is a delta at T0 for N = 0; it has no density")
    w = _omega(a)
    if s == 0.5:
        # Gamma chain: G^{N,0}_{0,N}( w T/T0 | - ; a - 1 ).
        pref = mp.mpf(w) / mp.mpf(T0)
        for x in a:
            pref /= mp.gamma(mp.mpf(x))
        g = meijerg(
            [[], []], [[mp.mpf(x) - 1 for x in a], []], mp.mpf(w) * T / T0,
            underflow_above=UNDERFLOW_Z,
        )
    else:
        # Inverse-gamma chain: G^{0,N}_{N,0}( T/(w T0) | -a - 1 ; - ).
        pref = 1 / (mp.mpf(T0) * mp.mpf(w))
        for x in a:
            pref /= mp.gamma(mp.mpf(x) + 1)
        # This orientation dies at *small* argument: f ~ exp(-omega T0 / T),
        # an essential singularity at T = 0. At the threshold the density is
        # already exp(-100).
        g = meijerg(
            [[-mp.mpf(x) - 1 for x in a], []], [[], []], mp.mpf(T) / (w * T0),
            underflow_below=UNDERFLOW_Z_SMALL,
        )
    return pref * g


def f_N(T, a, s, T0: float = 1.0):
    """The reservoir temperature law ``f_N(T_N)``, Eq. (33), corrected.

    Normalizes to 1 and has mean ``T0`` for every ``N``, ``a`` and both classes
    -- the two facts the sum rule Eq. (27) rests on.
    """
    a = _as_a(a)
    return vectorize(lambda t: _f_N_scalar(t, a, s, T0))(T)


def sample_TN(a, s, rng, T0: float = 1.0, size: int = 1):
    """Ancestral sampling of ``T_N`` down the chain, Eq. (30).

    ``T_0 = T0`` is fixed; each level is drawn from its conditional given the
    level above.  The empirical law of the result is Eq. (33) -- this is the
    Monte-Carlo check that the corrected branch assignment is the right one.
    """
    a = _as_a(a)
    T = np.full(int(size), float(T0))
    for a_i in a:
        T = sample_conditional(T, a_i, s, rng)
    return T


# --------------------------------------------------------------------------
# The dynamics, Eq. (29)
# --------------------------------------------------------------------------


def simulate_sde(
    a,
    s,
    gamma,
    rng,
    T0: float = 1.0,
    dt: float = 1e-3,
    n_steps: int = 200_000,
    burn_in: float = 0.1,
    n_paths: int = 1,
):
    """Euler-Maruyama on the hierarchy ``dT_i = -g_i (T_i - T_{i-1}) dt +
    k_i T_i^s T_{i-1}^{1-s} dW_i``, Eq. (29).

    ``gamma`` is the vector of relaxation rates ``g_i``; the noise amplitudes
    follow from ``a_i = 2 g_i / k_i**2``.  Eq. (29) presumes ``g_i >> g_{i-1}``
    so that each level sees the one above it as frozen -- that separation of
    scales is what makes the stationary law the *product* of conditionals in
    Eq. (30), and it is the caller's job to supply a ``gamma`` that respects it.

    Returns the post-burn-in trajectory of the innermost level, ``T_N``, whose
    stationary histogram is :func:`f_N`.  Convergence is O(sqrt(dt)) and the
    ``s = 1/2`` (CIR) case is only nonnegative in the exact process, so the
    diffusion is evaluated on the positive part -- the usual full-truncation
    scheme.
    """
    a = _as_a(a)
    N = a.size
    if N == 0:
        return np.full((n_paths, 1), float(T0))
    gamma = np.broadcast_to(np.asarray(gamma, dtype=float), (N,)).astype(float)
    kappa = kappa_from_a(a, gamma)
    if np.max(gamma) * dt > 0.1:
        raise ValueError(
            f"dt = {dt} is too coarse for gamma_max = {gamma.max()}; "
            f"need gamma_max * dt << 1"
        )

    T = np.full((n_paths, N), float(T0))
    n_burn = int(burn_in * n_steps)
    out = np.empty((n_paths, n_steps - n_burn), dtype=float)
    sqrt_dt = np.sqrt(dt)
    for step in range(n_steps):
        prev = np.concatenate([np.full((n_paths, 1), float(T0)), T[:, :-1]], axis=1)
        pos = np.maximum(T, 0.0)
        diff = kappa * pos**s * prev ** (1.0 - s)
        dW = rng.standard_normal((n_paths, N)) * sqrt_dt
        T = T - gamma * (T - prev) * dt + diff * dW
        T = np.maximum(T, 0.0)
        if step >= n_burn:
            out[:, step - n_burn] = T[:, -1]
    return out
