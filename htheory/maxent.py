"""Closing the loop: maximize ``S[p]`` and get ``p_i = B(E_i)/Z`` back.

This is the paper's central claim, Section 3.  The chain
``B -> s~(E) -> s(B) -> S[p]`` reconstructs an entropy from a prescribed
statistical factor; the claim is that maximizing that entropy under the
*standard* constraints -- normalization and mean energy, Eq. (8), with no escort
probabilities anywhere -- returns the distribution you started from, Eq. (6).

Nothing here is a new derivation; it is the regression test with the highest
value in the package, because it fails if any single link is wrong: the branch
assignment of Eq. (33), the ``omega``-versus-``a`` bookkeeping of Eq. (38), the
gauge ``alpha = -1`` of Eq. (16), or the ``1/Z`` scaling of Eq. (10).

The gradient is exact and cheap.  Differentiating Eq. (10),

    d/dp_i sum_j Z^{-1} s(Z p_j) = s'(Z p_i) = alpha + beta0 E(Z p_i),

by Eq. (12) -- so the derivative of the entropy *is* the inverse function
``E(B)``, and no numerical differentiation is needed.  At the target
``p_i = B(E_i)/Z`` we have ``E(Z p_i) = E_i``, so the gradient is affine in the
energies and the stationarity condition of Eq. (8) holds exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .background import _as_a
from .statfactor import B
from .entropy import make_interpolants

__all__ = ["target_distribution", "MaxEntResult", "maximize", "stationarity_residual"]


def target_distribution(energies, a, s, beta0: float = 1.0):
    """``p_i = B(E_i)/Z`` with ``Z = sum_i B(E_i)``, Eqs. (6) and (28).

    Returns ``(p, Z, U)`` with ``U = sum_i p_i E_i`` the mean energy that the
    maximization will be constrained to.
    """
    energies = np.asarray(energies, dtype=float)
    b = np.asarray(B(energies, a, s, beta0), dtype=float)
    Z = float(b.sum())
    p = b / Z
    return p, Z, float(np.dot(p, energies))


@dataclass
class MaxEntResult:
    """Outcome of :func:`maximize`."""

    p: np.ndarray
    p_target: np.ndarray
    Z: float
    U: float
    entropy: float
    entropy_target: float
    max_abs_error: float
    success: bool
    message: str = ""
    scipy_result: object = field(default=None, repr=False)


def maximize(
    energies,
    a,
    s,
    beta0: float = 1.0,
    x0=None,
    p_min: float = 1e-14,
    tol: float = 1e-12,
    maxiter: int = 500,
    n_table: int = 800,
):
    """Maximize ``S[p] = sum_i Z^{-1} s(Z p_i)``, Eq. (10), under Eq. (8).

    Constraints are exactly the two standard ones: ``sum_i p_i = 1`` and
    ``sum_i p_i E_i = U``, with ``U`` taken from :func:`target_distribution`.
    ``Z`` is held at the partition function of the target, as Eq. (10) requires
    -- it is a fixed property of the spectrum, not a variable of the problem.

    ``E(B)`` comes from the tabulated interpolant of
    :func:`~htheory.entropy.make_interpolants`; a per-iteration bisection would
    cost thousands of Meijer G evaluations and buy about six digits of nothing.

    SLSQP is run twice, the second time from the first answer: the entropy is
    very flat near its maximum (the gradient is affine in the energies there),
    so a single pass stops early on the objective long before the iterate has
    settled.
    """
    energies = np.asarray(energies, dtype=float)
    a = _as_a(a)
    p_target, Z, U = target_distribution(energies, a, s, beta0)

    # 800 knots, not the module default of 300: the recovered probabilities
    # inherit the table's accuracy directly, and 300 leaves the tail of the
    # power-law class good to only about three digits.
    interp = make_interpolants(
        a, s, beta0, n=n_table, B_min=min(1e-13, 0.01 * Z * p_target.min())
    )
    # Bounds are the *domain* of Eq. (10), not a hint at the answer. The entropy
    # density s(x) is built from E(x) by Eq. (13), and E is defined only for
    # x = B in (0, 1]: energy is bounded below by zero and B(0) = 1, Eq. (7).
    # So Z p_i must stay in (0, 1]. Without the upper bound the optimizer walks
    # into Z p_i > 1, where the tabulated s saturates instead of continuing, and
    # reports a "maximum" above the true one -- entropy for free. The lower
    # bound keeps the same search inside the tabulated range at the other end.
    p_min = max(p_min, interp.B_min / Z)
    p_max = 1.0 / Z

    def entropy(p):
        """``S[p] = sum_i Z^{-1} s(Z p_i)``, Eq. (10), with the linear term dropped."""
        x = np.clip(Z * np.asarray(p, dtype=float), interp.B_min, None)
        return float(np.sum(interp.s_of_B(x)) / Z)

    def neg_grad(p):
        # s'(B) = alpha + beta0 E(B) with alpha = -1, Eqs. (12) and (16).
        x = np.clip(Z * np.asarray(p, dtype=float), interp.B_min, None)
        return -(beta0 * np.asarray(interp.E_of_B(x), dtype=float) - 1.0)

    if x0 is not None:
        x0 = np.clip(np.asarray(x0, dtype=float), p_min, p_max)
    if x0 is None:
        # Deliberately not the answer: a tilted-uniform start that satisfies
        # both constraints without encoding the shape of B.
        x0 = _feasible_start(energies, U)

    constraints = (
        {"type": "eq", "fun": lambda p: np.sum(p) - 1.0, "jac": lambda p: np.ones_like(p)},
        {"type": "eq", "fun": lambda p: np.dot(p, energies) - U, "jac": lambda p: energies},
    )
    res = None
    for _ in range(2):
        res = minimize(
            lambda p: -entropy(p),
            x0,
            jac=neg_grad,
            bounds=[(p_min, p_max)] * energies.size,
            constraints=constraints,
            method="SLSQP",
            options={"maxiter": maxiter, "ftol": tol},
        )
        x0 = np.clip(res.x, p_min, p_max)
        x0 = x0 / x0.sum()
    p = np.clip(res.x, p_min, p_max)
    p = p / p.sum()
    return MaxEntResult(
        p=p,
        p_target=p_target,
        Z=Z,
        U=U,
        entropy=entropy(p),
        entropy_target=entropy(p_target),
        max_abs_error=float(np.max(np.abs(p - p_target))),
        success=bool(res.success),
        message=str(res.message),
        scipy_result=res,
    )


def stationarity_residual(p, energies, a, s, beta0: float = 1.0, Z: float | None = None):
    """Spread of ``s'(Z p_i) - beta0 E_i`` across states.

    Eq. (12) says this is the constant ``alpha`` at the maximizer, so the spread
    is zero there.  Reported as ``max - min``: a direct, optimizer-free check of
    the paper's claim, and the one that localizes a broken link fastest.
    """
    p = np.asarray(p, dtype=float)
    energies = np.asarray(energies, dtype=float)
    if Z is None:
        Z = float(np.sum(np.asarray(B(energies, a, s, beta0), dtype=float)))
    from .entropy import E_of_B as exact_E_of_B

    s_prime = beta0 * np.asarray(exact_E_of_B(Z * p, a, s, beta0), dtype=float) - 1.0
    resid = s_prime - beta0 * energies
    return float(resid.max() - resid.min())


def _feasible_start(energies, U):
    """A point on the simplex with the right mean energy, and nothing else.

    Uniform weights tilted linearly toward ``U``; if that leaves the simplex the
    tilt is scaled back until it does not.
    """
    n = energies.size
    p = np.full(n, 1.0 / n)
    centred = energies - energies.mean()
    denom = float(np.dot(centred, centred)) / n
    if denom <= 0:
        return p
    lam = (U - energies.mean()) / denom
    for _ in range(60):
        trial = p + lam * centred / n
        if trial.min() > 0:
            return trial / trial.sum()
        lam *= 0.5
    return p
