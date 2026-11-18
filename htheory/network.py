"""Hierarchical random networks: the mixed-Poisson application.

The network reading of the cascade takes energy as vertex degree and the
Boltzmann factor as a Poisson kernel. Concretely: each vertex
``v`` belongs to a nested sequence of communities, and its propensity to
connect is the product of per-level contributions, ``W_v = prod_i
X_{i,c_i(v)}`` -- exactly the multiplicative cascade of :func:`background.sample_TN`,
relabelled. Degrees are then mixed Poisson, ``D_v | W_v ~ Poisson(lambda W_v)``
(Chung-Lu / Norros-Reittu), so this module reuses :func:`background.sample_TN`
verbatim for ``W`` rather than reimplementing the cascade.

The identifiability chain (article draft §05, "a cadeia de identificabilidade
fecha"):

1. Mixed-Poisson laws are identifiable from the degree distribution alone
   (Feller 1943) -- the law of ``W`` is recoverable from ``D``.
2. By Theorem 3 (:mod:`htheory.identifiability`), ``(N, {a_i})`` is
   identifiable from the law of ``W``.
3. Hence the cascade depth is estimable from the degree sequence alone,
   without observing which vertex belongs to which community.

The estimator is the factorial-moment analogue of the log-cumulant one:
Poisson factorial moments satisfy ``E[D^{(k)}] = lambda^k E[W^k]`` exactly (a
mixed-Poisson identity, not an approximation), and ``E[W^k]`` has the closed
form :func:`log_factorial_moment_theory` -- linear in ``N`` for fixed ``a``,
so the same ratio trick as :func:`identifiability.fit_homogeneous` pins down
``a`` first and ``N`` from the scale.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special

from .background import sample_TN

__all__ = [
    "sample_degrees",
    "log_factorial_moment_theory",
    "empirical_log_factorial_moments",
    "HomogeneousNetworkFit",
    "fit_homogeneous_network",
    "DegreeDepthTestResult",
    "test_depth_degrees",
]


def sample_degrees(a, s, lam: float, rng, size: int, T0: float = 1.0):
    """Draw a mixed-Poisson degree sequence, ``D_v | W_v ~ Poisson(lambda W_v)``.

    ``W`` is drawn from :func:`background.sample_TN` with ``T0 = 1``, so
    ``E[W] = 1`` and ``E[D] = lambda`` exactly -- the estimator below uses
    this to pin down ``lambda`` from the sample mean before touching the
    higher factorial moments.
    """
    W = sample_TN(a, s, rng, T0=T0, size=size)
    return rng.poisson(lam * W)


def _level_log_moment(a: np.ndarray, s, k: int) -> np.ndarray:
    """Per-level ``log E[X_i^k]`` for the single multiplicative factor at each
    level (the ``s = 1/2``/``s = 1`` branch of :func:`background.sample_conditional`,
    with the parent temperature normalized away since each factor has mean 1)."""
    if s == 0.5:
        return special.gammaln(a + k) - special.gammaln(a) - k * np.log(a)
    if s == 1:
        if np.any(a + 1 - k <= 0):
            raise ValueError(
                f"k = {k} requires a_i > {k - 1} for all levels in the s = 1 "
                f"branch (E[W^k] diverges otherwise -- a heavy-tail effect, "
                f"not a bug); got a = {a}"
            )
        return k * np.log(a) + special.gammaln(a + 1 - k) - special.gammaln(a + 1)
    raise ValueError(f"s must be 1/2 or 1, got {s}")


def log_factorial_moment_theory(a, s, k: int) -> float:
    """``log E[W^k]``, summed over independent levels.

    Homogeneous case (``a_i = a`` for all ``i``) reduces to ``N * [log
    Gamma(a+k) - log Gamma(a) - k log a]`` for ``s = 1/2``, the formula in
    the article draft §05 -- linear in ``N``, which is what makes the depth
    a slope rather than a nonlinear fit.
    """
    a = np.atleast_1d(np.asarray(a, dtype=float))
    return float(np.sum(_level_log_moment(a, s, k))) if a.size else 0.0


def empirical_log_factorial_moments(D, k_max: int = 3):
    """``(lambda_hat, [L_2, .., L_k_max])`` with ``L_k = log Ehat[W^k]``.

    ``lambda_hat`` is the sample mean of ``D`` (exact since ``E[W] = 1``).
    ``Ehat[W^k]`` is the empirical ``k``-th factorial moment of ``D`` divided
    by ``lambda_hat^k`` -- ``E[D^{(k)}] = lambda^k E[W^k]`` is an exact
    mixed-Poisson identity, not an approximation. The falling-factorial
    product is exactly 0 whenever ``D < k`` (one factor is ``D - D``), so no
    special-casing of small degrees is needed.
    """
    D = np.asarray(D, dtype=float)
    lam_hat = float(D.mean())
    out = np.empty(k_max - 1)
    for idx, k in enumerate(range(2, k_max + 1)):
        fac = np.ones_like(D)
        for j in range(k):
            fac = fac * (D - j)
        out[idx] = np.log(fac.mean() / lam_hat**k)
    return lam_hat, out


@dataclass(frozen=True)
class HomogeneousNetworkFit:
    """Point estimate of ``(N, a)`` under the homogeneous model ``a_i = a``,
    from the degree sequence alone."""

    a: float
    N: float
    N_round: int
    lam: float


def fit_homogeneous_network(
    D, s, k_max: int = 3, a_bracket=None
) -> HomogeneousNetworkFit:
    """Recover ``(N, a, lambda)`` from a degree sequence, homogeneous model.

    Mirrors :func:`identifiability.fit_homogeneous`: the ratio ``L_k /
    L_2`` depends only on ``a`` (the common factor of ``N`` cancels), so a
    1-D root find on ``a`` comes first, then ``N = L_2 / [single-level log
    E[W^2]](a)`` from the scale. Uses ``k_max = 3`` by default (the ratio
    ``L_3 / L_2``); larger ``k_max`` would let a least-squares fit average
    over more moments; but note the ``s = 1`` branch requires ``a > k - 1``
    for ``E[W^k]`` to exist, so a small true ``a`` bounds how large ``k_max``
    can usefully go there.
    """
    lam_hat, L = empirical_log_factorial_moments(D, k_max=max(k_max, 3))
    L2, L3 = L[0], L[1]
    target = L3 / L2

    def ratio(a):
        c2 = _level_log_moment(np.array([a]), s, 2)[0]
        c3 = _level_log_moment(np.array([a]), s, 3)[0]
        return c3 / c2

    if a_bracket is None:
        # ratio() evaluates E[W^3], which needs a > 2 in the s = 1 branch
        # (see _level_log_moment); s = 1/2 has no such constraint.
        a_bracket = (2.0 + 1e-6, 1e3) if s == 1 else (1e-3, 1e3)
    lo, hi = a_bracket
    f = lambda a: ratio(a) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError(
            f"L3/L2 = {target} not bracketed by a in {a_bracket}; "
            "widen a_bracket or check the sample size"
        )
    a_hat = optimize.brentq(f, lo, hi, xtol=1e-10)
    c2_hat = _level_log_moment(np.array([a_hat]), s, 2)[0]
    N_hat = L2 / c2_hat
    return HomogeneousNetworkFit(
        a=a_hat, N=N_hat, N_round=max(0, round(N_hat)), lam=lam_hat
    )


@dataclass(frozen=True)
class DegreeDepthTestResult:
    """Parametric-bootstrap test of ``H0: depth = N0`` against a degree
    sequence, homogeneous model -- the factorial-moment analogue of
    :func:`identifiability.test_depth`.

    :func:`fit_homogeneous_network`'s point estimate ``N_round`` had no
    formal p-value attached to it; this closes that gap the same way
    :func:`identifiability.test_depth` does for the propensity sequence ``T``
    -- ``N`` sits at the boundary of its domain (a non-negative integer), so
    the usual chi-square asymptotics for a likelihood-ratio test do not apply
    (Self & Liang 1987), hence bootstrap rather than closed form.
    """

    statistic: float
    p_value: float
    a_hat: float
    lam_hat: float
    n_boot: int


def test_depth_degrees(
    D, N0: int, s, rng, n_boot: int = 2000, k_max: int = 3
) -> DegreeDepthTestResult:
    """Test ``H0: N = N0`` (homogeneous ``a``) against the degree sequence ``D``.

    Mirrors :func:`identifiability.test_depth` with ``L_2``/``L_3`` (the log
    factorial-moment ratios of :func:`empirical_log_factorial_moments`) in
    place of the log-cumulants ``kappa2``/``kappa3``: fits ``a`` under
    ``H0`` from ``L_2`` alone (``N0 * [single-level log E[W^2]](a) =
    L2_hat``), so the test statistic -- observed ``L_3`` minus the one
    predicted at ``(N0, a_hat)`` -- is not tautologically zero. Its null
    distribution is a parametric bootstrap: simulate ``n_boot`` degree
    sequences of depth ``N0``, shape ``a_hat`` and rate ``lam_hat`` (all held
    at their fitted values) via :func:`sample_degrees`, and refit ``a_hat``
    (hence the predicted ``L_3``) on each exactly as done on the real sample.

    The ``s = 1`` branch requires ``a > k - 1`` for ``E[W^k]`` to exist (see
    :func:`_level_log_moment`); the search bracket is restricted accordingly,
    matching :func:`fit_homogeneous_network`'s own default. A degree sequence
    whose ``L_3/L_2`` ratio cannot be matched by any ``a`` in that bracket
    (as happens for the power-law branch on the email-Eu-core data, see
    ``results/analyze_email_eu_core.json``) is a genuine model-selection signal, not something
    this test papers over -- it will raise the same ``ValueError`` as
    :func:`fit_homogeneous_network` rather than returning a fabricated p-value.
    """
    D = np.asarray(D, dtype=float)
    n = D.size
    lam_hat, L = empirical_log_factorial_moments(D, k_max=max(k_max, 3))
    L2_obs, L3_obs = L[0], L[1]

    a_lo, a_hi = (2.0 + 1e-6, 1e3) if s == 1 else (1e-6, 1e3)
    f2 = lambda a: _level_log_moment(np.array([a]), s, 2)[0]
    f3 = lambda a: _level_log_moment(np.array([a]), s, 3)[0]
    # f2 is strictly decreasing in a (E[W^2] shrinks toward 1 as a -> inf), so
    # val_hi_a < val_lo_a; clip a noisy L2 into that range rather than letting
    # brentq raise, mirroring identifiability.test_depth's own clipping.
    val_at_lo, val_at_hi = f2(a_lo), f2(a_hi)

    def discrepancy(L2, L3):
        target = np.clip(L2 / N0, val_at_hi * 1.0000001, val_at_lo * 0.9999999)
        a_hat = optimize.brentq(lambda a: f2(a) - target, a_lo, a_hi, xtol=1e-10)
        L3_pred = N0 * f3(a_hat)
        return float(L3 - L3_pred), a_hat

    stat_obs, a_hat_obs = discrepancy(L2_obs, L3_obs)

    boot_stats = np.empty(n_boot)
    a_vec_null = np.full(N0, a_hat_obs)
    for b in range(n_boot):
        D_boot = sample_degrees(a_vec_null, s, lam_hat, rng, size=n)
        lam_b, L_b = empirical_log_factorial_moments(D_boot, k_max=max(k_max, 3))
        boot_stats[b], _ = discrepancy(L_b[0], L_b[1])

    p_value = float(np.mean(np.abs(boot_stats) >= abs(stat_obs)))
    return DegreeDepthTestResult(
        statistic=stat_obs, p_value=p_value, a_hat=a_hat_obs, lam_hat=lam_hat, n_boot=n_boot
    )
