"""Identifiability of the cascade depth ``N`` and shapes ``a`` from ``T_N`` alone.

**Theorem 3** (``docs/theorems.md``): the cascade
is additive in log-temperature, and its log-cumulants of order ``r >= 2``
identify the multiset ``{a_i}`` -- hence ``N`` -- without observing any
intermediate level.

Derivation (recorded here because it is the load-bearing fact for everything
below, and does not appear elsewhere in the package). Write the increment at
level ``i`` as ``D_i = log T_i - log T_{i-1}``.

* ``s = 1/2``: :func:`background.sample_conditional` draws ``T_i = T_{i-1} *
  Y_i / a_i`` with ``Y_i ~ Gamma(a_i, 1)``, so ``D_i = log Y_i - log a_i`` --
  a function of ``Y_i`` alone, independent of ``T_{i-1}``.
* ``s = 1``: it draws ``T_i = a_i T_{i-1} / Y_i`` with ``Y_i ~ Gamma(a_i + 1,
  1)``, so ``D_i = log a_i - log Y_i``, again independent of ``T_{i-1}``.

In both branches the ``D_i`` are therefore mutually independent (each drawn
from a fresh ``Y_i``) and independent of the deterministic ``T_0``. Cumulants
add over independent summands, so ``kappa_r(log T_N) = sum_i kappa_r(D_i)``
for every order, and a constant shift (``log T_0``, ``-log a_i`` or
``log a_i``) contributes only to ``kappa_1``.

For ``Y ~ Gamma(alpha, rate=beta)``, the cumulant generating function of
``log Y`` is ``K(t) = log E[Y^t] = log Gamma(alpha + t) - log Gamma(alpha) - t
log beta``, so ``K^{(r)}(t) = psi^{(r-1)}(alpha + t)`` for ``r >= 2`` (the
``-t log beta`` term dies after one derivative) -- the rate drops out
entirely. Hence, for ``r >= 2``:

    s = 1/2:  kappa_r(D_i) = psi^{(r-1)}(a_i)
    s = 1:    kappa_r(D_i) = (-1)^r psi^{(r-1)}(a_i + 1)

and summing over levels gives :func:`log_cumulants_theory`. This is the same
identity behind the method of log-cumulants in SAR image analysis (Nicolas
2002): population cumulants of order >= 2 of a log-amplitude are scale-free,
so they see through the unknown reservoir scale ``T_0`` straight to the shape
parameters -- and, because they add over levels, straight to the depth ``N``
that the original paper's Eq. (33)/Fig. 1-2 argued was invisible.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special, stats

from .background import sample_TN

__all__ = [
    "log_cumulants_theory",
    "empirical_log_cumulants",
    "fit_homogeneous",
    "DepthTestResult",
    "test_depth",
    "power_curve",
]


def _shift(s) -> float:
    if s == 0.5:
        return 0.0
    if s == 1:
        return 1.0
    raise ValueError(f"s must be 1/2 or 1, got {s}")


def log_cumulants_theory(a, s, r_max: int = 4) -> np.ndarray:
    """Theoretical ``kappa_r(log T_N)`` for ``r = 2..r_max``, Theorem 3.

    Independent of ``T_0`` by construction (see module docstring); relies only
    on the shape vector ``a`` and the branch ``s``. ``a`` may be empty
    (``N = 0``), for which all cumulants are 0 -- ``T_N = T_0`` deterministically.
    """
    a = np.atleast_1d(np.asarray(a, dtype=float))
    off = _shift(s)
    out = np.zeros(r_max - 1)
    for idx, r in enumerate(range(2, r_max + 1)):
        vals = special.polygamma(r - 1, a + off) if a.size else np.array([])
        total = float(np.sum(vals))
        out[idx] = total if s == 0.5 else ((-1) ** r) * total
    return out


def empirical_log_cumulants(T, r_max: int = 4) -> np.ndarray:
    """Unbiased sample cumulants ``k_r`` of ``log T`` for ``r = 2..r_max``.

    Uses :func:`scipy.stats.kstat`, the standard unbiased k-statistic
    estimator, so this is directly comparable to :func:`log_cumulants_theory`
    without a finite-sample bias correction.
    """
    x = np.log(np.asarray(T, dtype=float))
    return np.array([stats.kstat(x, r) for r in range(2, r_max + 1)])


def _ratio_theory(a: float, s) -> float:
    """``kappa_3 / kappa_2`` for a single homogeneous level, as a function of ``a``.

    Delegates to :func:`log_cumulants_theory` rather than re-deriving the
    branch sign here, since ``kappa_3 / kappa_2`` for ``N`` homogeneous levels
    equals the single-level ratio regardless of ``N`` (the factor of ``N``
    cancels), and getting the ``s = 1`` sign flip of ``kappa_3`` right in two
    places independently is exactly how this went wrong once already.

    Monotone in ``a`` (both ``psi'`` and ``|psi''|`` are monotone decreasing on
    ``(0, inf)`` and their ratio inherits that), which is what makes
    :func:`fit_homogeneous` well posed as a scalar root find.
    """
    k2, k3 = log_cumulants_theory([a], s, r_max=3)
    return float(k3 / k2)


@dataclass(frozen=True)
class HomogeneousFit:
    """Point estimate of ``(N, a)`` under the homogeneous model ``a_i = a``."""

    a: float
    N: float
    N_round: int


def fit_homogeneous(kappa2: float, kappa3: float, s, a_bracket=(1e-4, 1e4)) -> HomogeneousFit:
    """Method-of-log-cumulants estimate of ``(a, N)`` assuming all ``a_i = a``.

    Two equations, two unknowns: the ratio ``kappa3 / kappa2`` depends on
    ``a`` alone (the ``N`` factor cancels), so ``a`` is recovered by a 1-D
    root find, then ``N = kappa2 / kappa_2^{(1)}(a)`` from the scale. ``N`` is
    returned continuous (a real-valued generalized depth) and separately
    rounded to the nearest non-negative integer; the article's model-selection
    machinery (:func:`test_depth`) works with a fixed candidate ``N`` and does
    not need the rounding.
    """
    target = kappa3 / kappa2
    lo, hi = a_bracket
    f = lambda a: _ratio_theory(a, s) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError(
            f"kappa3/kappa2 = {target} not bracketed by a in {a_bracket}; "
            "widen a_bracket or check the sample size"
        )
    a_hat = optimize.brentq(f, lo, hi, xtol=1e-10)
    off = _shift(s)
    N_hat = kappa2 / special.polygamma(1, a_hat + off)
    return HomogeneousFit(a=a_hat, N=N_hat, N_round=max(0, round(N_hat)))


@dataclass(frozen=True)
class DepthTestResult:
    """Parametric-bootstrap test of ``H0: depth = N0`` against a deeper cascade.

    ``N`` is a positive integer sitting at the boundary of its own domain, so
    the usual chi-square asymptotics for a likelihood-ratio test do not apply
    (Self & Liang 1987) -- hence the bootstrap rather than a closed-form
    p-value, per the article outline (§06, item 4).
    """

    statistic: float
    p_value: float
    a_hat: float
    n_boot: int


def test_depth(
    T,
    N0: int,
    s,
    rng,
    n_boot: int = 2000,
    T0: float = 1.0,
) -> DepthTestResult:
    """Test ``H0: N = N0`` (homogeneous ``a``) against the sample ``T``.

    Fits ``a`` under ``H0`` from ``kappa2`` alone (``a_hat`` solving
    ``N0 * psi'(a) = kappa2_hat``, the branch-appropriate scale identity), so
    the test statistic -- the discrepancy between the *observed* ``kappa3``
    and the one predicted at ``(N0, a_hat)`` -- is not tautologically zero.
    Its null distribution is approximated by simulating ``n_boot`` cascades of
    depth ``N0`` and shape ``a_hat`` and refitting ``a_hat`` (hence the
    predicted ``kappa3``) on each, exactly mirroring what happens on the
    observed sample; the p-value is the bootstrap tail probability of a
    discrepancy at least as large as observed.
    """
    T = np.asarray(T, dtype=float)
    n = T.size
    off = _shift(s)
    kappa2_obs, kappa3_obs = empirical_log_cumulants(T, r_max=3)

    a_lo, a_hi = 1e-6, 1e6
    psi_hi, psi_lo = special.polygamma(1, a_hi + off), special.polygamma(1, a_lo + off)

    def discrepancy(k2, k3):
        # psi'(a) is strictly decreasing, so N0*psi'(a) covers (N0*psi_hi,
        # N0*psi_lo) on the bracket; clip a noisy k2 (small n, or a bootstrap
        # draw landing near the boundary) into that range rather than letting
        # brentq raise, since the discrepancy is still well defined there --
        # it just means the null a_hat saturates at the edge of the bracket.
        target = np.clip(k2 / N0, psi_hi * 1.0000001, psi_lo * 0.9999999)
        a_hat = optimize.brentq(
            lambda a: special.polygamma(1, a + off) - target, a_lo, a_hi, xtol=1e-10
        )
        k3_pred = N0 * special.polygamma(2, a_hat + off)
        k3_pred = k3_pred if s == 0.5 else -k3_pred
        return float(k3 - k3_pred), a_hat

    stat_obs, a_hat_obs = discrepancy(kappa2_obs, kappa3_obs)

    boot_stats = np.empty(n_boot)
    a_vec_null = np.full(N0, a_hat_obs)
    for b in range(n_boot):
        T_boot = sample_TN(a_vec_null, s, rng, T0=T0, size=n)
        k2_b, k3_b = empirical_log_cumulants(T_boot, r_max=3)
        boot_stats[b], _ = discrepancy(k2_b, k3_b)

    p_value = float(np.mean(np.abs(boot_stats) >= abs(stat_obs)))
    return DepthTestResult(statistic=stat_obs, p_value=p_value, a_hat=a_hat_obs, n_boot=n_boot)


def power_curve(
    a_true,
    s,
    N0: int,
    rng,
    sample_sizes=(100, 300, 1000, 5000),
    n_reps: int = 200,
    n_boot: int = 400,
    alpha: float = 0.05,
    T0: float = 1.0,
):
    """Simulated power of :func:`test_depth` at ``H0: N = N0`` against ``a_true``.

    ``a_true`` is the true (homogeneous or not) shape vector generating the
    data -- its length is the true depth. Returns a dict ``{n: power}``. This
    is the numerical study of the claim, made in the original paper, that ``N``
    cannot be inferred from the foreground distribution: reject rate under a *true*
    deeper cascade, at a shallow null.
    """
    powers = {}
    for n in sample_sizes:
        rejects = 0
        for _ in range(n_reps):
            T = sample_TN(a_true, s, rng, T0=T0, size=n)
            result = test_depth(T, N0, s, rng, n_boot=n_boot, T0=T0)
            rejects += result.p_value < alpha
        powers[n] = rejects / n_reps
    return powers
