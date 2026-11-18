"""Cascade-depth identifiability test on real turbulence intermittency data.

Unlike the two network datasets (`analyze_email_eu_core.py`,
`analyze_openflights.py`), this is H-theory's *native* domain, not a
repurposing: the paper's own introduction motivates the hierarchy by
"statistical fluctuations in systems with coupled spatio-temporal scales,
akin to those observed in fully developed turbulence," and the cascade math
(Kolmogorov 1962's refined similarity hypothesis; Mandelbrot 1974's
multiplicative cascade for intermittency) is literally turbulence theory.

Data: hot-wire anemometry velocity time series, wind-tunnel wake behind a
50mm cylinder at Re=5.17e4, 1000 Hz, 2**17 samples (~131s). Source: "Turbulent
Flow Past Cylinders" (Mendeley Data, doi via
https://data.mendeley.com/datasets/xy8ggx9hn6/1), CC BY 4.0, no registration
required -- downloaded directly via Mendeley's public-files URL pattern.
Local copy: data/turbulence-mendeley-cylinder50mm-16ms.txt (two columns:
time in seconds, velocity in m/s).

Method: at a chosen lag tau (in samples), form the velocity increment
dv_tau(t) = v(t+tau) - v(t) and use its square (a local energy/dissipation
proxy, in the spirit of Kolmogorov's refined similarity hypothesis) as the
positive-valued sequence fed to `identifiability.fit_homogeneous`/
`test_depth` -- exactly the role `T_N` plays in the synthetic cascade. This
is a single-scale cut, not a full multifractal spectrum; the loop over
several tau below is a sanity check for scale-dependence, not a substitute
for one.

Serial dependence -- read this before quoting any p-value from here
--------------------------------------------------------------------
Squared increments of a turbulence signal are serially dependent (overlapping
increments at lag tau, plus the long memory of the dissipation field). The
depth test's null distribution comes from a parametric bootstrap that draws
*i.i.d.* samples of the same nominal size n; the data are not i.i.d., so the
full-sample p-values are anti-conservative (the effective sample size is
smaller than n). Three things are reported to make that visible:

1. the integral time scale of the squared increments at each lag (the
   autocorrelation function integrated, trapezoid rule, up to its first zero
   crossing; in samples and in ms);
2. the same depth test on a *thinned* series that keeps every k-th squared
   increment, k = max(lag, ceil(integral time scale in samples)) -- thinning
   reduces, but does not remove, the dependence;
3. 95% percentile intervals for the continuous estimates (N_hat, a_hat) from
   a moving-block bootstrap with block length max(lag, 5 x integral time
   scale) samples (for the thinned series the same block expressed in thinned
   steps), which respects the dependence the test's own bootstrap ignores.

Squared increments exactly equal to 0 are dropped (log 0); the count is
recorded.

Run from the repo root:

    PYTHONPATH=. python3 scripts/analyze_turbulence.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from htheory import identifiability  # noqa: E402
from scripts import common, data_loaders  # noqa: E402

SEED = 20260826
SYSTEM = "Turbulence"
LAGS = (2, 5, 10, 20, 50, 100)

DEPENDENCE_NOTE = (
    "The depth test's parametric bootstrap draws i.i.d. samples of size n; the "
    "squared increments are serially dependent, so full-sample p-values are "
    "anti-conservative. Thinned-series p-values and moving-block-bootstrap "
    "confidence intervals are reported as sensitivity checks."
)


def integral_time_scale(x):
    """Integral of the autocorrelation function of ``x`` (trapezoid rule, unit
    sample spacing) from lag 0 to its first zero crossing. Returns
    ``(T_int in samples, lag of the first zero crossing)``."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = x.size
    nfft = 1 << int(math.ceil(math.log2(2 * n)))
    f = np.fft.rfft(x, nfft)
    acov = np.fft.irfft(f * np.conj(f), nfft)[:n] / (n - np.arange(n))
    acf = acov / acov[0]
    neg = np.nonzero(acf <= 0)[0]
    j0 = int(neg[0]) if neg.size else n - 1
    rho = acf[: j0 + 1].copy()
    if neg.size:
        # linear interpolation of the crossing point between j0 - 1 and j0
        frac = rho[j0 - 1] / (rho[j0 - 1] - rho[j0])
        t_int = float(np.sum(0.5 * (rho[:-2] + rho[1:-1])) + 0.5 * rho[j0 - 1] * frac)
    else:
        t_int = float(np.sum(0.5 * (rho[:-1] + rho[1:])))
    return t_int, j0


def analyze_series(eps, lag, dt, variant, block_length, rng_test, rng_ci, **extra):
    """``eps``: squared increments in time order (zeros still included)."""
    pos = eps[eps > 0]
    k2, k3 = identifiability.empirical_log_cumulants(pos, r_max=3)
    mean, var = float(pos.mean()), float(pos.var())
    name = f"lag {lag * dt * 1000:g} ms ({variant})"
    rec = {"system": SYSTEM, "series": name, "variant": variant, "lag_samples": lag,
           "lag_ms": lag * dt * 1000, "n": int(pos.size),
           "n_zero_dropped": int(eps.size - pos.size), "mean": mean, "variance": var,
           "vmr": var / mean, "kappa2": float(k2), "kappa3": float(k3),
           "block_length": int(block_length)}
    rec.update(extra)
    rec["classes"] = {
        common.CLASS_KEY[s]: common.analyze_positive_class(
            pos, s, rng_test, rng_ci, block_length=block_length, x_for_ci=eps)
        for s in common.CLASSES
    }
    rec["selected_class"] = common.select_class(rec["classes"])
    common.print_series(rec)
    return rec


def main():
    t, v, dt = data_loaders.load_velocity()
    print(f"Loaded {v.size} samples, dt={dt*1000:.3f} ms, "
          f"mean velocity={v.mean():.3f} m/s, std={v.std():.3f} m/s")
    rng_test = np.random.default_rng(SEED)
    rng_ci = np.random.default_rng([SEED, 1])
    series, scales = [], []
    for lag in LAGS:
        eps = (v[lag:] - v[:-lag]) ** 2
        t_int, j0 = integral_time_scale(eps)
        block = int(max(lag, math.ceil(5 * t_int)))
        k = int(max(lag, math.ceil(t_int)))
        block_thinned = int(max(1, math.ceil(block / k)))
        scale = {"lag_samples": lag, "lag_ms": lag * dt * 1000,
                 "integral_time_scale_samples": t_int,
                 "integral_time_scale_ms": t_int * dt * 1000,
                 "acf_first_zero_crossing_samples": j0,
                 "block_length_samples": block, "thinning_factor": k,
                 "block_length_thinned_steps": block_thinned,
                 "n_zero_squared_increments": int(np.sum(eps == 0))}
        scales.append(scale)
        print(f"\n== lag {lag} samples ({lag*dt*1000:g} ms): integral time scale "
              f"{t_int:.2f} samples, ACF first zero at {j0}, block length {block}, "
              f"thinning k={k} ==")
        common_extra = {"integral_time_scale_samples": t_int}
        series.append(analyze_series(eps, lag, dt, "full", block, rng_test, rng_ci,
                                     thinning_factor=1,
                                     **common_extra))
        series.append(analyze_series(eps[::k], lag, dt, "thinned", block_thinned,
                                     rng_test, rng_ci, thinning_factor=k,
                                     **common_extra))
    common.write_json("analyze_turbulence", {
        "metadata": common.metadata(
            "scripts/analyze_turbulence.py", SEED,
            rng_streams={"depth_test": "default_rng(SEED)",
                         "bootstrap_ci": "default_rng([SEED, 1])"},
            n_boot_test=common.N_BOOT_TEST, n_boot_ci=common.N_BOOT_CI,
            sampling_interval_ms=dt * 1000, n_samples=int(v.size),
            note=common.SELECTION_NOTE, dependence_note=DEPENDENCE_NOTE),
        "time_scales": scales,
        "series": series,
    })


if __name__ == "__main__":
    main()
