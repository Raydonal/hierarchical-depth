"""Shared machinery for the analysis scripts: fitting both classes, the depth
test, nonparametric-bootstrap confidence intervals, and JSON output with
provenance metadata. No numerics of ``htheory`` are reimplemented here.

A caveat every consumer of these results must carry along: the depth test is
a test of ``H0: N = N_round`` where ``N_round`` is the rounded point estimate
*from the same data*. Its p-value is therefore a goodness-of-fit diagnostic
at the selected depth, not a confirmatory test. The bootstrap percentile
interval for the continuous estimate ``N_hat`` is reported next to it for
that reason.
"""

from __future__ import annotations

import json
import math
import platform
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from htheory import identifiability, network  # noqa: E402

RESULTS = REPO_ROOT / "results"
CLASSES = (0.5, 1.0)
CLASS_KEY = {0.5: "0.5", 1.0: "1"}
CLASS_NAME = {0.5: "gamma (stretched-exponential)", 1.0: "inverse-gamma (power-law)"}
N_BOOT_TEST = 1000     # parametric bootstrap of every empirical depth test
N_BOOT_CI = 1000       # nonparametric bootstrap for the confidence intervals

SELECTION_NOTE = (
    "The depth test is a test of H0: N = N_round with N_round selected from "
    "the same data; read its p-value as a goodness-of-fit diagnostic at the "
    "selected depth, not as a confirmatory test."
)


# --------------------------------------------------------------------------
# Provenance / JSON
# --------------------------------------------------------------------------


def metadata(script: str, seed, **extra) -> dict:
    import matplotlib
    import mpmath
    import scipy

    meta = {
        "script": script,
        "rng_seed": seed,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "mpmath": mpmath.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    meta.update(extra)
    return meta


def _clean(obj):
    """JSON-safe copy: numpy scalars -> Python, non-finite floats -> None."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        obj = float(obj)
        return obj if math.isfinite(obj) else None
    return obj


def write_json(name: str, payload: dict) -> Path:
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / f"{name}.json"
    path.write_text(json.dumps(_clean(payload), indent=2) + "\n")
    print(f"wrote {path.relative_to(REPO_ROOT)}")
    return path


def read_json(name: str) -> dict:
    return json.loads((RESULTS / f"{name}.json").read_text())


# --------------------------------------------------------------------------
# Bootstrap confidence intervals
# --------------------------------------------------------------------------


def _percentile_ci(values, n_boot, method, extra=None):
    values = np.asarray(values, dtype=float).reshape(-1, 2)
    n_ok = values.shape[0]
    out = {
        "ci_method": method,
        "ci_n_boot": n_boot,
        "ci_failure_fraction": 1.0 - n_ok / n_boot,
        "ci_level": 0.95,
        "ci_a_hat": None,
        "ci_N_hat": None,
    }
    if n_ok >= 20:
        out["ci_a_hat"] = np.percentile(values[:, 0], [2.5, 97.5]).tolist()
        out["ci_N_hat"] = np.percentile(values[:, 1], [2.5, 97.5]).tolist()
    if extra:
        out.update(extra)
    return out


def _safe_fit(fit_fn):
    """Run a fit; return (a, N) or None if it raises or is not finite/positive."""
    try:
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            fit = fit_fn()
    except (ValueError, FloatingPointError, ZeroDivisionError):
        return None
    if not (math.isfinite(fit.a) and math.isfinite(fit.N)):
        return None
    return fit.a, fit.N


def bootstrap_ci_degrees(D, s, rng, n_boot=N_BOOT_CI):
    """Resample vertices with replacement, refit, 95% percentile interval."""
    D = np.asarray(D, dtype=float)
    n = D.size
    vals = []
    for _ in range(n_boot):
        Db = D[rng.integers(0, n, size=n)]
        res = _safe_fit(lambda: network.fit_homogeneous_network(Db, s))
        if res is not None:
            vals.append(res)
    return _percentile_ci(vals, n_boot, "nonparametric bootstrap, vertices resampled "
                                        "with replacement, percentile interval")


def bootstrap_ci_positive(x, s, rng, block_length=1, n_boot=N_BOOT_CI):
    """Moving-block bootstrap (``block_length = 1`` is the i.i.d. bootstrap) of
    the log-cumulant fit on a positive series ``x`` (zeros dropped after
    resampling, as on the original series)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    L = int(max(1, min(block_length, n)))
    n_blocks = int(math.ceil(n / L))
    offsets = np.arange(L)
    vals = []
    for _ in range(n_boot):
        starts = rng.integers(0, n - L + 1, size=n_blocks)
        idx = (starts[:, None] + offsets[None, :]).ravel()[:n]
        xb = x[idx]
        xb = xb[xb > 0]

        def fit():
            k2, k3 = identifiability.empirical_log_cumulants(xb, r_max=3)
            return identifiability.fit_homogeneous(k2, k3, s)

        res = _safe_fit(fit)
        if res is not None:
            vals.append(res)
    method = ("moving-block bootstrap, percentile interval" if L > 1 else
              "nonparametric i.i.d. bootstrap, percentile interval")
    return _percentile_ci(vals, n_boot, method, {"ci_block_length": L})


# --------------------------------------------------------------------------
# Per-class analysis of one series
# --------------------------------------------------------------------------


def _failed(reason: str) -> dict:
    return {"status": "failed", "reason": reason, "a_hat": None, "N_hat": None,
            "N_round": None, "lambda_hat": None, "test_status": "not run",
            "N0": None, "statistic": None, "p_value": None, "a_hat_H0": None,
            "n_boot": None, "ci_a_hat": None, "ci_N_hat": None, "ci_n_boot": None,
            "ci_failure_fraction": None, "ci_method": None}


def analyze_degree_class(D, s, rng_test, rng_ci, n_boot=N_BOOT_TEST):
    D = np.asarray(D, dtype=float)
    try:
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            fit = network.fit_homogeneous_network(D, s)
        if not (math.isfinite(fit.a) and math.isfinite(fit.N)):
            raise ValueError("non-finite estimate (a factorial moment is zero or "
                             "the sample is under-dispersed)")
    except ValueError as exc:
        return _failed(str(exc))
    out = {"status": "ok", "reason": None, "a_hat": fit.a, "N_hat": fit.N,
           "N_round": fit.N_round, "lambda_hat": fit.lam}
    if fit.N_round >= 1:
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            res = network.test_depth_degrees(D, fit.N_round, s, rng_test, n_boot=n_boot)
        out.update({"test_status": "ok", "N0": fit.N_round, "statistic": res.statistic,
                    "p_value": res.p_value, "a_hat_H0": res.a_hat, "n_boot": res.n_boot})
    else:
        out.update({"test_status": "not run: N_round = 0 (the test needs N0 >= 1)",
                    "N0": 0, "statistic": None, "p_value": None, "a_hat_H0": None,
                    "n_boot": None})
    out.update(bootstrap_ci_degrees(D, s, rng_ci))
    return out


def analyze_positive_class(x, s, rng_test, rng_ci, block_length=1, x_for_ci=None,
                           n_boot=N_BOOT_TEST):
    """``x``: positive series (zeros already dropped). ``x_for_ci``: the series
    in time order *before* dropping zeros, so blocks are contiguous in time."""
    x = np.asarray(x, dtype=float)
    k2, k3 = identifiability.empirical_log_cumulants(x, r_max=3)
    try:
        fit = identifiability.fit_homogeneous(k2, k3, s)
        if not (math.isfinite(fit.a) and math.isfinite(fit.N)):
            raise ValueError("non-finite estimate")
    except ValueError as exc:
        return _failed(str(exc))
    out = {"status": "ok", "reason": None, "a_hat": fit.a, "N_hat": fit.N,
           "N_round": fit.N_round, "lambda_hat": None}
    if fit.N_round >= 1:
        res = identifiability.test_depth(x, fit.N_round, s, rng_test, n_boot=n_boot)
        out.update({"test_status": "ok", "N0": fit.N_round, "statistic": res.statistic,
                    "p_value": res.p_value, "a_hat_H0": res.a_hat, "n_boot": res.n_boot})
    else:
        out.update({"test_status": "not run: N_round = 0 (the test needs N0 >= 1)",
                    "N0": 0, "statistic": None, "p_value": None, "a_hat_H0": None,
                    "n_boot": None})
    out.update(bootstrap_ci_positive(x if x_for_ci is None else x_for_ci, s, rng_ci,
                                     block_length=block_length))
    return out


def select_class(classes: dict):
    """Class shown in the fit figures: the only admissible one, or, when both
    classes admit a fit, the one whose depth test has the larger p-value."""
    ok = [(k, c) for k, c in classes.items() if c["status"] == "ok"]
    if not ok:
        return None
    ok.sort(key=lambda kc: -(kc[1]["p_value"] if kc[1]["p_value"] is not None else -1.0))
    return ok[0][0]


def analyze_degree_series(system, name, D, rng_test, rng_ci, **extra):
    D = np.asarray(D)
    mean, var = float(D.mean()), float(D.var())
    rec = {"system": system, "series": name, "n": int(D.size), "mean": mean,
           "variance": var, "vmr": var / mean}
    rec.update(extra)
    rec["classes"] = {CLASS_KEY[s]: analyze_degree_class(D, s, rng_test, rng_ci)
                      for s in CLASSES}
    rec["selected_class"] = select_class(rec["classes"])
    print_series(rec)
    return rec


def _fmt_ci(ci):
    return "[--]" if ci is None else f"[{ci[0]:.3g}, {ci[1]:.3g}]"


def print_series(rec):
    print(f"\n-- {rec['system']} / {rec['series']} --  n={rec['n']}, "
          f"mean={rec['mean']:.4g}, var={rec['variance']:.4g}, VMR={rec['vmr']:.4g}")
    for key, c in rec["classes"].items():
        if c["status"] != "ok":
            print(f"  s={key}: FAILED: {c['reason']}")
            continue
        p = "n/a" if c["p_value"] is None else f"{c['p_value']:.4f}"
        print(f"  s={key}: a_hat={c['a_hat']:.4f} {_fmt_ci(c['ci_a_hat'])}  "
              f"N_hat={c['N_hat']:.4f} {_fmt_ci(c['ci_N_hat'])}  N_round={c['N_round']}  "
              f"p(H0:N={c['N0']})={p}  CI failure fraction={c['ci_failure_fraction']:.3f}")
