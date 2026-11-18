"""One-time oracle dump: golden values from the trusted Python `htheory` package.

Run this from the repo root with the project's Python environment (the one
`requirements.txt` describes), NOT as part of any R build/check/test:

    python3 -m htheoryR.data-raw.generate_golden_values   # not how you run it
    PYTHONPATH=. python3 htheoryR/data-raw/generate_golden_values.py

It writes `htheoryR/tests/testthat/fixtures/golden_htheory.csv`, which the R
port's `testthat` suite reads to check its own values against Python's
`mpmath`-backed ones without ever depending on Python at R build/check time.

Every expected value is stored as a **full-precision decimal string**, not a
rounded double, precisely so the R-side comparison (via `Rmpfr`) isn't capped
at double precision -- that would hide the 16th-25th digits this whole
exercise exists to validate. See GOLDEN_PROVENANCE.md for the exact htheory
commit/library versions this was generated against.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import mpmath as mp  # noqa: E402

from htheory import background, meijer, statfactor  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent.parent / "tests" / "testthat" / "fixtures" / "golden_htheory.csv"

# tests/conftest.py::PARAM_SETS, kept in lockstep with the Python suite.
A_VECTORS = [
    ([], "N=0"),
    ([1.0], "N=1 a=1"),
    ([0.5], "N=1 a=1/2"),
    ([2.0], "N=1 a=2"),
    ([1.0, 1.0], "N=2 a=1"),
    ([1.5, 2.0, 0.8], "N=3 unequal"),
]
CLASSES = [0.5, 1.0]

# Near-zero, interior, and tail values, in both the T (background) and E
# (statfactor) roles -- the two functions use the same argument range.
X_GRID = [1e-9, 1e-6, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 50.0, 200.0, 1e4]

DPS = 20  # matches meijer.DEFAULT_DPS, the runtime precision being validated


def mpf_str(x: mp.mpf) -> str:
    return mp.nstr(x, 30, strip_zeros=False)


def safe_call(fn, *args, **kwargs):
    """Return (value_str, ok). `ok=False` records a documented failure mode
    (e.g. MeijerGFailure/domain error) rather than silently skipping the row.
    """
    try:
        with mp.workdps(DPS):
            val = fn(*args, **kwargs)
        return mpf_str(mp.mpf(val)), True
    except Exception as exc:  # noqa: BLE001 -- deliberately broad, see docstring
        return f"ERROR: {type(exc).__name__}: {exc}", False


def raw_meijerg_rows(a, label, s):
    """Exercise `meijer.meijerg` directly, at the exact parameter shapes
    `f_N`/`B` build internally -- validates the lowest layer independent of
    any background/statfactor plumbing bug.
    """
    rows = []
    a_arr = background._as_a(a)
    omega = background._omega(a_arr)
    for x in X_GRID:
        if s == 0.5:
            a_s = [[], []]
            b_upper = [float(v) for v in (a_arr - 1)]
            b_s = [b_upper, []]
            z = x * omega
            underflow_above = meijer.UNDERFLOW_Z
            underflow_below = None
        else:
            a_upper = [float(v) for v in (-a_arr - 1)]
            a_s = [a_upper, []]
            b_s = [[], []]
            z = x / omega if omega else x
            underflow_above = None
            underflow_below = meijer.UNDERFLOW_Z_SMALL
        val_str, ok = safe_call(
            meijer.meijerg, a_s, b_s, z,
            dps=DPS, underflow_above=underflow_above, underflow_below=underflow_below,
        )
        rows.append(dict(
            func="meijerg_gamma_chain" if s == 0.5 else "meijerg_invgamma_chain",
            a_json=json.dumps(a), s=s, x=x, beta0="", T0="",
            expected_value_str=val_str, dps_used=DPS, ok=ok,
        ))
    return rows


def public_function_rows(a, label, s):
    rows = []
    beta0, T0 = 1.0, 1.0
    for x in X_GRID:
        # f_N(T) -- background law
        val_str, ok = safe_call(background._f_N_scalar, x, background._as_a(a), s, T0)
        rows.append(dict(func="f_N", a_json=json.dumps(a), s=s, x=x, beta0="", T0=T0,
                          expected_value_str=val_str, dps_used=DPS, ok=ok))

        # B(E), B_tilde(E), B_prime(E), s_tilde(E) -- statfactor/entropy layer
        val_str, ok = safe_call(statfactor.B, x, a, s, beta0)
        rows.append(dict(func="B", a_json=json.dumps(a), s=s, x=x, beta0=beta0, T0="",
                          expected_value_str=val_str, dps_used=DPS, ok=ok))

        val_str, ok = safe_call(statfactor.B_tilde, x, a, s, beta0)
        rows.append(dict(func="B_tilde", a_json=json.dumps(a), s=s, x=x, beta0=beta0, T0="",
                          expected_value_str=val_str, dps_used=DPS, ok=ok))

        val_str, ok = safe_call(statfactor.B_prime, x, a, s, beta0)
        rows.append(dict(func="B_prime", a_json=json.dumps(a), s=s, x=x, beta0=beta0, T0="",
                          expected_value_str=val_str, dps_used=DPS, ok=ok))
    return rows


def closed_form_rows(a, label, s):
    rows = []
    beta0 = 1.0
    for x in X_GRID:
        if s == 0.5 and len(a) == 1:
            val_str, ok = safe_call(statfactor.B_K_distribution, x, a[0], beta0)
            rows.append(dict(func="B_K_distribution", a_json=json.dumps(a), s=s, x=x,
                              beta0=beta0, T0="", expected_value_str=val_str, dps_used=DPS, ok=ok))
        if s == 0.5 and a == [0.5]:
            val_str, ok = safe_call(statfactor.B_stretched_exact, x, beta0)
            rows.append(dict(func="B_stretched_exact", a_json=json.dumps(a), s=s, x=x,
                              beta0=beta0, T0="", expected_value_str=val_str, dps_used=DPS, ok=ok))
        if s == 1.0 and len(a) == 1:
            val_str, ok = safe_call(statfactor.B_q_exponential, x, a[0], beta0)
            rows.append(dict(func="B_q_exponential", a_json=json.dumps(a), s=s, x=x,
                              beta0=beta0, T0="", expected_value_str=val_str, dps_used=DPS, ok=ok))
            val_str, ok = safe_call(statfactor.B_tilde_q_exponential, x, a[0], beta0)
            rows.append(dict(func="B_tilde_q_exponential", a_json=json.dumps(a), s=s, x=x,
                              beta0=beta0, T0="", expected_value_str=val_str, dps_used=DPS, ok=ok))
    return rows


def main():
    fieldnames = ["func", "a_json", "s", "x", "beta0", "T0", "expected_value_str", "dps_used", "ok"]
    rows = []
    for a, label in A_VECTORS:
        for s in CLASSES:
            rows.extend(raw_meijerg_rows(a, label, s))
            rows.extend(public_function_rows(a, label, s))
            rows.extend(closed_form_rows(a, label, s))

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    n_ok = sum(1 for r in rows if r["ok"])
    print(f"Wrote {len(rows)} rows ({n_ok} ok, {len(rows) - n_ok} documented errors) to {OUT_PATH}")


if __name__ == "__main__":
    main()
