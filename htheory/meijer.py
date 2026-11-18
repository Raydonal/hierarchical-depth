"""Robust, cached evaluation of the Meijer G function.

Every closed form in the paper (Eqs. 33, 34, 40, 41) is a Meijer G. ``mpmath``
is the only Python library that evaluates it for arbitrary parameters, but it
raises rather than returning a value when its internal hypergeometric expansion
cannot reach the requested precision -- which happens routinely in the tails,
where the true value underflows to zero. Every G evaluation in this package
therefore goes through :func:`meijerg`.

Convention (mpmath's, which matches the paper's):

    G^{m,n}_{p,q}( z | a_1..a_p ; b_1..b_q )
        a_s = [[a_1..a_n], [a_{n+1}..a_p]]
        b_s = [[b_1..b_m], [b_{m+1}..b_q]]
"""

from __future__ import annotations

import functools

import mpmath as mp
import numpy as np
from mpmath.libmp.libhyper import NoConvergence

__all__ = [
    "meijerg",
    "vectorize",
    "with_dps",
    "MeijerGFailure",
    "UNDERFLOW_Z",
    "UNDERFLOW_Z_SMALL",
]

DEFAULT_DPS = 20

# Where a failure to evaluate is credibly an underflow rather than a *bug*.
# Returning 0 unconditionally is how a G worth 0.9976 once became 0.0 in the
# small-argument end of the power-law class -- silently, because 0 is a
# perfectly plausible value for these functions. These are the tripwires:
# outside them the value really has decayed, inside them something is wrong
# and the caller must hear about it.
#
# Both directions are needed, because which end decays depends on the G, not
# on the argument. The two orientations that appear in the paper are mirror
# images: G^{N,0}_{0,N}(z | - ; a-1) carries the e^{-z} of the gamma chain and
# dies at large z, while its inverse G^{0,N}_{N,0}(z | -a-1 ; -) carries the
# e^{-1/z} of the inverse-gamma chain and dies at small z.
UNDERFLOW_Z = 100.0
UNDERFLOW_Z_SMALL = 1.0 / UNDERFLOW_Z


class MeijerGFailure(RuntimeError):
    """mpmath could not evaluate a G function where the value is not negligible."""


def with_dps(dps: int):
    """Decorator: evaluate the wrapped function at ``dps`` decimal digits."""

    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            with mp.workdps(dps):
                return fn(*args, **kwargs)

        return wrapper

    return deco


# mpmath signals "I cannot reach the requested accuracy" in several ways:
# ``NoConvergence`` from the series summation, ``ValueError`` from
# ``hypercomb`` when cancellation between terms eats the whole working
# precision, and ``ZeroDivisionError``/``OverflowError`` at extreme arguments.
# All of them mean the same thing here.
_NONCONVERGENCE = (NoConvergence, ValueError, ZeroDivisionError, OverflowError)

# Cap on the series length.  mpmath's default (10**6) turns every failure into
# a multi-minute stall, and failures are routine -- they happen at exactly the
# large arguments where the series is slowest and the answer is nearest zero.
# 10**4 terms is far more than any converging case here needs, and makes the
# non-converging ones fail in milliseconds.
MAXTERMS = 10**4


# Sentinel distinguishing "mpmath gave up" from "mpmath returned zero".
_FAILED = "<failed>"


@functools.lru_cache(maxsize=200_000)
def _meijerg_cached(a_upper, a_lower, b_upper, b_lower, z_str, dps):
    with mp.workdps(dps):
        try:
            return str(
                mp.meijerg(
                    [list(a_upper), list(a_lower)],
                    [list(b_upper), list(b_lower)],
                    mp.mpf(z_str),
                    maxterms=MAXTERMS,
                )
            )
        except _NONCONVERGENCE:
            return _FAILED


def meijerg(
    a_s,
    b_s,
    z,
    dps: int = DEFAULT_DPS,
    underflow_above: float | None = None,
    underflow_below: float | None = None,
) -> mp.mpf:
    """Evaluate ``G^{m,n}_{p,q}(z)``.

    Results are cached on the exact parameter/argument string, so repeated
    evaluation on the same grid (the plotting and inversion routines do this
    heavily) costs one mpmath call per distinct point.

    ``underflow_above`` is what a *failure* means. mpmath raises rather than
    returning a value when its hypergeometric expansion cannot reach the
    requested precision, and that happens for two unrelated reasons: the true
    value has underflowed (large ``z``), or the parameters are degenerate and
    the expansion cancels itself to nothing (small ``z``, coincident
    parameters -- routine in this package, where ``a_i = a`` for all levels).
    Only the first justifies returning 0.

    So a failure returns 0 only when ``z >= underflow_above`` or
    ``z <= underflow_below``; otherwise it raises :class:`MeijerGFailure`.
    Passing neither (the default) means the caller has no underflow regime and
    every failure is an error.

    Which end decays is a property of the G, not of the argument: the
    gamma-chain form dies at large ``z``, its inverse at small ``z``. A caller
    that passes a threshold must be able to say why the value is negligible
    past it. One that cannot should not pass one -- see
    :func:`htheory.statfactor.B`, which covers the small-argument end of the
    power-law class with an exact series rather than a guess.
    """
    z = mp.mpf(z)
    if z <= 0:
        raise ValueError(f"Meijer G argument must be positive, got {z}")
    key = (
        tuple(str(mp.mpf(x)) for x in a_s[0]),
        tuple(str(mp.mpf(x)) for x in a_s[1]),
        tuple(str(mp.mpf(x)) for x in b_s[0]),
        tuple(str(mp.mpf(x)) for x in b_s[1]),
        mp.nstr(z, dps + 5, strip_zeros=False),
        dps,
    )
    out = _meijerg_cached(*key)
    if out == _FAILED:
        if underflow_above is not None and z >= underflow_above:
            return mp.mpf(0)
        if underflow_below is not None and z <= underflow_below:
            return mp.mpf(0)
        raise MeijerGFailure(
            f"mpmath could not evaluate G at z = {mp.nstr(z, 8)} with "
            f"a_s = {a_s}, b_s = {b_s}. That is inside the range where the "
            f"value is not negligible (underflow_above={underflow_above}, "
            f"underflow_below={underflow_below}), so 0 would be wrong."
        )
    return mp.mpf(out)


def vectorize(scalar_fn):
    """Lift a scalar mpmath function to accept and return numpy arrays.

    Scalar input gives a Python float back; array input gives a float64 array of
    the same shape.
    """

    @functools.wraps(scalar_fn)
    def wrapper(x, *args, **kwargs):
        arr = np.asarray(x, dtype=float)
        if arr.ndim == 0:
            return float(scalar_fn(float(arr), *args, **kwargs))
        flat = np.array(
            [float(scalar_fn(float(v), *args, **kwargs)) for v in arr.ravel()],
            dtype=float,
        )
        return flat.reshape(arr.shape)

    return wrapper
