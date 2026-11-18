"""The paper's central claim: maximizing S[p] returns p_i = B(E_i)/Z.

The whole chain is under test here -- the branch assignment of Eq. (33), the
omega bookkeeping of Eq. (38), the gauge alpha = -1 of Eq. (16), the 1/Z scaling
of Eq. (10).  If any link is wrong, the maximizer moves.
"""

from __future__ import annotations

import numpy as np
import pytest

from htheory import maxent as mx
from htheory import statfactor as sf

# A spectrum that starts strictly above zero, so that B(E_i) < 1 for every
# state and the target sits in the interior of the domain Z p_i in (0, 1].
SPECTRUM = np.linspace(0.2, 9.0, 10)


def test_target_distribution_is_eq6(a_vec, s_class):
    p, Z, U = mx.target_distribution(SPECTRUM, a_vec, s_class)
    b = np.asarray(sf.B(SPECTRUM, a_vec, s_class), dtype=float)
    assert Z == pytest.approx(float(b.sum()))
    assert p == pytest.approx(b / Z)
    assert p.sum() == pytest.approx(1.0)
    assert U == pytest.approx(float(np.dot(p, SPECTRUM)))


def test_stationarity_holds_at_the_target(a_vec, s_class):
    """Eq. (12): s'(B_i) - beta0 E_i is the same constant alpha for every state.

    Optimizer-free, and the sharpest form of the claim -- it uses the exact
    bisection inverse, not the tabulated one, so it is limited only by the
    Meijer G evaluations themselves.
    """
    p, Z, _ = mx.target_distribution(SPECTRUM, a_vec, s_class)
    assert mx.stationarity_residual(p, SPECTRUM, a_vec, s_class, Z=Z) < 1e-9


def test_stationarity_fails_for_a_perturbed_distribution(a_vec, s_class):
    """The converse: the residual is not small for just any p on the simplex."""
    p, Z, _ = mx.target_distribution(SPECTRUM, a_vec, s_class)
    q = p * np.linspace(0.8, 1.2, p.size)
    q /= q.sum()
    assert mx.stationarity_residual(q, SPECTRUM, a_vec, s_class, Z=Z) > 1e-3


def test_maximization_recovers_the_distribution(a_vec, s_class):
    """The closure, run through an actual constrained optimizer.

    Accuracy is bounded by the tabulated E(B)/s(B) of
    :func:`~htheory.entropy.make_interpolants`, roughly six digits, not by the
    theory -- which :func:`test_stationarity_holds_at_the_target` pins to 1e-13.
    """
    res = mx.maximize(SPECTRUM, a_vec, s_class)
    assert res.success
    assert res.max_abs_error < 2e-4
    assert res.p.sum() == pytest.approx(1.0)
    assert float(np.dot(res.p, SPECTRUM)) == pytest.approx(res.U, rel=1e-6)


def test_maximum_is_a_maximum(a_vec, s_class):
    """Nothing feasible beats the target under the *exact* functional.

    The optimizer's own iterate can score a hair above it, because the entropy
    it sees is interpolated; the exact functional is the arbiter.
    """
    from htheory.entropy import S_simplex

    p_target, Z, U = mx.target_distribution(SPECTRUM, a_vec, s_class)
    best = S_simplex(p_target, a_vec, s_class, Z=Z)
    rng = np.random.default_rng(7)
    for _ in range(5):
        # Random feasible perturbation: keep sum p = 1 and sum p E = U.
        d = rng.normal(size=p_target.size)
        d -= d.mean()
        d -= np.dot(d, SPECTRUM - SPECTRUM.mean()) / np.dot(
            SPECTRUM - SPECTRUM.mean(), SPECTRUM - SPECTRUM.mean()
        ) * (SPECTRUM - SPECTRUM.mean())
        scale = 0.05 * p_target.min() / np.max(np.abs(d))
        q = p_target + scale * d
        assert q.min() > 0
        assert S_simplex(q, a_vec, s_class, Z=Z) < best


def test_boltzmann_case_gives_the_exponential():
    """N = 0: the construction must reproduce ordinary MaxEnt."""
    res = mx.maximize(SPECTRUM, [], 0.5)
    want = np.exp(-SPECTRUM)
    want /= want.sum()
    assert res.p == pytest.approx(want, abs=1e-3)


@pytest.mark.parametrize("a", [0.6, 2.0])
def test_tsallis_case_gives_the_q_exponential(a):
    """N = 1, s = 1: maximizing the Tsallis entropy under the *standard*
    constraints returns the q_stat-exponential of Eq. (66), with
    q_stat = 2 - q > 1 while the entropy index is q = a/(a+1) < 1."""
    res = mx.maximize(SPECTRUM, [a], 1.0)
    want = sf.B_q_exponential(SPECTRUM, a)
    want = want / want.sum()
    assert res.p == pytest.approx(want, abs=1e-4)


def test_result_reports_both_entropies(a_vec, s_class):
    res = mx.maximize(SPECTRUM, a_vec, s_class)
    assert np.isfinite(res.entropy)
    assert np.isfinite(res.entropy_target)
    assert res.entropy == pytest.approx(res.entropy_target, abs=1e-4)


def test_a_supplied_start_is_clipped_into_the_domain(a_vec, s_class):
    """A caller-supplied x0 outside Z p in (0, 1] must be pulled in, not
    silently used: s(B) has no meaning for B > 1, since E >= 0 and B(0) = 1."""
    x0 = np.full(SPECTRUM.size, 1.0 / SPECTRUM.size)
    res = mx.maximize(SPECTRUM, a_vec, s_class, x0=x0)
    assert res.success
    assert np.all(res.p * res.Z <= 1.0 + 1e-9)
