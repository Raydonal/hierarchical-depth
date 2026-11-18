# Theorems 1–3: statements and proofs

This note gives the formal statements and proofs of Theorems 1–3. It replaces the original paper's
page of Meijer-G identities (Theorem 1), its case-by-case derivations of the
entropy generator (Theorem 2), and closes the estimation gap the original
paper leaves open (Theorem 3). Notation: `B(E) = E[e^{-E/T}]`
is the generalized statistical factor, `T` the fluctuating reservoir
temperature, `β0 = 1/T0` the base inverse temperature, and — for a cascade —
`T = T0 ∏_{i=1}^N X_i` with the `X_i` the per-level multiplicative factors of
`background.sample_conditional`.

All three theorems have a numerical witness in the `htheory` package:
Theorem 1 in `htheory/statfactor.py:B_tilde` and the docstring of
`htheory/identifiability.py`; Theorem 2 in `htheory/entropy.py:s_tilde`, which
implements exactly the closed form derived below; Theorem 3 directly in
`htheory/identifiability.py`, with `tests/test_identifiability.py` and
`tests/test_network.py` as the empirical checks.

---

## Theorem 1 — the size-biasing identity behind Eq. (38)

**Setup.** Let `T > 0` be a random variable with `E[T] = 1/β0 < ∞`, and let

  B(E) = E[e^{-E/T}],  E ≥ 0.

By Bernstein–Widder, `B` completely monotone with `B(0) = 1` is *equivalent*
to `B` being of this form for some positive `T` — so nothing below assumes
more than "`B` is a Boltzmann average of a finite-mean temperature."

**Definitions.** The *size-biased* (length-biased) version `T*` of `T` is the
random variable with density `f*(t) = t f(t) / E[T]` (well defined whenever
`E[T] < ∞`; this is the standard construction from renewal theory, e.g.
Arratia–Goldstein–Kochman). Define the tail-integral partner

  B̃(E) := β0 ∫_E^∞ B(E') dE'.

**Theorem 1.** `B̃(E) = E[e^{-E/T*}]`. If moreover `T = T0 ∏_{i=1}^N X_i` with
`T0 > 0` constant and `X_1, …, X_N` mutually independent, then

  T* =_d T0 ∏_{i=1}^N X_i*,

i.e. size-biasing the cascade size-biases *each level independently*, leaving
the others (and `T0`) untouched. For the two families in
`background.sample_conditional`:

* `X_i ~ Gamma(a_i, rate = a_i)` (the `s = 1/2` branch) has
  `X_i* ~ Gamma(a_i + 1, rate = a_i)` — shape up by 1, rate unchanged.
* `X_i ~ InvGamma(a_i + 1, scale = a_i)` (the `s = 1` branch) has
  `X_i* ~ InvGamma(a_i, scale = a_i)` — shape down by 1, scale unchanged.

which is exactly the shift `a → a + 3 - 4s` of Eq. (38) (`+1` at `s = 1/2`,
`-1` at `s = 1`) acting on the shape parameter that appears inside the
Meijer-G kernel of Eq. (33), while `ω = ∏ a_i` — built entirely from the
*rates/scales*, which size-biasing never touches — stays fixed. That is the
fact the original paper's Appendix C does not explain and this proof does.

### Proof

**Step 1 (the tail integral is a size-biased expectation).** Since the
integrand is nonnegative, Fubini/Tonelli lets the `E'`-integral and the
expectation over `T` be swapped:

  ∫_E^∞ B(E') dE' = ∫_E^∞ E[e^{-E'/T}] dE' = E[∫_E^∞ e^{-E'/T} dE'] = E[T e^{-E/T}],

using `∫_E^∞ e^{-E'/T} dE' = T e^{-E/T}` for each fixed `T > 0`. Hence

  B̃(E) = β0 E[T e^{-E/T}] = (1/E[T]) E[T e^{-E/T}] = E[e^{-E/T*}]

directly from the definition `E[g(T*)] = E[T g(T)] / E[T]` applied to
`g(t) = e^{-E/t}`. This proves the first claim, for *any* finite-mean `T` —
no special functions, no independence assumption yet.

**Step 2 (size-biasing a product of independent factors).** Let
`Z = X_1 ⋯ X_N` with the `X_i` independent and each `E[X_i] < ∞`, so
`E[Z] = ∏ E[X_i] < ∞`. For any test function `g`,

  E[g(X_1* ⋯ X_N*)]
    = ∫⋯∫ g(x_1⋯x_N) ∏_i [x_i f_{X_i}(x_i) / E[X_i]] dx_1⋯dx_N
    = (1 / ∏_i E[X_i]) ∫⋯∫ (x_1⋯x_N) g(x_1⋯x_N) ∏_i f_{X_i}(x_i) dx_1⋯dx_N
    = (1 / E[Z]) E[Z g(Z)]
    = E[g(Z*)],

using independence to factor the joint density in the second line and
regrouping `∏_i x_i = x_1⋯x_N` into the single weight `z` in the third. So
`Z* =_d X_1^* ⋯ X_N^*` — size-biasing a product of independent factors is the
same in law as size-biasing every factor independently and multiplying. (The
argument is unchanged by inserting a constant `T0` as an extra "factor": a
constant is trivially its own size-biased version, since `f*(t) = t·δ(t-T0)/T0
= δ(t-T0)`.) Applying this with `Z = T0 X_1 ⋯ X_N = T` gives
`T* =_d T0 X_1^* ⋯ X_N^*`.

**Step 3 (size-biasing gamma and inverse-gamma).** For `X ~ Gamma(α, rate =
β)`, density `f(x) ∝ x^{α-1} e^{-βx}`, so `f*(x) ∝ x·x^{α-1}e^{-βx} =
x^α e^{-βx} ∝` the `Gamma(α+1, rate=β)` density: shape up by one, rate
untouched. For `X ~ InvGamma(α, scale=θ)` (`α>1`, so `E[X]<∞`), density
`f(x) ∝ x^{-α-1}e^{-θ/x}`, so `f*(x) ∝ x^{-α}e^{-θ/x} = x^{-(α-1)-1}e^{-θ/x}
∝` the `InvGamma(α-1, scale=θ)` density: shape down by one, scale untouched.
Substituting `α = a_i` (gamma branch) and `α = a_i + 1` (inverse-gamma
branch) reproduces the two bullet points above. ∎

**Remark.** This is also the two-line explanation, promised in the article
outline, for why `ω` in Eq. (38) is untouched: `ω` is built from the *rates*
(gamma branch) or *scales* (inverse-gamma branch) of the per-level factors,
and Step 3 shows size-biasing never moves those — only the shape.

---

## Theorem 2 — the entropy generator, without special functions

**Setup.** As in Theorem 1, `B(E) = E[e^{-E/T}]` for some `T` with `E[T] =
1/β0 < ∞`, and additionally `B` is *strictly* decreasing (true whenever `T`
is non-degenerate, and the content of Eq. (43) for the specific cascades of
this package). Assume the flat-measure sum rule (Eq. 27, `∫_0^∞ B(E) dE =
1/β0`) — automatic here, since `∫_0^∞ B = ∫_0^∞ E[e^{-E/T}] dE = E[T] = 1/β0`
by the same Fubini step as Theorem 1.

The inverse-MaxEnt construction (Naudts' φ-exponential families; the trace
form `S[p] = Σ Z^{-1} s(Z p_i)`, Eq. 10) asks for a generator `s` such that
maximizing `S[p]` subject to `Σ p_i = 1`, `Σ p_i E_i = U` is stationary
exactly at `p_i = B(E_i)/Z`. The stationarity (Lagrange) condition for
`S[p] = Σ Z^{-1}s(Zp_i)` is `s'(Z p_i) = α + β0 E_i` for a multiplier `α`; at
the target, `Z p_i = B(E_i)`, so writing `x = B(E)` and `E = E(x)` for the
(well-defined, since `B` is strictly monotone) inverse,

  s'(x) = α + β0 E(x).

The gauge `α = -1` (Eq. 16) fixes the free additive constant in `s` by
requiring `s(1) = 0` (a certain outcome, `E = 0` hence `x = B(0) = 1`, carries
no entropy) — see the derivation below for why this is the right normalizing
condition. Define `s̃(E) := s(B(E))`, the generator read as a function of
energy (Eq. 39).

**Theorem 2.** For any `B` as above,

  s̃(E) = (β0 E - 1) B(E) + B̃(E),

with `B̃` as in Theorem 1. Consequently `s̃` requires no special function
beyond whatever closed form `B` itself has — the cascade's Meijer G, in this
package, but the *formula relating `s̃` to `B`* is generic. Moreover, as a
function of `x = B(E)`, `s` is strictly concave, `s''(x) < 0`, if and only if
`B` is strictly decreasing.

### Proof

**Step 1 (integrate the stationarity condition).** From `s'(x) = -1 + β0
E(x)`, change variables to `E` via `x = B(E)`, `dx = B'(E) dE`:

  s̃'(E) = d/dE [s(B(E))] = s'(B(E)) B'(E) = (β0 E - 1) B'(E).

Integrate by parts, `u = β0E - 1`, `dv = B'(E)dE`:

  s̃(E) - s̃(0) = ∫_0^E (β0E' - 1) B'(E') dE'
              = [(β0E' - 1)B(E')]_0^E - β0∫_0^E B(E')dE'
              = (β0E - 1)B(E) - (-1)·B(0) - β0∫_0^E B(E')dE'
              = (β0E - 1)B(E) + 1 - β0∫_0^E B(E')dE',

using `B(0) = 1`. The gauge `s(1) = 0`, i.e. `s̃(0) = s(B(0)) = s(1) = 0`,
kills the boundary term on the left.

**Step 2 (rewrite the running integral via the sum rule).** By the sum rule,
`∫_0^E B dE' = 1/β0 - ∫_E^∞ B dE' = 1/β0 - B̃(E)/β0`, so

  -β0∫_0^E B(E')dE' = -1 + B̃(E).

Substituting into Step 1:

  s̃(E) = (β0E - 1)B(E) + 1 - 1 + B̃(E) = (β0E - 1)B(E) + B̃(E). ∎

**Step 3 (concavity).** Since `s'(x) = β0 E(x) - 1`,

  s''(x) = β0 E'(x) = β0 / B'(E(x))

by the inverse-function rule (`E'(x) = 1/B'(E)` at `E = E(x)`). As `β0 > 0`,
`s''(x)` has the sign of `B'(E(x))`. So `s''(x) < 0` everywhere on `(0,1]`
if and only if `B' < 0` everywhere — `B` strictly decreasing. ∎

**Remark.** This is the entire content of the original paper's §5: every
closed form there (the K-distribution entropy, the Tsallis case, etc.) is
this one formula specialized to a particular `B`, and concavity is not a
property that needs checking case by case — it *is* `B' < 0`, which for a
cascade of any depth and shape is Eq. (43), verified once for the general
Meijer-G form of `B` and inherited by every special case.

---

## Theorem 3 — identifiability of the depth `N`

**Setup.** Write the cascade in log-temperature: for `s = 1/2`,
`background.sample_conditional` draws `T_i = T_{i-1} · Y_i / a_i` with
`Y_i ~ Gamma(a_i, 1)`; for `s = 1`, it draws `T_i = a_i T_{i-1} / Y_i` with
`Y_i ~ Gamma(a_i+1, 1)`. In both cases the increment

  D_i := log T_i - log T_{i-1}

is a function of `Y_i` alone (plus a constant `± log a_i`) and hence
independent of `T_{i-1}` — and, since each `Y_i` is drawn independently, the
`D_i` are mutually independent. `T_0` is a deterministic constant.

**Theorem 3.** For `r ≥ 2`, the `r`-th cumulant of `log T_N` is

  s = 1/2:  κ_r(log T_N) = Σ_{i=1}^N ψ^{(r-1)}(a_i)
  s = 1:    κ_r(log T_N) = (-1)^r Σ_{i=1}^N ψ^{(r-1)}(a_i + 1)

where `ψ^{(r-1)}` is the polygamma function of order `r-1`, and these do
*not* depend on `T_0`. In the homogeneous case `a_i ≡ a`, the ratio
`κ_3/κ_2` depends on `a` alone (not on `N`), so `(N, a)` — and in particular
whether `N = 1` is consistent with the data — is recovered from `(κ_2,
κ_3)` of a single scalar observable, `log T_N`, with no access to any
intermediate level. This directly contradicts the original paper's
conclusion (its own conclusion section) that `N` "cannot, in general, be
inferred from the foreground distribution alone": that conclusion was read
off a visual comparison of curves, not derived as a statement about an
estimation problem, and Theorem 3 is precisely the estimation-theoretic
counterexample.

### Proof

**Step 1 (cumulants add over independent summands, and a constant only
shifts the mean).** `log T_N = log T_0 + Σ_{i=1}^N D_i` with the `D_i`
mutually independent and `T_0` constant. The cumulant generating function of
a sum of independent variables is the sum of their cumulant generating
functions, so `κ_r(log T_N) = κ_r(log T_0) + Σ_i κ_r(D_i)` for every `r`; and
`κ_r` of a constant is `0` for `r ≥ 2` (only `κ_1`, the mean, sees a
constant shift). Hence `κ_r(log T_N) = Σ_i κ_r(D_i)` for `r ≥ 2`.

**Step 2 (log-cumulants of a gamma variable are scale-free).** For `Y ~
Gamma(α, rate=β)`, the cumulant generating function of `log Y` is

  K(t) = log E[Y^t] = log Γ(α+t) - log Γ(α) - t log β

(using `E[Y^t] = Γ(α+t)β^{-t}/Γ(α)`, valid for `t > -α`). Then `K'(t) =
ψ(α+t) - log β` and, for `r ≥ 2`, `K^{(r)}(t) = ψ^{(r-1)}(α+t)` — the `-t log
β` term does not survive past the first derivative. So `κ_r(log Y) =
K^{(r)}(0) = ψ^{(r-1)}(α)` for `r ≥ 2`, independent of the rate `β`. This is
the identity behind the method of log-cumulants in SAR image analysis
(Nicolas 2002): population log-cumulants of order ≥ 2 see through an unknown
scale.

**Step 3 (apply Step 2 to each branch).**

* `s = 1/2`: `D_i = log Y_i - log a_i` with `Y_i ~ Gamma(a_i, 1)`
  (`β=1`). The additive constant `-log a_i` drops out for `r ≥ 2` (Step 1's
  constant-shift argument, applied at the level of a single summand), so
  `κ_r(D_i) = κ_r(log Y_i) = ψ^{(r-1)}(a_i)` by Step 2.
* `s = 1`: `D_i = log a_i - log Y_i` with `Y_i ~ Gamma(a_i+1, 1)`. Negation
  scales cumulants by `(-1)^r`: `κ_r(D_i) = (-1)^r κ_r(log Y_i) = (-1)^r
  ψ^{(r-1)}(a_i+1)`.

Summing over `i = 1, …, N` via Step 1 gives the theorem's formula. ∎

**Step 4 (the homogeneous ratio is a well-posed inverse problem).** For `a_i
≡ a`, `κ_2 = N k_2(a)` and `κ_3 = N k_3(a)` where `k_2(a) = ψ'(a)` and
`k_3(a) = ψ''(a)` (`s=1/2`; the `s=1` case carries the same argument with
`a → a+1` and an extra sign that cancels in the ratio). The factor of `N`
cancels in `r(a) := κ_3/κ_2 = k_3(a)/k_2(a)`, so `a` is recovered from `r(a)
= κ_3/κ_2` alone, and then `N = κ_2 / k_2(a)`. This is well posed provided
`r` is injective on the admissible range of `a`, which the endpoint behavior
supports: as `a → ∞`, `ψ'(a) ~ 1/a` and `ψ''(a) ~ -1/a^2`, so `r(a) → 0^-`;
as `a → 0^+`, `ψ'(a) ~ 1/a^2` and `ψ''(a) ~ -2/a^3`, so `r(a) → -∞`. `r`
sweeping the whole range `(-∞, 0)` monotonically is consistent with known
monotonicity results for ratios of polygamma functions (see e.g. Alzer's and
Guo–Qi's work on complete monotonicity of polygamma ratios), and is what
`htheory.identifiability.fit_homogeneous`'s bracketed root-find relies on
operationally; `tests/test_identifiability.py` checks the recovery is exact
(to solver tolerance) against the closed form and accurate against Monte
Carlo, across every `a` the test suite's fixtures cover. ∎ (numerically)

**Step 5 (from log-cumulants to a hypothesis test).** Because `N` is a
positive integer at the boundary of its own parameter space, the usual
chi-square asymptotics for a likelihood-ratio test of `H0: N=N0` do not apply
(Self & Liang 1987) — `htheory.identifiability.test_depth` therefore uses a
parametric bootstrap: fit `a` under `H0` from `κ_2` alone, predict `κ_3`
at `(N0, â)`, and calibrate the discrepancy's null distribution by
simulating `N0`-level cascades from the fitted `â`. `power_curve` in the
same module gives the resulting test power against a genuinely deeper
cascade as a function of sample size (Fig. evidence in the article draft).

**Remark (why this generalizes to the network application).** Feller (1943)
established that mixed-Poisson laws are identifiable from the mixture
distribution — here, the degree distribution. Composed with Theorem 3, this
gives a two-step identifiability chain: the law of the mixing variable `W`
is recoverable from the degree sequence, and `(N, {a_i})` is recoverable from
the law of `W` — so the hierarchical depth of a mixed-Poisson network is
estimable from degrees alone, never observing which vertex belongs to which
community. `htheory/network.py` implements the corresponding
factorial-moment estimator, `E[D^{(k)}] = λ^k E[W^k]` being an exact
mixed-Poisson identity (no approximation), with `log E[W^k]` given in closed
form and linear in `N` exactly as `κ_r` is above.
