# Package index

## Special functions

Arbitrary-precision evaluation of the two Meijer-G shapes the whole
package reduces to.

<!-- end list -->

  - `g_gamma_chain()` :
    
    `G^{N,0}_{0,N}(z | -- ; b_1,...,b_N)` – the gamma-chain shape (`s
    = 1/2`).

  - `g_invgamma_chain()` :
    
    `G^{0,N}_{N,0}(z | b_1,...,b_N ; --)` – the inverse-gamma-chain
    shape (`s = 1`), via the reflection identity `= G^{N,0}_{0,N}(1/z |
    -- ; 1-b)`.

## Reservoir background

The cascade’s temperature law and its dynamics.

<!-- end list -->

  - `f_N()` :
    
    The reservoir temperature law `f_N(T_N)`, Eq. (33), corrected.

  - `sample_TN()` :
    
    Ancestral sampling of `T_N` down the chain, Eq. (30).

  - `sample_conditional()` :
    
    Draw `T_i ~ f(T_i | T_{i-1})`, Eqs. (31)/(32).

  - `conditional_pdf()` :
    
    `f(T_i | T_{i-1})`, Eq. (31) for `s = 1/2` and Eq. (32) for `s = 1`.

  - `kappa_from_a()` :
    
    Invert `a_i = 2*gamma_i/kappa_i^2` for the noise amplitudes.

  - `simulate_sde()` :
    
    Euler-Maruyama on the hierarchy `dT_i = -gamma_i (T_i - T_{i-1}) dt
    + kappa_i T_i^s T_{i-1}^{1-s} dW_i`, Eq. (29).

## Statistical factor

`B_N^(s)(E; a)` and its tail integral.

<!-- end list -->

  - `B()` :
    
    `B_N^{(s)}(E; a)`, Eq. (34).

  - `B_tilde()` :
    
    `Btilde_N^{(s)}(E) = beta0 * int_E^inf B_N(E') dE'`, Eqs. (37)-(38).

  - `B_prime()` :
    
    `dB/dE`, Eq. (43): strictly negative, hence `s''(B) < 0`, Eq. (44).

  - `omega_of()` :
    
    `omega = prod(a_i)`, Eq. (33). The scalar Eq. (38) leaves alone.

  - `shift()` :
    
    The parameter shift of Eq. (38): `a -> a + 3 - 4s`.

## Entropy generator

<!-- end list -->

  - `s_tilde()` :
    
    `s~_N^{(s)}(E) = (beta0*E - 1)*B(E) + Btilde(E)`, Eq. (39).

  - `parametric_curve()` :
    
    The parametric curve `(B(E), s~(E))` – geometric in `1 + beta0*E` so
    both the `B ~ 1` knee and the tail are resolved. Outside the three
    solvable cases `E(B)` has no closed form, so `s(B)` exists only as
    this curve or through `E_of_B`.

  - `E_of_B()` :
    
    `E` such that `B(E) = Bval`, vectorized over `Bval`. Well-posed for
    every `(N, a)` and both classes because `B` is strictly decreasing
    from `B(0) = 1` to `B(inf) = 0`, Eqs. (7)/(43).

  - `make_interpolants()` :
    
    Build monotone-cubic interpolants of `E(B)` and `s~(B)`, amortizing
    the cost of repeated `B()`/`E_of_B()` evaluation (each one, worst
    case, an escalating-precision Meijer-G or quadrature call – see
    R/meijer.R and R/statfactor.R) over a table built once. Two curves
    are tabulated, both against `log(B)`, NOT `E` or `B` directly:

  - `make_E_of_B()` :
    
    Convenience wrapper returning just the tabulated `E_of_B` function.

  - `s_of_B()` :
    
    `s~(E_of_B(Bval))`, exact (not tabulated) composition.

  - `s_prime_of_B()` :
    
    `d(s~)/d(B) = beta0*E_of_B(B) - 1`, exact via the inverse function
    (Eq. 43 makes `s~` a function of `E` whose derivative w.r.t. `B`
    collapses to this, with no differentiation needed).

  - `S_simplex()` :
    
    `S[p] = sum_i s_of_B(Z*p_i)/Z + c*sum(p_i)`, Eq. (10). `Z` (the
    partition function) is mandatory – it is supplied externally, not
    recomputed here.

  - `s_boltzmann_shannon()` :
    
    `-B*log(B)`, Eq. (21): the `N=0` / plain-exponential special case,
    where `s~` collapses to this exactly.

  - `s_tsallis_kernel()` :
    
    `s(B) = -(a+1)*(B - B^q)`, Eq. (63) – the `N=1`, `s=1`
    exactly-solvable case, `q = a/(a+1)` the entropy index. Summed over
    the simplex this is exactly the Tsallis entropy, Eqs. (64)/(67).

  - `q_indices()` :
    
    The two Tsallis indices for the `s=1`, `N=1` case (Eq. 65) –
    distinct and must not be conflated: `q` is the *entropy* index,
    `q_stat` the *distribution* index.

  - `tsallis_entropy()` :
    
    The exact Tsallis entropy of a probability vector `p`, at entropy
    index `q = a/(a+1)`.

## Identifiability and estimation

Recovering the cascade depth from a single marginal.

<!-- end list -->

  - `log_cumulants_theory()` :
    
    Theoretical `kappa_r(log T_N)` for `r=2..r_max`. Independent of
    `T_0` by construction; relies only on the shape vector `a` and
    branch `s`. `a` may be empty (`N = 0`), for which all cumulants are
    0.

  - `empirical_log_cumulants()` :
    
    Unbiased sample cumulants `k_r` of `log T` for `r=2..r_max`, via the
    standard k-statistic (directly comparable to `log_cumulants_theory`
    with no finite-sample bias correction).

  - `fit_homogeneous()` :
    
    Point estimate of `(N, a)` under the homogeneous model `a_i = a`,
    from `(kappa2, kappa3)`. Two equations, two unknowns:
    `kappa3/kappa2` depends on `a` alone (the `N` factor cancels),
    recovered by a 1-D root find; then `N = kappa2 / psi'(a)` from the
    scale.

  - `test_depth()` :
    
    Parametric-bootstrap test of `H0: depth = N0` (homogeneous `a`)
    against the sample `T`. Fits `a` under `H0` from `kappa2` alone, so
    the test statistic – observed `kappa3` minus the value predicted at
    `(N0, a_hat)` – is not tautologically zero. Its null distribution is
    a parametric bootstrap: simulate `n_boot` depth-`N0` cascades at the
    fitted shape and refit on each. `N` sits at the boundary of its
    domain (a positive integer), so the usual chi-square asymptotics for
    a likelihood-ratio test do not apply (Self & Liang 1987), hence
    bootstrap rather than a closed-form p-value.

  - `power_curve()` :
    
    Simulated power of `test_depth` at `H0: N = N0` against a true
    cascade `a_true` (its length is the true depth). Returns a named
    list `n -> power`: reject rate under a true deeper cascade, at a
    shallow null.

## Networks

The mixed-Poisson application, depth from a degree sequence alone.

<!-- end list -->

  - `sample_degrees()` :
    
    Draw a mixed-Poisson degree sequence, `D_v | W_v ~
    Poisson(lambda*W_v)`. `W` is drawn from `sample_TN` with `T0 = 1`,
    so `E[W] = 1` and `E[D] = lambda` exactly.

  - `log_factorial_moment_theory()` :
    
    `log E[W^k]`, summed over independent levels. Homogeneous case (`a_i
    = a`) reduces to `N*[log Gamma(a+k) - log Gamma(a) - k*log(a)]` for
    `s=1/2`, linear in `N`.

  - `empirical_log_factorial_moments()` :
    
    `(lambda_hat, c(L_2,...,L_k_max))` with `L_k = log(Ehat[W^k])`.
    `lambda_hat` is the sample mean of `D` (exact since `E[W]=1`). The
    falling-factorial product is exactly 0 whenever `D < k`, so no
    special-casing of small degrees is needed.

  - `fit_homogeneous_network()` :
    
    Recover `(N, a, lambda)` from a degree sequence, homogeneous model.
    Mirrors `fit_homogeneous`: the ratio `L_3/L_2` depends on `a` alone
    (the common `N` factor cancels), a 1-D root find, then `N = L_2 /
    [single-level log E[W^2]](a)` from the scale.

  - `test_depth_degrees()` :
    
    Parametric-bootstrap test of `H0: N = N0` (homogeneous `a`) against
    a degree sequence `D` – the factorial-moment analogue of
    `test_depth`. `fit_homogeneous_network`'s point estimate has no
    formal p-value attached; this closes that gap the same way
    `test_depth` does for a propensity sequence. A degree sequence whose
    `L3/L2` ratio cannot be matched by any `a` in the bracket is a
    genuine model-selection signal (raises, like
    `fit_homogeneous_network`), not something papered over.

## Maximum entropy

Closing the loop – maximizing `S[p]` recovers `p = B(E)/Z`.

<!-- end list -->

  - `target_distribution()` :
    
    `p_i = B(E_i)/Z` with `Z = sum_i B(E_i)`. Returns `list(p, Z, U)`
    with `U = sum_i p_i E_i` the mean energy the maximization is
    constrained to.

  - `maximize()` :
    
    Maximize `S[p] = sum_i Z^{-1} s(Z p_i)` under `sum_i p_i = 1` and
    `sum_i p_i E_i = U` (`U` taken from `target_distribution`). `Z` is
    held at the partition function of the target, a fixed property of
    the spectrum, not a variable of the problem. `E(B)` comes from the
    tabulated interpolant of `make_interpolants` (n=800 knots, not the
    module default of 300: the recovered probabilities inherit the
    table's accuracy directly, and 300 leaves the power-law class's tail
    good to only about three digits). SLSQP-equivalent optimization
    (`nloptr`'s NLOPT\_LD\_SLSQP) is run twice, the second time
    warm-started from the first answer: the entropy is very flat near
    its maximum, so a single pass stops early.

  - `stationarity_residual()` :
    
    Spread of `s'(Z p_i) - beta0*E_i` across states: this is the
    constant `alpha` at the maximizer, so the spread is zero there.
    Reported as `max - min`: a direct, optimizer-free check, and the one
    that localizes a broken link fastest.

## Internal utilities

<!-- end list -->

  - `vectorize_scalar()` : Lift a scalar function to accept and return
    numeric vectors/arrays.
