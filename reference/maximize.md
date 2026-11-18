# Maximize `S[p] = sum_i Z^{-1} s(Z p_i)` under `sum_i p_i = 1` and `sum_i p_i E_i = U` (`U` taken from `target_distribution`). `Z` is held at the partition function of the target, a fixed property of the spectrum, not a variable of the problem. `E(B)` comes from the tabulated interpolant of `make_interpolants` (n=800 knots, not the module default of 300: the recovered probabilities inherit the table's accuracy directly, and 300 leaves the power-law class's tail good to only about three digits). SLSQP-equivalent optimization (`nloptr`'s NLOPT\_LD\_SLSQP) is run twice, the second time warm-started from the first answer: the entropy is very flat near its maximum, so a single pass stops early.

Maximize `S[p] = sum_i Z^{-1} s(Z p_i)` under `sum_i p_i = 1` and `sum_i
p_i E_i = U` (`U` taken from `target_distribution`). `Z` is held at the
partition function of the target, a fixed property of the spectrum, not
a variable of the problem. `E(B)` comes from the tabulated interpolant
of `make_interpolants` (n=800 knots, not the module default of 300: the
recovered probabilities inherit the table's accuracy directly, and 300
leaves the power-law class's tail good to only about three digits).
SLSQP-equivalent optimization (`nloptr`'s NLOPT\_LD\_SLSQP) is run
twice, the second time warm-started from the first answer: the entropy
is very flat near its maximum, so a single pass stops early.

## Usage

``` r
maximize(
  energies,
  a,
  s,
  beta0 = 1,
  x0 = NULL,
  p_min = 1e-14,
  tol = 1e-12,
  maxiter = 500L,
  n_table = 800L
)
```

## Arguments

  - energies:
    
    Numeric vector of energy levels `E_i`.

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - beta0:
    
    Inverse base temperature, `1/T0`.

  - x0:
    
    Optional starting point for the optimizer; defaults to a feasible
    tilted-uniform start.

  - p\_min:
    
    Lower bound on each `p_i` (avoids the tabulated interpolant's lower
    edge).

  - tol:
    
    Optimizer relative tolerance.

  - maxiter:
    
    Maximum optimizer iterations.

  - n\_table:
    
    Number of knots in the tabulated interpolant.

## Value

A list with components `p`, `p_target`, `Z`, `U`, `entropy`,
`entropy_target`, `max_abs_error`, `success`, `message`,
`nloptr_result`.
