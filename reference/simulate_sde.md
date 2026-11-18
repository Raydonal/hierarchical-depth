# Euler-Maruyama on the hierarchy `dT_i = -gamma_i (T_i - T_{i-1}) dt + kappa_i T_i^s T_{i-1}^{1-s} dW_i`, Eq. (29).

`gamma` is the vector of relaxation rates `gamma_i`; the noise
amplitudes follow from `a_i = 2*gamma_i/kappa_i^2`. Eq. (29) presumes
`gamma_i >> gamma_{i-1}` (scale separation), the caller's
responsibility. Returns the post-burn-in trajectory of the innermost
level `T_N` (matrix, `n_paths` x kept steps), whose stationary histogram
is `f_N`. Full-truncation scheme (the `s = 1/2` CIR case is only
guaranteed nonnegative in the exact process).

## Usage

``` r
simulate_sde(
  a,
  s,
  gamma,
  T0 = 1,
  dt = 0.001,
  n_steps = 200000L,
  burn_in = 0.1,
  n_paths = 1L
)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - gamma:
    
    Vector of relaxation rates (Eq. 29).

  - T0:
    
    Base (outermost) reservoir temperature.

  - dt:
    
    Euler-Maruyama time step.

  - n\_steps:
    
    Number of Euler-Maruyama steps to simulate.

  - burn\_in:
    
    Fraction of steps discarded as burn-in.

  - n\_paths:
    
    Number of independent paths to simulate.

## Value

A `n_paths` by `(n_steps - burn-in steps)` matrix of the post-burn-in
`T_N` trajectory.
