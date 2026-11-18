# `log E[W^k]`, summed over independent levels. Homogeneous case (`a_i = a`) reduces to `N*[log Gamma(a+k) - log Gamma(a) - k*log(a)]` for `s=1/2`, linear in `N`.

`log E[W^k]`, summed over independent levels. Homogeneous case (`a_i =
a`) reduces to `N*[log Gamma(a+k) - log Gamma(a) - k*log(a)]` for
`s=1/2`, linear in `N`.

## Usage

``` r
log_factorial_moment_theory(a, s, k)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - k:
    
    Moment order.

## Value

See details.
