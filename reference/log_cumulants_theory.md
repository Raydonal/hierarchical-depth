# Theoretical `kappa_r(log T_N)` for `r=2..r_max`. Independent of `T_0` by construction; relies only on the shape vector `a` and branch `s`. `a` may be empty (`N = 0`), for which all cumulants are 0.

Theoretical `kappa_r(log T_N)` for `r=2..r_max`. Independent of `T_0` by
construction; relies only on the shape vector `a` and branch `s`. `a`
may be empty (`N = 0`), for which all cumulants are 0.

## Usage

``` r
log_cumulants_theory(a, s, r_max = 4L)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - r\_max:
    
    Highest cumulant order to compute (cumulants of order `2..r_max` are
    returned).

## Value

See details.
