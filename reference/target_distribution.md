# `p_i = B(E_i)/Z` with `Z = sum_i B(E_i)`. Returns `list(p, Z, U)` with `U = sum_i p_i E_i` the mean energy the maximization is constrained to.

`p_i = B(E_i)/Z` with `Z = sum_i B(E_i)`. Returns `list(p, Z, U)` with
`U = sum_i p_i E_i` the mean energy the maximization is constrained to.

## Usage

``` r
target_distribution(energies, a, s, beta0 = 1)
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

## Value

See details.
