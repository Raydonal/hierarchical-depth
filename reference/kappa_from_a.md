# Invert `a_i = 2*gamma_i/kappa_i^2` for the noise amplitudes.

Invert `a_i = 2*gamma_i/kappa_i^2` for the noise amplitudes.

## Usage

``` r
kappa_from_a(a, gamma)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - gamma:
    
    Vector of relaxation rates (Eq. 29).

## Value

See details.
