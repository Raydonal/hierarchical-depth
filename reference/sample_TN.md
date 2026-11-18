# Ancestral sampling of `T_N` down the chain, Eq. (30).

`T_0 = T0` is fixed; each level is drawn from its conditional given the
level above. The empirical law of the result is Eq. (33) – this is the
Monte-Carlo check that the corrected branch assignment is the right one.

## Usage

``` r
sample_TN(a, s, T0 = 1, size = 1)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - T0:
    
    Base (outermost) reservoir temperature.

  - size:
    
    Number of draws.

## Value

See details.
