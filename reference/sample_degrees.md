# Draw a mixed-Poisson degree sequence, `D_v | W_v ~ Poisson(lambda*W_v)`. `W` is drawn from `sample_TN` with `T0 = 1`, so `E[W] = 1` and `E[D] = lambda` exactly.

Draw a mixed-Poisson degree sequence, `D_v | W_v ~ Poisson(lambda*W_v)`.
`W` is drawn from `sample_TN` with `T0 = 1`, so `E[W] = 1` and `E[D] =
lambda` exactly.

## Usage

``` r
sample_degrees(a, s, lam, size, T0 = 1)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - lam:
    
    Poisson rate parameter.

  - size:
    
    Number of draws.

  - T0:
    
    Base (outermost) reservoir temperature.

## Value

See details.
