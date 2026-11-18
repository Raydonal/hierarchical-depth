# Spread of `s'(Z p_i) - beta0*E_i` across states: this is the constant `alpha` at the maximizer, so the spread is zero there. Reported as `max - min`: a direct, optimizer-free check, and the one that localizes a broken link fastest.

Spread of `s'(Z p_i) - beta0*E_i` across states: this is the constant
`alpha` at the maximizer, so the spread is zero there. Reported as `max
- min`: a direct, optimizer-free check, and the one that localizes a
broken link fastest.

## Usage

``` r
stationarity_residual(p, energies, a, s, beta0 = 1, Z = NULL)
```

## Arguments

  - p:
    
    Probability vector on the simplex.

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

  - Z:
    
    Partition function (normalizing constant).

## Value

See details.
