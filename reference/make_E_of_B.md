# Convenience wrapper returning just the tabulated `E_of_B` function.

Convenience wrapper returning just the tabulated `E_of_B` function.

## Usage

``` r
make_E_of_B(a, s, beta0 = 1, n = 300L, B_min = 1e-12)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - beta0:
    
    Inverse base temperature, `1/T0`.

  - n:
    
    Number of grid points (or, in `sample_TN`/`power_curve`, a sample
    size).

  - B\_min:
    
    Smallest `B` value to tabulate down to.

## Value

See details.
