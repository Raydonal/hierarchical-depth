# The parametric curve `(B(E), s~(E))` – geometric in `1 + beta0*E` so both the `B ~ 1` knee and the tail are resolved. Outside the three solvable cases `E(B)` has no closed form, so `s(B)` exists only as this curve or through `E_of_B`.

The parametric curve `(B(E), s~(E))` – geometric in `1 + beta0*E` so
both the `B ~ 1` knee and the tail are resolved. Outside the three
solvable cases `E(B)` has no closed form, so `s(B)` exists only as this
curve or through `E_of_B`.

## Usage

``` r
parametric_curve(a, s, beta0 = 1, E = NULL, n = 400L, E_max = 60)
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

  - E:
    
    Energy: a nonnegative scalar or numeric vector.

  - n:
    
    Number of grid points (or, in `sample_TN`/`power_curve`, a sample
    size).

  - E\_max:
    
    See details.

## Value

See details.
