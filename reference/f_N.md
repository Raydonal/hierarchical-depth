# The reservoir temperature law `f_N(T_N)`, Eq. (33), corrected.

Normalizes to 1 and has mean `T0` for every `N`, `a` and both classes –
the two facts the sum rule Eq. (27) rests on.

## Usage

``` r
f_N(T, a, s, T0 = 1)
```

## Arguments

  - T:
    
    Sample of the reservoir temperature (or, in `conditional_pdf`, the
    value(s) at which to evaluate the density).

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - T0:
    
    Base (outermost) reservoir temperature.

## Value

See details.
