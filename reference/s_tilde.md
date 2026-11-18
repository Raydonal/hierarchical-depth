# `s~_N^{(s)}(E) = (beta0*E - 1)*B(E) + Btilde(E)`, Eq. (39).

`s~_N^{(s)}(E) = (beta0*E - 1)*B(E) + Btilde(E)`, Eq. (39).

## Usage

``` r
s_tilde(E, a, s, beta0 = 1)
```

## Arguments

  - E:
    
    Energy: a nonnegative scalar or numeric vector.

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
