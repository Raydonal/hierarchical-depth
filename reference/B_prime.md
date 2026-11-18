# `dB/dE`, Eq. (43): strictly negative, hence `s''(B) < 0`, Eq. (44).

`dB/dE`, Eq. (43): strictly negative, hence `s''(B) < 0`, Eq. (44).

## Usage

``` r
B_prime(E, a, s, beta0 = 1, omega = NULL)
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

  - omega:
    
    Optional override for `omega = prod(a_i)`, decoupled from `a`
    (needed by `B_tilde`'s shifted vector).

## Value

See details.
