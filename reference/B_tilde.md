# `Btilde_N^{(s)}(E) = beta0 * int_E^inf B_N(E') dE'`, Eqs. (37)-(38).

The same function with `a -> a + 3 - 4s` and `omega` untouched: a
size-biasing of the reservoir temperature, and size-biasing moves a
gamma's shape while leaving its rate alone – which is exactly why
`omega` does not move.

## Usage

``` r
B_tilde(E, a, s, beta0 = 1)
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
