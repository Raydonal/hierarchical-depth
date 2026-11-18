# `E` such that `B(E) = Bval`, vectorized over `Bval`. Well-posed for every `(N, a)` and both classes because `B` is strictly decreasing from `B(0) = 1` to `B(inf) = 0`, Eqs. (7)/(43).

`E` such that `B(E) = Bval`, vectorized over `Bval`. Well-posed for
every `(N, a)` and both classes because `B` is strictly decreasing from
`B(0) = 1` to `B(inf) = 0`, Eqs. (7)/(43).

## Usage

``` r
E_of_B(Bval, a, s, beta0 = 1, xtol = 1e-12)
```

## Arguments

  - Bval:
    
    Value(s) of the statistical factor `B`, in `(0, 1]`.

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - beta0:
    
    Inverse base temperature, `1/T0`.

  - xtol:
    
    Root-finding tolerance.

## Value

See details.
