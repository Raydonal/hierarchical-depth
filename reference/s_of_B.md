# `s~(E_of_B(Bval))`, exact (not tabulated) composition.

`s~(E_of_B(Bval))`, exact (not tabulated) composition.

## Usage

``` r
s_of_B(Bval, a, s, beta0 = 1)
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

## Value

See details.
