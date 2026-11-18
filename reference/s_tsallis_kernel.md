# `s(B) = -(a+1)*(B - B^q)`, Eq. (63) – the `N=1`, `s=1` exactly-solvable case, `q = a/(a+1)` the entropy index. Summed over the simplex this is exactly the Tsallis entropy, Eqs. (64)/(67).

`s(B) = -(a+1)*(B - B^q)`, Eq. (63) – the `N=1`, `s=1` exactly-solvable
case, `q = a/(a+1)` the entropy index. Summed over the simplex this is
exactly the Tsallis entropy, Eqs. (64)/(67).

## Usage

``` r
s_tsallis_kernel(Bval, a)
```

## Arguments

  - Bval:
    
    Value(s) of the statistical factor `B`, in `(0, 1]`.

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

## Value

See details.
