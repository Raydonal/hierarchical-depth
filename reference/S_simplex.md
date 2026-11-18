# `S[p] = sum_i s_of_B(Z*p_i)/Z + c*sum(p_i)`, Eq. (10). `Z` (the partition function) is mandatory – it is supplied externally, not recomputed here.

`S[p] = sum_i s_of_B(Z*p_i)/Z + c*sum(p_i)`, Eq. (10). `Z` (the
partition function) is mandatory – it is supplied externally, not
recomputed here.

## Usage

``` r
S_simplex(p, a, s, beta0 = 1, Z, c = 0)
```

## Arguments

  - p:
    
    Probability vector on the simplex.

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

  - c:
    
    Additive constant term.

## Value

See details.
