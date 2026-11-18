# The exact Tsallis entropy of a probability vector `p`, at entropy index `q = a/(a+1)`.

The exact Tsallis entropy of a probability vector `p`, at entropy index
`q = a/(a+1)`.

## Usage

``` r
tsallis_entropy(p, a)
```

## Arguments

  - p:
    
    Probability vector on the simplex.

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

## Value

See details.
