# Draw `T_i ~ f(T_i | T_{i-1})`, Eqs. (31)/(32).

Unlike the Python original (which threads an explicit `rng` object
through every call for reproducibility), this uses R's global RNG state
directly (idiomatic R – callers control reproducibility with
`set.seed()`). Bit- identical streams with the Python version are
unreachable regardless (R's Mersenne-Twister differs from numpy's
PCG64), so nothing is lost.

## Usage

``` r
sample_conditional(T_prev, a_i, s, size = length(T_prev))
```

## Arguments

  - T\_prev:
    
    Parent-level temperature, scalar or vector.

  - a\_i:
    
    Scalar shape parameter for this level.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - size:
    
    Number of draws.

## Value

See details.
