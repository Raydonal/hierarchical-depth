# `f(T_i | T_{i-1})`, Eq. (31) for `s = 1/2` and Eq. (32) for `s = 1`.

Both have mean `T_{i-1}`: the linear drift of Eq. (29) drives every
level to the temperature of the one above it.

## Usage

``` r
conditional_pdf(T, T_prev, a_i, s)
```

## Arguments

  - T:
    
    Sample of the reservoir temperature (or, in `conditional_pdf`, the
    value(s) at which to evaluate the density).

  - T\_prev:
    
    Parent-level temperature, scalar or vector.

  - a\_i:
    
    Scalar shape parameter for this level.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

## Value

See details.
