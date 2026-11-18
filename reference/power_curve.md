# Simulated power of `test_depth` at `H0: N = N0` against a true cascade `a_true` (its length is the true depth). Returns a named list `n -> power`: reject rate under a true deeper cascade, at a shallow null.

Simulated power of `test_depth` at `H0: N = N0` against a true cascade
`a_true` (its length is the true depth). Returns a named list `n ->
power`: reject rate under a true deeper cascade, at a shallow null.

## Usage

``` r
power_curve(
  a_true,
  s,
  N0,
  sample_sizes = c(100, 300, 1000, 5000),
  n_reps = 200L,
  n_boot = 400L,
  alpha = 0.05,
  T0 = 1
)
```

## Arguments

  - a\_true:
    
    True (possibly heterogeneous) shape vector generating the data; its
    length is the true depth.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - N0:
    
    The depth under the null hypothesis `H0: N = N0`.

  - sample\_sizes:
    
    Numeric vector of sample sizes to evaluate power at.

  - n\_reps:
    
    Number of Monte Carlo replications per sample size.

  - n\_boot:
    
    Number of bootstrap replications.

  - alpha:
    
    Nominal test significance level.

  - T0:
    
    Base (outermost) reservoir temperature.

## Value

A named list mapping each sample size (as a character string) to its
estimated power.
