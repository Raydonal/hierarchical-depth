# Parametric-bootstrap test of `H0: N = N0` (homogeneous `a`) against a degree sequence `D` – the factorial-moment analogue of `test_depth`. `fit_homogeneous_network`'s point estimate has no formal p-value attached; this closes that gap the same way `test_depth` does for a propensity sequence. A degree sequence whose `L3/L2` ratio cannot be matched by any `a` in the bracket is a genuine model-selection signal (raises, like `fit_homogeneous_network`), not something papered over.

Parametric-bootstrap test of `H0: N = N0` (homogeneous `a`) against a
degree sequence `D` – the factorial-moment analogue of `test_depth`.
`fit_homogeneous_network`'s point estimate has no formal p-value
attached; this closes that gap the same way `test_depth` does for a
propensity sequence. A degree sequence whose `L3/L2` ratio cannot be
matched by any `a` in the bracket is a genuine model-selection signal
(raises, like `fit_homogeneous_network`), not something papered over.

## Usage

``` r
test_depth_degrees(D, N0, s, n_boot = 2000L, k_max = 3L)
```

## Arguments

  - D:
    
    Integer (or numeric) vector, an observed degree sequence.

  - N0:
    
    The depth under the null hypothesis `H0: N = N0`.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - n\_boot:
    
    Number of bootstrap replications.

  - k\_max:
    
    Highest factorial-moment order to compute.

## Value

See details.
