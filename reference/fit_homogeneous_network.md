# Recover `(N, a, lambda)` from a degree sequence, homogeneous model. Mirrors `fit_homogeneous`: the ratio `L_3/L_2` depends on `a` alone (the common `N` factor cancels), a 1-D root find, then `N = L_2 / [single-level log E[W^2]](a)` from the scale.

Recover `(N, a, lambda)` from a degree sequence, homogeneous model.
Mirrors `fit_homogeneous`: the ratio `L_3/L_2` depends on `a` alone (the
common `N` factor cancels), a 1-D root find, then `N = L_2 /
[single-level log E[W^2]](a)` from the scale.

## Usage

``` r
fit_homogeneous_network(D, s, k_max = 3L, a_bracket = NULL)
```

## Arguments

  - D:
    
    Integer (or numeric) vector, an observed degree sequence.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - k\_max:
    
    Highest factorial-moment order to compute.

  - a\_bracket:
    
    Two-element vector, the search bracket for `a` in the root-find.

## Value

See details.
