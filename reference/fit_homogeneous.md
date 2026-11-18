# Point estimate of `(N, a)` under the homogeneous model `a_i = a`, from `(kappa2, kappa3)`. Two equations, two unknowns: `kappa3/kappa2` depends on `a` alone (the `N` factor cancels), recovered by a 1-D root find; then `N = kappa2 / psi'(a)` from the scale.

Point estimate of `(N, a)` under the homogeneous model `a_i = a`, from
`(kappa2, kappa3)`. Two equations, two unknowns: `kappa3/kappa2` depends
on `a` alone (the `N` factor cancels), recovered by a 1-D root find;
then `N = kappa2 / psi'(a)` from the scale.

## Usage

``` r
fit_homogeneous(kappa2, kappa3, s, a_bracket = c(1e-04, 10000))
```

## Arguments

  - kappa2:
    
    Observed (or theoretical) second log-cumulant.

  - kappa3:
    
    Observed (or theoretical) third log-cumulant.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - a\_bracket:
    
    Two-element vector, the search bracket for `a` in the root-find.

## Value

See details.
