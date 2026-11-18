# Parametric-bootstrap test of `H0: depth = N0` (homogeneous `a`) against the sample `T`. Fits `a` under `H0` from `kappa2` alone, so the test statistic – observed `kappa3` minus the value predicted at `(N0, a_hat)` – is not tautologically zero. Its null distribution is a parametric bootstrap: simulate `n_boot` depth-`N0` cascades at the fitted shape and refit on each. `N` sits at the boundary of its domain (a positive integer), so the usual chi-square asymptotics for a likelihood-ratio test do not apply (Self & Liang 1987), hence bootstrap rather than a closed-form p-value.

Parametric-bootstrap test of `H0: depth = N0` (homogeneous `a`) against
the sample `T`. Fits `a` under `H0` from `kappa2` alone, so the test
statistic – observed `kappa3` minus the value predicted at `(N0, a_hat)`
– is not tautologically zero. Its null distribution is a parametric
bootstrap: simulate `n_boot` depth-`N0` cascades at the fitted shape and
refit on each. `N` sits at the boundary of its domain (a positive
integer), so the usual chi-square asymptotics for a likelihood-ratio
test do not apply (Self & Liang 1987), hence bootstrap rather than a
closed-form p-value.

## Usage

``` r
test_depth(T, N0, s, n_boot = 2000L, T0 = 1)
```

## Arguments

  - T:
    
    Sample of the reservoir temperature (or, in `conditional_pdf`, the
    value(s) at which to evaluate the density).

  - N0:
    
    The depth under the null hypothesis `H0: N = N0`.

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - n\_boot:
    
    Number of bootstrap replications.

  - T0:
    
    Base (outermost) reservoir temperature.

## Value

See details.
