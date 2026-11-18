# Unbiased sample cumulants `k_r` of `log T` for `r=2..r_max`, via the standard k-statistic (directly comparable to `log_cumulants_theory` with no finite-sample bias correction).

Unbiased sample cumulants `k_r` of `log T` for `r=2..r_max`, via the
standard k-statistic (directly comparable to `log_cumulants_theory` with
no finite-sample bias correction).

## Usage

``` r
empirical_log_cumulants(T, r_max = 4L)
```

## Arguments

  - T:
    
    Sample of the reservoir temperature (or, in `conditional_pdf`, the
    value(s) at which to evaluate the density).

  - r\_max:
    
    Highest cumulant order to compute (cumulants of order `2..r_max` are
    returned).

## Value

See details.
