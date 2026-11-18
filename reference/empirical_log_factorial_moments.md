# `(lambda_hat, c(L_2,...,L_k_max))` with `L_k = log(Ehat[W^k])`. `lambda_hat` is the sample mean of `D` (exact since `E[W]=1`). The falling-factorial product is exactly 0 whenever `D < k`, so no special-casing of small degrees is needed.

`(lambda_hat, c(L_2,...,L_k_max))` with `L_k = log(Ehat[W^k])`.
`lambda_hat` is the sample mean of `D` (exact since `E[W]=1`). The
falling-factorial product is exactly 0 whenever `D < k`, so no
special-casing of small degrees is needed.

## Usage

``` r
empirical_log_factorial_moments(D, k_max = 3L)
```

## Arguments

  - D:
    
    Integer (or numeric) vector, an observed degree sequence.

  - k\_max:
    
    Highest factorial-moment order to compute.

## Value

See details.
