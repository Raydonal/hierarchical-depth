# Build monotone-cubic interpolants of `E(B)` and `s~(B)`, amortizing the cost of repeated `B()`/`E_of_B()` evaluation (each one, worst case, an escalating-precision Meijer-G or quadrature call – see R/meijer.R and R/statfactor.R) over a table built once. Two curves are tabulated, both against `log(B)`, NOT `E` or `B` directly:

  - `log(1 + beta0*E)` vs `log(B)`, inverted to give `E_of_B`;

  - the *bounded ratio* `u(B) = s~(B) / (B*(1+beta0*E))` vs `log(B)`,
    not `s~` itself – `s~ -> 0` and `B -> 0` together in the tail, so
    splining `s~` directly loses exactly the relative accuracy there
    that a downstream `maxent` optimizer needs; `u(B)` stays `O(1)`
    across the whole tail by construction, which is what actually gets
    splined well.

## Usage

``` r
make_interpolants(a, s, beta0 = 1, n = 300L, B_min = 1e-12)
```

## Arguments

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - beta0:
    
    Inverse base temperature, `1/T0`.

  - n:
    
    Number of grid points (or, in `sample_TN`/`power_curve`, a sample
    size).

  - B\_min:
    
    Smallest `B` value to tabulate down to.

## Value

See details.

## Details

Returns `list(E_of_B = function(B), s_of_B = function(B))`. Queries
outside the tabulated `log(B)` range are clamped to the nearest
tabulated endpoint rather than extrapolated.
