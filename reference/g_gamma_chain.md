# `G^{N,0}_{0,N}(z | -- ; b_1,...,b_N)` – the gamma-chain shape (`s = 1/2`).

Dies at *large* `z` (carries `e^{-z}`): pass `underflow_above`
(mirroring `meijer.UNDERFLOW_Z` in the Python original) to short-circuit
to exactly 0 past a threshold where the true value is confidently
negligible, rather than attempting the direct series there (see
`.g_core_raw`'s docs for why that regime is numerically treacherous, not
just slow).

## Usage

``` r
g_gamma_chain(z, b, prec_bits = DEFAULT_PREC_BITS, underflow_above = NULL)
```

## Arguments

  - z:
    
    Positive numeric scalar.

  - b:
    
    Numeric vector of length N (the lower/denominator parameters).

  - prec\_bits:
    
    Working precision in bits (default \~38 decimal digits).

  - underflow\_above:
    
    Optional threshold; `z` at or above it returns 0.

## Value

An `mpfr` number.
