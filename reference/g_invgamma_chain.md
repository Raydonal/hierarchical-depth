# `G^{0,N}_{N,0}(z | b_1,...,b_N ; --)` – the inverse-gamma-chain shape (`s = 1`), via the reflection identity `= G^{N,0}_{0,N}(1/z | -- ; 1-b)`.

Dies at *small* `z` (carries `e^{-1/z}`); pass `underflow_below`
(mirroring `meijer.UNDERFLOW_Z_SMALL`) to short-circuit to 0 below a
threshold – translated internally to the `underflow_above` threshold
`1/underflow_below` on the reflected argument `1/z`.

## Usage

``` r
g_invgamma_chain(z, b, prec_bits = DEFAULT_PREC_BITS, underflow_below = NULL)
```

## Arguments

  - z:
    
    Positive numeric scalar.

  - b:
    
    Numeric vector of length N (the lower/denominator parameters).

  - prec\_bits:
    
    Working precision in bits (default \~38 decimal digits).

  - underflow\_below:
    
    Optional threshold; `z` at or below it returns 0.
