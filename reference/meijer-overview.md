# Robust, arbitrary-precision evaluation of the two Meijer-G shapes this package needs.

Every closed form in the paper is one of exactly two restricted Meijer-G
shapes (confirmed by inspecting every call site in the Python original,
`background._f_N_scalar` and `statfactor._B_scalar`):

## Details

``` 
  G^{N,0}_{0,N}(z | -- ; b_1,...,b_N)   (gamma chain,        s = 1/2)
  G^{0,N}_{N,0}(z | b_1,...,b_N ; --)   (inverse-gamma chain, s = 1)
```

Both are "numerator-empty"/"denominator-empty" cases of the general
Meijer G, which is why they admit a closed-form reduction to a finite
sum of generalized hypergeometric terms rather than needing a black-box
Mellin-Barnes routine (R has no arbitrary-precision *complex* arithmetic
– `Rmpfr` wraps MPFR, real only – which is the concrete reason contour
quadrature was rejected as the implementation strategy; see the plan).

Derivation (residue sum of the Mellin-Barnes integral, closing the
contour at the simple poles of each Gamma(b\_j + s)): for
pairwise-distinct b\_i,

``` 
  G^{N,0}_{0,N}(z|--;b) = sum_j  z^{b_j} * prod_{i!=j} Gamma(b_i - b_j)
                               * 0F_{N-1}(--; {1+b_j-b_i}_{i!=j}; (-1)^N z)
```

Verified by hand against two independent, oracle-free identities already
pinned in the Python test suite (`tests/test_meijer.py`):

``` 
  N=1: reduces to z^{b_1} e^{-z}  (matches G^{1,0}_{0,1}(z|--;0) = e^{-z})
  N=2: reduces to 2 z^{a/2} K_a(2 sqrt(z))  (the b=(a,0) Bessel-K identity)
```

Both reductions are worked out in full in docs/meijer\_derivation.md.

The mirror shape needs no separate derivation: substituting s -\> -s in
its own Mellin-Barnes integral gives

``` 
  G^{0,N}_{N,0}(z|b;--) = G^{N,0}_{0,N}(1/z|--;1-b)
```

directly (this is DLMF 16.19.1 specialized to p=m=0), so only one
routine (`.g_core`) is implemented and hardened.
