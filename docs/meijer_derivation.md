# Closed form for the two restricted Meijer-G shapes this project needs

Every Meijer G evaluation in `htheory`/`htheoryR` is one of exactly two shapes
(confirmed by inspecting every call site in `background._f_N_scalar` and
`statfactor._B_scalar`):

```
G^{N,0}_{0,N}(z | -- ; b_1,...,b_N)      gamma chain,          s = 1/2
G^{0,N}_{N,0}(z | b_1,...,b_N ; --)      inverse-gamma chain,  s = 1
```

i.e. the two "one block entirely empty" cases of the general Meijer G. The
Python implementation treats these as opaque calls into `mpmath.meijerg`
(with a documented, hand-tuned underflow contract for the tails). The R port
instead derives a closed form for the first shape and gets the second for
free by reflection, which is what this note works out and checks against
identities the Python test suite already trusts.

## The reduction

Start from the Mellin-Barnes definition, specialized to `m = N, n = p = 0,
q = N` (no numerator Gammas at all):

```
G^{N,0}_{0,N}(z|--;b) = (1/2*pi*i) * Int_L  [ Prod_{j=1}^N Gamma(b_j + s) ]  z^{-s}  ds
```

This shape (`q > p`, entire in `1/z`) is closed by taking the contour to the
left and summing residues at the poles of each `Gamma(b_j + s)`, located at
`s = -b_j - k`, `k = 0, 1, 2, ...`, with residue `(-1)^k / k!`. For
pairwise-distinct `b_i`:

```
G = sum_{j=1}^N  sum_{k=0}^inf  [ (-1)^k / k! * Prod_{i != j} Gamma(b_i - b_j - k) ]  z^{b_j + k}
```

Using `Gamma(x - k) = (-1)^k Gamma(x) / (1 - x)_k` (Pochhammer in the
denominator) with `x = b_i - b_j`:

```
Prod_{i!=j} Gamma(b_i - b_j - k) = (-1)^{k(N-1)} [Prod_{i!=j} Gamma(b_i-b_j)] / Prod_{i!=j} (1+b_j-b_i)_k
```

so the inner sum over `k` collapses to a generalized hypergeometric series:

```
sum_k  [ ((-1)^N z)^k / k! ]  /  Prod_{i!=j} (1+b_j-b_i)_k
  = 0F_{N-1}( -- ; {1 + b_j - b_i}_{i!=j} ; (-1)^N z )
```

giving the closed form used by `htheoryR::g_gamma_chain`:

```
G^{N,0}_{0,N}(z|--;b) = sum_{j=1}^N  z^{b_j} * [Prod_{i!=j} Gamma(b_i-b_j)]
                                    * 0F_{N-1}( -- ; {1+b_j-b_i}_{i!=j} ; (-1)^N z )
```

Each term is an entire function of `z` (`0F_{N-1}` with `N-1 >= 1` denominator
parameters converges everywhere; for `N=1` it degenerates to `exp`), matching
the known fact that this shape (all parameters in the denominator block) is
entire and decays super-exponentially at large `z` -- the "dies at large z"
behavior `meijer.py`'s `UNDERFLOW_Z` cutoff is built around.

**Only valid as written for pairwise-distinct `b_i`.** Two `b_i` differing by
an integer (including 0) make `Gamma(b_i - b_j)` a pole at exactly the `k`
where the matching Pochhammer factor is zero -- a removable `0 * Inf`, the
same degeneracy that already breaks `mpmath`'s own hypergeometric expansion
(see `statfactor._B_powerlaw_series`'s docstring in the Python original,
which exists as a fallback for exactly this case). `htheoryR` regularizes it
by perturbing the colliding parameters by a tiny amount and evaluating at
boosted precision (`.g_core` in `R/meijer.R`) rather than deriving the
combinatorially heavier log-term confluent limit by hand.

## The mirror shape comes for free

`G^{0,N}_{N,0}(z|b;--)` has Mellin-Barnes integrand `Prod_j Gamma(1-b_j-s)
z^{-s}`. Substituting `s -> -s` (which reverses contour orientation, and
reversing it back costs a second sign flip that cancels) gives

```
G^{0,N}_{N,0}(z|b;--) = (1/2*pi*i) * Int_L [Prod_j Gamma((1-b_j)+s)] (1/z)^{-s} ds
                       = G^{N,0}_{0,N}(1/z | -- ; 1-b)
```

directly -- no separate derivation or implementation needed (this is DLMF
16.19.1 specialized to `p = m = 0`). `htheoryR::g_invgamma_chain` is a
one-line wrapper around `g_gamma_chain`.

## Verification against known identities (`tests/test_meijer.py`)

**N = 1** (`b = (b_1)`): the sum has one term, the product over `i != j` is
empty (`= 1`), and `0F_0(--;--;w) = e^w`, so

```
G^{1,0}_{0,1}(z|--;b_1) = z^{b_1} * e^{(-1)^1 z} = z^{b_1} e^{-z}
```

At `b_1 = 0` this is `e^{-z}`, matching `test_g_1_0_0_1_is_the_exponential`.

**N = 2** (`b = (a, 0)`): working through both terms and using
`I_nu(2*sqrt(z)) = [z^{nu/2}/Gamma(nu+1)] * 0F_1(--;nu+1;z)` plus the
reflection formula `Gamma(nu)Gamma(1-nu) = pi/sin(pi*nu)` collapses the two
terms to

```
G^{2,0}_{0,2}(z|--;a,0) = z^{a/2} * (pi/sin(a*pi)) * [I_{-a}(2 sqrt z) - I_a(2 sqrt z)]
                        = 2 * z^{a/2} * K_a(2 sqrt z)
```

using `K_nu = (pi/2) * (I_{-nu} - I_nu) / sin(nu*pi)` -- exactly the identity
`test_meijer.py` pins for the K-distribution special case. Both derivations
are carried out in full in the git history / session log around
2026-09-05; this file is the durable record.

**N >= 3**: not hand-verified symbolically (the algebra is the same
mechanism, just more terms). Cross-checked instead, numerically, against
`tests/testthat/fixtures/golden_htheory.csv`'s `a = [1.5, 2.0, 0.8]` row
(non-degenerate `b`) and the `a = [1.0, 1.0]` row (degenerate `b`, exercising
the perturbation fallback) -- see `tests/testthat/test-meijer.R`.
