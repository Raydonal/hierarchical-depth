# Closing the loop: maximize `S[p]` and get `p_i = B(E_i)/Z` back.

This is the paper's central claim: the chain `B -> s~(E) -> s(B) ->
S[p]` reconstructs an entropy from a prescribed statistical factor, and
maximizing that entropy under the standard constraints (normalization
and mean energy, no escort probabilities) returns the distribution
started from. This closes every link at once: the branch assignment, the
`omega`-versus-`a` bookkeeping, the gauge `alpha=-1`, and the `1/Z`
scaling – the highest-value regression test in the package.

## Details

The gradient is exact and cheap: differentiating `S[p] = sum_j Z^{-1}
s(Z p_j)` gives `s'(Z p_i) = alpha + beta0*E(Z p_i)`, so the derivative
of the entropy *is* the inverse function `E(B)` – no numerical
differentiation needed.
