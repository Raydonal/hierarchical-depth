# Hierarchical random networks: the mixed-Poisson application.

Each vertex `v` belongs to a nested sequence of communities, and its
propensity to connect is the product of per-level contributions, `W_v =
prod_i X_{i,c_i(v)}` – exactly the multiplicative cascade of
`sample_TN`, relabelled. Degrees are then mixed Poisson, `D_v | W_v ~
Poisson(lambda*W_v)` (Chung-Lu / Norros-Reittu).

## Details

The identifiability chain: (1) mixed-Poisson laws are identifiable from
the degree distribution alone (Feller 1943) – the law of `W` is
recoverable from `D`; (2) by Theorem 3 (`R/identifiability.R`), `(N,
{a_i})` is identifiable from the law of `W`; (3) hence the cascade depth
is estimable from the degree sequence alone, without observing which
vertex belongs to which community.

The estimator is the factorial-moment analogue of the log-cumulant one:
`E[D^(k)] = lambda^k * E[W^k]` exactly (a mixed-Poisson identity, not an
approximation), and `E[W^k]` has a closed form linear in `N` for fixed
`a`, so the same ratio-then-scale trick as `fit_homogeneous` applies.
