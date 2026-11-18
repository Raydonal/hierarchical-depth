# Identifiability of the cascade depth `N` and shapes `a` from `T_N` alone.

The cascade is additive in log-temperature, and its log-cumulants of
order `r>=2` identify the multiset `{a_i}` – hence `N` – without
observing any intermediate level. See `docs/theorems.md` Theorem 3 (in
the parent project) for the full derivation; the short version: writing
the increment at level `i` as `D_i = log T_i - log T_{i-1}`, `D_i` is a
function of the fresh per-level draw `Y_i` alone (independent of
`T_{i-1}`), so the `D_i` are mutually independent and cumulants add.
