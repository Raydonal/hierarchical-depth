# `B_N^{(s)}(E; a)`, Eq. (34).

`omega` defaults to `prod(a_i)` but can be overridden – the whole point
of the signature: Eq. (38) shifts the vector `a` while holding `omega`
at its *unshifted* value, so `B_tilde` must pass the two independently.
Recomputing `omega` from the shifted vector is the single easiest
mistake to make here (see `docs/meijer_derivation.md` in the
repository).

## Usage

``` r
B(E, a, s, beta0 = 1, omega = NULL)
```

## Arguments

  - E:
    
    Energy: a nonnegative scalar or numeric vector.

  - a:
    
    Numeric vector of per-level shape parameters `(a_1, ..., a_N)`;
    length is the cascade depth `N` (empty for `N = 0`).

  - s:
    
    Universality class: `0.5` (gamma / stretched-exponential) or `1`
    (inverse-gamma / power-law).

  - beta0:
    
    Inverse base temperature, `1/T0`.

  - omega:
    
    Optional override for `omega = prod(a_i)`, decoupled from `a`
    (needed by `B_tilde`'s shifted vector).

## Value

See details.
