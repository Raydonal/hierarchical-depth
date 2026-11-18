#' The generalized statistical factor `B_N^{(s)}(E; a)`, Eq. (34).
#'
#' `B_N(E)` is the Boltzmann factor averaged over the innermost reservoir
#' temperature, Eq. (24): `B_N(E) = int_0^inf exp(-E/T_N) f_N(T_N) dT_N`.
#'
#' Two invariants anchor the package (checked in `tests/testthat/test-statfactor.R`):
#' the flat-energy sum rule, Eq. (27), and the boundary conditions `B(0) = 1`,
#' `B(inf) = 0`, `B'(E) < 0`, Eq. (7)/(43).
#' @keywords internal
#' @name statfactor-overview
NULL

#' `omega = prod(a_i)`, Eq. (33). The scalar Eq. (38) leaves alone.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @return See details.
#' @export
omega_of <- function(a) .omega(.as_a(a))

#' The parameter shift of Eq. (38): `a -> a + 3 - 4s`.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @return See details.
#' @export
shift <- function(s) .shift(s)

#' `B_N^{(1)}` from the inverse moments of `T_N`, for small `beta0*E`.
#'
#' For the inverse-gamma cascade `T_N = omega*T0/prod(Y_i)` with
#' `Y_i ~ Gamma(a_i+1, 1)`, the *inverse* moments are elementary:
#' `E[T_N^{-k}] = prod_i Gamma(a_i+1+k)/Gamma(a_i+1) / (omega*T0)^k`, so
#' expanding `exp(-E/T_N)` term by term gives a series in `z = beta0*E/omega`
#' whose coefficients grow like `(k!)^N` -- asymptotic, not convergent: terms
#' shrink to a minimum then diverge. Summed to that minimum it is exact at
#' small `z`, which is precisely the regime where the Meijer G's own
#' hypergeometric expansion is degenerate (`a_i` differing by integers).
#'
#' Returns `list(value, converged)`; `converged = FALSE` once terms start
#' growing again, and the caller must discard the value in favor of the
#' fallback.
#' @noRd
.B_powerlaw_series <- function(E, a, beta0, omega, prec_bits = DEFAULT_PREC_BITS,
                                tol = 1e-18, kmax = 80L) {
  z <- Rmpfr::mpfr(beta0, prec_bits) * Rmpfr::mpfr(E, prec_bits) / omega
  term <- Rmpfr::mpfr(1, prec_bits)
  total <- Rmpfr::mpfr(1, prec_bits)
  a_mp <- Rmpfr::mpfr(a, prec_bits)
  smallest <- NULL
  for (k in seq_len(kmax - 1)) {
    term <- term * (-z / k)
    for (idx in seq_along(a_mp)) term <- term * (a_mp[idx] + k)
    total <- total + term
    size <- abs(term)
    if (as.numeric(size) < tol * as.numeric(abs(total))) return(list(value = total, converged = TRUE))
    if (!is.null(smallest) && as.numeric(size) > smallest) return(list(value = total, converged = FALSE))
    smallest <- as.numeric(size)
  }
  list(value = total, converged = FALSE)
}

#' `f_N`'s scalar core, but taking `omega` explicitly instead of recomputing
#' it from `a` -- needed because `statfactor.B`'s `omega` can be decoupled
#' from `a` (Eq. 38's shifted vector), the same reason `.B` below takes an
#' `omega` override. Used only by `.B_invgamma_fallback`.
#' @noRd
.f_N_omega_scalar <- function(T, a, s, T0, omega, prec_bits = DEFAULT_PREC_BITS) {
  if (T <= 0) return(0)
  if (length(a) == 0) stop("f_N is a delta at T0 for N = 0; it has no density")
  if (isTRUE(all.equal(s, 0.5))) {
    pref <- Rmpfr::mpfr(omega / T0, prec_bits)
    for (x in a) pref <- pref / gamma(Rmpfr::mpfr(x, prec_bits))
    g <- g_gamma_chain(omega * T / T0, a - 1, prec_bits, underflow_above = UNDERFLOW_Z)
    return(pref * g)
  }
  pref <- Rmpfr::mpfr(1 / (T0 * omega), prec_bits)
  for (x in a) pref <- pref / gamma(Rmpfr::mpfr(x, prec_bits) + 1)
  g <- g_invgamma_chain(T / (omega * T0), -a - 1, prec_bits, underflow_below = UNDERFLOW_Z_SMALL)
  pref * g
}

#' Fallback for `B_N^{(1)}(E)` past the power series' reach: direct
#' integration of Eq. (24), `int_0^inf exp(-beta0*E/T) f_N(T) dT`, using the
#' already-verified `f_N` core rather than deriving the mixed Meijer-G shape
#' `G^{1,N}_{N,1}(z|-a;0)` the Python original falls back to (that shape has
#' both numerator *and* denominator Gamma factors -- unlike the two "one
#' block empty" shapes `.g_core` handles, it doesn't reduce the same way, and
#' re-deriving it by hand was judged not worth the risk of a second unchecked
#' derivation; see the plan). This is exact (it's the defining integral, Eq.
#' 24) and needs no extra derivation, at the cost of ordinary double-precision
#' quadrature rather than arbitrary precision -- adequate here because the
#' Python original's own `B_quadrature` cross-check is *also* plain
#' double-precision numpy/scipy, never mpmath.
#'
#' Substitutes `T = exp(u)` to map `(0, Inf)` to all of the real line, which
#' both removes the singular endpoint at `T = 0` and turns the double-
#' exponential-decaying tails into a smooth, `integrate()`-friendly integrand.
#' @noRd
.B_invgamma_fallback <- function(E, a, beta0, omega, T0 = 1.0) {
  integrand <- function(u) {
    Tval <- exp(u)
    dens <- vapply(Tval, function(t) {
      as.numeric(.f_N_omega_scalar(t, a, 1.0, T0, omega))
    }, numeric(1))
    exp(-beta0 * E / Tval) * dens * Tval
  }
  # f_N has mean T0 and a power-law tail for s = 1 (Eq. 35): integrate over a
  # wide log-T window (T from ~1e-30*T0 to ~1e30*T0) to capture both the
  # essential-singularity decay at T -> 0 and the polynomial tail at T -> Inf.
  # rel.tol = 1e-6 (not tighter): the eventual output is cast to double
  # anyway, and each integrand evaluation can itself be an expensive
  # degenerate-parameter Meijer G call (see R/meijer.R) -- a tight quadrature
  # tolerance multiplies that cost across many more nodes for a precision
  # gain the double-precision return type can't even represent.
  res <- stats::integrate(integrand, lower = log(1e-30 * T0), upper = log(1e30 * T0),
                           rel.tol = 1e-6, subdivisions = 200L)
  res$value
}

.B_scalar <- function(E, a, s, beta0, omega, prec_bits = DEFAULT_PREC_BITS) {
  if (E < 0) stop(sprintf("energy must be nonnegative, got %s", E))
  if (E == 0) return(1.0)  # Eq. (7): the z -> 0 limit of Eq. (34) is exactly 1.
  if (length(a) == 0) return(exp(-beta0 * E))  # N = 0: Eq. (36), no special-casing
                                                 # needed for s -- both series and
                                                 # quadrature fallbacks degenerate
                                                 # awkwardly at an empty reservoir
                                                 # chain, so short-circuit directly
                                                 # to the plain Boltzmann factor.
  if (isTRUE(all.equal(s, 0.5))) {
    # G^{N+1,0}_{0,N+1}( beta0*omega*E | -- ; a, 0 ) / Gamma(a)
    pref <- Rmpfr::mpfr(1, prec_bits)
    for (x in a) pref <- pref / gamma(Rmpfr::mpfr(x, prec_bits))
    b_upper <- c(a, 0)
    g <- g_gamma_chain(beta0 * omega * E, b_upper, prec_bits, underflow_above = UNDERFLOW_Z)
    return(as.numeric(pref * g))
  }
  # s = 1: small argument first via the exact series; past its reach, the
  # defining integral (see .B_invgamma_fallback) rather than a hand-derived
  # mixed Meijer-G shape.
  series <- .B_powerlaw_series(E, a, beta0, omega, prec_bits)
  if (series$converged) return(as.numeric(series$value))
  .B_invgamma_fallback(E, a, beta0, omega)
}

#' `B_N^{(s)}(E; a)`, Eq. (34).
#'
#' `omega` defaults to `prod(a_i)` but can be overridden -- the whole point
#' of the signature: Eq. (38) shifts the vector `a` while holding `omega` at
#' its *unshifted* value, so `B_tilde` must pass the two independently.
#' Recomputing `omega` from the shifted vector is the single easiest mistake
#' to make here (see `docs/meijer_derivation.md` in the repository).
#' @param E Energy: a nonnegative scalar or numeric vector.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param omega Optional override for `omega = prod(a_i)`, decoupled from `a` (needed by `B_tilde`'s shifted vector).
#' @return See details.
#' @export
B <- function(E, a, s, beta0 = 1.0, omega = NULL) {
  a <- .as_a(a)
  .B(E, a, s, beta0, omega)
}

.B <- function(E, a, s, beta0, omega) {
  if (is.null(omega)) omega <- .omega(a)
  vectorize_scalar(.B_scalar)(E, a, s, beta0, omega)
}

#' `Btilde_N^{(s)}(E) = beta0 * int_E^inf B_N(E') dE'`, Eqs. (37)-(38).
#'
#' The same function with `a -> a + 3 - 4s` and `omega` untouched: a
#' size-biasing of the reservoir temperature, and size-biasing moves a
#' gamma's shape while leaving its rate alone -- which is exactly why
#' `omega` does not move.
#' @param E Energy: a nonnegative scalar or numeric vector.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @return See details.
#' @export
B_tilde <- function(E, a, s, beta0 = 1.0) {
  a <- .as_a(a)
  .B(E, a + .shift(s), s, beta0, .omega(a))
}

#' `dB/dE`, Eq. (43): strictly negative, hence `s''(B) < 0`, Eq. (44).
#' @param E Energy: a nonnegative scalar or numeric vector.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param omega Optional override for `omega = prod(a_i)`, decoupled from `a` (needed by `B_tilde`'s shifted vector).
#' @return See details.
#' @export
B_prime <- function(E, a, s, beta0 = 1.0, omega = NULL) {
  a <- .as_a(a)
  if (is.null(omega)) omega <- .omega(a)
  d_scalar <- function(e) {
    h <- max(1e-6, 1e-6 * e)
    if (e - h < 0) {
      # Forward difference near E = 0 (a central stencil would dip below the
      # domain, where B is undefined).
      (as.numeric(.B_scalar(e + h, a, s, beta0, omega)) -
         as.numeric(.B_scalar(e, a, s, beta0, omega))) / h
    } else {
      (as.numeric(.B_scalar(e + h, a, s, beta0, omega)) -
         as.numeric(.B_scalar(e - h, a, s, beta0, omega))) / (2 * h)
    }
  }
  vectorize_scalar(d_scalar)(E)
}
