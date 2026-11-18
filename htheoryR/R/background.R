#' The fluctuating background: the reservoir temperature law and its dynamics.
#'
#' Three descriptions of the same object, in increasing order of directness:
#' `f_N` (closed form, Eq. 33, a Meijer G), `sample_TN` (ancestral sampling
#' down the chain of conditionals, Eqs. 31/32), and `simulate_sde`
#' (Euler-Maruyama on the SDE hierarchy, Eq. 29). They agree only if the
#' branch labels of Eq. (33) are read as corrected -- see `ERRATUM_EQ33`.
#' `tests/testthat/test-background.R` pins that.
#' @keywords internal
#' @name background-overview
NULL

ERRATUM_EQ33 <- paste(
  "Eq. (33) as printed has its two branch labels crossed: the row labelled",
  "s = 1/2 is the form that compounds to the power-law factor, and the row",
  "labelled s = 1 the one that compounds to the stretched exponential --",
  "inconsistent with Eqs. (31)/(32) and (34). The corrected assignment,",
  "implemented here, is",
  "",
  "  s = 1/2 (gamma chain)     f_N = w/(T0 G(a))     G^{N,0}_{0,N}( w T_N/T0 | - ; a-1  )",
  "  s = 1   (inverse-gamma)   f_N = 1/(T0 w G(a+1)) G^{0,N}_{N,0}( T_N/(w T0) | -a-1 ; - )",
  "",
  "verified to 8 digits against quadrature of Eq. (24) and against Monte-Carlo",
  "ancestral sampling.",
  sep = "\n"
)

#' Invert `a_i = 2*gamma_i/kappa_i^2` for the noise amplitudes.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param gamma Vector of relaxation rates (Eq. 29).
#' @return See details.
#' @export
kappa_from_a <- function(a, gamma) {
  a <- .as_a(a)
  gamma <- rep_len(gamma, length(a))
  sqrt(2 * gamma / a)
}

#' Inverse-gamma density, `shape = alpha`, `scale = beta`: no base-R
#' equivalent to `scipy.stats.invgamma`, so hand-rolled from its log-density.
#' @noRd
.dinvgamma <- function(x, shape, scale) {
  ifelse(x > 0,
         exp(shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x),
         0)
}

#' `f(T_i | T_{i-1})`, Eq. (31) for `s = 1/2` and Eq. (32) for `s = 1`.
#'
#' Both have mean `T_{i-1}`: the linear drift of Eq. (29) drives every level
#' to the temperature of the one above it.
#' @param T Sample of the reservoir temperature (or, in `conditional_pdf`, the value(s) at which to evaluate the density).
#' @param T_prev Parent-level temperature, scalar or vector.
#' @param a_i Scalar shape parameter for this level.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @return See details.
#' @export
conditional_pdf <- function(T, T_prev, a_i, s) {
  .check_s(s)
  if (isTRUE(all.equal(s, 0.5))) {
    return(stats::dgamma(T, shape = a_i, scale = T_prev / a_i))
  }
  .dinvgamma(T, shape = a_i + 1, scale = a_i * T_prev)
}

#' Draw `T_i ~ f(T_i | T_{i-1})`, Eqs. (31)/(32).
#'
#' Unlike the Python original (which threads an explicit `rng` object through
#' every call for reproducibility), this uses R's global RNG state directly
#' (idiomatic R -- callers control reproducibility with `set.seed()`). Bit-
#' identical streams with the Python version are unreachable regardless (R's
#' Mersenne-Twister differs from numpy's PCG64), so nothing is lost.
#' @param T_prev Parent-level temperature, scalar or vector.
#' @param a_i Scalar shape parameter for this level.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param size Number of draws.
#' @return See details.
#' @export
sample_conditional <- function(T_prev, a_i, s, size = length(T_prev)) {
  .check_s(s)
  if (isTRUE(all.equal(s, 0.5))) {
    return(stats::rgamma(size, shape = a_i, rate = 1) * T_prev / a_i)
  }
  a_i * T_prev / stats::rgamma(size, shape = a_i + 1, rate = 1)
}

#' The scalar core of `f_N`, Eq. (33), corrected branch assignment.
#' @noRd
.f_N_scalar <- function(T, a, s, T0, prec_bits = DEFAULT_PREC_BITS) {
  if (T <= 0) return(0)
  if (length(a) == 0) stop("f_N is a delta at T0 for N = 0; it has no density")
  w <- .omega(a)
  if (isTRUE(all.equal(s, 0.5))) {
    pref <- Rmpfr::mpfr(w / T0, prec_bits)
    for (x in a) pref <- pref / gamma(Rmpfr::mpfr(x, prec_bits))
    g <- g_gamma_chain(w * T / T0, a - 1, prec_bits, underflow_above = UNDERFLOW_Z)
    return(as.numeric(pref * g))
  }
  pref <- Rmpfr::mpfr(1 / (T0 * w), prec_bits)
  for (x in a) pref <- pref / gamma(Rmpfr::mpfr(x, prec_bits) + 1)
  # This orientation dies at *small* argument: f ~ exp(-omega*T0/T), an
  # essential singularity at T = 0. At the threshold the density is already
  # exp(-100).
  g <- g_invgamma_chain(T / (w * T0), -a - 1, prec_bits, underflow_below = UNDERFLOW_Z_SMALL)
  as.numeric(pref * g)
}

#' The reservoir temperature law `f_N(T_N)`, Eq. (33), corrected.
#'
#' Normalizes to 1 and has mean `T0` for every `N`, `a` and both classes --
#' the two facts the sum rule Eq. (27) rests on.
#' @param T Sample of the reservoir temperature (or, in `conditional_pdf`, the value(s) at which to evaluate the density).
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param T0 Base (outermost) reservoir temperature.
#' @return See details.
#' @export
f_N <- function(T, a, s, T0 = 1.0) {
  a <- .as_a(a)
  vectorize_scalar(.f_N_scalar)(T, a, s, T0)
}

#' Ancestral sampling of `T_N` down the chain, Eq. (30).
#'
#' `T_0 = T0` is fixed; each level is drawn from its conditional given the
#' level above. The empirical law of the result is Eq. (33) -- this is the
#' Monte-Carlo check that the corrected branch assignment is the right one.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param T0 Base (outermost) reservoir temperature.
#' @param size Number of draws.
#' @return See details.
#' @export
sample_TN <- function(a, s, T0 = 1.0, size = 1) {
  a <- .as_a(a)
  T <- rep(as.numeric(T0), size)
  for (a_i in a) T <- sample_conditional(T, a_i, s, size = size)
  T
}

#' Euler-Maruyama on the hierarchy
#' `dT_i = -gamma_i (T_i - T_{i-1}) dt + kappa_i T_i^s T_{i-1}^{1-s} dW_i`,
#' Eq. (29).
#'
#' `gamma` is the vector of relaxation rates `gamma_i`; the noise amplitudes
#' follow from `a_i = 2*gamma_i/kappa_i^2`. Eq. (29) presumes `gamma_i >>
#' gamma_{i-1}` (scale separation), the caller's responsibility. Returns the
#' post-burn-in trajectory of the innermost level `T_N` (matrix, `n_paths` x
#' kept steps), whose stationary histogram is `f_N`. Full-truncation scheme
#' (the `s = 1/2` CIR case is only guaranteed nonnegative in the exact
#' process).
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param gamma Vector of relaxation rates (Eq. 29).
#' @param T0 Base (outermost) reservoir temperature.
#' @param dt Euler-Maruyama time step.
#' @param n_steps Number of Euler-Maruyama steps to simulate.
#' @param burn_in Fraction of steps discarded as burn-in.
#' @param n_paths Number of independent paths to simulate.
#' @return A `n_paths` by `(n_steps - burn-in steps)` matrix of the post-burn-in `T_N` trajectory.
#' @export
simulate_sde <- function(a, s, gamma, T0 = 1.0, dt = 1e-3, n_steps = 200000L,
                          burn_in = 0.1, n_paths = 1L) {
  a <- .as_a(a)
  N <- length(a)
  if (N == 0) return(matrix(as.numeric(T0), nrow = n_paths, ncol = 1))
  gamma <- rep_len(gamma, N)
  kappa <- kappa_from_a(a, gamma)
  if (max(gamma) * dt > 0.1) {
    stop(sprintf("dt = %s is too coarse for gamma_max = %s; need gamma_max * dt << 1",
                  dt, max(gamma)))
  }

  Tmat <- matrix(as.numeric(T0), nrow = n_paths, ncol = N)
  n_burn <- as.integer(burn_in * n_steps)
  out <- matrix(NA_real_, nrow = n_paths, ncol = n_steps - n_burn)
  sqrt_dt <- sqrt(dt)
  s_arr <- matrix(s, nrow = n_paths, ncol = N, byrow = TRUE)
  gamma_row <- matrix(gamma, nrow = n_paths, ncol = N, byrow = TRUE)
  kappa_row <- matrix(kappa, nrow = n_paths, ncol = N, byrow = TRUE)

  for (step in seq_len(n_steps)) {
    prev <- cbind(rep(as.numeric(T0), n_paths), Tmat[, -N, drop = FALSE])
    pos <- pmax(Tmat, 0)
    diff <- kappa_row * pos^s_arr * prev^(1 - s_arr)
    dW <- matrix(stats::rnorm(n_paths * N), nrow = n_paths, ncol = N) * sqrt_dt
    Tmat <- Tmat - gamma_row * (Tmat - prev) * dt + diff * dW
    Tmat <- pmax(Tmat, 0)
    if (step > n_burn) out[, step - n_burn] <- Tmat[, N]
  }
  out
}
