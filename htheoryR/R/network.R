#' Hierarchical random networks: the mixed-Poisson application.
#'
#' Each vertex `v` belongs to a nested sequence of communities, and its
#' propensity to connect is the product of per-level contributions, `W_v =
#' prod_i X_{i,c_i(v)}` -- exactly the multiplicative cascade of
#' `sample_TN`, relabelled. Degrees are then mixed Poisson, `D_v | W_v ~
#' Poisson(lambda*W_v)` (Chung-Lu / Norros-Reittu).
#'
#' The identifiability chain: (1) mixed-Poisson laws are identifiable from
#' the degree distribution alone (Feller 1943) -- the law of `W` is
#' recoverable from `D`; (2) by Theorem 3 (`R/identifiability.R`),
#' `(N, {a_i})` is identifiable from the law of `W`; (3) hence the cascade
#' depth is estimable from the degree sequence alone, without observing
#' which vertex belongs to which community.
#'
#' The estimator is the factorial-moment analogue of the log-cumulant one:
#' `E[D^(k)] = lambda^k * E[W^k]` exactly (a mixed-Poisson identity, not an
#' approximation), and `E[W^k]` has a closed form linear in `N` for fixed
#' `a`, so the same ratio-then-scale trick as `fit_homogeneous` applies.
#' @keywords internal
#' @name network-overview
NULL

#' Draw a mixed-Poisson degree sequence, `D_v | W_v ~ Poisson(lambda*W_v)`.
#' `W` is drawn from `sample_TN` with `T0 = 1`, so `E[W] = 1` and `E[D] =
#' lambda` exactly.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param lam Poisson rate parameter.
#' @param size Number of draws.
#' @param T0 Base (outermost) reservoir temperature.
#' @return See details.
#' @export
sample_degrees <- function(a, s, lam, size, T0 = 1.0) {
  W <- sample_TN(a, s, T0 = T0, size = size)
  rpois(size, lam * W)
}

#' Per-level `log E[X_i^k]` for the single multiplicative factor at each
#' level (parent temperature normalized away, mean 1). `s = 1` requires
#' `a_i > k - 1` for all levels (`E[W^k]` diverges otherwise -- a
#' heavy-tail effect, not a bug).
#' @noRd
.level_log_moment <- function(a, s, k) {
  if (isTRUE(all.equal(s, 0.5))) {
    return(lgamma(a + k) - lgamma(a) - k * log(a))
  }
  if (isTRUE(all.equal(s, 1))) {
    if (any(a + 1 - k <= 0)) {
      stop(sprintf("k = %d requires a_i > %d for all levels in the s = 1 branch (E[W^k] diverges otherwise); got a = %s",
                    k, k - 1, paste(a, collapse = ",")))
    }
    return(k * log(a) + lgamma(a + 1 - k) - lgamma(a + 1))
  }
  stop(sprintf("s must be 1/2 or 1, got %s", s))
}

#' `log E[W^k]`, summed over independent levels. Homogeneous case (`a_i =
#' a`) reduces to `N*[log Gamma(a+k) - log Gamma(a) - k*log(a)]` for
#' `s=1/2`, linear in `N`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param k Moment order.
#' @return See details.
#' @export
log_factorial_moment_theory <- function(a, s, k) {
  a <- as.numeric(a)
  if (length(a) == 0) return(0.0)
  sum(.level_log_moment(a, s, k))
}

#' `(lambda_hat, c(L_2,...,L_k_max))` with `L_k = log(Ehat[W^k])`.
#' `lambda_hat` is the sample mean of `D` (exact since `E[W]=1`). The
#' falling-factorial product is exactly 0 whenever `D < k`, so no
#' special-casing of small degrees is needed.
#' @param D Integer (or numeric) vector, an observed degree sequence.
#' @param k_max Highest factorial-moment order to compute.
#' @return See details.
#' @export
empirical_log_factorial_moments <- function(D, k_max = 3L) {
  D <- as.numeric(D)
  lam_hat <- mean(D)
  out <- numeric(k_max - 1)
  for (idx in seq_len(k_max - 1)) {
    k <- idx + 1L
    fac <- rep(1, length(D))
    for (j in 0:(k - 1)) fac <- fac * (D - j)
    out[idx] <- log(mean(fac) / lam_hat^k)
  }
  list(lam_hat = lam_hat, L = out)
}

#' Recover `(N, a, lambda)` from a degree sequence, homogeneous model.
#' Mirrors `fit_homogeneous`: the ratio `L_3/L_2` depends on `a` alone (the
#' common `N` factor cancels), a 1-D root find, then `N = L_2 /
#' [single-level log E[W^2]](a)` from the scale.
#' @param D Integer (or numeric) vector, an observed degree sequence.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param k_max Highest factorial-moment order to compute.
#' @param a_bracket Two-element vector, the search bracket for `a` in the root-find.
#' @return See details.
#' @export
fit_homogeneous_network <- function(D, s, k_max = 3L, a_bracket = NULL) {
  km <- empirical_log_factorial_moments(D, k_max = max(k_max, 3L))
  L2 <- km$L[1]; L3 <- km$L[2]
  target <- L3 / L2

  ratio <- function(a) {
    c2 <- .level_log_moment(a, s, 2)
    c3 <- .level_log_moment(a, s, 3)
    c3 / c2
  }
  if (is.null(a_bracket)) {
    a_bracket <- if (isTRUE(all.equal(s, 1))) c(2.0 + 1e-6, 1e3) else c(1e-3, 1e3)
  }
  lo <- a_bracket[1]; hi <- a_bracket[2]
  f <- function(a) ratio(a) - target
  if (f(lo) * f(hi) > 0) {
    stop(sprintf("L3/L2 = %s not bracketed by a in (%s, %s); widen a_bracket or check the sample size",
                  target, lo, hi))
  }
  a_hat <- uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
  c2_hat <- .level_log_moment(a_hat, s, 2)
  N_hat <- L2 / c2_hat
  list(a = a_hat, N = N_hat, N_round = max(0L, round(N_hat)), lam = km$lam_hat)
}

#' Parametric-bootstrap test of `H0: N = N0` (homogeneous `a`) against a
#' degree sequence `D` -- the factorial-moment analogue of `test_depth`.
#' `fit_homogeneous_network`'s point estimate has no formal p-value
#' attached; this closes that gap the same way `test_depth` does for a
#' propensity sequence. A degree sequence whose `L3/L2` ratio cannot be
#' matched by any `a` in the bracket is a genuine model-selection signal
#' (raises, like `fit_homogeneous_network`), not something papered over.
#' @param D Integer (or numeric) vector, an observed degree sequence.
#' @param N0 The depth under the null hypothesis `H0: N = N0`.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param n_boot Number of bootstrap replications.
#' @param k_max Highest factorial-moment order to compute.
#' @return See details.
#' @export
test_depth_degrees <- function(D, N0, s, n_boot = 2000L, k_max = 3L) {
  D <- as.numeric(D)
  n <- length(D)
  km <- empirical_log_factorial_moments(D, k_max = max(k_max, 3L))
  L2_obs <- km$L[1]; L3_obs <- km$L[2]; lam_hat <- km$lam_hat

  if (isTRUE(all.equal(s, 1))) { a_lo <- 2.0 + 1e-6; a_hi <- 1e3 } else { a_lo <- 1e-6; a_hi <- 1e3 }
  f2 <- function(a) .level_log_moment(a, s, 2)
  f3 <- function(a) .level_log_moment(a, s, 3)
  val_at_lo <- f2(a_lo); val_at_hi <- f2(a_hi)

  discrepancy <- function(L2, L3) {
    target <- min(max(L2 / N0, val_at_hi * 1.0000001), val_at_lo * 0.9999999)
    a_hat <- uniroot(function(a) f2(a) - target, lower = a_lo, upper = a_hi, tol = 1e-10)$root
    L3_pred <- N0 * f3(a_hat)
    list(stat = L3 - L3_pred, a_hat = a_hat)
  }

  obs <- discrepancy(L2_obs, L3_obs)
  boot_stats <- numeric(n_boot)
  a_vec_null <- rep(obs$a_hat, N0)
  for (b in seq_len(n_boot)) {
    D_boot <- sample_degrees(a_vec_null, s, lam_hat, size = n)
    km_b <- empirical_log_factorial_moments(D_boot, k_max = max(k_max, 3L))
    boot_stats[b] <- discrepancy(km_b$L[1], km_b$L[2])$stat
  }
  p_value <- mean(abs(boot_stats) >= abs(obs$stat))
  list(statistic = obs$stat, p_value = p_value, a_hat = obs$a_hat, lam_hat = lam_hat, n_boot = n_boot)
}
