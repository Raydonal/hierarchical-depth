#' Identifiability of the cascade depth `N` and shapes `a` from `T_N` alone.
#'
#' The cascade is additive in log-temperature, and its log-cumulants of
#' order `r>=2` identify the multiset `{a_i}` -- hence `N` -- without
#' observing any intermediate level. See `docs/theorems.md` Theorem 3 (in
#' the parent project) for the full derivation; the short version: writing
#' the increment at level `i` as `D_i = log T_i - log T_{i-1}`, `D_i` is a
#' function of the fresh per-level draw `Y_i` alone (independent of
#' `T_{i-1}`), so the `D_i` are mutually independent and cumulants add.
#' @keywords internal
#' @name identifiability-overview
NULL

#' @noRd
.shift_id <- function(s) {
  if (isTRUE(all.equal(s, 0.5))) return(0.0)
  if (isTRUE(all.equal(s, 1))) return(1.0)
  stop(sprintf("s must be 1/2 or 1, got %s", s))
}

#' Theoretical `kappa_r(log T_N)` for `r=2..r_max`. Independent of `T_0`
#' by construction; relies only on the shape vector `a` and branch `s`. `a`
#' may be empty (`N = 0`), for which all cumulants are 0.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param r_max Highest cumulant order to compute (cumulants of order `2..r_max` are returned).
#' @return See details.
#' @export
log_cumulants_theory <- function(a, s, r_max = 4L) {
  a <- as.numeric(a)
  off <- .shift_id(s)
  out <- numeric(r_max - 1)
  for (idx in seq_len(r_max - 1)) {
    r <- idx + 1L
    total <- if (length(a)) sum(psigamma(a + off, deriv = r - 1)) else 0
    out[idx] <- if (isTRUE(all.equal(s, 0.5))) total else ((-1)^r) * total
  }
  out
}

#' The unbiased k-statistic of order `r`. R has no base `kstat`; the
#' formulas below are the standard unbiased estimators in terms of the
#' power sums (David & Kendall's formulas, the ones `scipy.stats.kstat`
#' itself implements) -- cross-checked against `scipy.stats.kstat` directly
#' on a fixed vector during development.
#' @keywords internal
#' @noRd
.kstat <- function(x, r) {
  n <- length(x)
  m <- function(p) sum((x - mean(x))^p) / n  # central sample moment
  m2 <- m(2); m3 <- m(3)
  if (r == 2) return(n / (n - 1) * m2)
  if (r == 3) return(n^2 / ((n - 1) * (n - 2)) * m3)
  if (r == 4) {
    m4 <- m(4)
    return((n^2 * ((n + 1) * m4 - 3 * (n - 1) * m2^2)) / ((n - 1) * (n - 2) * (n - 3)))
  }
  stop("k-statistic only implemented for r in 2:4")
}

#' Unbiased sample cumulants `k_r` of `log T` for `r=2..r_max`, via the
#' standard k-statistic (directly comparable to `log_cumulants_theory` with
#' no finite-sample bias correction).
#' @param T Sample of the reservoir temperature (or, in `conditional_pdf`, the value(s) at which to evaluate the density).
#' @param r_max Highest cumulant order to compute (cumulants of order `2..r_max` are returned).
#' @return See details.
#' @export
empirical_log_cumulants <- function(T, r_max = 4L) {
  x <- log(as.numeric(T))
  vapply(2:r_max, function(r) .kstat(x, r), numeric(1))
}

#' `kappa_3/kappa_2` for a single homogeneous level, as a function of `a`.
#' Delegates to `log_cumulants_theory` (rather than re-deriving the branch
#' sign here) since the ratio for `N` homogeneous levels equals the
#' single-level ratio regardless of `N` -- the factor of `N` cancels.
#' @noRd
.ratio_theory <- function(a, s) {
  k23 <- log_cumulants_theory(a, s, r_max = 3)
  k23[2] / k23[1]
}

#' Point estimate of `(N, a)` under the homogeneous model `a_i = a`, from
#' `(kappa2, kappa3)`. Two equations, two unknowns: `kappa3/kappa2` depends
#' on `a` alone (the `N` factor cancels), recovered by a 1-D root find; then
#' `N = kappa2 / psi'(a)` from the scale.
#' @param kappa2 Observed (or theoretical) second log-cumulant.
#' @param kappa3 Observed (or theoretical) third log-cumulant.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param a_bracket Two-element vector, the search bracket for `a` in the root-find.
#' @return See details.
#' @export
fit_homogeneous <- function(kappa2, kappa3, s, a_bracket = c(1e-4, 1e4)) {
  target <- kappa3 / kappa2
  f <- function(a) .ratio_theory(a, s) - target
  lo <- a_bracket[1]; hi <- a_bracket[2]
  if (f(lo) * f(hi) > 0) {
    stop(sprintf("kappa3/kappa2 = %s not bracketed by a in (%s, %s); widen a_bracket or check the sample size",
                  target, lo, hi))
  }
  a_hat <- uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
  off <- .shift_id(s)
  N_hat <- kappa2 / psigamma(a_hat + off, deriv = 1)
  list(a = a_hat, N = N_hat, N_round = max(0L, round(N_hat)))
}

#' Parametric-bootstrap test of `H0: depth = N0` (homogeneous `a`) against
#' the sample `T`. Fits `a` under `H0` from `kappa2` alone, so the test
#' statistic -- observed `kappa3` minus the value predicted at `(N0,
#' a_hat)` -- is not tautologically zero. Its null distribution is a
#' parametric bootstrap: simulate `n_boot` depth-`N0` cascades at the
#' fitted shape and refit on each. `N` sits at the boundary of its domain
#' (a positive integer), so the usual chi-square asymptotics for a
#' likelihood-ratio test do not apply (Self & Liang 1987), hence bootstrap
#' rather than a closed-form p-value.
#' @param T Sample of the reservoir temperature (or, in `conditional_pdf`, the value(s) at which to evaluate the density).
#' @param N0 The depth under the null hypothesis `H0: N = N0`.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param n_boot Number of bootstrap replications.
#' @param T0 Base (outermost) reservoir temperature.
#' @return See details.
#' @export
test_depth <- function(T, N0, s, n_boot = 2000L, T0 = 1.0) {
  T <- as.numeric(T)
  n <- length(T)
  off <- .shift_id(s)
  k23_obs <- empirical_log_cumulants(T, r_max = 3)

  a_lo <- 1e-6; a_hi <- 1e6
  psi_hi <- psigamma(a_hi + off, deriv = 1)
  psi_lo <- psigamma(a_lo + off, deriv = 1)

  discrepancy <- function(k2, k3) {
    target <- min(max(k2 / N0, psi_hi * 1.0000001), psi_lo * 0.9999999)
    a_hat <- uniroot(function(a) psigamma(a + off, deriv = 1) - target,
                      lower = a_lo, upper = a_hi, tol = 1e-10)$root
    k3_pred <- N0 * psigamma(a_hat + off, deriv = 2)
    k3_pred <- if (isTRUE(all.equal(s, 0.5))) k3_pred else -k3_pred
    list(stat = k3 - k3_pred, a_hat = a_hat)
  }

  obs <- discrepancy(k23_obs[1], k23_obs[2])
  boot_stats <- numeric(n_boot)
  a_vec_null <- rep(obs$a_hat, N0)
  for (b in seq_len(n_boot)) {
    T_boot <- sample_TN(a_vec_null, s, T0 = T0, size = n)
    k23_b <- empirical_log_cumulants(T_boot, r_max = 3)
    boot_stats[b] <- discrepancy(k23_b[1], k23_b[2])$stat
  }
  p_value <- mean(abs(boot_stats) >= abs(obs$stat))
  list(statistic = obs$stat, p_value = p_value, a_hat = obs$a_hat, n_boot = n_boot)
}

#' Simulated power of `test_depth` at `H0: N = N0` against a true cascade
#' `a_true` (its length is the true depth). Returns a named list `n ->
#' power`: reject rate under a true deeper cascade, at a shallow null.
#' @param a_true True (possibly heterogeneous) shape vector generating the data; its length is the true depth.
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param N0 The depth under the null hypothesis `H0: N = N0`.
#' @param sample_sizes Numeric vector of sample sizes to evaluate power at.
#' @param n_reps Number of Monte Carlo replications per sample size.
#' @param n_boot Number of bootstrap replications.
#' @param alpha Nominal test significance level.
#' @param T0 Base (outermost) reservoir temperature.
#' @return A named list mapping each sample size (as a character string) to its estimated power.
#' @export
power_curve <- function(a_true, s, N0, sample_sizes = c(100, 300, 1000, 5000),
                         n_reps = 200L, n_boot = 400L, alpha = 0.05, T0 = 1.0) {
  powers <- list()
  for (n in sample_sizes) {
    rejects <- 0L
    for (rep in seq_len(n_reps)) {
      T <- sample_TN(a_true, s, T0 = T0, size = n)
      result <- test_depth(T, N0, s, n_boot = n_boot, T0 = T0)
      if (result$p_value < alpha) rejects <- rejects + 1L
    }
    powers[[as.character(n)]] <- rejects / n_reps
  }
  powers
}
