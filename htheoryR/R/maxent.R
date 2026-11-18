#' Closing the loop: maximize `S[p]` and get `p_i = B(E_i)/Z` back.
#'
#' This is the paper's central claim: the chain `B -> s~(E) -> s(B) ->
#' S[p]` reconstructs an entropy from a prescribed statistical factor, and
#' maximizing that entropy under the standard constraints (normalization
#' and mean energy, no escort probabilities) returns the distribution
#' started from. This closes every link at once: the branch assignment,
#' the `omega`-versus-`a` bookkeeping, the gauge `alpha=-1`, and the `1/Z`
#' scaling -- the highest-value regression test in the package.
#'
#' The gradient is exact and cheap: differentiating `S[p] = sum_j Z^{-1}
#' s(Z p_j)` gives `s'(Z p_i) = alpha + beta0*E(Z p_i)`, so the derivative
#' of the entropy *is* the inverse function `E(B)` -- no numerical
#' differentiation needed.
#' @keywords internal
#' @name maxent-overview
NULL

#' `p_i = B(E_i)/Z` with `Z = sum_i B(E_i)`. Returns `list(p, Z, U)` with
#' `U = sum_i p_i E_i` the mean energy the maximization is constrained to.
#' @param energies Numeric vector of energy levels `E_i`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @return See details.
#' @export
target_distribution <- function(energies, a, s, beta0 = 1.0) {
  energies <- as.numeric(energies)
  b <- as.numeric(B(energies, a, s, beta0))
  Z <- sum(b)
  p <- b / Z
  list(p = p, Z = Z, U = sum(p * energies))
}

#' Maximize `S[p] = sum_i Z^{-1} s(Z p_i)` under `sum_i p_i = 1` and
#' `sum_i p_i E_i = U` (`U` taken from `target_distribution`). `Z` is held
#' at the partition function of the target, a fixed property of the
#' spectrum, not a variable of the problem. `E(B)` comes from the tabulated
#' interpolant of `make_interpolants` (n=800 knots, not the module default
#' of 300: the recovered probabilities inherit the table's accuracy
#' directly, and 300 leaves the power-law class's tail good to only about
#' three digits). SLSQP-equivalent optimization (`nloptr`'s NLOPT_LD_SLSQP)
#' is run twice, the second time warm-started from the first answer: the
#' entropy is very flat near its maximum, so a single pass stops early.
#' @param energies Numeric vector of energy levels `E_i`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param x0 Optional starting point for the optimizer; defaults to a feasible tilted-uniform start.
#' @param p_min Lower bound on each `p_i` (avoids the tabulated interpolant's lower edge).
#' @param tol Optimizer relative tolerance.
#' @param maxiter Maximum optimizer iterations.
#' @param n_table Number of knots in the tabulated interpolant.
#' @return A list with components `p`, `p_target`, `Z`, `U`, `entropy`, `entropy_target`, `max_abs_error`, `success`, `message`, `nloptr_result`.
#' @export
maximize <- function(energies, a, s, beta0 = 1.0, x0 = NULL, p_min = 1e-14,
                      tol = 1e-12, maxiter = 500L, n_table = 800L) {
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("maximize() requires the 'nloptr' package")
  }
  energies <- as.numeric(energies)
  a <- .as_a(a)
  target <- target_distribution(energies, a, s, beta0)
  p_target <- target$p; Z <- target$Z; U <- target$U

  interp <- make_interpolants(a, s, beta0, n = n_table,
                               B_min = min(1e-13, 0.01 * Z * min(p_target)))
  # Bounds are the *domain* of the entropy trace form, not a hint at the
  # answer: E is defined only for x = B in (0, 1], so Z*p_i must stay
  # there. Without the upper bound the optimizer can walk into Z*p_i > 1,
  # where the tabulated s saturates instead of continuing, and reports a
  # spurious "maximum" above the true one.
  p_min <- max(p_min, interp$B_min / Z)
  p_max <- 1.0 / Z

  entropy_fn <- function(p) {
    x <- pmax(Z * p, interp$B_min)
    sum(interp$s_of_B(x)) / Z
  }
  neg_grad <- function(p) {
    x <- pmax(Z * p, interp$B_min)
    -(beta0 * interp$E_of_B(x) - 1.0)
  }

  n_e <- length(energies)
  if (is.null(x0)) {
    x0 <- .feasible_start(energies, U)
  } else {
    x0 <- pmin(pmax(x0, p_min), p_max)
  }

  # nloptr minimizes; neg_grad(p) is already d(-S)/dp, exactly the gradient
  # the objective -entropy_fn(p) needs.
  eval_f <- function(p) list(objective = -entropy_fn(p), gradient = neg_grad(p))
  eq_fn <- function(p) c(sum(p) - 1.0, sum(p * energies) - U)
  eq_jac <- function(p) rbind(rep(1, n_e), energies)

  res <- NULL
  for (iter in 1:2) {
    res <- nloptr::nloptr(
      x0 = x0,
      eval_f = eval_f,
      lb = rep(p_min, n_e), ub = rep(p_max, n_e),
      eval_g_eq = eq_fn, eval_jac_g_eq = eq_jac,
      opts = list(algorithm = "NLOPT_LD_SLSQP", maxeval = maxiter, xtol_rel = tol)
    )
    x0 <- pmin(pmax(res$solution, p_min), p_max)
    x0 <- x0 / sum(x0)
  }
  p <- pmin(pmax(res$solution, p_min), p_max)
  p <- p / sum(p)

  list(
    p = p, p_target = p_target, Z = Z, U = U,
    entropy = entropy_fn(p), entropy_target = entropy_fn(p_target),
    max_abs_error = max(abs(p - p_target)),
    success = res$status > 0, message = res$message, nloptr_result = res
  )
}

#' Spread of `s'(Z p_i) - beta0*E_i` across states: this is the constant
#' `alpha` at the maximizer, so the spread is zero there. Reported as
#' `max - min`: a direct, optimizer-free check, and the one that localizes
#' a broken link fastest.
#' @param p Probability vector on the simplex.
#' @param energies Numeric vector of energy levels `E_i`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param Z Partition function (normalizing constant).
#' @return See details.
#' @export
stationarity_residual <- function(p, energies, a, s, beta0 = 1.0, Z = NULL) {
  p <- as.numeric(p); energies <- as.numeric(energies)
  if (is.null(Z)) Z <- sum(as.numeric(B(energies, a, s, beta0)))
  s_prime <- beta0 * E_of_B(Z * p, a, s, beta0) - 1.0
  resid <- s_prime - beta0 * energies
  max(resid) - min(resid)
}

#' A point on the simplex with the right mean energy, and nothing else:
#' uniform weights tilted linearly toward `U`, scaled back if that would
#' leave the simplex.
#' @noRd
.feasible_start <- function(energies, U) {
  n <- length(energies)
  p <- rep(1 / n, n)
  centred <- energies - mean(energies)
  denom <- sum(centred^2) / n
  if (denom <= 0) return(p)
  lam <- (U - mean(energies)) / denom
  for (i in seq_len(60)) {
    trial <- p + lam * centred / n
    if (min(trial) > 0) return(trial / sum(trial))
    lam <- lam * 0.5
  }
  p
}
