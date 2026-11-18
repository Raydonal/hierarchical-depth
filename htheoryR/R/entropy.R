#' The entropy generator `s~_N^{(s)}(E)`, Eq. (39), the parametric curve
#' `(B(E), s~(E))`, numerical inversion `E(B)`, and the simplex functional
#' `S[p] = sum_i Z^{-1} s(Z p_i)`, Eq. (10).
#' @keywords internal
#' @name entropy-overview
NULL

#' `s~_N^{(s)}(E) = (beta0*E - 1)*B(E) + Btilde(E)`, Eq. (39).
#' @param E Energy: a nonnegative scalar or numeric vector.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @return See details.
#' @export
s_tilde <- function(E, a, s, beta0 = 1.0) {
  (beta0 * E - 1) * B(E, a, s, beta0) + B_tilde(E, a, s, beta0)
}

#' The parametric curve `(B(E), s~(E))` -- geometric in `1 + beta0*E` so both
#' the `B ~ 1` knee and the tail are resolved. Outside the three solvable
#' cases `E(B)` has no closed form, so `s(B)` exists only as this curve or
#' through `E_of_B`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param E Energy: a nonnegative scalar or numeric vector.
#' @param n Number of grid points (or, in `sample_TN`/`power_curve`, a sample size).
#' @param E_max See details.
#' @return See details.
#' @export
parametric_curve <- function(a, s, beta0 = 1.0, E = NULL, n = 400L, E_max = 60.0) {
  if (is.null(E)) E <- (exp(seq(0, log1p(beta0 * E_max), length.out = n)) - 1) / beta0
  list(B = B(E, a, s, beta0), s_tilde = s_tilde(E, a, s, beta0), E = E)
}

#' Bracket-then-bisect `E` such that `B(E; a, s, beta0) = Bval`, `0 < Bval <= 1`.
#' `B` is strictly decreasing (Eq. 43), so this is well-posed.
#' @noRd
.E_of_B_scalar <- function(Bval, a, s, beta0, xtol = 1e-12) {
  if (!(Bval > 0 && Bval <= 1)) stop(sprintf("B must lie in (0, 1], got %s", Bval))
  if (Bval == 1) return(0)
  hi <- 1 / beta0
  omega <- .omega(a)
  repeat {
    if (as.numeric(.B(hi, a, s, beta0, omega)) < Bval) break
    hi <- hi * 2
    if (hi > 2^200) stop(sprintf("could not bracket E(B) for B = %s", Bval))
  }
  f <- function(e) as.numeric(.B(e, a, s, beta0, omega)) - Bval
  stats::uniroot(f, lower = 0, upper = hi, tol = xtol)$root
}

#' `E` such that `B(E) = Bval`, vectorized over `Bval`. Well-posed for every
#' `(N, a)` and both classes because `B` is strictly decreasing from
#' `B(0) = 1` to `B(inf) = 0`, Eqs. (7)/(43).
#' @param Bval Value(s) of the statistical factor `B`, in `(0, 1]`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param xtol Root-finding tolerance.
#' @return See details.
#' @export
E_of_B <- function(Bval, a, s, beta0 = 1.0, xtol = 1e-12) {
  a <- .as_a(a)
  vectorize_scalar(.E_of_B_scalar)(Bval, a, s, beta0, xtol)
}

#' Build monotone-cubic interpolants of `E(B)` and `s~(B)`, amortizing the
#' cost of repeated `B()`/`E_of_B()` evaluation (each one, worst case, an
#' escalating-precision Meijer-G or quadrature call -- see R/meijer.R and
#' R/statfactor.R) over a table built once. Two curves are tabulated, both
#' against `log(B)`, NOT `E` or `B` directly:
#'
#' - `log(1 + beta0*E)` vs `log(B)`, inverted to give `E_of_B`;
#' - the *bounded ratio* `u(B) = s~(B) / (B*(1+beta0*E))` vs `log(B)`, not
#'   `s~` itself -- `s~ -> 0` and `B -> 0` together in the tail, so splining
#'   `s~` directly loses exactly the relative accuracy there that a
#'   downstream `maxent` optimizer needs; `u(B)` stays `O(1)` across the
#'   whole tail by construction, which is what actually gets splined well.
#'
#' Returns `list(E_of_B = function(B), s_of_B = function(B))`. Queries
#' outside the tabulated `log(B)` range are clamped to the nearest tabulated
#' endpoint rather than extrapolated.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param n Number of grid points (or, in `sample_TN`/`power_curve`, a sample size).
#' @param B_min Smallest `B` value to tabulate down to.
#' @return See details.
#' @export
make_interpolants <- function(a, s, beta0 = 1.0, n = 300L, B_min = 1e-12) {
  a <- .as_a(a)
  omega <- .omega(a)

  # Bracket E_max such that B(E_max) <= B_min, doubling from E = 1.
  E_max <- 1
  while (as.numeric(.B(E_max, a, s, beta0, omega)) > B_min) {
    E_max <- E_max * 2
    if (E_max > 1e18) stop("make_interpolants: could not bracket E_max")
  }

  # Log-uniform across the *whole* range, not just the top end: a single
  # near-zero anchor point followed by a jump straight to `E_max/n` (tried
  # first) leaves a huge, unresolved gap in B in-between whenever `E_max` is
  # large (e.g. E_max ~ 170 with n = 100 jumps straight from E = 1e-8 to
  # E = 1.7) -- caught via E_of_B(B=0.9) coming back 10x wrong.
  E_grid <- exp(seq(log(E_max * 1e-8), log(E_max), length.out = n))
  B_grid <- as.numeric(vapply(E_grid, function(e) as.numeric(.B(e, a, s, beta0, omega)), numeric(1)))
  s_grid <- s_tilde(E_grid, a, s, beta0)

  logB <- log(pmax(B_grid, .Machine$double.xmin))
  log1E <- log1p(beta0 * E_grid)
  u_grid <- s_grid / (B_grid * (1 + beta0 * E_grid))

  ord <- order(logB)  # B decreases as E increases; sort ascending in logB
  logB <- logB[ord]; log1E <- log1E[ord]; u_grid <- u_grid[ord]
  keep <- !duplicated(logB)
  logB <- logB[keep]; log1E <- log1E[keep]; u_grid <- u_grid[keep]

  spline_log1E <- stats::splinefun(logB, log1E, method = "monoH.FC")
  spline_u <- stats::splinefun(logB, u_grid, method = "monoH.FC")
  rng <- range(logB)

  clamp_logB <- function(Bval) pmin(pmax(log(pmax(Bval, .Machine$double.xmin)), rng[1]), rng[2])

  list(
    E_of_B = function(Bval) {
      lb <- clamp_logB(Bval)
      (exp(spline_log1E(lb)) - 1) / beta0
    },
    s_of_B = function(Bval) {
      lb <- clamp_logB(Bval)
      Eval <- (exp(spline_log1E(lb)) - 1) / beta0
      spline_u(lb) * Bval * (1 + beta0 * Eval)
    },
    B_min = exp(rng[1])  # smallest B actually tabulated; maxent.R clips queries here
  )
}

#' Convenience wrapper returning just the tabulated `E_of_B` function.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param n Number of grid points (or, in `sample_TN`/`power_curve`, a sample size).
#' @param B_min Smallest `B` value to tabulate down to.
#' @return See details.
#' @export
make_E_of_B <- function(a, s, beta0 = 1.0, n = 300L, B_min = 1e-12) {
  make_interpolants(a, s, beta0, n, B_min)$E_of_B
}

#' `s~(E_of_B(Bval))`, exact (not tabulated) composition.
#' @param Bval Value(s) of the statistical factor `B`, in `(0, 1]`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @return See details.
#' @export
s_of_B <- function(Bval, a, s, beta0 = 1.0) {
  Eval <- E_of_B(Bval, a, s, beta0)
  s_tilde(Eval, a, s, beta0)
}

#' `d(s~)/d(B) = beta0*E_of_B(B) - 1`, exact via the inverse function (Eq. 43
#' makes `s~` a function of `E` whose derivative w.r.t. `B` collapses to
#' this, with no differentiation needed).
#' @param Bval Value(s) of the statistical factor `B`, in `(0, 1]`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @return See details.
#' @export
s_prime_of_B <- function(Bval, a, s, beta0 = 1.0) {
  beta0 * E_of_B(Bval, a, s, beta0) - 1
}

#' `S[p] = sum_i s_of_B(Z*p_i)/Z + c*sum(p_i)`, Eq. (10). `Z` (the partition
#' function) is mandatory -- it is supplied externally, not recomputed here.
#' @param p Probability vector on the simplex.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @param s Universality class: `0.5` (gamma / stretched-exponential) or `1` (inverse-gamma / power-law).
#' @param beta0 Inverse base temperature, `1/T0`.
#' @param Z Partition function (normalizing constant).
#' @param c Additive constant term.
#' @return See details.
#' @export
S_simplex <- function(p, a, s, beta0 = 1.0, Z, c = 0.0) {
  sum(s_of_B(Z * p, a, s, beta0)) / Z + c * sum(p)
}

#' `-B*log(B)`, Eq. (21): the `N=0` / plain-exponential special case, where
#' `s~` collapses to this exactly.
#' @param Bval Value(s) of the statistical factor `B`, in `(0, 1]`.
#' @return See details.
#' @export
s_boltzmann_shannon <- function(Bval) ifelse(Bval > 0, -Bval * log(pmax(Bval, .Machine$double.xmin)), 0)

#' `s(B) = -(a+1)*(B - B^q)`, Eq. (63) -- the `N=1`, `s=1` exactly-solvable
#' case, `q = a/(a+1)` the entropy index. Summed over the simplex this is
#' exactly the Tsallis entropy, Eqs. (64)/(67).
#' @param Bval Value(s) of the statistical factor `B`, in `(0, 1]`.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @return See details.
#' @export
s_tsallis_kernel <- function(Bval, a) {
  q <- a / (a + 1)
  -(a + 1) * (Bval - Bval^q)
}

#' The two Tsallis indices for the `s=1`, `N=1` case (Eq. 65) -- distinct and
#' must not be conflated: `q` is the *entropy* index, `q_stat` the
#' *distribution* index.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @return See details.
#' @export
q_indices <- function(a) {
  q <- a / (a + 1)
  q_stat <- 2 - q
  list(q = q, q_stat = q_stat)
}

#' The exact Tsallis entropy of a probability vector `p`, at entropy index
#' `q = a/(a+1)`.
#' @param p Probability vector on the simplex.
#' @param a Numeric vector of per-level shape parameters `(a_1, ..., a_N)`; length is the cascade depth `N` (empty for `N = 0`).
#' @return See details.
#' @export
tsallis_entropy <- function(p, a) {
  q <- a / (a + 1)
  (1 - sum(p^q)) / (q - 1)
}
