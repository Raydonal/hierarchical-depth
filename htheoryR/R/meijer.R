#' Robust, arbitrary-precision evaluation of the two Meijer-G shapes this
#' package needs.
#'
#' Every closed form in the paper is one of exactly two restricted Meijer-G
#' shapes (confirmed by inspecting every call site in the Python original,
#' `background._f_N_scalar` and `statfactor._B_scalar`):
#'
#' ```
#'   G^{N,0}_{0,N}(z | -- ; b_1,...,b_N)   (gamma chain,        s = 1/2)
#'   G^{0,N}_{N,0}(z | b_1,...,b_N ; --)   (inverse-gamma chain, s = 1)
#' ```
#'
#' Both are "numerator-empty"/"denominator-empty" cases of the general
#' Meijer G, which is why they admit a closed-form reduction to a finite sum
#' of generalized hypergeometric terms rather than needing a black-box
#' Mellin-Barnes routine (R has no arbitrary-precision *complex* arithmetic
#' -- `Rmpfr` wraps MPFR, real only -- which is the concrete reason contour
#' quadrature was rejected as the implementation strategy; see the plan).
#'
#' Derivation (residue sum of the Mellin-Barnes integral, closing the
#' contour at the simple poles of each Gamma(b_j + s)): for pairwise-distinct
#' b_i,
#'
#' ```
#'   G^{N,0}_{0,N}(z|--;b) = sum_j  z^{b_j} * prod_{i!=j} Gamma(b_i - b_j)
#'                                * 0F_{N-1}(--; {1+b_j-b_i}_{i!=j}; (-1)^N z)
#' ```
#'
#' Verified by hand against two independent, oracle-free identities already
#' pinned in the Python test suite (`tests/test_meijer.py`):
#' ```
#'   N=1: reduces to z^{b_1} e^{-z}  (matches G^{1,0}_{0,1}(z|--;0) = e^{-z})
#'   N=2: reduces to 2 z^{a/2} K_a(2 sqrt(z))  (the b=(a,0) Bessel-K identity)
#' ```
#' Both reductions are worked out in full in docs/meijer_derivation.md.
#'
#' The mirror shape needs no separate derivation: substituting s -> -s in its
#' own Mellin-Barnes integral gives
#' ```
#'   G^{0,N}_{N,0}(z|b;--) = G^{N,0}_{0,N}(1/z|--;1-b)
#' ```
#' directly (this is DLMF 16.19.1 specialized to p=m=0), so only one routine
#' (`.g_core`) is implemented and hardened.
#' @keywords internal
#' @name meijer-overview
NULL

DEFAULT_PREC_BITS <- 128L  # ~= 38 decimal digits, comfortably above mpmath's
                            # 20-digit runtime default (meijer.DEFAULT_DPS)

# Mirrors meijer.py's UNDERFLOW_Z / UNDERFLOW_Z_SMALL exactly (same values,
# same rationale): the gamma-chain shape carries e^{-z} and dies at large z;
# its mirror carries e^{-1/z} and dies at small z. Past these thresholds the
# true value is confidently zero at any working precision that matters here,
# and -- critically -- attempting the direct series *inside* that regime is
# not just wasted work but numerically treacherous (see `.g_core_raw`'s
# large-|w| comment): the terms grow to astronomic, non-converged magnitudes
# long before kmax is reached, rather than merely being slow.
UNDERFLOW_Z <- 100.0
UNDERFLOW_Z_SMALL <- 1.0 / UNDERFLOW_Z

.pochhammer <- function(c_val, k) {
  # (c)_k = c(c+1)...(c+k-1), k >= 1. Vectorized over k via cumulative product
  # is avoided here because callers need it inside a running series anyway;
  # this is only used for a fresh (non-incremental) evaluation.
  if (k == 0) return(c_val * 0 + 1)
  out <- c_val
  if (k > 1) for (i in 1:(k - 1)) out <- out * (c_val + i)
  out
}

#' Generalized hypergeometric `0F(n-1)`, `n-1 >= 1`.
#'
#' Evaluated by direct series summation with a running term
#' (`term_k = term_(k-1) * w / (k * prod_l(c_l+k-1))`).
#'
#' Entire in `w`, order 1/(n-1): for large |w| the denominator only grows
#' like k^n for large k, so the term ratio |w|/k^n doesn't drop below 1 until
#' k ~ |w|^(1/n) -- for the extreme end of this package's grid (the
#' inverse-gamma branch reflects small T into large 1/T, so |w| can reach the
#' 1e9s even though no physical energy/temperature grid point exceeds 1e4)
#' that is tens of thousands of terms, not a few hundred. `kmax` is sized
#' adaptively from `|w|` and `length(c_vec)` (capped by `kmax_cap` so a
#' single pathological call can't run unbounded); the second return value
#' records whether the series actually converged within that budget, so the
#' caller can apply the same "confidently negligible past a threshold, else
#' this is a real failure" contract `meijer.meijerg` uses in the Python
#' original, rather than silently returning an unconverged partial sum.
#' @noRd
.hyp0F <- function(w, c_vec, prec_bits, extra_tol_digits = 15, kmax_cap = 6000L) {
  w_num <- abs(as.numeric(w))
  n_c <- length(c_vec)
  kmax <- min(kmax_cap, max(4000L, ceiling(3 * w_num^(1 / (n_c + 1))) + 50L))

  w <- Rmpfr::mpfr(w, prec_bits)
  one <- Rmpfr::mpfr(1, prec_bits)
  term <- one
  total <- one
  max_abs_term <- abs(one)
  tol_digits <- floor(prec_bits * log10(2)) - extra_tol_digits
  tol <- Rmpfr::mpfr(10, prec_bits)^(-tol_digits)
  converged <- FALSE
  streak <- 0L
  streak_needed <- 5L  # guards against a spurious near-zero crossing of the
                        # running partial sum in an alternating series (large
                        # negative w) being mistaken for convergence -- a
                        # single-iteration check on abs(term) < tol*abs(total)
                        # can trip right as `total` passes through such a
                        # crossing, long before the series has actually
                        # reached its decaying tail.
  k <- 0L
  repeat {
    k <- k + 1L
    denom <- Rmpfr::mpfr(k, prec_bits)
    for (idx in seq_along(c_vec)) denom <- denom * (c_vec[idx] + (k - 1))
    term <- term * w / denom
    total <- total + term
    if (abs(term) > max_abs_term) max_abs_term <- abs(term)
    if (abs(term) < tol * abs(total)) {
      streak <- streak + 1L
      if (streak >= streak_needed) { converged <- TRUE; break }
    } else {
      streak <- 0L
    }
    if (k >= kmax) break
  }
  # A large negative w (odd-N branches reflect small T into a large, negative
  # (-1)^N * z) makes this an *alternating* series: the terms climb to a peak
  # around k ~ |w|^(1/(n_c+1)) before decaying, and if that peak dwarfs the
  # converged total, the accumulated rounding from summing it at `prec_bits`
  # can dominate the answer even though the term-ratio test above reports
  # "converged". Downgrade that verdict whenever the peak-to-total ratio
  # would have eaten more precision than was available.
  if (converged && as.numeric(max_abs_term) > 0) {
    if (as.numeric(abs(total)) == 0) {
      converged <- FALSE
    } else {
      cancellation_bits <- as.numeric(log2(max_abs_term)) - as.numeric(log2(abs(total)))
      if (cancellation_bits > prec_bits - 20) converged <- FALSE
    }
  }
  list(value = total, converged = converged, k = k)
}

#' Detect near-degenerate pairwise differences among `b` (any two entries
#' within `tol` of differing by an integer, including 0) -- the case where
#' `Gamma(b_i - b_j)` blows up and the matching `0F_{N-1}` term has a
#' cancelling zero (a removable 0*Inf), which is exactly the case mpmath's
#' own series also cannot handle (see statfactor._B_powerlaw_series's
#' docstring in the Python original).
#' @noRd
.degenerate_pairs <- function(b, tol = 1e-6) {
  n <- length(b)
  if (n < 2) return(FALSE)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d <- b[i] - b[j]
      if (abs(d - round(d)) < tol) return(TRUE)
    }
  }
  FALSE
}

#' Core routine for the gamma-chain Meijer-G shape, `z > 0`. Verified correct at a
#' *given* precision does not mean the terms weren't 50+ orders of magnitude
#' bigger than the answer (any nearly-degenerate `b`, or simply N >= 3 with a
#' large enough `z` that the alternating `0F_{N-1}` series climbs to a huge
#' peak before collapsing back down) -- and no single cheap heuristic on the
#' computed total reliably flags that (tried: comparing the peak term to the
#' computed total; it fails when the *computed* total is itself already pure
#' rounding noise rather than a visibly-too-small residual). So this uses the
#' standard robust fix instead: compute at increasing precision until two
#' consecutive attempts agree, doubling from `prec_bits` up to `max_bits`.
#' Degenerate `b` additionally needs the perturbation in `.g_core_raw` to
#' even approach the pole at every precision tried.
#'
#' `underflow_above` mirrors `meijer.meijerg`'s contract: it is consulted
#' only once escalation has failed up to `max_bits`, never as a shortcut --
#' the N=1 path never needs it at all (closed form, exact at any magnitude).
#' @noRd
.g_core <- function(z, b, prec_bits = DEFAULT_PREC_BITS, underflow_above = NULL,
                     max_bits = 2048L, agree_tol_digits = 12) {
  n <- length(b)
  if (n == 0) stop("g_gamma_chain: b must have at least one entry")
  if (z <= 0) stop(sprintf("Meijer G argument must be positive, got %s", z))
  if (n == 1) return(.g_core_raw(z, b, prec_bits)$value)  # exact; no escalation needed

  # Cheap short-circuit for the *unambiguous* tail, e.g. far quadrature nodes
  # in statfactor's B_invgamma_fallback (T ranging over 1e-30..1e30 of T0):
  # z an order of magnitude past the threshold is nowhere near the boundary
  # cases (like z = 200 vs. threshold 100) that turned out non-negligible
  # during development -- skip the expensive escalation loop entirely rather
  # than confirm the obvious the hard way.
  if (!is.null(underflow_above) && z >= 10 * underflow_above) {
    return(Rmpfr::mpfr(0, prec_bits))
  }

  # Degenerate b gets a *fixed* perturbation, sized once against `max_bits`
  # (not the current attempt's precision): comparing across precision levels
  # to detect insufficient precision only works if every attempt is
  # approximating the *same* nearby point. Tying delta's size to each
  # attempt's own (smaller) prec_bits -- tried first -- made every escalation
  # level perturb to a *different* point, so results didn't converge, they
  # just changed incomparably (caught via a b = (2,2,0) case returning three
  # wildly different, suspiciously round-looking values across 128/512/2048
  # bits during development -- a real bug, not a precision shortfall). So
  # degenerate b instead gets exactly one attempt, at `max_bits` directly,
  # self-diagnosed via `.g_core_raw`'s own convergence/cancellation check.
  if (.degenerate_pairs(b)) {
    # A single attempt directly at max_bits, self-diagnosed via
    # `.g_core_raw`'s own convergence/cancellation check -- NOT a
    # cross-precision comparison. Tried comparing a cheaper first rung
    # (512 bits) against max_bits first; rejected it after finding a case
    # (b = (2,2,0), z = 0.7) where the 512-bit attempt's *own* diagnostic
    # falsely reported success with a value that was flatly wrong (-1
    # instead of the true ~0.677, confirmed independently against mpmath) --
    # so a lower-precision attempt cannot be trusted even as a comparison
    # baseline here. Going straight to max_bits and trusting its diagnostic
    # (verified correct on that same case) is both simpler and cheaper than
    # a ladder that has to distrust its own first rung.
    delta <- Rmpfr::mpfr(10, max_bits)^(-77)
    attempt <- .g_core_raw(z, b, max_bits, perturb_delta = delta)
    # Extra, package-specific sanity check on top of the internal diagnostic
    # (which the 512-bit case above showed can itself be fooled): every `b`
    # this package ever builds (`a - 1` for the gamma chain, `a` with a `0`
    # appended for `statfactor.B`, both with `a > 0`) makes this shape a
    # positive multiple of a genuine probability density's Laplace
    # transform, so a negative or non-finite result is unambiguously wrong,
    # not just imprecise -- treat it as a convergence failure too.
    valid <- !attempt$any_nonconverged && is.finite(as.numeric(attempt$value)) &&
      as.numeric(attempt$value) > 0
    if (valid) return(Rmpfr::roundMpfr(attempt$value, prec_bits))
    if (!is.null(underflow_above) && z >= underflow_above) return(Rmpfr::mpfr(0, prec_bits))
    stop(sprintf(
      paste("g_gamma_chain: degenerate b = %s did not resolve at max_bits =",
            "%d even with a fixed perturbation, and no underflow_above",
            "threshold confirmed the true value is negligible at z = %s."),
      paste(format(b, digits = 6), collapse = ","), max_bits, format(z, digits = 6)
    ))
  }

  # Non-degenerate b: plain precision escalation, comparing consecutive
  # attempts for agreement (handles large-|w| cancellation unrelated to any
  # parameter degeneracy, e.g. N = 3 at z in the tens of thousands).
  agree_tol <- Rmpfr::mpfr(10, max_bits)^(-agree_tol_digits)
  first <- .g_core_raw(z, b, prec_bits)
  if (!first$any_nonconverged) return(first$value)
  bits <- prec_bits
  prev <- NULL
  repeat {
    cur <- .g_core_raw(z, b, bits)
    if (!cur$any_nonconverged && !is.null(prev)) {
      denom <- max(abs(cur$value), Rmpfr::mpfr(1e-300, bits))
      if (as.numeric(abs(cur$value - prev) / denom) < as.numeric(agree_tol)) {
        return(Rmpfr::roundMpfr(cur$value, prec_bits))
      }
    }
    if (bits >= max_bits) {
      if (!is.null(underflow_above) && z >= underflow_above) {
        return(Rmpfr::mpfr(0, prec_bits))
      }
      stop(sprintf(
        paste("g_gamma_chain: value at z = %s, b = %s did not stabilize by",
              "%d bits of precision (N = %d), and no underflow_above",
              "threshold confirmed the true value is negligible there."),
        format(z, digits = 6), paste(format(b, digits = 6), collapse = ","),
        max_bits, n
      ))
    }
    prev <- cur$value
    bits <- min(max_bits, bits * 4L)
  }
}

#' One attempt at `G^{N,0}_{0,N}(z|--;b)` at a fixed `prec_bits`, handling the
#' degenerate-`b` perturbation internally (perturbation size fixed relative
#' to `prec_bits` so escalating precision also sharpens the regularization).
#' Returns `list(value, any_nonconverged)`; never applies the underflow
#' threshold itself -- that is `.g_core`'s job, after escalation is exhausted.
#' @noRd
.g_core_raw <- function(z, b, prec_bits, perturb_delta = NULL) {
  n <- length(b)
  b_mp <- Rmpfr::mpfr(b, prec_bits)
  if (!is.null(perturb_delta)) {
    # Offsets by index (1, 2, ..., n)/n so no two entries move by the same
    # amount -- breaks not just an exactly-repeated pair but also any
    # secondary integer-spacing among the *other* pairs (e.g. b = (2,2,0)
    # has both b1 = b2 and b1 - b3 = b2 - b3 = 2, both integer).
    offsets <- Rmpfr::mpfr(seq_len(n), prec_bits) / n
    b_mp <- b_mp + Rmpfr::mpfr(perturb_delta, prec_bits) * offsets
  }
  z_mp <- Rmpfr::mpfr(z, prec_bits)
  sign_pow <- if (n %% 2 == 0) 1 else -1
  w <- sign_pow * z_mp

  if (n == 1) {
    return(list(value = z_mp^b_mp[1] * exp(w), any_nonconverged = FALSE))
  }

  total <- Rmpfr::mpfr(0, prec_bits)
  any_nonconverged <- FALSE
  for (j in seq_len(n)) {
    others <- setdiff(seq_len(n), j)
    prefactor <- Rmpfr::mpfr(1, prec_bits)
    for (i in others) prefactor <- prefactor * gamma(b_mp[i] - b_mp[j])
    c_vec <- 1 + b_mp[j] - b_mp[others]
    hyp <- .hyp0F(w, c_vec, prec_bits)
    if (!hyp$converged) any_nonconverged <- TRUE
    total <- total + z_mp^b_mp[j] * prefactor * hyp$value
  }
  list(value = total, any_nonconverged = any_nonconverged)
}

#' `G^{N,0}_{0,N}(z | -- ; b_1,...,b_N)` -- the gamma-chain shape (`s = 1/2`).
#'
#' Dies at *large* `z` (carries `e^{-z}`): pass `underflow_above` (mirroring
#' `meijer.UNDERFLOW_Z` in the Python original) to short-circuit to exactly 0
#' past a threshold where the true value is confidently negligible, rather
#' than attempting the direct series there (see `.g_core_raw`'s docs for why
#' that regime is numerically treacherous, not just slow).
#'
#' @param z Positive numeric scalar.
#' @param b Numeric vector of length N (the lower/denominator parameters).
#' @param prec_bits Working precision in bits (default ~38 decimal digits).
#' @param underflow_above Optional threshold; `z` at or above it returns 0.
#' @return An `mpfr` number.
#' @export
g_gamma_chain <- function(z, b, prec_bits = DEFAULT_PREC_BITS, underflow_above = NULL) {
  .g_core(z, b, prec_bits, underflow_above = underflow_above)
}

#' `G^{0,N}_{N,0}(z | b_1,...,b_N ; --)` -- the inverse-gamma-chain shape
#' (`s = 1`), via the reflection identity `= G^{N,0}_{0,N}(1/z | -- ; 1-b)`.
#'
#' Dies at *small* `z` (carries `e^{-1/z}`); pass `underflow_below` (mirroring
#' `meijer.UNDERFLOW_Z_SMALL`) to short-circuit to 0 below a threshold --
#' translated internally to the `underflow_above` threshold `1/underflow_below`
#' on the reflected argument `1/z`.
#'
#' @inheritParams g_gamma_chain
#' @param underflow_below Optional threshold; `z` at or below it returns 0.
#' @export
g_invgamma_chain <- function(z, b, prec_bits = DEFAULT_PREC_BITS, underflow_below = NULL) {
  if (z <= 0) stop(sprintf("Meijer G argument must be positive, got %s", z))
  underflow_above <- if (is.null(underflow_below)) NULL else 1 / underflow_below
  .g_core(1 / z, 1 - b, prec_bits, underflow_above = underflow_above)
}
