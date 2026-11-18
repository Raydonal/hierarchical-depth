#' Coerce the shape parameter to a numeric vector (N = 0 gives length 0).
#'
#' `require_positive` is on for the physical parameter `a_i = 2 gamma_i /
#' kappa_i^2`, which is positive by construction. It is off only for the
#' *shifted* vector `a + 3 - 4s` of Eq. (38), where the `s = 1` shift is `-1`
#' and drives entries to zero or below whenever `a_i <= 1` -- still a
#' perfectly admissible Meijer G parameter, just no longer a reservoir shape.
#'
#' Mirrors `htheory.background._as_a` (Python) verbatim, including the error
#' message content, so ported tests can assert on the same substrings.
#' @keywords internal
#' @noRd
.as_a <- function(a, require_positive = TRUE) {
  arr <- as.numeric(a)
  if (!is.null(dim(a)) && length(dim(a)) > 1) {
    stop(sprintf("a must be a scalar or 1-D sequence, got shape %s",
                  paste(dim(a), collapse = "x")))
  }
  if (require_positive && length(arr) && any(arr <= 0)) {
    stop(sprintf("a must be strictly positive, got %s", paste(arr, collapse = ", ")))
  }
  arr
}

#' `omega = prod(a_i)`, Eq. (33). The scalar that Eq. (38) leaves alone.
#' @keywords internal
#' @noRd
.omega <- function(a) {
  if (length(a) == 0) return(1)
  prod(a)
}

#' Validate the universality-class selector. Mirrors the Python
#' `if s == 0.5 / elif s == 1 / else raise` style verbatim (no enum).
#' @keywords internal
#' @noRd
.check_s <- function(s) {
  if (!isTRUE(all.equal(s, 0.5)) && !isTRUE(all.equal(s, 1))) {
    stop(sprintf("s must be 0.5 or 1, got %s", s))
  }
  invisible(s)
}

#' The parameter shift of Eq. (38): `a -> a + 3 - 4s`. `+1` for `s = 1/2`,
#' `-1` for `s = 1`.
#' @keywords internal
#' @noRd
.shift <- function(s) {
  .check_s(s)
  if (isTRUE(all.equal(s, 0.5))) 1 else -1
}

#' Lift a scalar function to accept and return numeric vectors/arrays.
#'
#' Scalar input gives a scalar `double` back; array input gives a `double`
#' array of the same shape. This is the direct analogue of Python's
#' `meijer.vectorize`, built on `vapply` rather than `base::Vectorize()` so
#' the scalar-in/scalar-out, array-in/same-shape-array-out contract is
#' enforced by the type signature, not merely likely.
#'
#' @param scalar_fn A function taking a scalar (plus optional extra
#'   arguments) and returning a scalar.
#' @return A function with the same extra arguments as `scalar_fn`, taking
#'   `x` as a scalar or array and returning a value of the same shape.
#' @export
vectorize_scalar <- function(scalar_fn) {
  function(x, ...) {
    orig_dim <- dim(x)
    flat <- vapply(as.vector(x), function(v) as.double(scalar_fn(v, ...)), numeric(1))
    if (length(flat) == 1) return(flat[[1]])
    if (!is.null(orig_dim)) dim(flat) <- orig_dim
    flat
  }
}


#' Stringify a numeric value at a given number of significant decimal digits,
#' for use as an exact-match cache key (mirrors `mp.nstr(z, dps + 5, ...)`).
#' @keywords internal
#' @noRd
.num_key <- function(x, digits = 25) {
  formatC(x, digits = digits, format = "e")
}
