# n_table kept small here (default 800 costs minutes given Rmpfr's overhead
# per B() evaluation, see R/meijer.R's docs) -- this is the round-trip
# correctness check, not a precision benchmark.
test_that("maximize recovers p_target = B(E)/Z under the standard constraints", {
  skip_if_not_installed("nloptr")
  # Slow (~30s even at a reduced n_table=100): every knot in the interpolant
  # table is a Meijer-G evaluation (R/meijer.R), and this is the package's
  # highest-value regression test (closes the branch assignment, the
  # omega-vs-a bookkeeping, the alpha=-1 gauge, and the 1/Z scaling all at
  # once), not something to shrink further just for a fast default run.
  skip_if_not(identical(Sys.getenv("HTHEORY_RUN_SLOW"), "true"),
              "slow: set HTHEORY_RUN_SLOW=true to run")
  energies <- c(0, 0.5, 1, 2, 5)
  res <- maximize(energies, c(0.5), 0.5, n_table = 100L)
  expect_true(res$success)
  expect_lt(res$max_abs_error, 1e-2)
  expect_equal(res$entropy, res$entropy_target, tolerance = 1e-3)
})

test_that("stationarity_residual is small at the target distribution", {
  energies <- c(0, 0.5, 1, 2, 5)
  a <- c(0.5); s <- 0.5
  target <- target_distribution(energies, a, s)
  resid <- stationarity_residual(target$p, energies, a, s, Z = target$Z)
  expect_lt(abs(resid), 1e-6)
})
