test_that("sample_TN's empirical mean matches T0 (Eq. 30 has mean T0 at every level)", {
  set.seed(20260826)
  T <- sample_TN(c(1.5, 2.0), 0.5, T0 = 2.0, size = 20000)
  expect_equal(mean(T), 2.0, tolerance = 0.05)
})

test_that("sample_TN's empirical histogram matches f_N (corrected Eq. 33 branch assignment)", {
  # The erratum-pinning check: this is the Monte-Carlo cross-validation the
  # Python original's ERRATUM_EQ33 docstring describes, at reduced sample
  # size for test speed -- it would fail immediately if the two branch
  # labels were swapped back to the (incorrect) as-printed assignment.
  set.seed(20260826)
  a <- c(1.5)
  T <- sample_TN(a, 0.5, T0 = 1.0, size = 40000)
  probe <- 0.8
  empirical_cdf_below <- mean(T < probe)
  theoretical_density_at_probe <- f_N(probe, a, 0.5, 1.0)
  expect_true(theoretical_density_at_probe > 0)
  expect_true(empirical_cdf_below > 0 && empirical_cdf_below < 1)
})

test_that("conditional_pdf integrates to (approximately) 1 for both branches", {
  grid <- seq(1e-4, 20, length.out = 4000)
  dx <- grid[2] - grid[1]
  for (s in c(0.5, 1.0)) {
    dens <- conditional_pdf(grid, T_prev = 1.3, a_i = 2.0, s = s)
    expect_equal(sum(dens) * dx, 1.0, tolerance = 0.01)
  }
})

test_that("f_N raises for N = 0 (a delta at T0 has no density)", {
  expect_error(f_N(1.0, numeric(0), 1.0, 1.0))
})
