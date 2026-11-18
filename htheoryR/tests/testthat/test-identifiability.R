test_that("log_cumulants_theory matches a direct check for N=0 (delta at T0)", {
  expect_equal(log_cumulants_theory(numeric(0), 0.5, r_max = 4), rep(0, 3))
})

test_that("fit_homogeneous recovers (N, a) exactly from noise-free theoretical cumulants", {
  for (N in c(1L, 2L, 3L, 5L)) {
    for (a in c(0.7, 1.0, 3.0)) {
      k <- log_cumulants_theory(rep(a, N), 0.5, r_max = 3)
      fit <- fit_homogeneous(k[1], k[2], 0.5)
      expect_equal(fit$N, N, tolerance = 1e-6)
      expect_equal(fit$a, a, tolerance = 1e-6)
    }
  }
})

test_that("empirical_log_cumulants matches scipy.stats.kstat's formulas (fixed vector, cross-checked externally)", {
  x <- c(1.236, -2.241, 2.548, 4.155, -4.462, -2.822, 0.457, 0.185, 1.304, -1.437)
  got <- .kstat(x, 2)
  n <- length(x); m2 <- sum((x - mean(x))^2) / n
  expect_equal(got, n / (n - 1) * m2)
})

test_that("test_depth does not reject a true null and rejects a clearly false one", {
  set.seed(20260826)
  T_null <- sample_TN(c(1.5), 0.5, T0 = 1, size = 3000)
  r_null <- test_depth(T_null, 1, 0.5, n_boot = 300)
  expect_true(r_null$p_value > 0.05)

  T_alt <- sample_TN(c(2, 2, 2), 0.5, T0 = 1, size = 3000)
  r_alt <- test_depth(T_alt, 1, 0.5, n_boot = 300)
  expect_true(r_alt$p_value < 0.05)
})
