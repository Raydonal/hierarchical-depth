test_that("fit_homogeneous_network recovers a known depth from a large synthetic sample", {
  set.seed(7)
  D <- sample_degrees(rep(3, 3), 0.5, lam = 15, size = 400000)
  fit <- fit_homogeneous_network(D, 0.5)
  expect_equal(fit$N_round, 3L)
  expect_equal(fit$a, 3, tolerance = 0.2)
})

test_that("fit_homogeneous_network distinguishes N=1 from N=3 at matched mean", {
  set.seed(8)
  D1 <- sample_degrees(rep(3, 1), 0.5, lam = 15, size = 200000)
  D3 <- sample_degrees(rep(3, 3), 0.5, lam = 15, size = 200000)
  fit1 <- fit_homogeneous_network(D1, 0.5)
  fit3 <- fit_homogeneous_network(D3, 0.5)
  expect_equal(fit1$N_round, 1L)
  expect_equal(fit3$N_round, 3L)
})

test_that("test_depth_degrees does not reject the true depth and rejects a false one", {
  set.seed(9)
  D_true <- sample_degrees(rep(2, 1), 0.5, lam = 20, size = 5000)
  r_true <- test_depth_degrees(D_true, 1, 0.5, n_boot = 300)
  expect_true(r_true$p_value > 0.05)

  D_false <- sample_degrees(rep(2, 4), 0.5, lam = 20, size = 5000)
  r_false <- test_depth_degrees(D_false, 1, 0.5, n_boot = 300)
  expect_true(r_false$p_value < 0.05)
})
