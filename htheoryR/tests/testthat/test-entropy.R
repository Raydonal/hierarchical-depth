# All tests here use single-element, non-integer `a` (e.g. c(0.5)) to avoid
# the degenerate-parameter path (R/meijer.R), which is correct but expensive
# (see its comments) -- irrelevant to what these tests check.

test_that("s_tilde(0) = 0 (Eq. 39 boundary: B(0)=1, Btilde(0)=1)", {
  expect_equal(s_tilde(0, c(0.5), 0.5), 0, tolerance = 1e-12)
  expect_equal(s_tilde(0, c(0.7), 1.0), 0, tolerance = 1e-12)
})

test_that("E_of_B inverts B", {
  for (Eval in c(0.3, 1.3, 5.0)) {
    Bval <- B(Eval, c(0.5), 0.5)
    expect_equal(E_of_B(Bval, c(0.5), 0.5), Eval, tolerance = 1e-8)
  }
})

test_that("s_of_B matches the exact (non-tabulated) s_tilde composition", {
  Eval <- 1.3
  Bval <- B(Eval, c(0.5), 0.5)
  expect_equal(s_of_B(Bval, c(0.5), 0.5), s_tilde(Eval, c(0.5), 0.5), tolerance = 1e-10)
})

test_that("make_interpolants reproduces exact E_of_B/s_of_B to ~1e-4 across the B range", {
  interp <- make_interpolants(c(0.5), 0.5, n = 60)
  for (Bval in c(0.9, 0.5, 0.1, 0.01, 1e-4)) {
    expect_equal(interp$E_of_B(Bval), E_of_B(Bval, c(0.5), 0.5), tolerance = 1e-3)
    expect_equal(interp$s_of_B(Bval), s_of_B(Bval, c(0.5), 0.5), tolerance = 1e-3)
  }
})

test_that("q_indices and tsallis_entropy match the N=1 s=1 exact special case (Eq. 65-67)", {
  a <- 1.3
  idx <- q_indices(a)
  expect_equal(idx$q, a / (a + 1))
  expect_equal(idx$q_stat, 2 - idx$q)
  expect_true(idx$q < 1)
  expect_true(idx$q_stat > 1)

  p <- c(0.2, 0.3, 0.5)
  # tsallis_entropy(p, a) should equal sum_i s_tsallis_kernel(p_i, a) at Z=1
  # (Eq. 64/67) -- the closed form s_tilde reduces to at N=1, s=1.
  expect_equal(tsallis_entropy(p, a), sum(s_tsallis_kernel(p, a)), tolerance = 1e-10)
})

test_that("s_boltzmann_shannon is -B*log(B)", {
  expect_equal(s_boltzmann_shannon(0.3), -0.3 * log(0.3))
})
