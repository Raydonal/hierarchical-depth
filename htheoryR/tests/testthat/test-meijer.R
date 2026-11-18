test_that("G^{1,0}_{0,1}(z|--;0) is the exponential", {
  z <- 2.3
  expect_equal(as.numeric(g_gamma_chain(z, c(0))), exp(-z), tolerance = 1e-12)
})

test_that("G^{2,0}_{0,2}(z|--;a,0) is the K-distribution form 2*z^(a/2)*K_a(2*sqrt(z))", {
  a <- 0.7; z <- 3.1
  got <- as.numeric(g_gamma_chain(z, c(a, 0)))
  expected <- 2 * z^(a / 2) * besselK(2 * sqrt(z), a)
  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("degenerate b = (0,0) does not error and is internally consistent", {
  g <- g_gamma_chain(1.7, c(0, 0))
  expect_true(is.finite(as.numeric(g)))
  expect_true(as.numeric(g) > 0)
})

test_that("g_invgamma_chain is the reflection of g_gamma_chain", {
  b <- c(1.3, 0.4); z <- 0.6
  direct <- g_invgamma_chain(z, b)
  via_reflection <- g_gamma_chain(1 / z, 1 - b)
  expect_equal(as.numeric(direct), as.numeric(via_reflection), tolerance = 1e-12)
})

test_that("f_N matches the Python mpmath oracle across the golden grid", {
  csv <- golden_csv()
  rows <- csv[csv$func == "f_N" & csv$ok == "True", ]
  skip_if(nrow(rows) == 0, "golden_htheory.csv not generated")

  n_checked <- 0
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    if (KNOWN_SLOW_EXTREME(row)) next
    a <- parse_a_json(row$a_json)
    expected <- as.numeric(row$expected_value_str)
    got <- f_N(row$x, a, row$s, row$T0)
    if (abs(expected) < 1e-290) next  # both sides underflow double; not meaningful
    n_checked <- n_checked + 1
    expect_equal(got, expected,
                 tolerance = 1e-6,
                 label = sprintf("f_N(x=%s, a=%s, s=%s)", row$x, row$a_json, row$s))
  }
  expect_gt(n_checked, 50)  # sanity: the grid actually exercised something
})

test_that("f_N has no density at N = 0 (matches the Python original's ValueError)", {
  expect_error(f_N(1.0, numeric(0), 0.5, 1.0))
})

test_that("B and B_tilde match the Python mpmath oracle across the golden grid", {
  # Slow (~4.5 min): several golden rows have a degenerate b (e.g. a = (1,1)),
  # and B's s=1 fallback for those runs a numerical quadrature whose every
  # node can itself be an expensive degenerate Meijer-G evaluation. Gated the
  # same way the Python suite gates its own mpmath-heavy cross-validation
  # (`pytest -m "not slow"`) -- run with
  # `HTHEORY_RUN_SLOW=true` before trusting a change to R/meijer.R or
  # R/statfactor.R.
  skip_if_not(identical(Sys.getenv("HTHEORY_RUN_SLOW"), "true"),
              "slow: set HTHEORY_RUN_SLOW=true to run")
  csv <- golden_csv()
  for (func_name in c("B", "B_tilde")) {
    rows <- csv[csv$func == func_name & csv$ok == "True", ]
    skip_if(nrow(rows) == 0, "golden_htheory.csv not generated")
    n_checked <- 0
    for (i in seq_len(nrow(rows))) {
      row <- rows[i, ]
      if (row$s == 1 && row$x %in% c(1e-9, 1e-6)) next  # see KNOWN_SLOW_EXTREME
      a <- parse_a_json(row$a_json)
      expected <- as.numeric(row$expected_value_str)
      got <- if (func_name == "B") B(row$x, a, row$s, row$beta0) else B_tilde(row$x, a, row$s, row$beta0)
      if (abs(expected) < 1e-290) next
      n_checked <- n_checked + 1
      expect_equal(got, expected,
                   tolerance = 1e-6,
                   label = sprintf("%s(x=%s, a=%s, s=%s)", func_name, row$x, row$a_json, row$s))
    }
    expect_gt(n_checked, 30)
  }
})

test_that("B(0) = 1 and B is monotonically decreasing (Eq. 7 boundary conditions)", {
  # N = 0 and N = 1 only in the fast default run: any `a` containing an
  # integer makes statfactor.B's s=0.5 shape degenerate (it always appends a
  # literal 0, see R/statfactor.R), and the s=1 branch's quadrature fallback
  # is expensive per-node -- both cases are slow one at a time, and this test
  # calls B() vectorized over 6 E values, multiplying that cost. Covered at
  # a = (1.5, 2, 0.8) in the slow, opt-in pass below instead.
  for (a in list(numeric(0), c(1.0))) {
    for (s in c(0.5, 1.0)) {
      expect_equal(B(0, a, s), 1.0)
      Es <- c(0.1, 0.5, 1, 2, 5, 10)
      vals <- B(Es, a, s)
      expect_true(all(diff(vals) < 0))
    }
  }
})

test_that("B(0) = 1 and B is monotonically decreasing for N = 3, unequal a (slow)", {
  skip_if_not(identical(Sys.getenv("HTHEORY_RUN_SLOW"), "true"),
              "slow: set HTHEORY_RUN_SLOW=true to run")
  a <- c(1.5, 2.0, 0.8)
  for (s in c(0.5, 1.0)) {
    expect_equal(B(0, a, s), 1.0)
    Es <- c(0.1, 0.5, 1, 2, 5, 10)
    vals <- B(Es, a, s)
    expect_true(all(diff(vals) < 0))
  }
})

test_that("B_tilde shifts the vector a but leaves omega untouched (Eq. 38)", {
  a <- c(1.0, 1.0)
  s <- 0.5
  direct <- B_tilde(0.7, a, s)
  manual <- B(0.7, a + shift(s), s, omega = omega_of(a))
  expect_equal(direct, manual, tolerance = 1e-10)
})
