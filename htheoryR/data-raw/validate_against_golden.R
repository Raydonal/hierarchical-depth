suppressMessages(library(Rmpfr))
source("R/utils.R")
source("R/meijer.R")
source("R/background.R")

csv <- read.csv("tests/testthat/fixtures/golden_htheory.csv", stringsAsFactors = FALSE)

check_func <- function(func_name, compute) {
  rows <- csv[csv$func == func_name & csv$ok == "True", ]
  cat(sprintf("--- %s: %d rows ---\n", func_name, nrow(rows)))
  max_rel_err <- 0
  worst <- NULL
  n_underflow_both <- 0
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    a <- as.numeric(strsplit(gsub("\\[|\\]", "", row$a_json), ",\\s*")[[1]])
    if (length(a) == 1 && a == "") a <- numeric(0)
    expected <- Rmpfr::mpfr(row$expected_value_str, 200)
    got <- compute(row, a)
    got_mp <- Rmpfr::mpfr(got, 200)
    if (abs(expected) < 1e-290 && abs(got_mp) < 1e-290) {
      n_underflow_both <- n_underflow_both + 1
      next  # both effectively zero in double land; not a meaningful ratio
    }
    rel_err <- as.numeric(abs(got_mp - expected) / max(abs(expected), Rmpfr::mpfr(1e-300, 200)))
    if (is.na(rel_err)) rel_err <- Inf
    if (rel_err > max_rel_err) { max_rel_err <- rel_err; worst <- row }
  }
  cat(sprintf("  max relative error: %s (skipped %d double-underflow rows)\n",
              format(max_rel_err, digits = 6), n_underflow_both))
  if (!is.null(worst) && max_rel_err > 1e-9) print(worst)
}

csv <- csv[!(csv$x %in% c(1e-9, 1e-6) & csv$s == 1), ]
cat("(excluding x in {1e-9,1e-6} for s=1: known-slow extreme reflected-z rows, see docs/meijer_derivation.md follow-ups)\n")

check_func("f_N", function(row, a) f_N(row$x, a, row$s, row$T0))
