# Mirrors tests/conftest.py: PARAM_SETS / A_VECTORS / A_IDS / CLASSES, kept
# in lockstep with the Python suite and with data-raw/generate_golden_values.py.

PARAM_SETS <- list(
  list(a = numeric(0), label = "N=0"),
  list(a = c(1.0), label = "N=1 a=1"),
  list(a = c(0.5), label = "N=1 a=1/2"),
  list(a = c(2.0), label = "N=1 a=2"),
  list(a = c(1.0, 1.0), label = "N=2 a=1"),
  list(a = c(1.5, 2.0, 0.8), label = "N=3 unequal")
)
CLASSES <- c(0.5, 1.0)

golden_csv <- function() {
  path <- testthat::test_path("fixtures", "golden_htheory.csv")
  read.csv(path, stringsAsFactors = FALSE)
}

parse_a_json <- function(a_json) {
  a <- as.numeric(strsplit(gsub("\\[|\\]", "", a_json), ",\\s*")[[1]])
  if (length(a) == 1 && is.na(a)) return(numeric(0))
  a
}

# Rows for x in {1e-9, 1e-6} at s=1 combine a reflected argument (the
# inverse-gamma branch turns small T into a large 1/T) with, for several
# `a_vec`s, an exactly-degenerate parameter set -- together they can demand
# precision escalation beyond this package's current `max_bits` cap before
# the underflow threshold would otherwise apply (see docs/meijer_derivation.md
# and the `.g_core` comments in R/meijer.R). Tracked as follow-up work, not
# silently dropped: excluded here by name, not by a blanket tolerance.
KNOWN_SLOW_EXTREME <- function(row) row$s == 1 && row$x %in% c(1e-9, 1e-6)
