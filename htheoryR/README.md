# htheoryR <img src="man/figures/logo.png" align="right" height="139" alt="htheoryR logo" />

Native R implementation of the H-theory construction in **"Hierarchical
entropies"** (de Moraes, Brum, Macêdo & Vasconcelos, *Physica A* **699**
(2026) 131899), plus the identifiability theory, estimators, and real-data
applications that extend it in an accompanying manuscript (in preparation).
The package lives in the `htheoryR/` folder of the repository
<https://github.com/Raydonal/hierarchical-depth>, next to the Python package
and the pipeline that reproduces the manuscript.

This is a from-scratch R port, not a wrapper: the Meijer-G special-function
core (`R/meijer.R`) is a closed-form reduction derived and verified
independently (`docs/meijer_derivation.md` in the parent repository), backed
by `Rmpfr` for arbitrary precision, and cross-checked against the trusted
Python package (`../htheory`) via a golden-value oracle
(`data-raw/GOLDEN_PROVENANCE.md`).

## What it reproduces

- **The original paper's central claim** — maximizing the entropy `S[p]`
  under the standard constraints (normalization and mean energy only)
  returns exactly `p_i = B(E_i)/Z` — via `maximize()` and
  `stationarity_residual()`.
- **The identifiability, estimation, and inference results** developed on
  top of it: the cascade depth `N` is identifiable from the log-cumulants of
  a single marginal (`log_cumulants_theory()`, `fit_homogeneous()`), with a
  calibrated parametric-bootstrap test (`test_depth()`) and measured power
  (`power_curve()`).
- **The network application** — the same identifiability chain, through
  Feller's (1943) mixed-Poisson identifiability, recovers a hierarchical
  depth from a plain degree sequence alone (`fit_homogeneous_network()`,
  `test_depth_degrees()`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("Raydonal/hierarchical-depth", subdir = "htheoryR")
```

`Rmpfr` needs GMP and MPFR development headers at compile time (Debian/Ubuntu:
`sudo apt install libgmp-dev libmpfr-dev`); `nloptr` is used for the
maximum-entropy optimizer.

## A minimal reproduction

The paper's own worked example (a stretched-exponential cascade, depth
`N = 1`, shape `a = 0.5`) round-tripped through the whole chain:

```r
library(htheoryR)

energies <- c(0, 0.5, 1, 2, 5, 10)
a <- c(0.5); s <- 0.5

# B(E) -> the target distribution p_i = B(E_i)/Z
target <- target_distribution(energies, a, s)

# Maximize S[p] under normalization + mean energy alone, and recover p_target
fit <- maximize(energies, a, s, n_table = 200L)
max(abs(fit$p - fit$p_target))  # small: the round trip closes
```

Identifiability, applied to a synthetic three-level homogeneous cascade --
the depth is recovered from the log-cumulants of the innermost temperature
alone, with no access to the intermediate levels (`fit_homogeneous` targets
exactly this homogeneous case, `a_i = a`; a heterogeneous cascade returns an
effective homogeneous depth rather than the literal per-level shapes):

```r
set.seed(1)
T <- sample_TN(c(2, 2, 2), s = 0.5, size = 200000)
k <- empirical_log_cumulants(T, r_max = 3)
fit_homogeneous(k[1], k[2], s = 0.5)  # recovers N_round = 3, a close to 2

# Is a shallower depth (N = 1) rejected by the data?
test_depth(T, N0 = 1, s = 0.5, n_boot = 500)  # p_value = 0: rejected
```

The network application, applied to a synthetic mixed-Poisson degree
sequence -- the depth is recovered from the degree sequence alone, without
ever observing which vertex belongs to which community:

```r
set.seed(2)
D <- sample_degrees(a = c(3, 3, 3), s = 0.5, lam = 15, size = 200000)
fit_homogeneous_network(D, s = 0.5)          # point estimate
test_depth_degrees(D, N0 = 1, s = 0.5, n_boot = 500)  # rejects a shallow null
```

## Real-data results from the manuscript

The four case studies of the manuscript (`email-Eu-core`, OpenFlights,
wind-tunnel turbulence, and the *C. elegans* connectome) were run with the
Python package's equivalent functions (`scripts/analyze_*.py` in the
repository, results in `results/*.json`) rather than this package, since this
port's `identifiability`/`network` modules were completed after those
analyses. Re-running them with `htheoryR`, to confirm that the two
implementations agree on the reported numbers, is open work.

## Performance note

Every Meijer-G evaluation in `R/meijer.R` costs one or more `Rmpfr`
operations, which are slower per call than `mpmath`'s compiled series
(roughly 0.1-0.5s per call for non-degenerate parameters, more for
degenerate ones needing escalated precision). Functions that build a
tabulated interpolant over hundreds of points (`make_interpolants()`,
`maximize()`) can take tens of seconds to minutes as a result — this is a
known, accepted characteristic of the current implementation, not a bug;
see the comments in `R/meijer.R` for the numerical reasoning behind the
precision-escalation strategy. The test suite marks the slowest checks with
`Sys.getenv("HTHEORY_RUN_SLOW") == "true"`, matching the parent Python
project's own `pytest -m slow` convention — run

```r
Sys.setenv(HTHEORY_RUN_SLOW = "true")
testthat::test_dir("tests/testthat")
```

before trusting a change to `R/meijer.R`, `R/statfactor.R`, or `R/maxent.R`.

## Relationship to the Python package

`../htheory` (Python, `mpmath`/`numpy`/`scipy`) is the original, trusted
reference implementation and remains the oracle this package's golden-value
tests are checked against (`data-raw/generate_golden_values.py`,
`data-raw/GOLDEN_PROVENANCE.md`). The two are kept independent on purpose —
`htheoryR` does not call into Python at runtime (no `reticulate`) — so that
numerical agreement between them is a genuine cross-implementation check,
not a shared-code assumption.
