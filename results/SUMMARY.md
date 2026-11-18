# Summary of the regenerated results

Produced by `make results tables figures` (fixed seeds; see `REPRODUCING.md`).
Wall times on this machine: `make results` 655 s (10.9 min), `make tables`
0.8 s, `make figures` 67 s, `make test-fast` 120 s (393 passed, 2 skipped,
38 deselected).

Conventions: 1000 parametric-bootstrap draws per depth test; 95% percentile
intervals from 1000 nonparametric resamples (vertices for networks, moving
blocks for turbulence). "CI fail" is the fraction of resamples on which the fit
was not admissible (the interval is over the successful ones). The depth test
is of `H0: N = N_round` with `N_round` chosen from the same data, so `p` is a
goodness-of-fit diagnostic, not a confirmatory test.

## 1. Empirical results (every series, both classes)

Dispersion column: variance-to-mean ratio for networks, kappa_2 of the
log-series for turbulence.

| System | Series | n | VMR / kappa_2 | s | a_hat [95% CI] | N_hat [95% CI] | N_round | p | CI fail |
|---|---|---|---|---|---|---|---|---|---|
| email-Eu-core | out-degree | 1005 | 43.6 | 1/2 | 0.768 [0.396, 1.27] | 1.18 [0.731, 1.75] | 1 | 0.391 | 0.000 |
| email-Eu-core | out-degree | 1005 | 43.6 | 1 | no admissible fit | — | — | — | — |
| email-Eu-core | in-degree | 1005 | 30.7 | 1/2 | 0.802 [0.555, 1.13] | 0.956 [0.727, 1.24] | 1 | 0.825 | 0.000 |
| email-Eu-core | in-degree | 1005 | 30.7 | 1 | no admissible fit | — | — | — | — |
| email-Eu-core | total-degree | 1005 | 71.0 | 1/2 | 0.825 [0.495, 1.29] | 1.09 [0.743, 1.58] | 1 | 0.645 | 0.000 |
| email-Eu-core | total-degree | 1005 | 71.0 | 1 | no admissible fit | — | — | — | — |
| OpenFlights | out-degree | 3315 | 144 | 1/2 | 0.216 [0.126, 0.333] | 1.21 [0.930, 1.53] | 1 | 0.194 | 0.000 |
| OpenFlights | out-degree | 3315 | 144 | 1 | no admissible fit | — | — | — | — |
| OpenFlights | in-degree | 3321 | 144 | 1/2 | 0.215 [0.127, 0.334] | 1.20 [0.938, 1.54] | 1 | 0.193 | 0.000 |
| OpenFlights | in-degree | 3321 | 144 | 1 | no admissible fit | — | — | — | — |
| OpenFlights | total-degree | 3330 | 288 | 1/2 | 0.212 [0.125, 0.317] | 1.20 [0.920, 1.50] | 1 | 0.193 | 0.000 |
| OpenFlights | total-degree | 3330 | 288 | 1 | no admissible fit | — | — | — | — |
| C. elegans | out-degree (chemical, weighted) | 253 | 16.7 | 1/2 | 8.06 [1.03, 61.0] | 4.13 [0.570, 28.4] | 4 | 0.989 | 0.047 |
| C. elegans | out-degree | 253 | 16.7 | 1 | no admissible fit | — | — | — | — |
| C. elegans | in-degree (chemical, weighted) | 268 | 42.5 | 1/2 | 1.09 [0.409, 1.81] | 1.55 [0.667, 2.23] | 2 | 0.770 | 0.000 |
| C. elegans | in-degree | 268 | 42.5 | 1 | no admissible fit | — | — | — | — |
| C. elegans | total-degree (out + in + gap) | 279 | 51.6 | 1/2 | no admissible fit | — | — | — | — |
| C. elegans | total-degree | 279 | 51.6 | 1 | 7.55 [3.92, 134] | 4.77 [1.85, 98.8] | 5 | 0.995 | 0.264 |
| Turbulence | lag 2 ms, full | 131035 | 4.72 | 1/2 | 0.595 [0.583, 0.607] | 1.28 [1.24, 1.31] | 1 | <0.001 | 0.000 |
| Turbulence | lag 2 ms, thinned k=2 | 65519 | 4.71 | 1/2 | 0.591 [0.575, 0.609] | 1.26 [1.22, 1.31] | 1 | <0.001 | 0.000 |
| Turbulence | lag 5 ms, full | 131043 | 3.78 | 1/2 | 0.503 [0.494, 0.513] | 0.775 [0.757, 0.794] | 1 | <0.001 | 0.000 |
| Turbulence | lag 5 ms, thinned k=5 | 26209 | 3.72 | 1/2 | 0.513 [0.489, 0.536] | 0.788 [0.741, 0.833] | 1 | <0.001 | 0.000 |
| Turbulence | lag 10 ms, full | 131027 | 4.66 | 1/2 | 0.536 [0.527, 0.547] | 1.06 [1.04, 1.09] | 1 | 0.001 | 0.000 |
| Turbulence | lag 10 ms, thinned k=10 | 13104 | 4.66 | 1/2 | 0.537 [0.506, 0.571] | 1.06 [0.986, 1.15] | 1 | 0.254 | 0.000 |
| Turbulence | lag 20 ms, full | 131029 | 3.83 | 1/2 | 0.491 [0.482, 0.500] | 0.751 [0.733, 0.769] | 1 | <0.001 | 0.000 |
| Turbulence | lag 20 ms, thinned k=20 | 6553 | 3.83 | 1/2 | 0.474 [0.436, 0.518] | 0.708 [0.631, 0.797] | 1 | 0.002 | 0.000 |
| Turbulence | lag 50 ms, full | 130989 | 4.66 | 1/2 | 0.522 [0.513, 0.532] | 1.02 [0.989, 1.05] | 1 | 0.354 | 0.000 |
| Turbulence | lag 50 ms, thinned k=50 | 2620 | 4.64 | 1/2 | 0.499 [0.435, 0.582] | 0.937 [0.780, 1.15] | 1 | 0.578 | 0.000 |
| Turbulence | lag 100 ms, full | 130945 | 4.71 | 1/2 | 0.519 [0.510, 0.529] | 1.02 [0.990, 1.05] | 1 | 0.385 | 0.000 |
| Turbulence | lag 100 ms, thinned k=100 | 1309 | 4.40 | 1/2 | 0.604 [0.538, 0.689] | 1.22 [1.04, 1.43] | 1 | 0.202 | 0.000 |

Turbulence, class s = 1: no admissible fit at any lag, full or thinned
(kappa_3/kappa_2 is negative, outside the range of the inverse-gamma class).

### How this differs from the old numbers

| Item | Old | New |
|---|---|---|
| email-Eu-core, out-degree | a_hat = 0.77, N_round = 1, p = 0.39 | unchanged (0.768, 1, 0.39); now with CI N_hat in [0.73, 1.75]; in- and total-degree are also tested (N_round = 1, p = 0.82 and 0.65) |
| OpenFlights | a_hat about 0.21–0.22, N_round = 1, p = 0.16–0.20 | unchanged within Monte Carlo error (0.212–0.216, 1, p = 0.19); CI N_hat about [0.92, 1.54] |
| Turbulence, full sample | N_round = 1 at all lags; p < 1e-4 for lags <= 20 ms; p = 0.36 / 0.37 at 50 / 100 ms | N_round = 1 at all lags; p < 0.001 at 2, 5, 20 ms; p = 0.001 at 10 ms; p = 0.35 / 0.39 at 50 / 100 ms. Same picture. |
| Turbulence, thinned (new) | — | rejection survives at 2, 5 ms (p < 0.001) and 20 ms (p = 0.002) but **not at 10 ms (p = 0.25)**; 50 / 100 ms: p = 0.58 / 0.20 |
| C. elegans, out-degree | power-law class, N_round = 4, p >= 0.998 | **class changes**: only s = 1/2 fits; N_round = 4, p = 0.99, but CI N_hat = [0.57, 28.4] |
| C. elegans, in-degree | power-law class, N_round = 26 | only s = 1/2 fits; N_round = 2, p = 0.77, CI N_hat = [0.67, 2.23]. The depth of 26 was an artefact of the NMJ pseudo-hub. |
| C. elegans, total-degree | power-law class, N_round = 7 | only s = 1 fits; N_round = 5, p = 0.995, CI N_hat = [1.85, 98.8], 26% of resamples have no admissible fit |
| C. elegans topology | 283 nodes, 2403 edges, density 0.060, mean degree 17.0, clustering 0.62, assortativity −0.15 | 279 nodes, 2287 edges, density 0.0590, mean degree 16.4, transitivity 0.213, assortativity −0.093 |
| Transitivity, all three networks | email 0.80, OpenFlights 0.75, C. elegans 0.62 | **0.267, 0.249, 0.213** — see the bug below |

C. elegans loader counts: 6417 rows; 2575 S/Sp rows kept (2194 directed pairs,
6394 contacts); 2658 R/Rp rows dropped; 153 NMJ rows dropped; 3 self-loops
dropped (all gap junctions); 514 gap-junction pairs (887 contacts), with 0
disagreements between the two listed directions. 281 neurons; with degree 0
and dropped: 28 (out), 13 (in), 2 (total).

### Things found on the way that the manuscript must take into account

1. **Transitivity was three times too large for all three networks.**
   `network_topology._clustering_and_assortativity` counts each triangle once
   per corner and then multiplied by 3 again. Fixed (in `scripts/`, not
   `htheory/`). The corrected email value, 0.267, matches the standard figure
   for that data set. Any text quoting the old values (0.80, 0.75, 0.62) is
   wrong.
2. **C. elegans conclusions are weak after the fix.** The three series no
   longer agree on a class, the intervals for N_hat span more than an order of
   magnitude for out- and total-degree, and p close to 1 with n about 260 means
   the test has no power here, not that the depth is confirmed.
3. **Turbulence: the fit is what a Gaussian increment would give.** At every
   lag a_hat is about 0.5 and N_hat about 1, i.e. the squared increment is close
   to a chi-square with one degree of freedom (a gamma with shape 1/2) — the
   law of the square of a Gaussian variable. "N = 1" here is not evidence of a
   one-level cascade specifically.
4. **OpenFlights fit misses the low-degree end** (Fig. 6b): the fitted model
   puts about half its mass at degree 0, while the analysed series has zeros
   removed (15 of 3330 airports). Moments match, the body and tail match, the
   first decade does not.
5. The size of the test at n = 5000 is 0.09 with Monte Carlo s.e. 0.03 (100
   replicates per cell); the power table cannot distinguish that from 0.05.
   More replicates would be needed to claim the nominal level holds.
6. Zero-degree policy differs by data set (kept for email, dropped for
   OpenFlights and C. elegans, as in the earlier scripts); it is recorded per
   series in the JSON (`zero_policy`, `n_zero`).

## 2. Integral time scale and thinning factors (turbulence)

Sampling interval 1 ms, 131 072 samples. Integral time scale = autocorrelation
of the squared increments integrated (trapezoid) to its first zero crossing.

| Lag (ms = samples) | Integral time scale (samples = ms) | ACF first zero (samples) | Thinning k = max(lag, ceil(T_int)) | Block length = max(lag, ceil(5 T_int)) | Block length, thinned series (steps) | Zeros dropped |
|---|---|---|---|---|---|---|
| 2 | 0.86 | 2 | 2 | 5 | 3 | 35 |
| 5 | 1.19 | 3 | 5 | 6 | 2 | 24 |
| 10 | 1.32 | 3 | 10 | 10 | 1 | 35 |
| 20 | 1.27 | 3 | 20 | 20 | 1 | 23 |
| 50 | 1.36 | 3 | 50 | 50 | 1 | 33 |
| 100 | 1.35 | 3 | 100 | 100 | 1 | 27 |

The integral time scale is about one sample at every lag, so the lag itself
always sets both k and the block length. Two cautions: (i) the ACF integrated
only to its first zero crossing does not see the correlation that overlapping
increments have at separation equal to the lag; (ii) after thinning by k = lag,
neighbouring values still share one velocity sample, so the thinned series is
less dependent, not independent. The test's own bootstrap is i.i.d. in both
cases (stated in the script docstring and in the JSON `dependence_note`).

## 3. Figures

All thirteen files were generated and inspected; none is known to be wrong.
Points that depart from the brief or that a reader should know:

- **Fonts.** "CMU Serif" and "Computer Modern Roman" are not installed here;
  `figstyle.py` adds "Latin Modern Roman" (the TeX Live build of Computer
  Modern) to the serif list, ahead of DejaVu Serif, and that is what the files
  use.
- **Fig. 6a–c model.** The envelope is simulated at the tested null — depth
  N_round, rate lambda_hat, and the shape refitted under H0 (`a_hat_H0` in the
  JSON) — not at the unconstrained a_hat, because that is the model the test
  simulates from and it matches the sample variance at integer depth. When both
  classes fit, the one with the larger p is shown (never happens in the current
  data: every series has exactly one admissible class).
- **Fig. 6c** carries three series with three envelopes in a half-width panel;
  readable, but dense. The bands run to the bottom of the axis where the lower
  envelope is zero.
- **Fig. 6d.** p = 0 out of 1000 draws is drawn at the resolution limit 0.001
  with a downward arrow; the legend sits above the axes.
- **Fig. 5c.** Chemical synapses (black) and gap junctions (orange) are
  separated by colour only; the two differ strongly in lightness, so the
  distinction survives greyscale, but there is no second encoding at
  half-pixel cell size. Gap junctions are drawn on top where both exist.
- **Fig. 3a** plots |rho(a)| on log–log axes (−rho for s = 1/2, +rho for
  s = 1), since the two classes have opposite signs and very different ranges.
- **Fig. 4b** uses a seed chosen for legibility of the illustration.
- **Fig. 2** is 2.55 in tall (brief: about 2.4) to hold the shared legend.
- **tab_empirical.tex** uses `\tabcolsep` 3 pt, not 4 pt, and two stacked
  tabulars (networks; turbulence with a narrow k column, k = 1 meaning the
  full sample, under the group headings "full sample" and "thinned"). With
  4 pt and three-digit intervals it overflowed 390 pt; checked by compiling in
  `elsarticle` `preprint,12pt` with no overfull box.

Not done / left for the owner: `make article` was not run (the manuscript is
owned by another agent and does not yet `\input` the generated tables or
reference the new figure names). The command list in `CLAUDE.md` still points
to `figures/article_numbers.json` and the old figure scripts. The old figure
files and the stale `figures/article_numbers.json` /
`figures/network_topology.json` were left untouched.
