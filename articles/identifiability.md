# Recovering the depth of a cascade

Every number and plot on this page is computed by `htheoryR` when the
site is built. The page runs in a few seconds.

## The cascade

A cascade of depth \(N\) starts at \(T_{0}\) and draws each temperature
from a gamma law (class \(s = 1/2\)) or an inverse-gamma law (class
\(s = 1\)) whose mean is the temperature of the level above, with shape
\(a_{i}\). Every level keeps the mean,
\(E\lbrack T_{N}\rbrack = T_{0}\), but adds dispersion. `sample_TN()`
draws \(T_{N}\) directly.

``` r
set.seed(1)
logT <- lapply(1:3, function(N) log(sample_TN(rep(2, N), s = 0.5, size = 1e5)))
dens <- lapply(logT, density, n = 512)
par_default()
plot(NULL, xlim = c(-5, 2), ylim = c(0, max(sapply(dens, function(d) max(d$y)))),
     xlab = expression(log~T[N]), ylab = "density")
for (N in 1:3) lines(dens[[N]], col = pal[N], lwd = 2, lty = N)
legend("topleft", legend = paste("N =", 1:3), col = pal[1:3], lty = 1:3,
       lwd = 2, bty = "n")
```

![Density of log T\_N for depths N = 1, 2, 3 (class s = 1/2, a = 2, 10^5
draws each). Deeper cascades are more dispersed, with a longer left
tail.](identifiability_files/figure-html/cascade-1.png)

Density of log T\_N for depths N = 1, 2, 3 (class s = 1/2, a = 2, 10^5
draws each). Deeper cascades are more dispersed, with a longer left
tail.

## Why the log-cumulants identify the depth

The cascade is multiplicative, so \(\log T_{N}\) is a sum of \(N\)
independent terms. Its cumulants of order two and higher are sums over
the levels. For equal shapes, every level contributes the same amount:

\[\kappa_{r}(\log T_{N}) = N\,\kappa_{r}^{(1)}(a),\qquad r \geq 2.\]

The ratio \(\rho_{s}(a) = \kappa_{3}/\kappa_{2}\) therefore does not
depend on \(N\). It is negative for \(s = 1/2\) and positive for
\(s = 1\), so its sign gives the class and its value gives the shape
\(a\). Once \(a\) is known, \(\kappa_{2}\) gives \(N\).
`log_cumulants_theory()` returns the closed forms, which are polygamma
functions of \(a\).

``` r
a <- exp(seq(log(0.05), log(50), length.out = 200))
rho <- function(a, s) { k <- log_cumulants_theory(a, s, r_max = 3); k[2] / k[1] }
r_half <- sapply(a, rho, s = 0.5)
r_one  <- sapply(a, rho, s = 1)
c(range_s_half = range(r_half), range_s_one = range(r_one))
#> range_s_half1 range_s_half2  range_s_one1  range_s_one2 
#>  -39.85259934   -0.02020065    0.01980069    1.37576136
par_default()
plot(a, abs(r_half), log = "xy", type = "l", lwd = 2, col = pal[1],
     xlab = "shape a", ylab = expression("|" * kappa[3] / kappa[2] * "|"),
     ylim = range(abs(c(r_half, r_one))))
lines(a, r_one, lwd = 2, col = pal[2], lty = 2)
abline(h = 2 * 1.2020569 / (pi^2 / 6), lty = 3, col = "grey40")
legend("bottomleft", lwd = 2, col = pal[1:2], lty = 1:2, bty = "n",
       legend = c("s = 1/2 (ratio < 0)", "s = 1 (ratio > 0)"))
```

![Absolute log-cumulant ratio |kappa\_3 / kappa\_2| as a function of the
shape a. The two classes have opposite signs, and for s = 1 the ratio
stays below 2 zeta(3) / zeta(2), about 1.46 (dotted
line).](identifiability_files/figure-html/ratio-1.png)

Absolute log-cumulant ratio |kappa\_3 / kappa\_2| as a function of the
shape a. The two classes have opposite signs, and for s = 1 the ratio
stays below 2 zeta(3) / zeta(2), about 1.46 (dotted line).

## Estimating the depth from a sample

`empirical_log_cumulants()` computes \({\widehat{\kappa}}_{2}\) and
\({\widehat{\kappa}}_{3}\) from a sample of \(T_{N}\), and
`fit_homogeneous()` inverts the closed forms. Here the data come from a
cascade with \(N = 3\) and \(a = 2\); the intermediate levels are never
observed.

``` r
set.seed(2)
T <- sample_TN(c(2, 2, 2), s = 0.5, size = 1e5)
k_hat <- empirical_log_cumulants(T, r_max = 3)
k_true <- log_cumulants_theory(c(2, 2, 2), s = 0.5, r_max = 3)
data.frame(cumulant = c("kappa_2", "kappa_3"), sample = k_hat, theory = k_true)
#>   cumulant    sample    theory
#> 1  kappa_2  1.946289  1.934802
#> 2  kappa_3 -1.215907 -1.212341
str(fit_homogeneous(k_hat[1], k_hat[2], s = 0.5))
#> List of 3
#>  $ a      : num 2.01
#>  $ N      : num 3.03
#>  $ N_round: num 3
```

The estimator is consistent, but its spread depends on the sample size.
The next plot repeats the fit on 300 independent samples at each size.

``` r
set.seed(3)
N_hat <- sapply(c(1000, 10000), function(n) replicate(300, {
  k <- empirical_log_cumulants(sample_TN(c(2, 2, 2), s = 0.5, size = n), r_max = 3)
  fit <- tryCatch(fit_homogeneous(k[1], k[2], s = 0.5), error = function(e) NULL)
  if (is.null(fit)) NA else fit$N
}))
colnames(N_hat) <- c("n = 1000", "n = 10000")
apply(N_hat, 2, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
#>     n = 1000 n = 10000
#> 5%  2.155289  2.695967
#> 50% 3.074596  3.006167
#> 95% 4.515377  3.344705
par_default()
boxplot(N_hat, col = paste0(pal[1:2], "55"), border = pal[1:2],
        ylab = expression(hat(N)), outline = FALSE)
abline(h = 3, lty = 2)
```

![Sampling distribution of the depth estimate for a true depth N = 3 (s
= 1/2, a = 2), over 300 samples of each size. The dashed line marks the
true depth.](identifiability_files/figure-html/sampling-1.png)

Sampling distribution of the depth estimate for a true depth N = 3 (s =
1/2, a = 2), over 300 samples of each size. The dashed line marks the
true depth.

## Testing a fixed depth

`test_depth()` is a parametric-bootstrap test of \(H_{0}:N = N_{0}\).
With \(n = 1000\) draws from a two-level cascade, the test of
\(N_{0} = 1\) rejects; with draws from a one-level cascade, it does not.
The manuscript measures the size and power of this test in a simulation
study (see [Results of the
manuscript](https://raydonal.github.io/hierarchical-depth/articles/results.html#simulation-size-and-power-of-the-depth-test)).

``` r
set.seed(4)
two_levels <- sample_TN(c(2, 2), s = 0.5, size = 1000)
one_level  <- sample_TN(2, s = 0.5, size = 1000)
rbind(true_N_2 = unlist(test_depth(two_levels, N0 = 1, s = 0.5, n_boot = 200)),
      true_N_1 = unlist(test_depth(one_level,  N0 = 1, s = 0.5, n_boot = 200)))
#>            statistic p_value    a_hat n_boot
#> true_N_2 0.648695491   0.005 1.185367    200
#> true_N_1 0.004459086   0.935 2.060048    200
```

## The same argument for networks

In the hierarchical mixed-Poisson model, each vertex has a propensity
\(W_{v}\) that is the product of factors attached to its nested
communities, and its degree is Poisson with mean \(\lambda W_{v}\).
Mixed-Poisson laws are identifiable from the degrees, and the factorial
moments
\(E\lbrack D^{(k)}\rbrack = \lambda^{k}\, E\lbrack W^{k}\rbrack\) carry
the log-cumulants of \(W\). The depth can therefore be read from the
degree sequence alone.

``` r
set.seed(5)
D <- lapply(0:3, function(N) sample_degrees(rep(2, N), s = 0.5, lam = 20, size = 1e5))
par_default()
plot(NULL, xlim = c(1, max(unlist(D))), ylim = c(1e-5, 1), log = "xy",
     xlab = "degree k", ylab = expression(P(D >= k)))
for (i in 1:4) {
  k <- sort(unique(D[[i]]))
  ccdf <- 1 - ecdf(D[[i]])(k - 1)
  lines(pmax(k, 1), ccdf, col = pal[i], lwd = 2, lty = i)
}
legend("bottomleft", bty = "n", col = pal[1:4], lty = 1:4, lwd = 2,
       legend = sprintf("N = %d, VMR = %.1f", 0:3,
                        sapply(D, function(d) var(d) / mean(d))))
```

![Complementary cumulative degree distribution P(D \>= k) of 10^5
simulated vertices for depths N = 0 to 3 (s = 1/2, a = 2, lambda = 20).
At equal lambda and a, depth spreads the
degrees.](identifiability_files/figure-html/degrees-1.png)

Complementary cumulative degree distribution P(D \>= k) of 10^5
simulated vertices for depths N = 0 to 3 (s = 1/2, a = 2, lambda = 20).
At equal lambda and a, depth spreads the degrees.

``` r
set.seed(6)
D3 <- sample_degrees(c(3, 3, 3), s = 0.5, lam = 15, size = 1e5)
net_fit <- fit_homogeneous_network(D3, s = 0.5)
str(net_fit)
#> List of 4
#>  $ a      : num 2.59
#>  $ N      : num 2.65
#>  $ N_round: num 3
#>  $ lam    : num 15.1
```

The rounded depth is 3, against a true depth \(N = 3\), and no community
membership was observed. The point estimates (\(\widehat{N} = 2.65\),
\(\widehat{a} = 2.59\) against \(a = 3\)) are less precise than those
obtained from \(T_{N}\) itself above, because the degrees observe \(W\)
only through Poisson noise. The page *Results of the manuscript* applies
the same estimator to real networks and to turbulence data.
