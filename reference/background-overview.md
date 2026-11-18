# The fluctuating background: the reservoir temperature law and its dynamics.

Three descriptions of the same object, in increasing order of
directness: `f_N` (closed form, Eq. 33, a Meijer G), `sample_TN`
(ancestral sampling down the chain of conditionals, Eqs. 31/32), and
`simulate_sde` (Euler-Maruyama on the SDE hierarchy, Eq. 29). They agree
only if the branch labels of Eq. (33) are read as corrected – see
`ERRATUM_EQ33`. `tests/testthat/test-background.R` pins that.
