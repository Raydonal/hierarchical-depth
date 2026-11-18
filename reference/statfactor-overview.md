# The generalized statistical factor `B_N^{(s)}(E; a)`, Eq. (34).

`B_N(E)` is the Boltzmann factor averaged over the innermost reservoir
temperature, Eq. (24): `B_N(E) = int_0^inf exp(-E/T_N) f_N(T_N) dT_N`.

## Details

Two invariants anchor the package (checked in
`tests/testthat/test-statfactor.R`): the flat-energy sum rule, Eq. (27),
and the boundary conditions `B(0) = 1`, `B(inf) = 0`, `B'(E) < 0`, Eq.
(7)/(43).
