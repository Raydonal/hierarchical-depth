# hierarchical-depth

Code, data and results for the identifiability of the hierarchical depth in
superstatistical cascades (H-theory), with applications to networks and
turbulence. The accompanying manuscript is in preparation.

The starting point is the construction of **"Hierarchical entropies"** (de Moraes, Brum, Macêdo
and Vasconcelos, *Physica A* **699** (2026) 131899,
[doi:10.1016/j.physa.2026.131899](https://doi.org/10.1016/j.physa.2026.131899)): a system coupled
to `N` nested reservoirs whose temperatures fluctuate, with gamma (`s = 1/2`) or inverse-gamma
(`s = 1`) conditionals. Averaging the Boltzmann factor over the innermost temperature gives a
statistical factor `B_N^(s)(E; a)` in closed form as a Meijer G-function, and the maximum-entropy
principle under the standard constraints returns the entropy that produces it.

The cascade is additive in `log T`, so the log-cumulants of the innermost temperature are sums of
polygamma functions of the shape parameters. With equal shapes, the sign of `kappa_3 / kappa_2`
identifies the class, its value the shape, and `kappa_2` the depth `N`. Through the
identifiability of mixed-Poisson laws, the same holds for the degree sequence of a network.

## Contents

| Path | What it is |
|---|---|
| `htheory/` | Python package: Meijer G (`meijer`), reservoir law and samplers (`background`), statistical factor (`statfactor`), entropy (`entropy`), maximum entropy (`maxent`), log-cumulant estimator and depth test (`identifiability`), mixed-Poisson networks (`network`) |
| `htheoryR/` | Native R port of the same package ([documentation site](https://raydonal.github.io/hierarchical-depth/)) |
| `tests/` | `pytest` suite; tests marked `slow` cross-validate closed forms against quadrature and Monte Carlo |
| `scripts/` | Analyses of the four data sets, tables and figures |
| `data/` | Raw data (sources and licences in `REPRODUCING.md`) |
| `results/` | JSON outputs of every analysis, with seeds and package versions |
| `figures/` | Figures generated from `results/` |
| `docs/` | Proofs of the theorems (`theorems.md`) and the Meijer G reduction used by the R port (`meijer_derivation.md`) |

## Installation

Python (tested with Python 3.10, exact versions in `requirements-lock.txt`):

```bash
pip install -r requirements.txt
```

R package:

```r
# install.packages("remotes")
remotes::install_github("Raydonal/hierarchical-depth", subdir = "htheoryR")
```

`Rmpfr` needs the GMP and MPFR development headers (Debian/Ubuntu:
`sudo apt install libgmp-dev libmpfr-dev`).

## Reproducing the results

From the repository root:

```bash
make results     # analyses -> results/*.json (about 13 min, fixed seeds)
make tables      # LaTeX tables from results/*.json
make figures     # figures/*.pdf and .png previews
make test-fast   # pytest without the slow cross-validation
python3 -m pytest tests/ -q   # full suite, a few minutes
```

`REPRODUCING.md` lists the data sources, licences, seeds and run times.

## Data

| Data set | Source | Licence |
|---|---|---|
| `email-Eu-core` | [SNAP](https://snap.stanford.edu/data/email-Eu-core.html) | Research use, citation requested by SNAP |
| OpenFlights airports and routes | [OpenFlights](https://github.com/jpatokal/openflights) | ODbL |
| *C. elegans* connectome | [WormAtlas](https://www.wormatlas.org/) `NeuronConnect.xls` | Public research resource; cite WormAtlas |
| Wind-tunnel turbulence | [Mendeley Data](https://data.mendeley.com/datasets/xy8ggx9hn6/1) | CC BY 4.0 |

The data files keep their original licences; the MIT licence of this repository covers the code.

## Licence

Code: MIT (see `LICENSE`). Author: Raydonal Ospina
([ORCID 0000-0002-9884-9090](https://orcid.org/0000-0002-9884-9090)), Universidade Federal da
Bahia and Universidade Federal de Pernambuco.
