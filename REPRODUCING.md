# Reproducing the numbers, tables and figures

Everything the manuscript quotes is regenerated from the raw files in `data/`
by one command, run from the repository root:

```bash
make all          # results -> tables -> figures -> article
```

## Environment

Python 3.10.12 with the exact package versions in `requirements-lock.txt`
(mpmath 1.3.0, numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9, pytest 9.1.1):

```bash
pip3 install --user -r requirements-lock.txt
```

No other Python package is needed (in particular, not networkx). `make article`
needs a TeX Live installation with `latexmk` and `elsarticle`. The figures use
the font "Latin Modern Roman" when it is installed (it ships with TeX Live) so
that figure text matches the LaTeX body; otherwise matplotlib falls back to
DejaVu Serif and only the look of the text changes.

Another interpreter can be selected with `make PYTHON=.venv/bin/python3 all`.
Every target sets `PYTHONPATH=.`.

## Targets

| Target | What it does | Output | Wall time observed |
|---|---|---|---|
| `make results` | runs the six analysis scripts | `results/*.json` | about 11 min (655 s; `compute_article_numbers.py` about 8 min, `analyze_turbulence.py` about 2.5 min, the rest seconds) |
| `make tables` | `scripts/make_tables.py` | `article/generated/tab_{topology,power,empirical}.tex` | < 1 s |
| `make figures` | `scripts/make_theory_figures.py`, `scripts/make_data_figures.py` | `figures/fig1_…` to `fig6d_…` (`.pdf` + 300 dpi `.png`) | about 1 min (67 s; the entropy curves of Fig. 2 are mpmath-bound) |
| `make article` | `cd article && latexmk -pdf main.tex` | `article/main.pdf` | — |
| `make test-fast` | `pytest tests/ -q -m "not slow"` | — | 2 min (393 passed, 2 skipped, 38 deselected) |
| `make all` | `results tables figures article` | | |

Times are for one core of the machine this was last run on; nothing is
parallelised.

## Pipeline

```
data/*  --scripts/data_loaders.py-->  scripts/analyze_*.py, network_topology.py,
                                      compute_article_numbers.py
                                          |
                                          v
                                   results/*.json
                                    |           |
                 scripts/make_tables.py     scripts/make_*_figures.py
                                    |           |
                   article/generated/*.tex   figures/fig*.pdf
```

* `scripts/data_loaders.py` is the only place where raw files are parsed.
* `scripts/common.py` holds the shared fit / depth-test / bootstrap-interval
  code and writes the JSON files with provenance metadata (script, seed,
  package versions, UTC timestamp).
* `scripts/figstyle.py` is the single figure style (sizes, fonts, palette).
* Tables contain no hand-typed number: `make_tables.py` reads only
  `results/*.json`.

## Random seeds

| Script | Seed | Streams |
|---|---|---|
| `analyze_email_eu_core.py` | 20260905 | depth test: `default_rng(seed)`; bootstrap intervals: `default_rng([seed, 1])` |
| `analyze_openflights.py` | 20260905 | same |
| `analyze_celegans.py` | 20260826 | same |
| `analyze_turbulence.py` | 20260826 | same |
| `compute_article_numbers.py` | 20260826 | one stream |
| `network_topology.py` | — | deterministic |
| `make_theory_figures.py` | 0 (Fig. 1a), 12 (Fig. 4b), 20260826 (Fig. 4c) | |
| `make_data_figures.py` | 0 (route sample, Fig. 5a), 20260905 (envelopes, Fig. 6a–c) | |

Bootstrap sizes: 1000 parametric-bootstrap draws for every empirical depth
test, 1000 resamples for every confidence interval, 500 simulated sequences
for each envelope in Fig. 6. The size/power study keeps its own design (100
replicates x 200 bootstrap draws per cell; Monte Carlo standard errors are
stored next to each rate).

With the same package versions the JSON files are reproduced exactly except
for the `timestamp_utc` field.

## Data sources

| File(s) in `data/` | Source | Licence / terms | Retrieved |
|---|---|---|---|
| `email-Eu-core.txt(.gz)`, `email-Eu-core-department-labels.txt(.gz)` | SNAP, Stanford Large Network Dataset Collection, `https://snap.stanford.edu/data/email-Eu-core.html` (1005 nodes, 25 571 directed edges, 42 departments) | The script docstring records no licence; SNAP distributes the collection for research use and asks that the data set and collection be cited (Leskovec & Krevl 2014; Yin et al. 2017; Leskovec, Kleinberg & Faloutsos 2007). Confirm the terms before redistribution. | local files dated 2026-08-26 |
| `openflights-airports.dat`, `openflights-routes.dat` | OpenFlights Airport and Route Databases, `https://github.com/jpatokal/openflights` | Open Database License (ODbL): free to use with attribution, no registration | 2026-09-05 (terms checked that day) |
| `celegans-NeuronConnect.csv` | WormAtlas `NeuronConnect.xls` (`wormatlas.org/images/NeuronConnect.xls`), adult hermaphrodite connectome, converted to CSV with LibreOffice | described in the script docstring as a public-domain research resource, no registration; cite WormAtlas and the underlying connectome papers | local file dated 2026-09-05 |
| `turbulence-mendeley-cylinder50mm-16ms.txt` | "Turbulent Flow Past Cylinders", Mendeley Data, `https://data.mendeley.com/datasets/xy8ggx9hn6/1`; hot-wire velocity in the wake of a 50 mm cylinder, Re = 5.17e4, 1000 Hz, 2^17 samples | CC BY 4.0 | local file dated 2026-09-05 |

The continent of each airport is not in the OpenFlights data; it comes from the
hand-built `COUNTRY_TO_CONTINENT` table in `scripts/analyze_openflights.py`.

## Data-processing decisions that affect results

* **C. elegans.** Chemical network from the `S`/`Sp` rows only (the `R`/`Rp`
  rows are mirror images of the same synapses); gap junctions (`EJ`) as
  undirected pairs counted once; the 153 `NMJ` rows, which all point to one
  pseudo-node, and self-loops are dropped. Weights are the contact counts
  `Nbr`. Series: weighted chemical out-degree, in-degree, and total degree =
  out + in + gap-junction degree; neurons with degree 0 in a series are
  dropped from it (counts in `results/analyze_celegans.json`).
* **OpenFlights.** Degree = number of rows of the route table (a pair served
  by several airlines counts several times); airports with degree 0 in a
  series are dropped.
* **email-Eu-core.** All 1005 vertices are kept, including zeros.
* **Turbulence.** Squared velocity increments at lags 2–100 ms; exact zeros
  are dropped. The depth test's bootstrap is i.i.d. while the data are
  serially dependent; thinned-series tests and moving-block-bootstrap
  intervals are reported alongside (see the script docstring).
* **Depth test.** It tests `H0: N = N_round` with `N_round` chosen from the
  same data, so its p-value is a goodness-of-fit diagnostic, not a
  confirmatory test; percentile intervals for the continuous `N_hat` are
  reported with it.

## Legacy scripts

`scripts/make_article_figures.py` and `scripts/make_network_figures.py`
produced the earlier figure files (`fig_cascade_and_B`, `fig_identifiability`,
`fig_network_schematic`, `fig_cross_dataset_summary`, `fig_real_networks`).
They are superseded and refuse to run without `--legacy`; their old outputs
are left in `figures/` untouched. `figures/article_numbers.json` and
`figures/network_topology.json` are stale copies of what now lives in
`results/`.
