"""Cascade-depth identifiability test on the C. elegans connectome.

Data: `NeuronConnect.xls` from WormAtlas (wormatlas.org/images/NeuronConnect.xls,
the Cook/White-lineage adult hermaphrodite connectome, public domain research
resource, no registration required), converted to CSV with LibreOffice
(`data/celegans-NeuronConnect.csv`). Columns: Neuron 1, Neuron 2, Type, Nbr.
`Type` is one of EJ (gap junction, electrical), S/Sp (chemical synapse,
sending), R/Rp (chemical synapse, receiving), NMJ (neuromuscular junction);
`Nbr` is the number of anatomical contacts for that (pre, post, type) triple.

Processing (all of it in `scripts/data_loaders.load_celegans`, shared with the
topology and figure scripts):

* chemical network: directed, from the S and Sp rows only (Neuron 1 =
  presynaptic, Neuron 2 = postsynaptic), weight = Nbr. The R/Rp rows are the
  mirror images of the same synapses and are dropped -- summing over all rows,
  as an earlier version of this script did, mixes sends and receives.
* gap junctions: undirected, each unordered pair once, weight = Nbr.
* the 153 NMJ rows all point to one pseudo-node called "NMJ", which is not a
  neuron; they are dropped. Self-loops are dropped.

Degree series fed to the estimator: weighted chemical out-degree, weighted
chemical in-degree, and weighted total degree = out + in + gap-junction
degree. Only neurons with positive degree in the series at hand are kept; the
number dropped is recorded.

No functional hierarchy label (sensory/inter/motor) is attached, so this is a
*flat* data set: a real degree sequence with no multi-level ground truth to
check the estimated depth against.

Run from the repo root:

    PYTHONPATH=. python3 scripts/analyze_celegans.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts import common, data_loaders  # noqa: E402

SEED = 20260826
SYSTEM = "C. elegans"


def main():
    data = data_loaders.load_celegans(verbose=True)
    rng_test = np.random.default_rng(SEED)
    rng_ci = np.random.default_rng([SEED, 1])
    series = []
    for name, D in (("out-degree", data.out_deg), ("in-degree", data.in_deg),
                    ("total-degree", data.total_deg)):
        n_zero = int(np.sum(D == 0))
        print(f"{name}: {D.size} neurons, {n_zero} with degree 0 dropped")
        series.append(common.analyze_degree_series(
            SYSTEM, name, D[D > 0], rng_test, rng_ci,
            n_zero=n_zero, zero_policy="dropped", weighted=True))
    common.write_json("analyze_celegans", {
        "metadata": common.metadata(
            "scripts/analyze_celegans.py", SEED,
            rng_streams={"depth_test": "default_rng(SEED)",
                         "bootstrap_ci": "default_rng([SEED, 1])"},
            n_boot_test=common.N_BOOT_TEST, n_boot_ci=common.N_BOOT_CI,
            note=common.SELECTION_NOTE),
        "loader": data.info,
        "series": series,
    })


if __name__ == "__main__":
    main()
