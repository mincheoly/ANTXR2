"""Fibroblast interferon lean: subtype-composition check
(planning_summaries/program_coexpression.md, "Fibroblast interferon lean").

The landscape's pooled fibroblasts are grouped donor x ontology class, and the
"fibroblast" class pools three author subtypes (LP fibroblast ADAMDEC1,
Fibroblast CPM, Fibroblast C3) that differ in ANTXR2 and interferon level.
This reruns the landscape's prep + score for that population with groups =
donor x author_cell_type, reusing the main run's weights and bootstrap, and
tabulates the per-subtype mean levels.

    LANDSCAPE_OUT/fib_author_subtype/
        subtype_means.csv   per state x author subtype: donors, cells,
                            donor-mean log-normalised ANTXR2, its detection,
                            and hallmark IFN-alpha / IFN-gamma scores (mean of
                            per-gene z over the hallmark genes)
        cells.csv, landscape.csv, genes/, boot/, prep/, score_meta.json
                            as program_landscape.py, author-subtype groups

Needs program_landscape.py all (weights.csv, programs.csv, boot_counts.npz)
and the IBDverse h5ad (author_cell_type is read from its obs).

Usage:
    python program_landscape_fib_subtype.py
"""
import os
import shutil
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fib_cell_extract as F                       # noqa: E402
import program_config                              # noqa: E402
from program_config import (                       # noqa: E402
    LANDSCAPE_CELLS, LANDSCAPE_FIB, LANDSCAPE_FIB_TYPES, LANDSCAPE_OUT, LANDSCAPE_STATES,
    SOURCES, pop_slug,
)

MAIN = LANDSCAPE_OUT
SUB = MAIN / "fib_author_subtype"
SUB.mkdir(parents=True, exist_ok=True)
for f in ("weights.csv", "programs.csv", "boot_counts.npz"):
    shutil.copy(MAIN / f, SUB / f)
program_config.LANDSCAPE_OUT = SUB                 # before program_landscape binds it
import program_landscape as PL                     # noqa: E402
from program_scores import load_cells, load_library, size_factors  # noqa: E402

IFN = ["hallmark:INTERFERON_ALPHA_RESPONSE", "hallmark:INTERFERON_GAMMA_RESPONSE"]


def author_labels():
    """author_cell_type of the pooled fibroblasts, in the npz's row order
    (program_extract's selection, h5ad order)."""
    with h5py.File(SOURCES["ibdverse"]["h5ad"], "r") as f:
        o = pd.DataFrame({c: F.rc(f["obs"], c)
                          for c in ["tissue", "cell_type", "donor_id", "author_cell_type"]})
    o = o[(o.tissue == "ileum") & o.cell_type.isin(LANDSCAPE_FIB_TYPES)]
    return o.author_cell_type.astype(str).values, o.donor_id.astype(str).values


def subtype_means(X, obs, genes, lib):
    sf = size_factors(X, genes)
    L = (sp.diags(np.median(sf) / sf) @ X).tocsr()
    L.data = np.log1p(L.data)
    t = L[:, genes.get_loc("ANTXR2")].toarray().ravel()
    d = obs.assign(state=obs.disease.map({v: k for k, v in LANDSCAPE_STATES.items()}),
                   ANTXR2=t, ANTXR2_det=t > 0)
    cols = ["ANTXR2", "ANTXR2_det"]
    for p in IFN:
        M = L[:, [genes.get_loc(g) for g in lib[lib.program_id == p].gene if g in genes]].toarray()
        d[p] = ((M - M.mean(0)) / (M.std(0) + 1e-9)).mean(1)
        cols.append(p)
    by_donor = d.groupby(["state", "fine", "donor"])
    dm = by_donor[cols].mean().join(by_donor.size().rename("cells"))
    out = dm.groupby(level=[0, 1]).agg(n_donors=("cells", "size"), cells=("cells", "sum"),
                                       **{c: (c, "mean") for c in cols})
    out.index.names = ["state", "author_cell_type"]
    return out


def main():
    author, donor = author_labels()
    _load = PL.load_cells

    def load_with_author(path):
        X, obs, genes = _load(path)
        assert (obs.donor.values == donor).all(), "npz rows not in h5ad selection order"
        obs["fine"] = author
        return X, obs, genes

    PL.load_cells = load_with_author
    lib = load_library()
    X, obs, genes = load_with_author(LANDSCAPE_CELLS / f"cells_{pop_slug(LANDSCAPE_FIB)}.npz")
    m = subtype_means(X, obs, genes, lib)
    m.to_csv(SUB / "subtype_means.csv")
    print(m.round(3).to_string())
    PL.prep([LANDSCAPE_FIB], lib)
    PL.score([LANDSCAPE_FIB], False)


if __name__ == "__main__":
    main()
