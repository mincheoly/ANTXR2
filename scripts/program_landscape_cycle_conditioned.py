"""Proliferation tilt, step 2: does the co-expression tilt survive accounting
for dividing cells? (planning_summaries/program_coexpression.md,
"Proliferation tilt").

Reruns the landscape's prep + score for the epithelial and lymphocyte
populations with the main run's weights and donor bootstrap, in one of two
modes, then the lineage consensus:

    bins   groups = donor x within-donor tertile of a cell-cycle score
           (Tirosh S + G2M; program_landscape_cycling), so every correlation
           is within cells of similar cycling -- conditioning, as the
           maturation bins did for enterocytes. Group floor unchanged, so
           donors with < 3 x MIN_GROUP_CELLS cells lose bins.
    drop   cells with any MKI67 count removed, groups as in the main run.
    compare  the proliferation / tilt programs side by side: main, bins, drop
           (per-cell excess z and lineage consensus z) -> cycle_compare.csv

Outputs in LANDSCAPE_OUT/cycle_<mode>/ (same layout as program_landscape.py,
plus report/consensus_<lineage>.csv) and LANDSCAPE_OUT/cycle_compare.csv.

Usage:
    python program_landscape_cycle_conditioned.py {bins,drop,compare}
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import program_config                              # noqa: E402
from program_config import LANDSCAPE_LINEAGES, LANDSCAPE_OUT  # noqa: E402
from program_landscape_cycling import G2M_GENES, S_GENES  # noqa: E402

MAIN = LANDSCAPE_OUT
POPS = LANDSCAPE_LINEAGES["epithelial"] + LANDSCAPE_LINEAGES["lymphocyte"]
N_BINS = 3
PROGRAMS = ["hallmark:E2F_TARGETS", "hallmark:G2M_CHECKPOINT", "hallmark:MYC_TARGETS_V1",
            "hallmark:MYC_TARGETS_V2", "hallmark:DNA_REPAIR", "hallmark:MTORC1_SIGNALING",
            "collectri:E2F1", "collectri:E2F2", "collectri:E2F4", "collectri:E2F6",
            "cytosig:TGFB3", "progeny:Hypoxia", "collectri:ZNF331"]


def cycle_bins(X, obs, genes):
    """Within-donor tertile of the S + G2M score (log1p CP10k, per-gene z over
    the population, genes detected in >= 1% of cells)."""
    idx = [genes.get_loc(g) for g in S_GENES + G2M_GENES if g in genes]
    L = (sp.diags(1e4 / obs.total_counts.values.astype(float)) @ X[:, idx]).toarray()
    L = np.log1p(L[:, (L > 0).mean(0) >= 0.01])
    s = pd.Series(((L - L.mean(0)) / (L.std(0) + 1e-9)).mean(1))
    by = s.groupby(obs.donor.values)
    b = (by.rank(method="first") - 1) * N_BINS // by.transform("size")
    return "cy" + (b.astype(int) + 1).astype(str)


def run(mode):
    sub = MAIN / f"cycle_{mode}"
    sub.mkdir(parents=True, exist_ok=True)
    for f in ("weights.csv", "programs.csv", "boot_counts.npz"):
        shutil.copy(MAIN / f, sub / f)
    program_config.LANDSCAPE_OUT = sub             # before the landscape modules bind it
    import program_landscape as PL
    _load = PL.load_cells

    def load(path):
        X, obs, genes = _load(path)
        if mode == "bins":
            obs["fine"] = obs.fine + "|" + cycle_bins(X, obs, genes).values
            return X, obs, genes
        mki = genes.get_loc("MKI67")
        keep = np.where(X[:, mki].toarray().ravel() == 0)[0]
        return X[keep], obs.iloc[keep].reset_index(drop=True), genes

    PL.load_cells = load
    from program_scores import load_library
    PL.prep(POPS, load_library())
    PL.score(POPS, False)
    import program_landscape_report as R
    R.consensus(pd.read_csv(sub / "cells.csv"))


def compare():
    rows, cons = [], []
    for name, d in [("main", MAIN), ("bins", MAIN / "cycle_bins"), ("drop", MAIN / "cycle_drop")]:
        L = pd.read_csv(d / "landscape.csv")
        L = L[(L.scheme == "pc1") & L.population.isin(POPS) & L.unit.isin(PROGRAMS)]
        rows.append(L[["population", "state", "unit", "excess", "excess_z"]].assign(run=name))
        for lin in ("epithelial", "lymphocyte"):
            c = pd.read_csv(d / "report" / f"consensus_{lin}.csv", index_col=0)
            c = c.loc[c.index.intersection(PROGRAMS), ["consensus_z", "rank_among_genes"]]
            cons.append(c.assign(lineage=lin, run=name))
    per = pd.concat(rows)
    cells = {name: pd.read_csv(d / "cells.csv") for name, d in
             [("main", MAIN), ("bins", MAIN / "cycle_bins"), ("drop", MAIN / "cycle_drop")]}
    lineage_of = {p: lin for lin, ps in LANDSCAPE_LINEAGES.items() for p in ps}
    per["lineage"] = per.population.map(lineage_of)
    # per-cell z only over cells scored in all three runs
    w = per.pivot_table(index=["lineage", "population", "state", "unit"], columns="run",
                        values="excess_z").dropna()
    summ = w.groupby(["lineage", "unit"]).agg(
        n_cells=("main", "size"), **{f"{r}_mean_z": (r, "mean") for r in ("main", "bins", "drop")},
        **{f"{r}_neg": (r, lambda v: (v < 0).mean()) for r in ("main", "bins", "drop")})
    con = pd.concat(cons).reset_index(names="unit").pivot_table(
        index=["lineage", "unit"], columns="run", values=["consensus_z", "rank_among_genes"])
    con.columns = [f"{r}_{v.replace('_among_genes', '')}" for v, r in con.columns]
    out = summ.join(con)
    out.to_csv(MAIN / "cycle_compare.csv")
    w.to_csv(MAIN / "cycle_compare_cells.csv")
    pd.set_option("display.width", 250)
    print(out.round(2).to_string())
    for name, c in cells.items():
        f = c[c.population.isin(POPS) & c.filled]
        print(f"{name}: {len(f)} filled cells, median donors {f.n_donors.median():.0f}, "
              f"median cells in groups {f.cells_in_groups.median():.0f}")


if __name__ == "__main__":
    mode = sys.argv[1]
    assert mode in ("bins", "drop", "compare"), mode
    compare() if mode == "compare" else run(mode)
