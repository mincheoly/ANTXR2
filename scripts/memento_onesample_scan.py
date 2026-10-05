"""Genome-wide one-sample memento test for one target gene, native inference.

Uses memento's ht_2d_moments(treatment=None) one-sample test (PR #83 on
yelabucsf/scrna-parameter-estimation): per pair, the cell-count-weighted mean
of the group correlations, with memento's cell bootstrap for SE and p-value.
Groups with |r| == 1 (memento's variance-placeholder pairs) are skipped by
memento itself.

shrinkage is the centering knob: on the pooled lymphatic population, with
pegged groups excluded, the median target-gene correlation is -0.027 / +0.004
/ +0.030 at shrinkage 0 / 0.25 / 0.5 (memento default), so 0.25 is used to
put "no relationship" near r = 0 for a test against zero. A residual
expression trend remains (rho(log raw mean, r) = +0.31), so results are
reported by expression quintile.

Same group definition and filters as coexpr_point_scan.py (donor x
cell_type, q, filter_mean_thresh, min_cell_count, min_perc_group).

Usage: python memento_onesample_scan.py <pop.npz> <out.csv> [--target ANTXR2]
           [--shrinkage 0.25] [--trim 0.1] [--boot 5000] [--backend gpu]
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from antxr2_coexpr import (FILTER_MEAN_THRESH, MIN_CELL_COUNT, MIN_PERC_GROUP,  # noqa: E402
                           Q_DEFAULT, bh, load_pop)

warnings.filterwarnings("ignore")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("out")
    ap.add_argument("--target", default="ANTXR2")
    ap.add_argument("--shrinkage", type=float, default=0.25)
    ap.add_argument("--trim", type=float, default=0.1)
    ap.add_argument("--q", type=float, default=Q_DEFAULT)
    ap.add_argument("--boot", type=int, default=5000)
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    ap.add_argument("--seed", type=int, default=5)
    a = ap.parse_args()
    import memento
    t0 = time.time()
    A = load_pop(a.npz, q=a.q)
    n = A.shape[0]
    Xc = A.X.tocsc()
    stats = pd.DataFrame({"raw_mean": np.asarray(Xc.sum(0)).ravel() / n,
                          "det_rate": np.diff(Xc.indptr) / n}, index=A.var_names)
    del Xc
    A = A[:, (stats.det_rate > 0).values].copy()
    memento.setup_memento(A, q_column="q", filter_mean_thresh=FILTER_MEAN_THRESH,
                          min_cell_count=MIN_CELL_COUNT, shrinkage=a.shrinkage,
                          trim_percent=a.trim)
    memento.create_groups(A, label_columns=["donor", "cell_type"])
    memento.compute_1d_moments(A, min_perc_group=MIN_PERC_GROUP, filter_genes=True)
    genes = [g for g in A.var_names if g != a.target]
    memento.compute_2d_moments(A, [(a.target, g) for g in genes])
    memento.ht_2d_moments(A, num_boot=a.boot, backend=a.backend, random_state=a.seed,
                          verbose=0, num_cpus=1)
    r = memento.get_2d_ht_result(A).rename(columns={"gene_2": "gene"})
    r = r[["gene", "corr_coef", "corr_se", "corr_pval"]].copy()
    r["fdr"] = bh(r.corr_pval.values)
    r = r.join(stats, on="gene")
    r.insert(0, "target", a.target)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    r.sort_values("corr_pval").to_csv(a.out, index=False)
    print(f"{os.path.basename(a.npz)}: {len(r)} pairs, {len(A.uns['memento']['groups'])} groups, "
          f"backend={a.backend}, shrinkage={a.shrinkage}, boot={a.boot}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
