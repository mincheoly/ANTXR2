"""Point-estimate genome-wide co-expression scan for one target gene per population.

Bootstrap-free counterpart of antxr2_coexpr.py, for analyses that rank genes by
correlation and never use per-gene p-values (pathway GSEA, anchor controls,
size-factor sensitivity). The statistic is the cell-count-weighted mean of the
per-group (donor x cell_type) memento point estimates. On the existing ANTXR2
scans this reproduces antxr2_coexpr.py's bootstrap corr_coef at Spearman
0.98-0.999 with no offset (0.88 for blood, where 18% of group estimates are
pegged at |r|=1).

Why not the bootstrap: memento's current release drops all-constant treatment
columns (_prefilter_constant_treatment), so the all-ones one-sample branch that
antxr2_coexpr.py was built on no longer exists on either backend.

Group estimates with |r| == 1 are dropped before averaging, as in memento's
ht_2d_moments (see the comment in scan()). Versions of this script before
2026-10-04 did not do this.

Same filters as antxr2_coexpr.py (q, filter_mean_thresh, min_cell_count,
min_perc_group); shrinkage / trim_percent are arguments so both the memento
defaults (0.5 / 0.1, used for the existing scans) and the project's earlier
corrected settings (0 / 0.5) can be run.

Usage:
    python coexpr_point_scan.py <pop.npz> <out.csv> [--targets ANTXR2 G1 ...]
        [--shrinkage 0.5] [--trim 0.1] [--q 0.15]
Writes one long CSV: target, gene, corr, n_groups_valid, raw_mean, det_rate.
With --groups-out PREFIX also writes PREFIX_corr.parquet (per-group
correlations), PREFIX_rawmean.parquet (genes x groups raw means) and
PREFIX_ncells.csv.
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from antxr2_coexpr import (FILTER_MEAN_THRESH, MIN_CELL_COUNT, MIN_PERC_GROUP,
                           Q_DEFAULT, load_pop)

warnings.filterwarnings("ignore")

PAIR_CHUNK = 2000


def scan(npz, targets, shrinkage, trim, q, groups_out=None):
    import memento
    A = load_pop(npz, q=q)
    Xc = A.X.tocsc()
    n = A.shape[0]
    stats = pd.DataFrame({"raw_mean": np.asarray(Xc.sum(0)).ravel() / n,
                          "det_rate": np.diff(Xc.indptr) / n}, index=A.var_names)
    del Xc
    # Genes with no counts in the population can never clear memento's mean
    # filter; dropping them first changes no result but roughly halves the
    # per-gene arrays memento allocates (Tabula Sapiens extracts carry 60,606
    # genes, and two concurrent epithelial runs exhausted a 13 GB host).
    A = A[:, (stats.det_rate > 0).values].copy()
    memento.setup_memento(A, q_column="q", filter_mean_thresh=FILTER_MEAN_THRESH,
                          min_cell_count=MIN_CELL_COUNT, shrinkage=shrinkage,
                          trim_percent=trim)
    memento.create_groups(A, label_columns=["donor", "cell_type"])
    if not list(A.uns["memento"]["groups"]):
        return pd.DataFrame()
    memento.compute_1d_moments(A, min_perc_group=MIN_PERC_GROUP, filter_genes=True)
    genes = list(A.var_names)
    targets = [t for t in targets if t in set(genes)]
    pairs = [(t, g) for t in targets for g in genes if g != t]
    # memento's covariance estimator slices data[:, idx1] with the target's
    # column index repeated once per pair, i.e. it materializes the target
    # column ~len(pairs) times. All pairs at once peaked at 12.4 GB on a
    # 25k-cell population; chunking bounds it without changing any estimate.
    chunks = []
    for i in range(0, len(pairs), PAIR_CHUNK):
        memento.compute_2d_moments(A, pairs[i:i + PAIR_CHUNK])
        m, cell_counts = memento.get_2d_moments(A, groupby=None)
        chunks.append(m)
    moments = pd.concat(chunks, ignore_index=True)
    cols = [c for c in moments.columns if c.startswith("sg^")]
    M = moments[cols].to_numpy(float, copy=True)  # writable under copy-on-write
    # memento's _corr_from_cov fills pairs whose variance is <= 0 in a group
    # with a 5.0 placeholder and then clips to exactly +/-1 (still true in the
    # current release). ht_2d_moments skips groups with |r| == 1; do the same
    # here, or one group where the target is undetected (e.g. a 282-cell
    # lymphatic group with ANTXR2 raw mean 0) pegs every pair at +1 and drags
    # the cell-weighted mean of all of them upward.
    M[np.abs(M) >= 1] = np.nan
    w = np.array([cell_counts[c] for c in cols], float)
    ok = np.isfinite(M)
    corr = np.nansum(np.where(ok, M, 0) * w, 1) / np.where(ok.any(1), (ok * w).sum(1), np.nan)
    out = pd.DataFrame({"target": moments["gene_1"].values, "gene": moments["gene_2"].values,
                        "corr": corr, "n_groups_valid": ok.sum(1)})
    if groups_out:
        # per-group point estimates plus each gene's raw mean within each
        # group, for donor-level tests and per-group expression adjustment
        os.makedirs(os.path.dirname(os.path.abspath(groups_out)), exist_ok=True)
        names = [c.split("^", 1)[1] for c in cols]
        g = pd.DataFrame(M.astype(np.float32), columns=names)
        g.insert(0, "gene", moments["gene_2"].values)
        g.insert(0, "target", moments["gene_1"].values)
        g.to_parquet(groups_out + "_corr.parquet", index=False)
        lab = A.obs["memento_group"].astype(str).str.split("^", n=1).str[1].values
        X = A.X.tocsr()
        rm = {n: np.asarray(X[lab == n].mean(0)).ravel() for n in names}
        pd.DataFrame(rm, index=A.var_names).to_parquet(groups_out + "_rawmean.parquet")
        pd.Series({n: int(cell_counts["sg^" + n]) for n in names}).to_csv(groups_out + "_ncells.csv")
    return out.join(stats, on="gene")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("out")
    ap.add_argument("--targets", nargs="+", default=["ANTXR2"])
    ap.add_argument("--shrinkage", type=float, default=0.5)
    ap.add_argument("--trim", type=float, default=0.1)
    ap.add_argument("--q", type=float, default=Q_DEFAULT)
    ap.add_argument("--groups-out", default=None,
                    help="prefix: also write per-group correlations, per-group raw means, group sizes")
    a = ap.parse_args()
    t0 = time.time()
    df = scan(a.npz, a.targets, a.shrinkage, a.trim, a.q, a.groups_out)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"{os.path.basename(a.npz)}: {df.target.nunique() if len(df) else 0} targets, "
          f"{len(df)} pairs, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
