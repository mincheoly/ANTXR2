"""Genome-wide ANTXR2 co-expression per population, memento one-sample test.

For one population (a CSR .npz written by fib_cell_extract.py) this pairs
ANTXR2 against every gene clearing memento's raw-mean filter and runs
memento's one-sample branch (all-ones treatment -> cell-count-weighted mean of
the bootstrapped correlations, p from _compute_asl on the bootstrap null).

Conventions carried from the project:
  q = 0.15                capture rate, the corrected point estimate
  filter_mean_thresh 0.07 unnormalized raw mean, the house reliability floor
  min_cell_count 100      donor x cell_type groups smaller than this are dropped
                          (the 250-cell floor was derived at q=0.10; 100 is the
                          equivalent at q=0.15)
  groups = donor x cell_type, so within-compartment subtype composition cannot
                          leak into the correlation

Writes <out>/<pop>.csv (one row per partner gene), <out>/<pop>_groups.parquet
(per-group point estimates, for pegged/undefined diagnostics) and
<out>/<pop>.json (run metadata).

Usage:
    python antxr2_coexpr.py <pop.npz> <outdir> [--boot N] [--cpus N] [--q F]
"""
import argparse
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")

TARGET = "ANTXR2"
Q_DEFAULT = 0.15
FILTER_MEAN_THRESH = 0.07
MIN_CELL_COUNT = 100
MIN_PERC_GROUP = 0.1
SHRINKAGE = 0.5
TRIM_PERCENT = 0.1


def load_pop(npz_path, q=Q_DEFAULT):
    import anndata as ad
    z = np.load(npz_path, allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.DataFrame({"donor": z["donor"].astype(str),
                        "cell_type": z["cell_type"].astype(str),
                        "total_counts": z["total_counts"]})
    obs.index = [f"c{i}" for i in range(X.shape[0])]
    A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=z["genes"].astype(str)))
    A.var_names_make_unique()
    A.obs["q"] = float(q)
    return A


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    ranked = pv[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def pick_celltypes(obs, k):
    """Top-k cell types by the number of cells sitting in >=MIN_CELL_COUNT
    donor x cell_type groups. Bootstrap cost is linear in the group count, so
    the compartment populations (up to ~20 subtypes x ~10 donors) are capped
    here; the retained subtypes are recorded in the run metadata."""
    g = obs.groupby(["cell_type", "donor"]).size()
    g = g[g >= MIN_CELL_COUNT]
    if g.empty:
        return sorted(obs.cell_type.unique())
    tot = g.groupby(level=0).sum().sort_values(ascending=False)
    return list(tot.index[:k])


def run(npz_path, outdir, num_boot=5000, num_cpus=1, q=Q_DEFAULT, seed=5,
        max_celltypes=None):
    import memento
    pop = os.path.splitext(os.path.basename(npz_path))[0]
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    A = load_pop(npz_path, q=q)
    kept_ct = None
    if max_celltypes and A.obs.cell_type.nunique() > max_celltypes:
        kept_ct = pick_celltypes(A.obs, max_celltypes)
        A = A[A.obs.cell_type.isin(kept_ct)].copy()
    n_cells, n_genes_all = A.shape

    # raw-count statistics on the unfiltered matrix, kept for post-hoc floors
    Xc = A.X.tocsc()
    raw_mean = np.asarray(Xc.sum(0)).ravel() / n_cells
    det = np.diff(Xc.indptr) / n_cells
    stats = pd.DataFrame({"raw_mean": raw_mean, "det_rate": det}, index=A.var_names)
    del Xc
    antxr2_mean = float(stats.loc[TARGET, "raw_mean"]) if TARGET in stats.index else float("nan")

    meta = dict(population=pop, n_cells=int(n_cells), n_genes_input=int(n_genes_all),
                n_donors=int(A.obs.donor.nunique()), q=float(q), num_boot=int(num_boot),
                filter_mean_thresh=FILTER_MEAN_THRESH, min_cell_count=MIN_CELL_COUNT,
                min_perc_group=MIN_PERC_GROUP, shrinkage=SHRINKAGE,
                trim_percent=TRIM_PERCENT, seed=seed, target=TARGET,
                antxr2_raw_mean=antxr2_mean, max_celltypes=max_celltypes,
                celltypes_used=kept_ct if kept_ct is not None
                else sorted(A.obs.cell_type.unique().tolist()))

    memento.setup_memento(A, q_column="q", filter_mean_thresh=FILTER_MEAN_THRESH,
                          min_cell_count=MIN_CELL_COUNT, shrinkage=SHRINKAGE,
                          trim_percent=TRIM_PERCENT)
    memento.create_groups(A, label_columns=["donor", "cell_type"])
    if not list(A.uns["memento"]["groups"]):
        # no donor x cell_type group clears MIN_CELL_COUNT; compute_1d_moments
        # would raise on the empty gene-filter stack
        meta.update(status="skipped_no_groups", n_groups=0, n_genes_tested=0)
        json.dump(meta, open(os.path.join(outdir, f"{pop}.json"), "w"), indent=1)
        print(f"SKIP {pop}: skipped_no_groups", flush=True)
        return meta
    memento.compute_1d_moments(A, min_perc_group=MIN_PERC_GROUP, filter_genes=True)

    groups = list(A.uns["memento"]["groups"])
    meta["n_groups"] = len(groups)
    meta["n_genes_tested"] = int(A.shape[1])
    meta["cells_in_groups"] = int(A.shape[0])
    if TARGET not in set(A.var_names) or len(groups) < 1:
        meta["status"] = "skipped_no_target" if TARGET not in set(A.var_names) else "skipped_no_groups"
        json.dump(meta, open(os.path.join(outdir, f"{pop}.json"), "w"), indent=1)
        print(f"SKIP {pop}: {meta['status']}", flush=True)
        return meta

    partners = [g for g in A.var_names if g != TARGET]
    memento.compute_2d_moments(A, [(TARGET, g) for g in partners])
    trt = pd.DataFrame(np.ones((len(groups), 1)), columns=["intercept"])
    memento.ht_2d_moments(A, treatment=trt, num_boot=num_boot, num_cpus=num_cpus,
                          verbose=0, random_state=seed)
    ht = memento.get_2d_ht_result(A).copy()
    ht = ht.rename(columns={"gene_2": "gene"})[["gene", "corr_coef", "corr_se", "corr_pval"]]
    ht["qval"] = bh(ht.corr_pval.values)
    ht = ht.join(stats, on="gene")
    ht.insert(0, "population", pop)
    ht.to_csv(os.path.join(outdir, f"{pop}.csv"), index=False)

    # per-group point estimates: pegged (|r|=1) and undefined rates by group size
    moments, cell_counts = memento.get_2d_moments(A, groupby=None)
    cols = [c for c in moments.columns if c.startswith("sg^")]
    M = moments[cols].to_numpy(dtype=np.float32)
    gdf = pd.DataFrame(M, columns=[c.split("^", 1)[1] for c in cols])
    gdf.insert(0, "gene", moments["gene_2"].values)
    gdf.to_parquet(os.path.join(outdir, f"{pop}_groups.parquet"), index=False)
    nc = np.array([cell_counts.get(c, np.nan) for c in cols], float)
    meta["group_cells"] = {c.split("^", 1)[1]: (None if np.isnan(v) else int(v))
                           for c, v in zip(cols, nc)}
    meta["pct_pegged"] = float(np.nanmean(np.abs(M) >= 0.999) * 100)
    meta["pct_undefined"] = float(np.mean(~np.isfinite(M)) * 100)
    meta["n_sig_q05"] = int((ht.qval < 0.05).sum())
    meta["status"] = "ok"
    meta["runtime_s"] = round(time.time() - t0, 1)
    json.dump(meta, open(os.path.join(outdir, f"{pop}.json"), "w"), indent=1)
    print(f"OK {pop}: {n_cells} cells, {len(groups)} groups, {A.shape[1]} genes, "
          f"{meta['n_sig_q05']} at q<0.05, {meta['runtime_s']}s", flush=True)
    return meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("outdir")
    ap.add_argument("--boot", type=int, default=5000)
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--q", type=float, default=Q_DEFAULT)
    ap.add_argument("--max-celltypes", type=int, default=None)
    a = ap.parse_args()
    run(a.npz, a.outdir, num_boot=a.boot, num_cpus=a.cpus, q=a.q,
        max_celltypes=a.max_celltypes)
