"""ANTXR2 and IFN-gamma in mouse colon organoids (user, 2026-11-20: "I would
still try the IFN-gamma organoid dataset"). GEO GSE255366: one organoid line,
untreated / IFN-gamma 24 h / 48 h, 10x (`/data/antxr2/raw/colon_organoid_ifng/`).

Mouse genes are mapped to human symbols by MGI one-to-one homology
(HOM_MouseHumanSequence.rpt), so the human program library, IBDverse weights and
gene spans apply. Cells are clustered (Leiden) on genes IFN-gamma does not induce and clusters are
labelled by lineage markers, so the strata cannot absorb the IFN response.
Clusters still split by sample (IFN-gamma remodels the organoid); Antxr2 is
confined to the stem / cycling clusters (1.3-2.0 counts / cell; ~0.02
elsewhere), so the test runs there, strata = sample x cluster.

1. mean    ANTXR2 per 10k by lineage x sample.
2. couple  ANTXR2's noise-corrected correlation with each program per sample,
           pooled within lineage strata (and per lineage), excess over 200
           matched anchors, cell bootstrap (program_htec_ifn.couple).

    -> OUT_DIR/organoid_ifng/{cells.csv, mean.csv, coupling_<lineage>.csv}

Usage:
    PYTHONPATH=<memento worktree> ANTXR2_LANDSCAPE_RUN=ref python program_organoid_ifng.py [--units ...] [--boot N]
"""
import argparse

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

from program_config import OUT_DIR, ROOT, TARGET
import program_htec_ifn as H

RAW = ROOT / "raw" / "colon_organoid_ifng"
HOM = ROOT / "annotation" / "HOM_MouseHumanSequence.rpt"
OUT = OUT_DIR / "organoid_ifng"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = {"GSM8070573_SP191_007": "untreated", "GSM8070574_SP191_008": "IFNg_24h",
           "GSM8070575_SP191_009": "IFNg_48h"}
MIN_UMI = 1000
MAX_MT = 0.25
MIN_GROUP = 20   # cells per sample x cluster stratum
LINEAGE = {   # mouse symbols; no IFN-responsive genes
    "stem": ["Lgr5", "Ascl2", "Olfm4", "Smoc2"],
    "cycling": ["Mki67", "Top2a", "Cenpf", "Ccnb1"],
    "colonocyte": ["Car1", "Slc26a3", "Krt20", "Aqp8", "Car4"],
    "goblet": ["Muc2", "Clca1", "Zg16", "Tff3", "Fcgbp"],
    "enteroendocrine": ["Chga", "Chgb", "Neurod1"],
    "tuft": ["Dclk1", "Pou2f3", "Trpm5"],
}


def one_to_one():
    h = pd.read_csv(HOM, sep="\t", usecols=["DB Class Key", "Common Organism Name", "Symbol"])
    h["sp"] = np.where(h["Common Organism Name"].str.startswith("mouse"), "mouse", "human")
    n = h.groupby(["DB Class Key", "sp"]).size().unstack(fill_value=0)
    keys = n.index[(n.get("mouse", 0) == 1) & (n.get("human", 0) == 1)]
    h = h[h["DB Class Key"].isin(keys)].pivot(index="DB Class Key", columns="sp", values="Symbol")
    return pd.Series(h.human.values, index=h.mouse.values)


def lineage_clusters(X, o, mouse):
    """Leiden clusters on genes IFN-gamma does not induce (pseudobulk log2FC <= 1
    in both treated samples vs untreated), labelled by lineage-marker means."""
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(X=X.astype(np.float32), obs=o.copy(), var=pd.DataFrame(index=pd.Index(mouse).astype(str)))
    a.var_names_make_unique()
    sc.pp.normalize_total(a, target_sum=1e4)
    pb = pd.DataFrame({s: np.asarray(a.X[(o["sample"] == s).values].mean(0)).ravel() for s in SAMPLES.values()},
                      index=a.var_names)
    lfc = np.log2((pb[["IFNg_24h", "IFNg_48h"]] + 0.01).div(pb.untreated + 0.01, axis=0))
    induced = (lfc > 1).any(axis=1).values
    sc.pp.log1p(a)
    b = a[:, ~induced & ~a.var_names.str.startswith("mt-")].copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000)
    b = b[:, b.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30, random_state=0)
    sc.pp.neighbors(b, random_state=0)
    sc.tl.leiden(b, resolution=1.0, random_state=0)
    cl = b.obs.leiden.astype(str).values
    sc_ = {}
    for k, gs in LINEAGE.items():
        gs = [g for g in gs if g in a.var_names]
        sc_[k] = np.asarray(a[:, gs].X.mean(1)).ravel()
    S = pd.DataFrame(sc_)
    Sz = (S - S.mean()) / S.std()
    M_ = Sz.groupby(cl).mean()
    lab = M_.idxmax(1).where(M_.max(1) > 0.25, "unassigned")
    comp = pd.crosstab(cl, o["sample"].values, normalize="index")
    print(pd.concat([M_.round(2), lab.rename("label"), comp.round(2),
                     pd.Series(cl).value_counts().rename("n")], axis=1).to_string(), flush=True)
    print(f"{induced.sum()} IFN-induced genes left out of clustering")
    return lab.reindex(cl).values, cl


def load():
    mats, obs = [], []
    for s, lab in SAMPLES.items():
        X = sio.mmread(RAW / f"{s}_matrix.mtx.gz").tocsr().T.tocsr()
        f = pd.read_csv(RAW / f"{s}_features.tsv.gz", sep="\t", header=None)
        mouse = f[1].values
        tot = np.asarray(X.sum(1)).ravel()
        mt = np.asarray(X[:, np.char.startswith(mouse.astype(str), "mt-")].sum(1)).ravel() / tot
        keep = (tot >= MIN_UMI) & (mt <= MAX_MT)
        mats.append(X[keep])
        obs.append(pd.DataFrame({"sample": lab, "total_counts": tot[keep], "pct_mt": mt[keep]}))
    X = sp.vstack(mats).tocsr()
    o = pd.concat(obs, ignore_index=True)
    o["lineage"], o["cluster"] = lineage_clusters(X, o, mouse)
    # map to human one-to-one orthologs
    o2o = one_to_one()
    m = np.array([g in o2o.index for g in mouse])
    human = o2o.reindex(mouse[m]).values
    Xh = X[:, np.where(m)[0]]
    dup = pd.Index(human).duplicated()
    return Xh[:, ~dup], o, np.asarray(human[~dup]), X, mouse


def cmd_mean(o, Xm, mouse):
    j = list(mouse).index("Antxr2")
    o = o.assign(umi=np.asarray(Xm[:, j].todense()).ravel())
    t = o.groupby(["lineage", "sample"]).agg(cells=("umi", "size"), umis=("umi", "sum"),
                                             depth=("total_counts", "sum")).reset_index()
    t["per10k"] = t.umis / t.depth * 1e4
    t["raw_mean"] = t.umis / t.cells
    t.to_csv(OUT / "mean.csv", index=False)
    print(t.pivot(index="lineage", columns="sample", values="cells").to_string())
    print(t.pivot(index="lineage", columns="sample", values="per10k").round(3).to_string())
    print(t.pivot(index="lineage", columns="sample", values="raw_mean").round(3).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["mean", "couple"])
    ap.add_argument("--lineages", nargs="*", default=["stem_cycling", "all"])
    ap.add_argument("--units", nargs="*")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"])
    a = ap.parse_args()
    X, o, genes, Xm, mouse = load()
    assert TARGET in genes
    o.to_csv(OUT / "cells.csv", index=False)
    if a.cmd == "mean":
        cmd_mean(o, Xm, mouse)
        return
    tag = "" if not a.units else "_subset"
    for lin in a.lineages:
        # strata: sample x Leiden cluster (clusters on non-induced genes)
        m = {"all": np.ones(len(o), bool), "stem_cycling": o.lineage.isin(["stem", "cycling"]).values}.get(
            lin, (o.lineage == lin).values)
        n = o[m].groupby(["sample", "cluster"]).size()
        big = o.set_index(["sample", "cluster"]).index.isin(n.index[n >= MIN_GROUP])
        m = m & big
        obs = pd.DataFrame({"donor": "org", "condition": o["sample"].values[m], "fine": o.cluster.values[m],
                            "memento_size_factor": o.total_counts.values[m] / 1e4})
        out = H.couple(X[m], obs, genes, list(SAMPLES.values()), a.boot, a.units, a.backend, lin)
        out.to_csv(OUT / f"coupling_{lin}{tag}.csv", index=False)
        print(f"-> {OUT / f'coupling_{lin}{tag}.csv'}", flush=True)


if __name__ == "__main__":
    main()
