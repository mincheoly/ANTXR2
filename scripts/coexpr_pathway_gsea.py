"""Pre-ranked GSEA on the per-population ANTXR2 one-sample co-expression scans.

Input: <population>.csv files from antxr2_coexpr.py (gene, corr_coef, corr_se,
qval, raw_mean, det_rate). Exploratory characterization only -- the per-gene
p-values from the one-sample bootstrap are within-group sampling p-values
(median 4 donors per population) and are not used anywhere here.

Two ranking statistics per population:

  raw    corr_coef as estimated.
  resid  corr_coef minus its expected value at the gene's expression level:
         a running median of corr_coef over genes sorted by log10(raw_mean)
         (window WINDOW genes). In 98/100 populations the scan shows
         Spearman(log raw_mean, corr_coef) > 0 (median +0.25), and the whole
         distribution is shifted positive (memento default shrinkage/trim), so
         an unadjusted ranking pushes well-expressed gene sets to the top
         regardless of ANTXR2. The residual removes both the global shift and
         the expression trend; it does not remove a cell type's dominant
         program travelling with any ANTXR2-like gene (that needs anchors).

Gene sets: MSigDB Hallmark 2020, Reactome 2022, KEGG 2021, GO BP 2023
(Enrichr libraries), size 15-500 within each population's tested universe.

Usage:
    python coexpr_pathway_gsea.py <outdir> --pops <csv> [<csv> ...]
        [--stat resid raw] [--perm 1000] [--workers 10] [--libs HALLMARK REACTOME KEGG]

Memory: each worker holds the gene-set dict and a full prerank run; with all
four libraries (~7.6k terms) 8 workers exhausted a 13 GB host. Without GOBP
(~2.2k terms) and 3 workers it fits comfortably.
"""
import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

WINDOW = 301
LIBS = {"HALLMARK": "MSigDB_Hallmark_2020", "REACTOME": "Reactome_2022",
        "KEGG": "KEGG_2021_Human", "GOBP": "GO_Biological_Process_2023"}
MIN_SIZE, MAX_SIZE = 15, 500


def load_genesets(cache, libs=tuple(LIBS)):
    import json
    if os.path.exists(cache):
        return json.load(open(cache))
    import gseapy as gp
    gs = {}
    for tag in libs:
        lib = LIBS[tag]
        for term, genes in gp.get_library(lib, organism="Human").items():
            gs[f"{tag}::{term}"] = sorted(set(g.upper() for g in genes))
    json.dump(gs, open(cache, "w"))
    return gs


def population_name(path):
    pop = os.path.splitext(os.path.basename(path))[0]
    # fibroblast scans (coexpr_out) carry no compartment prefix
    if not pop.split("_")[0] in ("endothelial", "epithelial", "immune", "fibroblast"):
        pop = "fibroblast_" + pop
    return pop


def ranking(path, stat):
    df = pd.read_csv(path).dropna(subset=["corr_coef", "raw_mean"])
    df["gene"] = df.gene.astype(str).str.upper()
    df = df[df.gene != "ANTXR2"].drop_duplicates("gene")
    df = df.sort_values("raw_mean").reset_index(drop=True)
    if stat == "raw":
        s = df.corr_coef
    else:
        trend = df.corr_coef.rolling(WINDOW, center=True, min_periods=WINDOW // 3).median()
        s = df.corr_coef - trend
    return pd.Series(s.values, index=df.gene.values).dropna().sort_values(ascending=False)


def run_one(args):
    path, stat, outdir, gs_cache, perm, seed = args
    import gseapy as gp
    pop = population_name(path)
    out = os.path.join(outdir, stat, f"{pop}.csv")
    if os.path.exists(out):
        return pop, stat, "cached", 0.0
    t0 = time.time()
    rnk = ranking(path, stat)
    gs = load_genesets(gs_cache)
    universe = set(rnk.index)
    gs = {k: [g for g in v if g in universe] for k, v in gs.items()}
    gs = {k: v for k, v in gs.items() if MIN_SIZE <= len(v) <= MAX_SIZE}
    res = gp.prerank(rnk=rnk, gene_sets=gs, permutation_num=perm, min_size=MIN_SIZE,
                     max_size=MAX_SIZE, threads=1, seed=seed, outdir=None, verbose=False)
    r = res.res2d.rename(columns={"Term": "term", "NES": "nes", "NOM p-val": "pval",
                                  "FDR q-val": "fdr", "Lead_genes": "lead_genes",
                                  "Tag %": "tag_pct"})
    r = r[["term", "nes", "pval", "fdr", "tag_pct", "lead_genes"]]
    r.insert(0, "population", pop)
    r.insert(1, "stat", stat)
    r["n_genes_ranked"] = len(rnk)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    r.to_csv(out, index=False)
    return pop, stat, "ok", time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--pops", nargs="+", required=True)
    ap.add_argument("--stat", nargs="+", default=["resid"])
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--libs", nargs="+", default=list(LIBS), choices=list(LIBS),
                    help="gene-set libraries; the cache in <outdir> is built once from these")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    cache = os.path.join(a.outdir, "genesets.json")
    load_genesets(cache, a.libs)
    jobs = [(p, s, a.outdir, cache, a.perm, a.seed) for s in a.stat for p in a.pops]
    with ProcessPoolExecutor(a.workers) as ex:
        for pop, stat, status, dt in ex.map(run_one, jobs):
            print(f"{status} {stat} {pop} {dt:.0f}s", flush=True)


if __name__ == "__main__":
    main()
