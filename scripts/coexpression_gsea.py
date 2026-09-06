"""Pre-ranked GSEA (Subramanian et al.) on ANTXR2's per-cell-type gene
ranking, via gseapy's prerank -- as opposed to coexpression_enrichment.py's
over-representation analysis (Enrichr, a hypergeometric test on a fixed
top-50 gene list).

This uses the ENTIRE tested gene universe per cell type (~3,700+ genes), not
just the top 50. That's the statistically appropriate use of GSEA: the
running-sum statistic needs the full ranking to detect a gene set skewed
toward one end, even if no individual member of that set makes a fixed top-N
cutoff. Positive NES means a gene set's members skew toward the positive end
of the ranking; negative NES means they skew toward the negative end.

Two ranking METRICS are supported (--metric):
  mean    (default, original) -- donor-averaged point-estimate correlation
          (`.uns['memento_correlations']['donor_averaged']`). Simple, but
          treats a noisy, imprecise estimate the same as a tight, confident
          one of equal magnitude.
  zscore  -- coef/se from coexpression_full_ht.py's full-dataset one-sample
          bootstrap test (ALL usable donors, not a discovery half -- see that
          script's docstring). This down-weights spurious/noisy correlations
          (large magnitude but large se) relative to precise ones, which a
          plain point-estimate ranking cannot distinguish. Per user request
          this session, run alongside (not instead of) the mean-based
          ranking, to compare whether pathway-level conclusions are
          metric-sensitive. Still uses every gene with a valid test --
          matching the "full ranking, not top-50" principle above.

Requires internet access to Enrichr (gseapy pulls the named gene-set
libraries from the same backend as coexpression_enrichment.py's ORA). Run in
the `antxr2` conda env:
`conda run -n antxr2 python scripts/coexpression_gsea.py [--metric mean|zscore]`.
"""
import argparse
import os

import anndata
import gseapy as gp
import numpy as np
import pandas as pd

from config import COEXPR_CELL_TYPES, COEXPR_FIGURES_DIR, COEXPR_OUTPUT_H5AD
from coexpression_full_ht import OUT_DIR as HT_DIR

GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
FDR_THRESH = 0.25  # GSEA's own convention (Subramanian et al.) -- looser than ORA's 0.05,
                    # since NES/FDR here is testing gene-SET skew, not per-gene significance
OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "gsea")
OUT_DIR_ZSCORE = os.path.join(COEXPR_FIGURES_DIR, "gsea_zscore")


def build_ranking(ct):
    a = anndata.read_h5ad(COEXPR_OUTPUT_H5AD)
    mc = a.uns["memento_correlations"]
    sym_map = mc["gene_symbol"]
    gene_order = mc["antxr2_gene_order"]
    da = mc["donor_averaged"][ct]
    scores = da["antxr2_vs_all_mean"]

    symbols = np.array([sym_map.get(g, g) for g in gene_order])
    df = pd.DataFrame({"gene": symbols, "score": scores})
    df = df.dropna(subset=["score"])
    # gseapy needs one score per gene symbol; a small number of Ensembl IDs can
    # share a symbol in principle (not the case in this dataset, verified
    # earlier, but guard anyway) -- keep the entry with the larger |score|.
    df["abs_score"] = df["score"].abs()
    df = df.sort_values("abs_score", ascending=False).drop_duplicates("gene", keep="first")
    rnk = df.set_index("gene")["score"].sort_values(ascending=False)
    return rnk


def build_ranking_zscore(ct):
    """z = coef/se from the full-dataset (all usable donors) one-sample bootstrap
    test -- see coexpression_full_ht.py. Uses every gene with a valid test, same
    "full ranking" principle as build_ranking, just a different per-gene score."""
    df = pd.read_csv(os.path.join(HT_DIR, f"{ct}_full_dataset_ht.csv"), comment="#")
    df = df.dropna(subset=["z"]).rename(columns={"gene_symbol": "gene", "z": "score"})
    df["abs_score"] = df["score"].abs()
    df = df.sort_values("abs_score", ascending=False).drop_duplicates("gene", keep="first")
    rnk = df.set_index("gene")["score"].sort_values(ascending=False)
    return rnk


def run_gsea(ct, metric="mean"):
    rnk = build_ranking(ct) if metric == "mean" else build_ranking_zscore(ct)
    print(f"\n  ranking ({metric}): {len(rnk)} genes, score range [{rnk.min():.3f}, {rnk.max():.3f}]")
    pre = gp.prerank(
        rnk=rnk, gene_sets=GENE_SETS, organism="human",
        min_size=5, max_size=1000, permutation_num=1000,
        outdir=None, seed=0, threads=4, no_plot=True,
    )
    return pre.res2d


def summarize(res, ct, top_n_terms=15):
    if res.empty:
        return res
    res = res.copy()
    # gseapy's prerank res2d has no separate library column -- Term is actually
    # "<gene_set_library>__<term name>"; split it back apart for reporting.
    split = res["Term"].str.split("__", n=1, expand=True)
    res["Gene_set"] = split[0]
    res["Term"] = split[1]
    res["FDR q-val"] = pd.to_numeric(res["FDR q-val"], errors="coerce")
    res["NES"] = pd.to_numeric(res["NES"], errors="coerce")
    sig = res[res["FDR q-val"] < FDR_THRESH].sort_values("FDR q-val").reset_index(drop=True)
    print(f"\n=== {ct}: {len(sig)} significant terms (FDR q < {FDR_THRESH}) "
          f"across {res['Gene_set'].nunique()} libraries ({len(res)} terms tested) ===")
    pos = sig[sig["NES"] > 0]
    neg = sig[sig["NES"] < 0]
    print(f"  {len(pos)} positively enriched (skew toward ANTXR2-positively-correlated genes), "
          f"{len(neg)} negatively enriched (skew toward ANTXR2-negatively-correlated genes)")
    for label, sub_all in [("POSITIVE (NES > 0)", pos), ("NEGATIVE (NES < 0)", neg)]:
        sub = sub_all.head(top_n_terms)
        if sub.empty:
            continue
        print(f"\n  -- {label}, top {len(sub)} of {len(sub_all)} --")
        for _, row in sub.iterrows():
            genes = row["Lead_genes"].split(";")
            gene_preview = ";".join(genes[:8]) + (f" (+{len(genes)-8} more)" if len(genes) > 8 else "")
            print(f"    [{row['Gene_set']:26s}] {row['Term']:55s} "
                  f"NES={row['NES']:+.2f}  FDR={row['FDR q-val']:.3f}  lead_genes={gene_preview}")
    return sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", choices=["mean", "zscore"], default="mean",
                     help="mean: original donor-averaged point estimate; zscore: coef/se from "
                          "the full-dataset one-sample bootstrap test (coexpression_full_ht.py)")
    args = ap.parse_args()
    out_dir = OUT_DIR if args.metric == "mean" else OUT_DIR_ZSCORE

    os.makedirs(out_dir, exist_ok=True)
    for ct in COEXPR_CELL_TYPES:
        print(f"\n\n######## {ct} ({args.metric}) ########")
        res = run_gsea(ct, metric=args.metric)
        res.to_csv(os.path.join(out_dir, f"{ct}_gsea_all.csv"), index=False)
        sig = summarize(res, ct)
        sig.to_csv(os.path.join(out_dir, f"{ct}_gsea_significant.csv"), index=False)
    print(f"\n\nwrote per-cell-type GSEA ({args.metric}) results to {out_dir}/")


if __name__ == "__main__":
    main()
