"""Remove contaminating / non-target cells from endothelial population matrices.

Motivation: in the endothelial ANTXR2 scans, ANTXR2's correlation with
interstitial fibroblast genes beats its matched anchors in 11/25 tissues, so
the endothelial ECM enrichment could come from endothelial-fibroblast
doublets or ambient fibroblast RNA (ANTXR2 is high in fibroblasts). This
writes filtered copies of the comp_cell_extract.py .npz files with cells
removed whose summed marker counts reach a threshold.

Filters (each off unless its threshold is > 0):
  fibroblast  DCN, LUM, COL1A1, COL1A2        --min-counts (default 2)
  lymphatic   PROX1, CCL21, plus any cell labelled lymphatic
                                              --lymphatic-min (default 0 = off)
  macrophage  C1QA, C1QB, CD163, TYROBP       --macrophage-min (default 0 = off)

The lymphatic filter exists because ANTXR2 is 2-15x higher in lymphatic than
blood-vessel endothelium (14 tissues, 39/40 donors), so within a mixed
"endothelial cell" label its correlations partly encode lymphatic identity.
PROX1+CCL21 >= 1 recovers 93-100% of labelled lymphatic cells in most
tissues; LYVE1/PDPN/MMRN1 are NOT used because sinusoidal and some blood
endothelium express them (they flag 57% of liver and 87% of spleen ECs).
Spleen (ambient CCL21) and stomach (all cells labelled lymphatic) cannot be
separated this way and should be excluded with --exclude.

Single ambient counts are common (up to 66% of cells carry >=1 fibroblast
count in fat and ovary), so the fibroblast threshold is >=2 by default; >=1
is the strict sensitivity version. Filter genes' own correlations are
mechanically truncated, so each contaminant must be read on held-out genes
(endothelial_source_check.py), and filter genes should be dropped from any
downstream ranking (ALL_FILTER_GENES).

Caveat: deeper cells carry more ambient counts, so marker filters
preferentially remove deep cells; ANTXR2 raw mean and detection
before/after are recorded.

Usage: python endothelial_decontam.py <comp_cells dir> <out dir> [--min-counts 2]
           [--lymphatic-min 1] [--macrophage-min 2] [--exclude ts_Spleen ts_Stomach]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

FILTER_GENES = ["DCN", "LUM", "COL1A1", "COL1A2"]
LYMPHATIC_GENES = ["PROX1", "CCL21"]
LYMPHATIC_LABEL = "endothelial cell of lymphatic vessel"
MACROPHAGE_GENES = ["C1QA", "C1QB", "CD163", "TYROBP"]
ALL_FILTER_GENES = FILTER_GENES + LYMPHATIC_GENES + MACROPHAGE_GENES


def marker_sum(X, genes, names):
    return np.asarray(X[:, [genes.index(g) for g in names]].sum(1)).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("outdir")
    ap.add_argument("--min-counts", type=int, default=2)
    ap.add_argument("--lymphatic-min", type=int, default=0)
    ap.add_argument("--macrophage-min", type=int, default=0)
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="population suffixes to skip, e.g. ts_Spleen")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(a.indir, "endothelial_*.npz"))):
        pop = os.path.basename(f)[:-4]
        if pop[len("endothelial_"):] in a.exclude:
            continue
        z = dict(np.load(f, allow_pickle=False))
        genes = list(z["genes"].astype(str))
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        drop_fib = marker_sum(X, genes, FILTER_GENES) >= a.min_counts
        drop_lym = np.zeros(X.shape[0], bool)
        if a.lymphatic_min > 0:
            drop_lym = ((marker_sum(X, genes, LYMPHATIC_GENES) >= a.lymphatic_min)
                        | (z["cell_type"].astype(str) == LYMPHATIC_LABEL))
        drop_mac = np.zeros(X.shape[0], bool)
        if a.macrophage_min > 0:
            drop_mac = marker_sum(X, genes, MACROPHAGE_GENES) >= a.macrophage_min
        keep = ~(drop_fib | drop_lym | drop_mac)
        j = genes.index("ANTXR2")
        a_before = X[:, j].toarray().ravel()
        umi_before = float(np.median(z["total_counts"]))
        Xk = X[keep]
        z.update(data=Xk.data, indices=Xk.indices, indptr=Xk.indptr,
                 shape=np.array(Xk.shape), donor=z["donor"][keep],
                 cell_type=z["cell_type"][keep], total_counts=z["total_counts"][keep])
        np.savez(os.path.join(a.outdir, os.path.basename(f)), **z)
        rows.append(dict(population=pop, n_before=len(keep), n_after=int(keep.sum()),
                         frac_removed=1 - keep.mean(), frac_fibroblast=drop_fib.mean(),
                         frac_lymphatic=drop_lym.mean(), frac_macrophage=drop_mac.mean(),
                         antxr2_mean_before=a_before.mean(), antxr2_mean_after=a_before[keep].mean(),
                         antxr2_det_before=(a_before > 0).mean(), antxr2_det_after=(a_before[keep] > 0).mean(),
                         med_umi_before=umi_before,
                         med_umi_after=float(np.median(z["total_counts"]))))
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(a.outdir, "decontam_report.csv"), index=False)
    pd.set_option("display.width", 250)
    print(r.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
