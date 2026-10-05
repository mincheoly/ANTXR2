"""ANTXR2 vs Wnt target output along progenitor -> differentiated axes, per tissue.

Descriptive, means-level only. For each tissue with a progenitor/differentiated
split in its epithelial extract (comp_cell_extract.py output), per donor:
pseudobulk CP10K of every gene in the progenitor and the differentiated
cell-type pool, then log2(diff / prog). Donor-paired: only donors with
>= MIN_CELLS in both pools contribute, and the per-gene effect is the median
over those donors.

ANTXR2 and the Wnt target set are then placed within the genome-wide log2FC
distribution of genes expressed in that tissue's pair (percentile), which is
the comparison the gut result was missing: many genes move with
differentiation, so "ANTXR2 goes up" means little without its rank.

Wnt targets are the general, non-tissue-specific ones (feedback genes plus
pan-tissue transcriptional targets); the intestinal effectors (ASCL2, OLFM4,
SMOC2, ...) are reported separately for the gut rows only.

Usage: python wnt_axis_cross_tissue.py <comp_cells dir> <outdir>
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

MIN_CELLS = 20
EXPR_FLOOR = 0.1       # CP10K, in either pool (donor-median), for the ranked universe
PSEUDO = 0.05          # CP10K pseudocount for log2FC

WNT_GENERAL = ["AXIN2", "NKD1", "LEF1", "TCF7", "NOTUM", "APCDD1", "SP5",
               "RNF43", "ZNRF3", "LGR5", "LGR6"]
WNT_GUT = ["ASCL2", "OLFM4", "SMOC2", "EPHB2", "EPHB3", "SOX9", "CD44", "MYC"]

# tissue -> (file, progenitor labels, differentiated labels)
AXES = {
    "Small_Intestine": ("epithelial_ts_Small_Intestine.npz",
        ["intestinal crypt stem cell of small intestine", "transit amplifying cell of small intestine"],
        ["enterocyte of epithelium proper of ileum", "enterocyte of epithelium proper of duodenum",
         "enterocyte of epithelium proper of jejunum", "enterocyte of epithelium proper of small intestine",
         "BEST4+ enterocyte"]),
    "Large_Intestine": ("epithelial_ts_Large_Intestine.npz",
        ["intestinal crypt stem cell of colon", "transit amplifying cell of colon"],
        ["enterocyte of epithelium of large intestine", "BEST4+ enterocyte"]),
    "Tongue": ("epithelial_ts_Tongue.npz",
        ["basal cell"], ["stratified squamous epithelial cell"]),
    "Trachea": ("epithelial_ts_Trachea.npz",
        ["basal cell"], ["multiciliated columnar cell of tracheobronchial tree",
                         "mucus secreting cell", "tracheal goblet cell"]),
    "Lung_airway": ("epithelial_ts_Lung.npz",
        ["basal cell"], ["club cell", "respiratory tract goblet cell",
                         "lung multiciliated epithelial cell"]),
    "Prostate": ("epithelial_ts_Prostate.npz",
        ["basal cell of prostate epithelium"], ["luminal cell of prostate epithelium"]),
    "Mammary": ("epithelial_ts_Mammary.npz",
        ["basal cell"], ["luminal epithelial cell of mammary gland"]),
    "Salivary_Gland": ("epithelial_ts_Salivary_Gland.npz",
        ["basal cell"], ["duct epithelial cell"]),
}


def load(path):
    z = np.load(path, allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return X, z["genes"].astype(str), z["donor"].astype(str), z["cell_type"].astype(str), z["total_counts"]


def pool_cp10k(X, lib, mask):
    return np.asarray(X[mask].sum(0)).ravel() / lib[mask].sum() * 1e4


def det_rate(X, mask, j):
    return float((X[mask][:, j].toarray() > 0).mean())


def run_axis(tissue, path, prog, diff):
    X, genes, donor, ct, lib = load(path)
    gidx = {g: i for i, g in enumerate(genes)}
    is_p, is_d = np.isin(ct, prog), np.isin(ct, diff)
    fcs, ps, ds, used = [], [], [], []
    for dn in np.unique(donor):
        mp, md = is_p & (donor == dn), is_d & (donor == dn)
        if mp.sum() < MIN_CELLS or md.sum() < MIN_CELLS:
            continue
        p, d = pool_cp10k(X, lib, mp), pool_cp10k(X, lib, md)
        fcs.append(np.log2((d + PSEUDO) / (p + PSEUDO)))
        ps.append(p); ds.append(d); used.append((dn, int(mp.sum()), int(md.sum())))
    if not used:
        return None, None
    fc = np.median(np.vstack(fcs), 0)
    p_med, d_med = np.median(np.vstack(ps), 0), np.median(np.vstack(ds), 0)
    n_same_sign = (np.sign(np.vstack(fcs)) == np.sign(fc)).sum(0)
    expressed = np.maximum(p_med, d_med) >= EXPR_FLOOR
    univ = fc[expressed]

    def pct(v):
        return float((univ < v).mean() * 100)

    rows = []
    for g in ["ANTXR2"] + WNT_GENERAL + (WNT_GUT if "Intestine" in tissue else []):
        if g not in gidx:
            continue
        j = gidx[g]
        rows.append(dict(tissue=tissue, gene=g,
                         set="target" if g == "ANTXR2" else ("wnt_general" if g in WNT_GENERAL else "wnt_gut"),
                         log2fc_diff_vs_prog=fc[j], pct_in_universe=pct(fc[j]),
                         donors_same_sign=int(n_same_sign[j]), n_donors=len(used),
                         cp10k_prog=p_med[j], cp10k_diff=d_med[j],
                         in_universe=bool(expressed[j]),
                         det_prog=det_rate(X, is_p, j), det_diff=det_rate(X, is_d, j)))
    genes_df = pd.DataFrame(rows)
    w = genes_df[(genes_df.set == "wnt_general") & genes_df.in_universe]
    a = genes_df[genes_df.gene == "ANTXR2"].iloc[0]
    summ = dict(tissue=tissue, n_donors=len(used),
                donors=";".join(f"{d}({a_}/{b})" for d, a_, b in used),
                n_universe=int(expressed.sum()),
                antxr2_log2fc=a.log2fc_diff_vs_prog, antxr2_pct=a.pct_in_universe,
                antxr2_donors_same_sign=int(a.donors_same_sign),
                antxr2_cp10k_prog=a.cp10k_prog, antxr2_cp10k_diff=a.cp10k_diff,
                antxr2_det_prog=a.det_prog, antxr2_det_diff=a.det_diff,
                antxr2_in_universe=bool(a.in_universe),
                wnt_n_expressed=len(w), wnt_median_log2fc=w.log2fc_diff_vs_prog.median(),
                wnt_n_prog_side=int((w.log2fc_diff_vs_prog < 0).sum()),
                wnt_median_pct=w.pct_in_universe.median())
    return summ, genes_df


def main(indir, outdir):
    os.makedirs(outdir, exist_ok=True)
    summ, genes = [], []
    for t, (f, prog, diff) in AXES.items():
        s, g = run_axis(t, os.path.join(indir, f), prog, diff)
        if s is None:
            print(f"{t}: no donor with >= {MIN_CELLS} cells in both pools")
            continue
        summ.append(s); genes.append(g)
    S, G = pd.DataFrame(summ), pd.concat(genes)
    S.to_csv(os.path.join(outdir, "wnt_axis_cross_tissue_summary.csv"), index=False)
    G.to_csv(os.path.join(outdir, "wnt_axis_cross_tissue_genes.csv"), index=False)
    pd.set_option("display.width", 250)
    print(S.drop(columns="donors").round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
