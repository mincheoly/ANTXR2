"""How good a lymphatic endothelial marker is ANTXR2, against canonical markers?

Per tissue with >= MIN_CELLS cells labelled "endothelial cell of lymphatic
vessel" in the endothelial extract (comp_cell_extract.py output):

  1. Cell-level separation of lymphatic (LEC) from blood-vessel (BEC)
     endothelium: AUROC of depth-normalized expression (counts / library *
     1e4), plus detection rate in each. AUROC ties at zero count are scored
     0.5, so a sparsely detected gene is capped well below 1 even when every
     detected cell is lymphatic -- that cap is the point of reporting it.
  2. Specificity beyond endothelium: pooled CP10K in LEC, BEC, and the
     fibroblast / immune / epithelial extracts of the same tissue.

Descriptive, cell-weighted pooling within tissue; donors are not separately
weighted here (the donor-paired LEC vs BEC fold-change is in the session
notes). Labels are the atlases' own.

Usage: python lymphatic_marker_strength.py <workspace dir> <out prefix>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

MIN_CELLS = 20
LEC = "endothelial cell of lymphatic vessel"
GENES = ["ANTXR2", "PROX1", "LYVE1", "CCL21", "TFF3", "MMRN1", "FLT4", "PDPN", "RELN"]


def load(path, genes):
    z = np.load(path, allow_pickle=False)
    names = list(z["genes"].astype(str))
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    cols = X[:, [names.index(g) for g in genes]].toarray().astype(float)
    return cols, z["total_counts"].astype(float), z["cell_type"].astype(str)


def cp10k(cols, lib, mask):
    return cols[mask].sum(0) / lib[mask].sum() * 1e4


def main(ws, out):
    sep, spec = [], []
    for f in sorted(glob.glob(os.path.join(ws, "comp_cells", "endothelial_*.npz"))):
        tissue = os.path.basename(f)[len("endothelial_"):-4]
        cols, lib, ct = load(f, GENES)
        lec = ct == LEC
        if lec.sum() < MIN_CELLS or (~lec).sum() < MIN_CELLS:
            continue
        norm = cols / lib[:, None] * 1e4
        for i, g in enumerate(GENES):
            auc = mannwhitneyu(norm[lec, i], norm[~lec, i]).statistic / (lec.sum() * (~lec).sum())
            sep.append(dict(tissue=tissue, gene=g, auroc=auc, n_lec=int(lec.sum()), n_bec=int((~lec).sum()),
                            det_lec=(cols[lec, i] > 0).mean(), det_bec=(cols[~lec, i] > 0).mean()))
        a = GENES.index("ANTXR2")
        row = dict(tissue=tissue, LEC=cp10k(cols, lib, lec)[a], BEC=cp10k(cols, lib, ~lec)[a])
        for comp, path in [("fibroblast", os.path.join(ws, "fib_cells", f"{tissue}.npz")),
                           ("immune", os.path.join(ws, "comp_cells", f"immune_{tissue}.npz")),
                           ("epithelial", os.path.join(ws, "comp_cells", f"epithelial_{tissue}.npz"))]:
            if os.path.exists(path):
                c, l, _ = load(path, ["ANTXR2"])
                row[comp] = cp10k(c, l, np.ones(len(l), bool))[0]
        spec.append(row)
    s, p = pd.DataFrame(sep), pd.DataFrame(spec)
    s.to_csv(out + "_separation.csv", index=False)
    p.to_csv(out + "_compartments.csv", index=False)
    pd.set_option("display.width", 250)
    print("AUROC, lymphatic vs blood endothelium (rows: tissue)")
    print(s.pivot(index="tissue", columns="gene", values="auroc")[GENES].round(2).to_string())
    print("\nmedian across tissues")
    print(s.groupby("gene")[["auroc", "det_lec", "det_bec"]].median().loc[GENES].round(3).to_string())
    print("\nANTXR2 CP10K by compartment, same tissue")
    print(p.round(2).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
