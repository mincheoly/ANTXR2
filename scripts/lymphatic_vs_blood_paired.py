"""Donor-paired lymphatic vs blood-vessel endothelium fold change, per tissue.

Means-level (pseudobulk) comparison, the project's strongest evidence class.
For every endothelial extract (comp_cell_extract.py output) with cells
labelled "endothelial cell of lymphatic vessel": per donor with >= MIN_CELLS
lymphatic and >= MIN_CELLS other endothelial cells, pseudobulk CP10K of every
gene in each pool and log2((LEC + PSEUDO) / (BEC + PSEUDO)). Per tissue: the
median over donors, the number of donors in which the gene is higher in LEC,
and the gene's percentile among genes expressed >= EXPR_FLOOR CP10K in either
pool (canonical lymphatic markers reported alongside as label checks).

Usage: python lymphatic_vs_blood_paired.py <comp_cells dir> <out.csv>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

LEC = "endothelial cell of lymphatic vessel"
MIN_CELLS = 20
PSEUDO = 0.05
EXPR_FLOOR = 0.1
REPORT = ["ANTXR2", "PROX1", "LYVE1", "MRC1", "CCL21", "TFF3"]


def main(indir, out):
    rows = []
    for f in sorted(glob.glob(os.path.join(indir, "endothelial_*.npz"))):
        z = np.load(f, allow_pickle=False)
        ct = z["cell_type"].astype(str)
        if (ct == LEC).sum() < MIN_CELLS:
            continue
        genes = z["genes"].astype(str)
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        lib, don = z["total_counts"], z["donor"].astype(str)
        fcs, pl, pb = [], [], []
        for d in np.unique(don):
            a, b = (don == d) & (ct == LEC), (don == d) & (ct != LEC)
            if a.sum() < MIN_CELLS or b.sum() < MIN_CELLS:
                continue
            ca = np.asarray(X[a].sum(0)).ravel() / lib[a].sum() * 1e4
            cb = np.asarray(X[b].sum(0)).ravel() / lib[b].sum() * 1e4
            fcs.append(np.log2((ca + PSEUDO) / (cb + PSEUDO)))
            pl.append(ca)
            pb.append(cb)
        if not fcs:
            continue
        F = np.vstack(fcs)
        fc, cl, cb = np.median(F, 0), np.median(pl, 0), np.median(pb, 0)
        expressed = np.maximum(cl, cb) >= EXPR_FLOOR
        row = dict(tissue=os.path.basename(f)[len("endothelial_"):-4], n_donors=len(fcs),
                   n_lec=int((ct == LEC).sum()))
        idx = {g: i for i, g in enumerate(genes)}
        for g in REPORT:
            i = idx[g]
            row[f"{g}_pct"] = float((fc[expressed] < fc[i]).mean() * 100)
        a = idx["ANTXR2"]
        row.update(antxr2_cp10k_lec=cl[a], antxr2_cp10k_bec=cb[a], antxr2_log2fc=fc[a],
                   antxr2_donors_up=int((F[:, a] > 0).sum()))
        rows.append(row)
    r = pd.DataFrame(rows)
    r.to_csv(out, index=False)
    pd.set_option("display.width", 250)
    print(r.round(2).to_string(index=False))
    print(f"\n{len(r)} tissues; ANTXR2 higher in LEC in {r.antxr2_donors_up.sum()}/{r.n_donors.sum()} "
          f"donors; log2FC range {r.antxr2_log2fc.min():.2f}-{r.antxr2_log2fc.max():.2f}; "
          f"percentile range {r.ANTXR2_pct.min():.1f}-{r.ANTXR2_pct.max():.1f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
