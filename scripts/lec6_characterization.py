"""LEC6 (ADAMTS4+) vs other postnatal lymphatic subtypes, gut atlas, donor-paired.

Input: gut_gut.npz written by lymphatic_extract.py (all gut-atlas LECs, all
genes; cell_type = "<author subtype>|<assay>").

1. Paired pseudobulk: for every donor x chemistry with >= MIN_CELLS LEC6 and
   >= MIN_CELLS cells of LEC1/LEC3/LEC5 (the other postnatal subtypes), CP10K
   in each pool, log2((LEC6 + PSEUDO) / (other + PSEUDO)); per gene the median
   over paired units and the number of units in which it is higher in LEC6.
   Genes expressed >= EXPR_FLOOR CP10K in either pool are ranked.
2. Contamination / stress checks in the same paired units: fraction of cells
   with >= 2 counts summed over fibroblast markers, fibroblast and stress gene
   CP10K, and ANTXR2 recomputed in cells with zero fibroblast-marker counts.

Usage: python lec6_characterization.py <lymphatic dir>
Writes LEC6_vs_otherLEC_paired.csv and LEC6_contamination_check.csv.
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

LEC6 = ["LEC6 (ADAMTS4+)"]
OTHER = ["LEC1 (ACKR4+)", "LEC3 (ADGRG3+)", "LEC5 (CLDN11+)"]
MIN_CELLS = 20
PSEUDO = 0.05
EXPR_FLOOR = 0.5
FIB = ["DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "PDGFRA"]
STRESS = ["FOS", "JUN", "EGR1", "HSPA1A"]


def main(lym):
    z = np.load(os.path.join(lym, "gut_gut.npz"), allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    genes = z["genes"].astype(str)
    gi = {g: i for i, g in enumerate(genes)}
    lib, don = z["total_counts"], z["donor"].astype(str)
    lab = pd.Series(z["cell_type"].astype(str))
    sub, assay = lab.str.split("|").str[0].values, lab.str.split("|").str[1].values
    fib_counts = np.asarray(X[:, [gi[g] for g in FIB]].sum(1)).ravel()
    a2 = X[:, gi["ANTXR2"]].toarray().ravel()

    def cp(m, cols=None):
        v = np.asarray(X[m].sum(0)).ravel() / lib[m].sum() * 1e4
        return v if cols is None else v[[gi[c] for c in cols]]

    units, checks = [], []
    for d in np.unique(don):
        for a in np.unique(assay):
            m6 = (don == d) & (assay == a) & np.isin(sub, LEC6)
            mo = (don == d) & (assay == a) & np.isin(sub, OTHER)
            if m6.sum() < MIN_CELLS or mo.sum() < MIN_CELLS:
                continue
            units.append((cp(m6), cp(mo)))
            for name, m in (("LEC6", m6), ("other", mo)):
                clean = m & (fib_counts < 1)
                checks.append(dict(donor=d, assay=a, pool=name, n_cells=int(m.sum()),
                                   frac_fib_ge2=float((fib_counts[m] >= 2).mean()),
                                   **{f"{g}_cp10k": v for g, v in zip(FIB[:3] + STRESS, cp(m, FIB[:3] + STRESS))},
                                   n_fib_free=int(clean.sum()),
                                   antxr2_cp10k_fib_free=float(a2[clean].sum() / lib[clean].sum() * 1e4),
                                   antxr2_det_fib_free=float((a2[clean] > 0).mean())))
    F = np.vstack([np.log2((u6 + PSEUDO) / (uo + PSEUDO)) for u6, uo in units])
    P6, PO = np.median([u[0] for u in units], 0), np.median([u[1] for u in units], 0)
    res = pd.DataFrame({"gene": genes, "log2fc": np.median(F, 0), "n_up": (F > 0).sum(0),
                        "n_units": len(units), "cp10k_LEC6": P6, "cp10k_other": PO})
    res = res[np.maximum(P6, PO) >= EXPR_FLOOR].drop_duplicates("gene")
    res["pct"] = res.log2fc.rank(pct=True) * 100
    res = res.sort_values("log2fc", ascending=False)
    res.to_csv(os.path.join(lym, "LEC6_vs_otherLEC_paired.csv"), index=False)
    pd.DataFrame(checks).to_csv(os.path.join(lym, "LEC6_contamination_check.csv"), index=False)

    pd.set_option("display.width", 250)
    up_all = res[res.n_up == len(units)].reset_index(drop=True)
    print(f"{len(units)} paired donor x chemistry units; {len(res)} genes ranked")
    print(res[res.gene == "ANTXR2"].round(3).to_string(index=False))
    print("ANTXR2 rank among genes up in every unit:", int(up_all.index[up_all.gene == "ANTXR2"][0]) + 1)
    print(up_all.head(30)[["gene", "log2fc", "cp10k_LEC6", "cp10k_other"]].round(2).to_string(index=False))
    print(pd.DataFrame(checks).round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
