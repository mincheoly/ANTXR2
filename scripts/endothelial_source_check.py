"""Which gene groups drive ANTXR2's endothelial ECM signal, before/after decontamination.

For each endothelial population scan (coexpr_point_scan.py output: ANTXR2 +
matched anchors vs all genes), each target's correlation is residualized on
expression exactly as in coexpr_pathway_gsea.py (running median over genes
sorted by raw mean). The mean residual over a gene group is compared between
ANTXR2 and its anchors: ANTXR2 is the max of 6 with probability 1/6 under
exchangeability, tested by one-sided binomial across populations.

Each contaminant (fibroblast, lymphatic, macrophage) is read on HELD-OUT
genes, never on the genes endothelial_decontam.py filters on, whose
correlations the filter truncates mechanically.

Usage: python endothelial_source_check.py <out.csv> <label>=<scan_dir> ...
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

WINDOW = 301
GROUPS = {
    "EC_scavenger": ["STAB1", "STAB2", "SCARF1", "MRC1", "LYVE1", "CD36"],
    "basement_membrane": ["COL4A1", "COL4A2", "COL15A1", "COL18A1", "LAMB1", "LAMA4", "LAMC1",
                          "HSPG2", "NID1", "NID2", "SPARC"],
    "fibroblast_heldout": ["COL3A1", "PDGFRA", "FBLN1", "MFAP5", "C7", "SFRP2", "MMP2", "PDGFRL"],
    "pericyte_SMC": ["PDGFRB", "RGS5", "ACTA2", "MYH11", "NOTCH3"],
    # held out of endothelial_decontam.py's lymphatic (PROX1, CCL21) and
    # macrophage (C1QA, C1QB, CD163, TYROBP) filters
    "lymphatic_heldout": ["TFF3", "RELN", "FLT4", "MMRN1", "PKHD1L1"],
    "macrophage_heldout": ["AIF1", "CD68", "FCER1G", "PTPRC", "LYZ"],
}


def residualize(d):
    out = []
    for _, g in d.dropna(subset=["corr"]).groupby("target"):
        g = g.sort_values("raw_mean")
        g = g.assign(resid=g["corr"] - g["corr"].rolling(WINDOW, center=True,
                                                         min_periods=WINDOW // 3).median())
        out.append(g)
    return pd.concat(out)


def per_population(scan_dir, label):
    rows = []
    for f in sorted(glob.glob(os.path.join(scan_dir, "endothelial_*.csv"))):
        if os.path.getsize(f) < 10:  # no donor x cell_type group cleared min_cell_count
            continue
        d = residualize(pd.read_csv(f))
        for grp, genes in GROUPS.items():
            x = d[d.gene.isin(genes)].groupby("target").resid.mean()
            if "ANTXR2" not in x or len(x) < 6:
                continue
            anc = x.drop("ANTXR2")
            rows.append(dict(condition=label, population=os.path.basename(f)[:-4], group=grp,
                             antxr2=x["ANTXR2"], anchor_median=anc.median(),
                             antxr2_max=bool((x["ANTXR2"] > anc).all())))
    return pd.DataFrame(rows)


def main():
    out = sys.argv[1]
    r = pd.concat([per_population(d, lab) for lab, d in (a.split("=", 1) for a in sys.argv[2:])])
    r.to_csv(out, index=False)
    s = r.groupby(["condition", "group"]).agg(n=("population", "size"), antxr2_med=("antxr2", "median"),
                                              anchor_med=("anchor_median", "median"),
                                              n_max=("antxr2_max", "sum")).reset_index()
    s["binom_p"] = [binomtest(int(k), int(n), 1 / 6, alternative="greater").pvalue
                    for k, n in zip(s.n_max, s.n)]
    s.to_csv(out.replace(".csv", "_summary.csv"), index=False)
    pd.set_option("display.width", 200)
    print(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
