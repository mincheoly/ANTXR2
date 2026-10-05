"""How much did dropping |r|==1 group estimates change each scan's ranking?

coexpr_point_scan.py versions before 2026-10-04 averaged memento's +/-1
placeholder group estimates (pairs whose variance is <= 0 in a group) into
the cell-weighted mean. Every downstream analysis consumed the
expression-adjusted ranking (residual from a running median over genes
sorted by raw mean; coexpr_pathway_gsea.ranking / endothelial_source_check),
so that is what is compared here, per population x target, old vs fixed.

Usage: python compare_fixed_scans.py <out.csv> <old_dir>=<fixed_dir> ...
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endothelial_source_check import residualize   # noqa: E402


def main(out, pairs):
    rows = []
    for spec in pairs:
        old_dir, new_dir = spec.split("=", 1)
        for f_new in sorted(glob.glob(os.path.join(new_dir, "*.csv"))):
            f_old = os.path.join(old_dir, os.path.basename(f_new))
            if not os.path.exists(f_old) or os.path.getsize(f_new) < 10 or os.path.getsize(f_old) < 10:
                continue
            o = residualize(pd.read_csv(f_old)).set_index(["target", "gene"])
            n = residualize(pd.read_csv(f_new)).set_index(["target", "gene"])
            raw_o = pd.read_csv(f_old).set_index(["target", "gene"])["corr"]
            raw_n = pd.read_csv(f_new).set_index(["target", "gene"])["corr"]
            for t in n.index.get_level_values(0).unique():
                a, b = o.loc[t, "resid"], n.loc[t, "resid"]
                j = a.index.intersection(b.index)
                m = a.loc[j].notna() & b.loc[j].notna()
                rows.append(dict(scan_set=os.path.basename(os.path.normpath(new_dir)),
                                 population=os.path.basename(f_new)[:-4], target=t,
                                 is_antxr2=t == "ANTXR2", n_genes=int(m.sum()),
                                 rho_adjusted=spearmanr(a.loc[j][m], b.loc[j][m])[0],
                                 median_raw_old=raw_o.loc[t].median(),
                                 median_raw_new=raw_n.loc[t].median()))
    r = pd.DataFrame(rows)
    r.to_csv(out, index=False)
    pd.set_option("display.width", 220)
    s = r.groupby(["scan_set", "is_antxr2"]).rho_adjusted.describe(percentiles=[.05, .5])
    print(s[["count", "min", "5%", "50%"]].round(3).to_string())
    bad = r[r.rho_adjusted < 0.95].sort_values("rho_adjusted")
    print(f"\npopulation x target rankings with rho < 0.95: {len(bad)} of {len(r)}")
    print(bad.head(30).round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
