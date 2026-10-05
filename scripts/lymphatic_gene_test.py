"""Per-gene test: is ANTXR2 co-expressed with gene g at all, within lymphatic EC?

Earlier tables compared ANTXR2's number of hits with control genes' numbers
of hits. That was the wrong null: control genes have real co-expression
partners of their own, so their hit counts are biology, not false positives,
and "is ANTXR2 more connected than a random gene" is not a prior worth
testing. This asks the per-gene question directly.

For each partner gene g:
  1. Per donor x subtype/tissue group (memento point estimates from
     coexpr_point_scan.py --groups-out), the ANTXR2-g correlation is
     expression-adjusted within that group: minus the running median
     (WINDOW genes) of ANTXR2's correlations with genes of similar raw mean
     in that same group. This removes the estimator's global positive shift
     and its detection-dependent trend, so "no relationship" sits near 0.
  2. Adjusted values are averaged within donor (donors contributing two
     subtypes count once): the donor is the replicate.
  3. One-sample t-test across donors (>= MIN_DONORS), BH across genes.

Output columns also give the per-source mean (gut groups hold lymphatic
subtype fixed; ts / site groups are donor x tissue and can mix subtypes) and
the LEC6-vs-other fold change, so composition-driven genes can be spotted.

Usage: python lymphatic_gene_test.py <lymphatic dir> [groups_subdir] [out_suffix]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphatic_gene_table import annotate          # noqa: E402

WINDOW = 301
MIN_DONORS = 8


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    o = np.argsort(p[ok])
    r = p[ok][o] * ok.sum() / (np.arange(ok.sum()) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    v = np.empty_like(r)
    v[o] = np.minimum(r, 1)
    q[ok] = v
    return q


def adjusted(corr, rawmean):
    out = {}
    for grp in corr.columns:
        c = corr[grp].dropna()
        m = rawmean[grp].reindex(c.index)
        c = c[m.sort_values().index]
        out[grp] = c - c.rolling(WINDOW, center=True, min_periods=WINDOW // 3).median()
    return pd.DataFrame(out)


def main(lym, groups_dir="scan_groups", suffix=""):
    pre = os.path.join(lym, groups_dir, "antxr2")
    corr = pd.read_parquet(pre + "_corr.parquet").set_index("gene").drop(columns="target")
    corr = corr.where(corr.abs() < 1)   # memento +/-1 placeholder groups
    rawmean = pd.read_parquet(pre + "_rawmean.parquet")
    adj = adjusted(corr, rawmean)
    donor = pd.Series({g: g.split("^", 1)[0] for g in adj.columns})
    source = donor.str.split(":").str[0]
    by_donor = adj.T.groupby(donor).mean().T
    n = by_donor.notna().sum(axis=1)
    t = ttest_1samp(by_donor, 0, axis=1, nan_policy="omit")
    res = pd.DataFrame({"mean_adj_r": by_donor.mean(axis=1), "n_donors": n,
                        "frac_donors_pos": (by_donor > 0).sum(axis=1) / n,
                        "t": t.statistic, "p": t.pvalue})
    res.loc[res.n_donors < MIN_DONORS, ["t", "p"]] = np.nan
    res["fdr"] = bh(res.p)
    for s in sorted(source.unique()):
        res[f"mean_{s}"] = adj.loc[:, source[source == s].index].mean(axis=1)
    lec6 = pd.read_csv(os.path.join(lym, "LEC6_vs_otherLEC_paired.csv")).set_index("gene")
    res["lec6_log2fc"] = lec6.log2fc.reindex(res.index)
    res = res.sort_values("p")
    res.to_csv(os.path.join(lym, f"antxr2_lymphatic_gene_test{suffix}.csv"))

    tested = res.p.notna()
    print(f"groups {adj.shape[1]}, donors {by_donor.shape[1]}, genes tested {tested.sum()}")
    print("p-value histogram (deciles):",
          np.histogram(res.p[tested], bins=np.linspace(0, 1, 11))[0].tolist())
    for th in (0.05, 0.1, 0.25):
        print(f"FDR<{th}: {(res.fdr < th).sum()}")
    hits = res[res.fdr < 0.1]
    if len(hits):
        ann = annotate(list(hits.index), os.path.join(lym, "mygene_cache.json"))
        hits = hits.assign(name=[ann.get(g, {}).get("name", "") for g in hits.index])
        pd.set_option("display.width", 260)
        pd.set_option("display.max_colwidth", 55)
        cols = ["mean_adj_r", "n_donors", "frac_donors_pos", "p", "fdr",
                "mean_gut", "mean_site", "mean_ts", "lec6_log2fc", "name"]
        print(hits[cols].round(4).to_string())
        hits.to_csv(os.path.join(lym, f"antxr2_lymphatic_gene_hits{suffix}.csv"))


if __name__ == "__main__":
    main(*sys.argv[1:])
