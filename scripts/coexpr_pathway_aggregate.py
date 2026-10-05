"""Aggregate per-population GSEA (coexpr_pathway_gsea.py) into recurrence tables.

Two axes, both treating populations as replicates (caveat: Tabula Sapiens
donors overlap across tissues, so tissues are not fully independent replicates,
and compartments within a tissue share donors by construction):

  across_tissue   per compartment x term: in how many tissues is the term
                  enriched (FDR < FDR) positive / negative; median NES;
                  two-sided sign test on NES signs across tissues.
  within_tissue   per tissue x term: how many of that tissue's compartments
                  enrich the term in the same direction. A term counts as
                  "tissue-shared" if it is significant in >= MIN_COMP
                  compartments of a tissue, all in the same direction; the
                  summary then counts in how many tissues that happens.

Usage: python coexpr_pathway_aggregate.py <gsea_dir> [--stat resid] [--fdr 0.05]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import binomtest

MIN_COMP = 3


def load(gsea_dir, stat):
    fs = glob.glob(os.path.join(gsea_dir, stat, "*.csv"))
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df["compartment"] = df.population.str.split("_").str[0]
    df["tissue"] = df.population.str.split("_", n=1).str[1].str.replace(
        r"^(ts|oral|syn|tach|tquad)_", "", regex=True)
    return df


def across_tissue(df, fdr):
    rows = []
    n_pop = df.groupby("compartment").population.nunique()
    for (comp, term), g in df.groupby(["compartment", "term"]):
        pos, neg = (g.nes > 0).sum(), (g.nes < 0).sum()
        sp, sn = ((g.fdr < fdr) & (g.nes > 0)).sum(), ((g.fdr < fdr) & (g.nes < 0)).sum()
        rows.append(dict(compartment=comp, term=term, n_tested=len(g), n_pops=n_pop[comp],
                         median_nes=g.nes.median(), n_pos=pos, n_neg=neg,
                         n_sig_pos=sp, n_sig_neg=sn,
                         sign_p=binomtest(int(pos), int(pos + neg)).pvalue if pos + neg else np.nan))
    out = pd.DataFrame(rows)
    out["frac_sig_same_dir"] = np.where(out.median_nes > 0, out.n_sig_pos, out.n_sig_neg) / out.n_pops
    out["sign_fdr"] = out.groupby("compartment").sign_p.transform(_bh)
    return out.sort_values(["compartment", "frac_sig_same_dir"], ascending=[True, False])


def within_tissue(df, fdr):
    rows = []
    for (tis, term), g in df.groupby(["tissue", "term"]):
        sp = g[(g.fdr < fdr) & (g.nes > 0)].compartment.tolist()
        sn = g[(g.fdr < fdr) & (g.nes < 0)].compartment.tolist()
        rows.append(dict(tissue=tis, term=term, n_comp=g.compartment.nunique(),
                         n_sig_pos=len(sp), n_sig_neg=len(sn),
                         comps_pos=",".join(sorted(sp)), comps_neg=",".join(sorted(sn)),
                         median_nes=g.nes.median()))
    w = pd.DataFrame(rows)
    w["shared_dir"] = np.select([(w.n_sig_pos >= MIN_COMP) & (w.n_sig_neg == 0),
                                 (w.n_sig_neg >= MIN_COMP) & (w.n_sig_pos == 0)],
                                ["pos", "neg"], "")
    n_tis = df.groupby("tissue").compartment.nunique()
    eligible = n_tis[n_tis >= MIN_COMP].index
    s = (w[w.tissue.isin(eligible) & (w.shared_dir != "")]
         .groupby(["term", "shared_dir"])
         .agg(n_tissues=("tissue", "nunique"), tissues=("tissue", lambda x: ",".join(sorted(x))))
         .reset_index().sort_values("n_tissues", ascending=False))
    s["n_eligible_tissues"] = len(eligible)
    return w, s


def _bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    r = p[o] * len(p) / (np.arange(len(p)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty_like(r)
    out[o] = np.minimum(r, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gsea_dir")
    ap.add_argument("--stat", default="resid")
    ap.add_argument("--fdr", type=float, default=0.05)
    a = ap.parse_args()
    df = load(a.gsea_dir, a.stat)
    at = across_tissue(df, a.fdr)
    w, s = within_tissue(df, a.fdr)
    at.to_csv(os.path.join(a.gsea_dir, f"across_tissue_{a.stat}.csv"), index=False)
    w.to_csv(os.path.join(a.gsea_dir, f"within_tissue_{a.stat}.csv"), index=False)
    s.to_csv(os.path.join(a.gsea_dir, f"tissue_shared_{a.stat}.csv"), index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    for comp, g in at.groupby("compartment"):
        print(f"\n== {comp} ({g.n_pops.iloc[0]} populations): top recurrent terms")
        print(g.head(15)[["term", "median_nes", "n_sig_pos", "n_sig_neg", "frac_sig_same_dir",
                          "sign_fdr"]].round(3).to_string(index=False))
    print(f"\n== terms shared across >= {MIN_COMP} compartments within a tissue")
    print(s.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
