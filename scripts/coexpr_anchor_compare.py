"""ANTXR2 vs its matched anchors: is a recurrent pathway specific to ANTXR2?

Two modes.

  split   Explode coexpr_point_scan.py output (one CSV per population holding
          ANTXR2 + anchors) into one ranking CSV per population x target,
          named <population>__<target>.csv, in the column layout
          coexpr_pathway_gsea.py reads. ANTXR2 is re-ranked from the same
          point-scan run as its anchors so both go through an identical
          pipeline.

  compare After coexpr_pathway_gsea.py has run on the split files. Per
          population x term, the rank of ANTXR2's NES among the 6 targets
          (ANTXR2 + 5 anchors). Under the null that ANTXR2 is exchangeable
          with its anchors, ANTXR2 is the maximum of the 6 with probability
          1/6 and the minimum with probability 1/6. The direction is NOT
          chosen from ANTXR2's own sign: "most extreme in ANTXR2's direction"
          has null probability ~1/3, and testing it against 1/6 inflated an
          earlier version of this table. Per compartment x term:
            n_antxr2_sig     populations with ANTXR2 FDR<0.05
            mean_anchor_sig  mean fraction of anchors with FDR<0.05 in the
                             same direction as ANTXR2 (how generic the term is)
            n_max / n_min    populations where ANTXR2 is the max / min of 6
            top_p            2 x min(binomial P(n_max | 1/6), P(n_min | 1/6)),
                             one-sided each, capped at 1; dir = the tail used

          --standardize divides each target's NES by its SD across all terms
          in that population before ranking. ANTXR2's NES spread is wider
          than its anchors' in immune and epithelial populations (widest of 6
          in 48% / 38% vs 1/6 expected), which lets it "win" many terms in
          both tails without any term-specific signal; standardizing asks
          whether a term is extreme *for that gene*.

  single  One population with many anchors (e.g. 30). Per term: z of
          ANTXR2's NES against the anchors' NES (mean, SD), two-sided normal
          p, BH across terms, plus ANTXR2's empirical rank among the anchors
          and the fraction of anchors that are themselves FDR<0.05 in
          ANTXR2's direction. With few populations the rank test above has
          no resolution (minimum p = 1/(n_anchors+1)), hence the z.

Usage:
    python coexpr_anchor_compare.py split <point_scan_dir> <ranks_dir>
    python coexpr_anchor_compare.py single <gsea_dir> <out.csv> [--stat resid]
    python coexpr_anchor_compare.py compare <gsea_dir> [--stat resid] [--standardize]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import binomtest

FDR = 0.05


def split(scan_dir, ranks_dir):
    os.makedirs(ranks_dir, exist_ok=True)
    n = 0
    for f in sorted(glob.glob(os.path.join(scan_dir, "*.csv"))):
        if os.path.getsize(f) < 10:  # no donor x cell_type group cleared min_cell_count
            continue
        pop = os.path.basename(f)[:-4]
        df = pd.read_csv(f)
        for t, g in df.groupby("target"):
            g.rename(columns={"corr": "corr_coef"})[["gene", "corr_coef", "raw_mean", "det_rate"]] \
                .to_csv(os.path.join(ranks_dir, f"{pop}__{t}.csv"), index=False)
            n += 1
    print(f"wrote {n} ranking files")


def compare(gsea_dir, stat, standardize=False):
    fs = glob.glob(os.path.join(gsea_dir, stat, "*__*.csv"))
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    if standardize:
        df["nes"] = df.nes / df.groupby("population").nes.transform("std")
    df[["pop", "target"]] = df.population.str.split("__", n=1, expand=True)
    df["compartment"] = df["pop"].str.split("_").str[0]
    a = df[df.target == "ANTXR2"].set_index(["pop", "term"])
    anc = df[df.target != "ANTXR2"]
    rows = []
    for (pop, term), g in anc.groupby(["pop", "term"]):
        if (pop, term) not in a.index:
            continue
        x = a.loc[(pop, term)]
        s = np.sign(x.nes)
        rows.append(dict(pop=pop, compartment=x.compartment, term=term, antxr2_nes=x.nes,
                         antxr2_fdr=x.fdr, n_anchors=len(g),
                         anchor_median_nes=g.nes.median(),
                         anchor_sig_same_dir=((g.fdr < FDR) & (np.sign(g.nes) == s)).mean(),
                         antxr2_max=bool((x.nes > g.nes).all()),
                         antxr2_min=bool((x.nes < g.nes).all())))
    r = pd.DataFrame(rows)
    tag = f"{stat}_std" if standardize else stat
    r.to_csv(os.path.join(gsea_dir, f"anchor_compare_by_population_{tag}.csv"), index=False)

    out = []
    for (comp, term), g in r.groupby(["compartment", "term"]):
        sig = g[g.antxr2_fdr < FDR]
        n_pop = g["pop"].nunique()
        n_max, n_min = int(g.antxr2_max.sum()), int(g.antxr2_min.sum())
        p_max = binomtest(n_max, n_pop, 1 / 6, alternative="greater").pvalue
        p_min = binomtest(n_min, n_pop, 1 / 6, alternative="greater").pvalue
        out.append(dict(compartment=comp, term=term, n_pops=n_pop,
                        median_antxr2_nes=g.antxr2_nes.median(),
                        median_anchor_nes=g.anchor_median_nes.median(),
                        n_antxr2_sig=len(sig),
                        mean_anchor_sig=g.anchor_sig_same_dir.mean(),
                        n_max=n_max, n_min=n_min,
                        dir="+" if p_max <= p_min else "-",
                        top_p=min(1.0, 2 * min(p_max, p_min))))
    o = pd.DataFrame(out)
    o["top_fdr"] = o.groupby("compartment").top_p.transform(_bh)
    o = o.sort_values(["compartment", "top_p"])
    o.to_csv(os.path.join(gsea_dir, f"anchor_compare_{tag}.csv"), index=False)
    return r, o


def single(gsea_dir, out, stat):
    from scipy.stats import norm
    fs = glob.glob(os.path.join(gsea_dir, stat, "*__*.csv"))
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df["target"] = df.population.str.split("__", n=1).str[1]
    a = df[df.target == "ANTXR2"].set_index("term")
    anc = df[df.target != "ANTXR2"]
    g = anc.groupby("term")
    r = pd.DataFrame({"antxr2_nes": a.nes, "antxr2_fdr": a.fdr,
                      "anchor_mean": g.nes.mean(), "anchor_sd": g.nes.std(),
                      "n_anchors": g.nes.size()}).dropna(subset=["antxr2_nes"])
    r["z"] = (r.antxr2_nes - r.anchor_mean) / r.anchor_sd
    r["p"] = 2 * norm.sf(r.z.abs())
    r["fdr"] = _bh(r.p.values)
    r["rank_among_anchors"] = [int((anc[anc.term == t].nes < r.at[t, "antxr2_nes"]).sum()) for t in r.index]
    sig_same = anc.assign(s=np.sign(anc.term.map(r.antxr2_nes)) == np.sign(anc.nes),
                          f=anc.fdr < FDR).groupby("term").apply(lambda x: (x.s & x.f).mean())
    r["anchor_frac_sig_same_dir"] = sig_same
    r = r.sort_values("p")
    r.to_csv(out)
    return r


def _bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    s = sub.add_parser("split")
    s.add_argument("scan_dir")
    s.add_argument("ranks_dir")
    sg = sub.add_parser("single")
    sg.add_argument("gsea_dir")
    sg.add_argument("out")
    sg.add_argument("--stat", default="resid")
    c = sub.add_parser("compare")
    c.add_argument("gsea_dir")
    c.add_argument("--stat", default="resid")
    c.add_argument("--standardize", action="store_true")
    a = ap.parse_args()
    if a.mode == "split":
        split(a.scan_dir, a.ranks_dir)
    elif a.mode == "single":
        r = single(a.gsea_dir, a.out, a.stat)
        pd.set_option("display.width", 250)
        pd.set_option("display.max_colwidth", 75)
        print(f"{(r.fdr < 0.05).sum()} of {len(r)} terms FDR<0.05 vs anchors")
        print(r.head(25).round(3).to_string())
    else:
        r, o = compare(a.gsea_dir, a.stat, a.standardize)
        pd.set_option("display.width", 250)
        pd.set_option("display.max_colwidth", 70)
        for comp, g in o.groupby("compartment"):
            print(f"\n== {comp}: terms where ANTXR2 beats its anchors most often")
            print(g.head(12).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
