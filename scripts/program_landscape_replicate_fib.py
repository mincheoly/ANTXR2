"""Cross-dataset replication, step 2: the IBDverse fibroblast co-expression
profile in GCA fibroblasts (planning_summaries/program_coexpression.md,
"2. Fibroblast co-expression: plan").

Reads the IBDverse fibroblast consensus (OUT_DIR/landscape/report) and the
GCA fibroblast cell (OUT_DIR/landscape_gca_cov50, after
ANTXR2_LANDSCAPE=gca ANTXR2_CELL_COVERAGE=0.5 program_landscape.py all
--pops fibroblast_pooled).

    F1  the IFN set: (a) mean GCA excess z vs N_DRAWS random same-size sets
        of GCA-scored library programs, one-sided greater; (b) ANTXR2's mean
        rho over the set ranked among its 50 anchors'
    F2  Spearman of excess z, IBDverse consensus mean z vs GCA, shared programs
    F3  the other IBDverse rank >= 50 programs, one by one

Writes <GCA run>/replication_fibroblast.json and .csv (the shared programs).

Usage:
    python program_landscape_replicate_fib.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import LANDSCAPE_FIB, OUT_DIR, SEED, pop_slug  # noqa: E402

IBD, GCA = OUT_DIR / "landscape", OUT_DIR / "landscape_gca_cov50"
CELL = f"{pop_slug(LANDSCAPE_FIB)}@normal"
N_DRAWS = 10_000
IFN = ["collectri:IRF1", "cytosig:IFN1", "cytosig:IFNG", "hallmark:INTERFERON_GAMMA_RESPONSE",
       "hallmark:INTERFERON_ALPHA_RESPONSE", "collectri:STAT1"]
OTHER = ["cytosig:IL27", "matrisome:NABA_ECM_REGULATORS", "cytosig:FGF2", "collectri:ZNF91",
         "collectri:FOXO4"]
FIXED = ("discovery:", "bracq:")


def set_rank(units):
    """Percentile of ANTXR2 among its anchors on mean rho over `units` (PC1)."""
    z = np.load(GCA / "boot" / f"{pop_slug(LANDSCAPE_FIB)}.npz")
    k = [i for i, (u, s) in enumerate(zip(z["units"], z["schemes"])) if s == "pc1" and u in units]
    m = np.nanmean(z["normal_point"][k][:, :1 + int(z["n_pop_anchors"])], 0)
    return float((m[1:] < m[0]).mean() * 100), float(m[0]), float(m[1:].mean())


def main():
    ibd = pd.read_csv(IBD / "report" / "consensus_fibroblast.csv", index_col=0)
    g = pd.read_csv(GCA / "landscape.csv")
    g = g[(g.cell == CELL) & (g.scheme == "pc1")].set_index("unit")
    lib = g[~g.index.str.startswith(FIXED)]
    res = dict(n_programs_gca=len(lib))

    S = [u for u in IFN if u in lib.index]
    if len(S) < 3:
        res["F1"] = dict(verdict="untestable", scored=S)
    else:
        obs = lib.excess_z[S].mean()
        rng = np.random.default_rng(SEED)
        z = lib.excess_z.dropna().values
        null = np.array([rng.choice(z, len(S), replace=False).mean() for _ in range(N_DRAWS)])
        pct, rho, arho = set_rank(S)
        res["F1"] = dict(scored=S, per_program=lib.loc[S, ["rho", "anchor_mean", "excess",
                                                           "excess_z", "coverage"]].round(3)
                         .to_dict("index"),
                         mean_z=float(obs), p=float((1 + (null >= obs).sum()) / (N_DRAWS + 1)),
                         set_rho=rho, set_anchor_rho=arho, anchor_pct=pct)
        a, b = res["F1"]["p"] < 0.05, pct >= 90
        res["F1"]["verdict"] = "replicated" if a and b else "partial" if a or b else "not replicated"

    both = ibd.join(lib[["excess_z", "coverage", "rho"]], how="inner")
    r, p = spearmanr(both.mean_z, both.excess_z)
    hi = both[both.coverage >= 0.8]
    res["F2"] = dict(n=len(both), spearman=float(r), p=float(p), n_cov80=len(hi),
                     spearman_cov80=float(spearmanr(hi.mean_z, hi.excess_z)[0]))
    res["F3"] = {u: (round(float(lib.excess_z[u]), 2) if u in lib.index else None) for u in OTHER}
    both.to_csv(GCA / "replication_fibroblast.csv")
    (GCA / "replication_fibroblast.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    pd.set_option("display.width", 200)
    s = both.sort_values("excess_z")[["mean_z", "consensus_z", "excess_z", "coverage"]].round(2)
    print(pd.concat([s.head(10), s.tail(10)]).to_string())


if __name__ == "__main__":
    main()
