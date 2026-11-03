"""Pre-registered replication tests for the enterocyte-leaning ANTXR2 program
difference (planning_summaries/program_coexpression.md, "Replication
pre-registration"), from two program_scores.py runs.

    1. discovery:IEG_core   excess < 0 with donor-bootstrap CI excluding 0
    2. DISCOVERY_SET        mean excess z, one-sided vs 10,000 random same-size
                            sets of replication-comparable library programs
    3. DISCOVERY_SET        how many have excess < 0
    4. whole profile        Spearman of excess z, programs comparable in both

Verdict: supported if 1 passes and 2 has p < 0.05; partial if one does; not
replicated otherwise. PC1 weights throughout, as in discovery.

Writes <replication run>/replication.json and replication.png.

Usage:
    python program_replicate.py <discovery run> <replication run>
"""
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
from scipy.stats import spearmanr                   # noqa: E402

from program_config import DISCOVERY_SET, OUT_DIR, SEED  # noqa: E402
from program_report import BLUE, GRID, INK2, MUTED, ORANGE, SURFACE, short  # noqa: E402

N_DRAWS = 10_000
CORE = "discovery:IEG_core"


def load(run):
    d = pd.read_csv(OUT_DIR / run / "differential.csv")
    return d[d.scheme == "pc1"].set_index("unit")


def library(d):
    """Screened library programs only (no fixed discovery / hypothesis modules)."""
    return d[~d.index.str.startswith(("discovery:", "bracq:"))]


def set_vs_random(lib, S, side):
    """Mean excess z of set S vs N_DRAWS random same-size sets of lib's
    programs; one-sided p in the direction `side` ("less" or "greater")."""
    obs = lib.loc[S, "excess_z"].mean()
    rng = np.random.default_rng(SEED)
    z = lib.excess_z.dropna().values
    null = np.array([rng.choice(z, len(S), replace=False).mean() for _ in range(N_DRAWS)])
    hits = (null <= obs) if side == "less" else (null >= obs)
    return obs, float((1 + hits.sum()) / (N_DRAWS + 1)), null


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("discovery")
    ap.add_argument("replication")
    a = ap.parse_args()
    disc, rep = load(a.discovery), load(a.replication)
    res = dict(discovery=a.discovery, replication=a.replication)

    # 1. the pre-specified module
    if CORE in rep.index:
        r = rep.loc[CORE]
        res["test1"] = dict(excess=r.excess, ci=[r.excess_ci_lo, r.excess_ci_hi],
                            excess_z=r.excess_z, rho_fib=r.rho_a, rho_ent=r.rho_b,
                            n_genes=int(r.n_genes),
                            passed=bool(r.excess < 0 and r.excess_ci_hi < 0))
    else:
        res["test1"] = dict(passed=False, note="IEG_core not comparable in replication")

    # 2 + 3. the discovery set against random same-size sets
    lib = library(rep)
    S = [u for u in DISCOVERY_SET if u in lib.index]
    obs, p, null = set_vs_random(lib, S, "less")
    res["test2"] = dict(n_in_replication=len(S), missing=sorted(set(DISCOVERY_SET) - set(S)),
                        mean_excess_z=obs, p_one_sided=p,
                        null_mean=float(null.mean()), null_sd=float(null.std()))
    res["test3"] = dict(n_negative=int((lib.loc[S, "excess"] < 0).sum()), n=len(S),
                        per_program={u: round(float(lib.loc[u, "excess_z"]), 2) for u in S})

    # 4. whole profile
    sh = disc.index.intersection(lib.index)
    rho, p = spearmanr(disc.loc[sh, "excess_z"], lib.loc[sh, "excess_z"])
    res["test4"] = dict(n_shared=len(sh), spearman=float(rho), p=float(p))

    t1, t2 = res["test1"]["passed"], res["test2"]["p_one_sided"] < 0.05
    res["verdict"] = "supported" if t1 and t2 else "partial" if t1 or t2 else "not replicated"
    out = OUT_DIR / a.replication
    json.dump(res, open(out / "replication.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=lambda x: round(float(x), 4)))

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    x, y = disc.loc[sh, "excess_z"], lib.loc[sh, "excess_z"]
    ax.axhline(0, color=GRID, lw=1)
    ax.axvline(0, color=GRID, lw=1)
    ax.scatter(x, y, s=10, color="#b5b4ae", lw=0, label=f"{len(sh)} shared programs")
    inS = [u for u in S if u in sh]
    ax.scatter(disc.loc[inS, "excess_z"], lib.loc[inS, "excess_z"], s=30, color=ORANGE,
               edgecolor=SURFACE, lw=1.2, zorder=3, label="discovery set S")
    for u in inS:
        ax.annotate(short(u), (disc.loc[u, "excess_z"], lib.loc[u, "excess_z"]), fontsize=6,
                    color=INK2, xytext=(4, 2), textcoords="offset points")
    if CORE in rep.index and CORE in disc.index:
        ax.scatter([disc.loc[CORE, "excess_z"]], [rep.loc[CORE, "excess_z"]], s=60, marker="D",
                   color=BLUE, edgecolor=SURFACE, lw=1.2, zorder=4, label="IEG_core module")
    ax.set_xlabel(f"excess z, discovery ({a.discovery})")
    ax.set_ylabel(f"excess z, replication ({a.replication})")
    ax.set_title(f"verdict: {res['verdict']}   Spearman {rho:.2f}", fontsize=8, color=INK2, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.12))
    fig.tight_layout()
    fig.savefig(out / "replication.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
