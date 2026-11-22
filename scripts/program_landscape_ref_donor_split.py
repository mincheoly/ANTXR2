"""Donor split-half check of the reference-scaled coupling profile inside one
population (default: fibroblasts). How reproducible is ANTXR2's program
percentile profile between two random halves of the donors of the SAME
dataset? That is the ceiling for any cross-dataset agreement.

Per population x state cell: the group (donor x subtype) sufficient statistics
are computed once with the main run's anchors and weights; each of N_SPLITS
random donor halvings then re-derives beta, and ANTXR2's percentile among its
anchors, in each half. r = Spearman between halves over the programs active in
the full cell (PC1 scheme, activity_sd >= 0.2). Spearman-Brown 2r / (1 + r)
estimates the reliability of the full-donor profile.

Writes LANDSCAPE_OUT/donor_split/{splits,summary}.csv and
{cell}_example.{csv,png} (the split with the median r).

Usage:
    ANTXR2_LANDSCAPE_RUN=ref [ANTXR2_LANDSCAPE=gca] python program_landscape_ref_donor_split.py \
        [--pops fibroblast_pooled]
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402
import scipy.sparse as sp         # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import LANDSCAPE_ATLAS, LANDSCAPE_OUT, SEED, TARGET, pop_slug  # noqa: E402
import program_landscape as PL                     # noqa: E402
import program_landscape_ref as REF                # noqa: E402
from program_scores import group_stats             # noqa: E402

OUT = LANDSCAPE_OUT / "donor_split"
OUT.mkdir(parents=True, exist_ok=True)
N_SPLITS = 200
MIN_ACTIVITY_SD = 0.2
SCHEME = "pc1"


def percentiles(S, g, sigma, member):
    """ANTXR2's anchor percentile per program from the groups in `g`."""
    n, s_sum, t_sum, t_sq, t_noise, ST = (a[g] for a in S)
    N = n.sum()
    vt = (t_sq - t_sum ** 2 / n[:, None] - t_noise).sum(0) / N
    cov = (ST - s_sum[:, :, None] * t_sum[:, None, :] / n[:, None, None]).sum(0) / N
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = np.where(member, np.nan, cov / (np.sqrt(np.where(vt > 0, vt, np.nan))[None, :]
                                               * sigma[:, None]))
        return (np.nansum(beta[:, 1:] < beta[:, :1], 1) / np.isfinite(beta[:, 1:]).sum(1) * 100,
                beta[:, 0])


def run(pop, rng):
    cells = REF.filled_cells()
    fc = cells[cells.population == pop]
    progs = REF.load_programs()
    units = list(progs)
    sigma = np.array([progs[u][f"sigma_{SCHEME}"] for u in units])
    anchors = json.load(open(LANDSCAPE_OUT / "score_meta.json"))[pop]["anchors"]
    targets = [TARGET] + anchors
    X, Y, sf, obs, genes, groups = PL.load_population(pop)
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    T = gpos[targets].values
    member = np.array([[g in set(progs[u]["genes_all"]) for g in targets] for u in units])
    full = pd.read_csv(LANDSCAPE_OUT / "coupling.csv")
    full = full[(full.scheme == SCHEME) & (full.population == pop)]
    splits, summ = [], []
    for c, state, scale in fc[["cell", "state", "scale"]].itertuples(index=False):
        gc = groups[groups.cell_type == c].reset_index(drop=True)
        rr = np.concatenate([gpos[progs[u]["genes"]].values for u in units])
        cc = np.concatenate([np.full(len(progs[u]["genes"]), k) for k, u in enumerate(units)])
        vv = np.concatenate([progs[u][SCHEME] / (progs[u]["s_ref"] * scale) for u in units])
        Wm = sp.csr_matrix((vv, (rr, cc)), shape=(len(genes), len(units)))
        (n, s_sum, _, _), (t_sum, t_sq, t_noise), ST = group_stats(X, Y, sf, gc, Wm, T)
        S = (n, s_sum, t_sum, t_sq, t_noise, ST)
        fa = full[full.cell == c].set_index("unit")
        active = fa.activity_sd.reindex(units).values >= MIN_ACTIVITY_SD
        pf, _ = percentiles(S, np.arange(len(gc)), sigma, member)
        assert np.allclose(pf[active], fa.pct.reindex(units).values[active], equal_nan=True), c
        donors = gc.donor.unique()
        res = []
        for i in range(N_SPLITS):
            a = set(rng.permutation(donors)[: len(donors) // 2])
            ina = gc.donor.isin(a).values
            pa, _ = percentiles(S, np.where(ina)[0], sigma, member)
            pb, _ = percentiles(S, np.where(~ina)[0], sigma, member)
            ok = active & np.isfinite(pa) & np.isfinite(pb)
            res.append((spearmanr(pa[ok], pb[ok])[0], spearmanr(pa[ok], pf[ok])[0], int(n[ina].sum()),
                        int(n[~ina].sum()), pa, pb, ok))
        r = np.array([x[0] for x in res])
        splits += [dict(cell=c, split=i, r_half=x[0], r_half_vs_full=x[1], cells_a=x[2], cells_b=x[3])
                   for i, x in enumerate(res)]
        med = float(np.median(r))
        summ.append(dict(population=pop, state=state, cell=c, donors=len(donors), cells=int(n.sum()),
                         programs=int(active.sum()), r_half_median=med, r_half_q10=np.quantile(r, .1),
                         r_half_q90=np.quantile(r, .9), share_positive=(r > 0).mean(),
                         reliability_full=2 * med / (1 + med)))
        _, _, _, _, pa, pb, ok = res[int(np.argsort(r)[len(r) // 2])]
        ex = pd.DataFrame({"half_a": pa[ok], "half_b": pb[ok]}, index=np.array(units)[ok])
        ex.to_csv(OUT / f"{c}_example.csv")
        plot(ex, c, summ[-1])
        print(pd.Series(summ[-1]).to_string(), flush=True)
    return splits, summ


def plot(ex, c, s):
    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=160)
    ax.axhline(50, color="#d0d0d0", lw=0.8, zorder=0)
    ax.axvline(50, color="#d0d0d0", lw=0.8, zorder=0)
    ax.plot([0, 100], [0, 100], color="#b0b0b0", lw=0.8, ls=(0, (4, 4)), zorder=0)
    ax.scatter(ex.half_a, ex.half_b, s=14, color="#3b6fb6", alpha=0.55, lw=0)
    slope, icpt = np.polyfit(ex.half_a, ex.half_b, 1)
    ax.plot([0, 100], [icpt, icpt + 100 * slope], color="#1f3f6e", lw=2)
    half = s["donors"] // 2
    ax.set(xlim=(0, 100), ylim=(0, 100), aspect="equal",
           xlabel=f"Donor half A ({half} donors): ANTXR2 percentile per program",
           ylabel=f"Donor half B ({s['donors'] - half} donors): ANTXR2 percentile per program")
    ax.set_title(f"{LANDSCAPE_ATLAS} {c}: two random halves of the donors\n{len(ex)} programs · median split shown, "
                 f"Spearman {s['r_half_median']:.2f} · {N_SPLITS} splits: {s['r_half_q10']:.2f} to "
                 f"{s['r_half_q90']:.2f} (10th–90th pct)", fontsize=9.5, loc="left")
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    ax.tick_params(labelsize=8.5)
    ax.xaxis.label.set_size(9)
    ax.yaxis.label.set_size(9)
    fig.tight_layout()
    fig.savefig(OUT / f"{c}_example.png")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pops", nargs="*", default=["fibroblast_pooled"])
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)
    splits, summ = [], []
    for pop in [p for p in PL.POPS if pop_slug(p) in a.pops]:
        s, m = run(pop, rng)
        splits += s
        summ += m
    out = pd.DataFrame(summ)
    for name, df in (("splits", pd.DataFrame(splits)), ("summary", out)):
        path = OUT / f"{name}.csv"
        if path.exists():                  # --pops reruns replace their own cells
            old = pd.read_csv(path)
            df = pd.concat([old[~old.cell.isin(out.cell)], df], ignore_index=True)
        df.to_csv(path, index=False)
    out = pd.read_csv(OUT / "summary.csv")
    pd.set_option("display.width", 250)
    print(out.round(3).to_string(index=False))
