"""Tables and figures for the between-population difference (population A -
population B, from the run's contrast) in ANTXR2 x program association, from
program_scores.py (OUT_DIR/<run>/).

The headline statistic is the EXCESS difference: ANTXR2's delta
(rho_a - rho_b) minus the mean delta of 50 genes matched to ANTXR2 in both
populations, and excess_z = excess / SD of the anchors' deltas. A program can
differ between cell types for any gene; only the excess is about ANTXR2.

Per run:
    top_programs.csv   comparable programs ranked by |excess_z| (PC1 weights),
                       with the equal-weight excess beside it and coverage
    forest.png         top programs: ANTXR2 delta with donor-bootstrap CI over
                       the anchors' delta band
    calibration.png    excess_z across all programs against N(0,1): how many
                       large values chance alone would give
    scatter.png        rho_a vs rho_b, |excess_z| > 2 highlighted
    sensitivity.png    excess, PC1 vs equal weights
    vs_<baseline>.png  (with --baseline) excess here vs in another run

Exploratory throughout: weights and associations come from the same cells,
programs share genes (so they are not independent), and nothing is
corrected for the number of programs screened.

Usage:
    python program_report.py <run> [--top 25] [--baseline <run>]
"""
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
from scipy.stats import norm                        # noqa: E402

from program_config import OUT_DIR                  # noqa: E402

SURFACE = "#fcfcfb"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#cfcec8"
BLUE, ORANGE = "#2a78d6", "#eb6834"                  # leans population A, leans B
Z = 2.0

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "font.size": 9, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
})


def short(pid):
    col, name = pid.split(":", 1)
    return f"{name.replace('HALLMARK_', '').replace('NABA_', '')} ({col})"


def load(run):
    d = OUT_DIR / run
    pops = json.load(open(d / "run.json"))["pops"]
    return d, pd.read_csv(d / "differential.csv"), pd.read_csv(d / "programs.csv"), pops


def top_table(D, P, pops):
    pc = D[D.scheme == "pc1"].set_index("unit")
    eq = D[D.scheme == "equal"].set_index("unit")
    t = pc.join(eq[["excess", "excess_z"]].add_suffix("_equal"))
    t = t.join(P.set_index("program_id")[["collection", "n_members", "n_measured", "pc1_ve",
                                          "cos_pc1_equal"]])
    t["excess_ci_excludes_0"] = (t.excess_ci_lo > 0) | (t.excess_ci_hi < 0)
    t["leans"] = np.where(t.excess > 0, pops[0], pops[1])
    t = t.reindex(t.excess_z.abs().sort_values(ascending=False).index)
    t.insert(0, "rank", np.arange(1, len(t) + 1))
    return t.reset_index().drop(columns="scheme")


def forest(t, pops, path, title):
    t = t.iloc[::-1]
    y = np.arange(len(t))
    fig, ax = plt.subplots(figsize=(6.8, 0.26 * len(t) + 2.0))
    ax.barh(y, 2 * 1.96 * t.anchor_delta_sd, left=t.anchor_delta_mean - 1.96 * t.anchor_delta_sd,
            height=0.55, color=BAND, label="matched anchors: mean ± 1.96 SD of Δρ")
    col = np.where(t.excess > 0, BLUE, ORANGE)
    ax.hlines(y, t.delta_ci_lo, t.delta_ci_hi, color=col, lw=1.5)
    ax.scatter(t.delta, y, s=26, color=col, edgecolor=SURFACE, lw=1.2, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(y, [f"{short(u)}  z={z:+.1f}" for u, z in zip(t.unit, t.excess_z)], fontsize=7)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(f"Δρ = ρ({pops[0]}) − ρ({pops[1]}), PC1 weights\n"
                  "ANTXR2 point with donor-bootstrap 95% CI")
    ax.scatter([], [], color=BLUE, s=20, label=f"ANTXR2 leans {pops[0]} vs anchors")
    ax.scatter([], [], color=ORANGE, s=20, label=f"ANTXR2 leans {pops[1]} vs anchors")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(-0.4, -0.1 * 20 / (len(t) + 5)))
    ax.set_title(title, fontsize=8, color=INK2, loc="right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def calibration(t, path, title):
    z = t.excess_z.dropna().values
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    bins = np.linspace(-4.5, 4.5, 37)
    ax.hist(z, bins=bins, color=BLUE, alpha=0.85, rwidth=0.9, label=f"{len(z)} programs")
    x = np.linspace(-4.5, 4.5, 200)
    ax.plot(x, len(z) * (bins[1] - bins[0]) * norm.pdf(x), color=INK2, lw=1.2,
            label="N(0,1): ANTXR2 as one more anchor")
    n_hi = int((np.abs(z) > Z).sum())
    ax.set_title(f"{title}\n|excess z| > {Z:g}: {n_hi} observed, "
                 f"{2 * norm.sf(Z) * len(z):.1f} expected by chance", fontsize=8, color=INK2, loc="left")
    ax.set_xlabel("excess z (ANTXR2 Δρ vs matched anchors' Δρ)")
    ax.set_ylabel("programs")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def scatter(t, pops, path, title):
    hi = t.excess_z.abs() > Z
    lim = max(0.05, np.nanmax(np.abs(t[["rho_a", "rho_b"]].values))) * 1.08
    fig, ax = plt.subplots(figsize=(5.4, 5.8))
    ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=1, zorder=1)
    ax.scatter(t.rho_a[~hi], t.rho_b[~hi], s=12, color="#b5b4ae", lw=0, zorder=2,
               label=f"|excess z| ≤ {Z:g}")
    for side, c, lab in [(t.excess > 0, BLUE, f"leans {pops[0]}"), (t.excess < 0, ORANGE, f"leans {pops[1]}")]:
        s = hi & side
        ax.scatter(t.rho_a[s], t.rho_b[s], s=26, color=c, edgecolor=SURFACE, lw=1.5, zorder=3,
                   label=f"|excess z| > {Z:g}, {lab}")
    # label up to five per side in stacked columns away from the cloud
    for side, x0, y0 in [(t.excess < 0, -0.97, 0.95), (t.excess > 0, 0.30, -0.45)]:
        lab = t[hi & side].head(5)
        for i, (_, r) in enumerate(lab.iterrows()):
            ax.annotate(short(r.unit), (r.rho_a, r.rho_b),
                        xytext=(x0 * lim, (y0 - i * 0.075) * lim), textcoords="data",
                        fontsize=6, color=INK2, va="center",
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"ANTXR2–program ρ in {pops[0]}")
    ax.set_ylabel(f"ANTXR2–program ρ in {pops[1]}")
    ax.set_title(title, fontsize=8, color=INK2, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.12))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def paired(x, y, xl, yl, path, title):
    ok = x.notna() & y.notna()
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    lim = np.nanmax(np.abs(np.r_[x[ok], y[ok]])) * 1.05
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=1)
    ax.scatter(x[ok], y[ok], s=10, color=BLUE, lw=0, alpha=0.75)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(f"{title}\nr = {r:.2f} across {ok.sum()} programs", fontsize=8, color=INK2, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--baseline", default=None)
    a = ap.parse_args()
    d, D, P, pops = load(a.run)
    t = top_table(D, P, pops)
    t.to_csv(d / "top_programs.csv", index=False)
    title = a.run.replace("_", " ")
    forest(t.head(a.top), pops, d / "forest.png", f"{title}: top {min(a.top, len(t))} by |excess z|")
    calibration(t, d / "calibration.png", title)
    scatter(t, pops, d / "scatter.png",
            f"{title}\n{pops[0]} vs {pops[1]}, comparable programs, PC1 weights")
    paired(t.excess_equal, t.excess, "excess Δρ, signed equal weights", "excess Δρ, PC1 weights",
           d / "sensitivity.png", f"{title}: weighting sensitivity")
    if a.baseline:
        b = top_table(*load(a.baseline)[1:]).set_index("unit")
        s = t.set_index("unit")
        paired(b.excess.reindex(s.index), s.excess, f"excess Δρ, {a.baseline}",
               f"excess Δρ, {a.run}", d / f"vs_{a.baseline}.png", "same programs, two runs")
    print(f"wrote tables and figures to {d}")


if __name__ == "__main__":
    main()
