"""Tables and figures for the ANTXR2 co-expression landscape
(program_landscape.py output in LANDSCAPE_OUT) -> LANDSCAPE_OUT/report/.

The value in every cell is excess z: ANTXR2's association with a program
minus the mean of 50 genes matched to ANTXR2 in that population, over the
anchors' SD (PC1 weights). Positive = ANTXR2 tracks the program more than
genes like it do.

    measurability.png  population x state: filled or why not, ANTXR2 raw mean,
                       reliability
    calibration.png    per cell: how many programs reach |z| > 2 for ANTXR2,
                       as a percentile of the same count for each anchor
                       (leave-one-out against the other anchors); the median
                       z; and the gene-length diagnostic (Spearman of
                       gene-level z with partner genomic span, ~0 once anchors
                       are span-matched)
    heatmap.png        the most variable programs x filled cells, clustered
    similarity.png     cell x cell correlation of z profiles, clustered
    wide_z.csv         programs x cells, excess z
    programs.csv       per program: cells scored, cells z > 2 / z < -2,
                       mean z, where
    genes_top.csv      per cell, top gene-level partners by z (both signs)
    genes_recurrent.csv  genes with |z| > 3 in >= 3 cells
    consensus_<lineage>.csv / consensus.png
                       lineage-wide profile: per program, ANTXR2's mean z over
                       the lineage's filled cells against the same mean for
                       each of 50 lineage anchors (one gene across all cells,
                       program_config.LANDSCAPE_LINEAGES); and per gene, how
                       many programs reach |consensus z| > 2 -- ANTXR2's
                       percentile there is the headline "is its landscape
                       more structured than a matched gene's"

Descriptive. Programs share genes and the matrix is screened as a whole, so
read every value against the calibration plot, not as a test.

Usage:
    python program_landscape_report.py [--top 60]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402

from scipy.stats import spearmanr                   # noqa: E402

from program_config import (                        # noqa: E402
    LANDSCAPE_FIB, LANDSCAPE_LINEAGES, LANDSCAPE_OUT, LANDSCAPE_POPS, LANDSCAPE_STATES,
    gene_spans, pop_slug,
)
from program_report import BLUE, GRID, INK2, MUTED, ORANGE, SURFACE, short  # noqa: E402

REPORT = LANDSCAPE_OUT / "report"
REPORT.mkdir(parents=True, exist_ok=True)
STATES = list(LANDSCAPE_STATES)
STATE_LABEL = {"normal": "normal", "cd_uninflamed": "CD uninfl.", "cd_inflamed": "CD infl."}
POPS = LANDSCAPE_POPS + [LANDSCAPE_FIB]
DIVERGING = LinearSegmentedColormap.from_list("div", [ORANGE, "#f6c4ab", "#f0efec", "#b3cdf0", BLUE])


def label(cell, cells):
    r = cells.set_index("cell").loc[cell]
    return f"{r.population}  ·  {STATE_LABEL[r.state]}"


def cluster_order(M):
    """Leaf order of average-linkage clustering on correlation distance;
    missing values count as 0."""
    if len(M) < 3:
        return np.arange(len(M))
    A = np.nan_to_num(M)
    A = A - A.mean(1, keepdims=True)
    nrm = np.linalg.norm(A, axis=1, keepdims=True)
    A = A / np.where(nrm > 0, nrm, 1)
    return leaves_list(linkage(A, "average", metric="euclidean"))


def measurability(cells, path):
    fig, ax = plt.subplots(figsize=(6.6, 0.28 * len(POPS) + 1.4))
    for i, p in enumerate(POPS):
        for j, s in enumerate(STATES):
            r = cells[(cells.population == p) & (cells.state == s)]
            if r.empty:
                continue
            r = r.iloc[0]
            if r.filled:
                ax.scatter(j, i, s=150, color=BLUE, edgecolor=SURFACE, lw=1.5, zorder=3)
                rel = r.get("antxr2_reliability", np.nan)
                ax.annotate(f"{r.antxr2_raw_mean:.2f} · rel {rel:.2f} · {int(r.n_donors)}d",
                            (j, i), xytext=(9, 0), textcoords="offset points", va="center",
                            fontsize=6, color=INK2)
            else:
                ax.scatter(j, i, s=60, facecolor=SURFACE, edgecolor=MUTED, lw=1, zorder=3)
                txt = (f"{r.antxr2_raw_mean:.2f} below floor" if r.reason == "ANTXR2 below floor"
                       else f"{int(r.n_donors)} donors")
                ax.annotate(txt, (j, i), xytext=(9, 0), textcoords="offset points", va="center",
                            fontsize=6, color=MUTED)
    ax.set_yticks(range(len(POPS)), POPS, fontsize=7)
    ax.set_xticks(range(len(STATES)), [STATE_LABEL[s] for s in STATES])
    ax.set_xlim(-0.4, len(STATES) - 0.1)
    ax.invert_yaxis()
    ax.grid(False)
    ax.set_title("ANTXR2 measurability: filled cells show raw mean · reliability · donors",
                 fontsize=8, color=INK2, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def anchor_calibration(cells):
    """Per cell: programs with |z| > 2 for ANTXR2 and for each anchor, each
    scored against the other anchors (ANTXR2 excluded from an anchor's
    reference). PC1 weights."""
    span = np.log10(gene_spans())
    rows = []
    for pop in cells[cells.filled].population.unique():
        z = np.load(LANDSCAPE_OUT / "boot" / f"{pop_slug(pop)}.npz")
        pc = z["schemes"] == "pc1"
        for state in cells[(cells.population == pop) & cells.filled].state:
            c = f"{pop_slug(pop)}@{state}"
            if f"{state}_point" not in z:
                continue
            P = z[f"{state}_point"][pc].astype(float)          # programs x targets
            P = P[np.isfinite(P[:, 0])]
            A = P[:, 1:]
            n_hi, med = [], []
            for t in range(P.shape[1]):
                ref = np.delete(A, t - 1, axis=1) if t > 0 else A
                with np.errstate(invalid="ignore"):
                    zt = (P[:, t] - np.nanmean(ref, 1)) / np.nanstd(ref, 1)
                n_hi.append(int(np.nansum(np.abs(zt) > 2)))
                med.append(np.nanmedian(zt))
            g = pd.read_csv(LANDSCAPE_OUT / "genes" / f"{c}.csv.gz")
            g = g[~g.is_anchor].assign(span=lambda d: d.gene.map(span)).dropna(subset=["span", "z"])
            rows.append(dict(cell=c, n_programs=len(P), antxr2_n_hi=n_hi[0],
                             anchor_n_hi_median=float(np.median(n_hi[1:])),
                             antxr2_n_hi_pct=float(np.mean(np.array(n_hi[1:]) < n_hi[0]) * 100),
                             antxr2_median_z=med[0],
                             anchor_median_z_sd=float(np.nanstd(med[1:])),
                             gene_z_vs_span=spearmanr(g.z, g.span)[0]))
    return pd.DataFrame(rows).set_index("cell")


def calibration(cells, path):
    t = anchor_calibration(cells)
    t = t.loc[[c for c in cells.cell if c in t.index]]
    y = np.arange(len(t))
    fig, axs = plt.subplots(1, 3, figsize=(10.4, 0.24 * len(t) + 1.8), sharey=True)
    axs[0].barh(y, t.antxr2_n_hi_pct, color=BLUE, height=0.6)
    axs[0].axvline(50, color=INK2, lw=1)
    axs[0].set_xlim(0, 100)
    axs[0].set_xlabel("programs with |z| > 2: ANTXR2's percentile\namong its 50 matched genes")
    axs[1].barh(y, t.antxr2_median_z, color=np.where(t.antxr2_median_z > 0, BLUE, ORANGE), height=0.6)
    axs[1].axvline(0, color=INK2, lw=1)
    axs[1].set_xlabel("median excess z across programs\n(cell-wide shift)")
    axs[2].barh(y, t.gene_z_vs_span, color=MUTED, height=0.6)
    axs[2].axvline(0, color=INK2, lw=1)
    axs[2].set_xlabel("gene-level z vs partner genomic span\n(Spearman; ~0 = no length axis left)")
    axs[0].set_yticks(y, [label(c, cells) for c in t.index], fontsize=6)
    axs[0].invert_yaxis()
    for a in axs:
        a.grid(axis="y", visible=False)
    fig.suptitle("Calibration per cell (PC1 weights)", fontsize=8, color=INK2, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return t


def heatmap(Z, cells, path, top):
    keep = Z.notna().mean(1) >= 0.75
    Zk = Z[keep]
    sel = Zk.abs().max(1).sort_values(ascending=False).index[:top]
    M = Zk.loc[sel]
    M = M.iloc[cluster_order(M.values), cluster_order(M.values.T)]
    fig, ax = plt.subplots(figsize=(0.24 * M.shape[1] + 4.2, 0.15 * len(M) + 3.2))
    im = ax.imshow(M.values, cmap=DIVERGING, vmin=-4, vmax=4, aspect="auto",
                   interpolation="nearest")
    na = np.isnan(M.values)
    ax.imshow(np.where(na, 1, np.nan), cmap=LinearSegmentedColormap.from_list("na", [GRID, GRID]),
              aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(M)), [short(u) for u in M.index], fontsize=5.5)
    ax.set_xticks(range(M.shape[1]), [label(c, cells) for c in M.columns], rotation=90, fontsize=5.5)
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("excess z (ANTXR2 vs matched genes)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(f"{len(M)} programs with the largest |excess z| (scored in ≥ 75% of filled cells); "
                 "grey = not scored", fontsize=8, color=INK2, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def similarity(Z, cells, path):
    C = Z.corr(min_periods=100)
    o = cluster_order(C.values)
    C = C.iloc[o, o]
    fig, ax = plt.subplots(figsize=(0.22 * len(C) + 3.6, 0.22 * len(C) + 2.6))
    im = ax.imshow(C.values, cmap=DIVERGING, vmin=-1, vmax=1, interpolation="nearest")
    lab = [label(c, cells) for c in C.index]
    ax.set_xticks(range(len(C)), lab, rotation=90, fontsize=5.5)
    ax.set_yticks(range(len(C)), lab, fontsize=5.5)
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("correlation of excess-z profiles over shared programs", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Which cells look alike from ANTXR2's side", fontsize=8, color=INK2, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return C


def cell_target_z(P, targets, nb):
    """programs x genes: z of each target vs the population anchors (columns
    1..nb-1), leaving the target itself out of its own reference."""
    A, pa = P[:, 1:nb], np.array(targets[1:nb])
    out = {}
    for t, g in enumerate(targets):
        if t >= nb or t == 0:
            ref = A[:, pa != g] if t > 0 else A
            with np.errstate(invalid="ignore"):
                out[g if t else "ANTXR2"] = (P[:, t] - np.nanmean(ref, 1)) / np.nanstd(ref, 1)
    return pd.DataFrame(out)


def consensus(cells, min_frac=0.75):
    res = {}
    for lin, members in LANDSCAPE_LINEAGES.items():
        fc = cells[cells.population.isin(members) & cells.filled & (cells.n_programs > 0)]
        if len(fc) < 3:
            continue
        stack, units = [], None
        for pop in fc.population.unique():
            z = np.load(LANDSCAPE_OUT / "boot" / f"{pop_slug(pop)}.npz")
            pc = z["schemes"] == "pc1"
            units = z["units"][pc]
            targets, nb = list(z["targets"].astype(str)), int(z["n_pop_anchors"]) + 1
            for st in fc[fc.population == pop].state:
                P = z[f"{st}_point"][pc].astype(float)
                d = cell_target_z(P, targets, nb)
                d[~np.isfinite(P[:, 0])] = np.nan
                stack.append(d.set_index(units))
        genes = sorted(set.intersection(*[set(d.columns) for d in stack]))
        M = np.stack([d[genes].values for d in stack])            # cells x programs x genes
        keep = np.isfinite(M[:, :, genes.index("ANTXR2")]).mean(0) >= min_frac
        with np.errstate(invalid="ignore"):
            C = np.nanmean(M[:, keep], 0)                         # programs x genes
        C = pd.DataFrame(C, index=units[keep], columns=genes)
        others = [g for g in genes if g != "ANTXR2"]
        Z = pd.DataFrame(index=C.index)
        for g in genes:                                           # leave-one-out standardisation
            ref = C[[o for o in others if o != g]]
            Z[g] = (C[g] - ref.mean(1)) / ref.std(1)
        n_hi = (Z.abs() > 2).sum()
        t = pd.DataFrame(dict(mean_z=C["ANTXR2"], null_mean=C[others].mean(1),
                              null_sd=C[others].std(1), consensus_z=Z["ANTXR2"],
                              rank_among_genes=C.rank(axis=1)["ANTXR2"].astype(int),
                              n_genes=len(genes)))
        t.sort_values("consensus_z").to_csv(REPORT / f"consensus_{lin}.csv")
        res[lin] = dict(t=t, n_cells=len(stack), n_hi=n_hi,
                        pct=float((n_hi[others] < n_hi["ANTXR2"]).mean() * 100))
    return res


def consensus_plot(res, path, n=12):
    lins = list(res)
    fig, axs = plt.subplots(1, len(lins), figsize=(4.6 * len(lins), 0.21 * 2 * n + 2.2), squeeze=False)
    for ax, lin in zip(axs[0], lins):
        t = res[lin]["t"].dropna(subset=["consensus_z"]).sort_values("consensus_z")
        t = pd.concat([t.head(n), t.tail(n)])
        y = np.arange(len(t))
        ax.barh(y, t.consensus_z, color=np.where(t.consensus_z > 0, BLUE, ORANGE), height=0.6)
        ax.axvspan(-2, 2, color=GRID, zorder=0)
        ax.axvline(0, color=INK2, lw=1)
        ax.set_yticks(y, [short(u) for u in t.index], fontsize=6)
        ax.grid(axis="y", visible=False)
        r = res[lin]
        ax.set_title(f"{lin}: {r['n_cells']} cells\nANTXR2: {int(r['n_hi']['ANTXR2'])} programs |z| > 2,\n"
                     f"percentile {r['pct']:.0f} among {len(r['n_hi']) - 1} matched genes",
                     fontsize=7.5, color=INK2, loc="left")
        ax.set_xlabel("lineage-wide z vs 50 matched genes\n(grey band: |z| ≤ 2)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def gene_tables(cells, n_top=25):
    top, hits = [], []
    for c in cells[cells.filled].cell:
        g = pd.read_csv(LANDSCAPE_OUT / "genes" / f"{c}.csv.gz")
        g = g[~g.is_anchor].dropna(subset=["z"])
        for side, asc in [("pos", False), ("neg", True)]:
            t = g.sort_values("z", ascending=asc).head(n_top)
            top.append(t.assign(cell=c, side=side))
        hits.append(g[g.z.abs() > 3].assign(cell=c))
    pd.concat(top).to_csv(REPORT / "genes_top.csv", index=False)
    H = pd.concat(hits)
    rec = (H.assign(pos=H.z > 0).groupby("gene")
           .agg(n_cells=("cell", "size"), n_pos=("pos", "sum"), mean_z=("z", "mean"),
                cells=("cell", lambda s: ";".join(s)))
           .query("n_cells >= 3").sort_values("n_cells", ascending=False))
    rec.to_csv(REPORT / "genes_recurrent.csv")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=60)
    a = ap.parse_args()
    cells = pd.read_csv(LANDSCAPE_OUT / "cells.csv")
    order = {p: i for i, p in enumerate(POPS)}
    cells = cells.assign(_o=cells.population.map(order), _s=cells.state.map(STATES.index)) \
        .sort_values(["_o", "_s"]).drop(columns=["_o", "_s"])
    L = pd.read_csv(LANDSCAPE_OUT / "landscape.csv")
    L = L[L.scheme == "pc1"]
    Z = L.pivot(index="unit", columns="cell", values="excess_z")
    Z = Z[[c for c in cells.cell if c in Z.columns]]
    Z.to_csv(REPORT / "wide_z.csv")

    measurability(cells, REPORT / "measurability.png")
    cal = calibration(cells, REPORT / "calibration.png")
    cal.to_csv(REPORT / "calibration.csv")
    heatmap(Z, cells, REPORT / "heatmap.png", a.top)
    C = similarity(Z, cells, REPORT / "similarity.png")
    C.to_csv(REPORT / "similarity.csv")

    P = pd.DataFrame(dict(n_cells=Z.notna().sum(1), n_pos=(Z > 2).sum(1), n_neg=(Z < -2).sum(1),
                          mean_z=Z.mean(1)))
    P["pos_cells"] = [";".join(Z.columns[Z.loc[u] > 2]) for u in Z.index]
    P["neg_cells"] = [";".join(Z.columns[Z.loc[u] < -2]) for u in Z.index]
    P.sort_values("mean_z", ascending=False).to_csv(REPORT / "programs.csv")
    rec = gene_tables(cells)
    res = consensus(cells)
    consensus_plot(res, REPORT / "consensus.png")
    for lin, r in res.items():
        print(f"consensus {lin}: {r['n_cells']} cells; ANTXR2 programs |z|>2 = {int(r['n_hi']['ANTXR2'])}, "
              f"matched genes median {r['n_hi'].drop('ANTXR2').median():.0f}, percentile {r['pct']:.0f}")
    print(f"{cells.filled.sum()} filled cells, {len(Z)} programs; report in {REPORT}")
    print(cal.round(2).to_string())
    print(rec.head(30).drop(columns="cells").round(2).to_string())


if __name__ == "__main__":
    main()
