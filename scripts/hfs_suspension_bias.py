"""Matched control for single-nucleus vs whole-cell bias in the COL6:ANTXR2 ratio.

Why this exists
---------------
The tendon atlases are single-nucleus and every other tissue in this arm is
whole-cell, so tendon's very low ``col6_per_antxr2`` could be either biology or
platform. An earlier version of this analysis asserted the bias was specific to
ANTXR2; that claim was wrong, because it compared tendon nuclei to *other
tissues'* cells, confounding platform with tissue, and never asked whether
ANTXR2 shifts more than genes in general.

This script does the comparison properly, using the two datasets that contain
both suspension types for the same tissue:

  kidney  interstitial fibroblasts, restricted to donors sampled both ways
  heart   cardiac fibroblasts, region-matched, plus a single-donor check

and reports the ANTXR2 shift as a *percentile of the genome-wide shift*, which
is the statistic that separates "this gene is special" from "everything moved".

What it shows: the direction reproduces — nuclei recover relatively more ANTXR2
and relatively less collagen — but ANTXR2 sits in the upper-middle of the
genome-wide distribution, not the tail. So the ratio is biased from both ends,
and tendon's value is understated by a bounded but dataset-dependent factor,
which is why tendon is reported as unplaceable rather than low.

Outputs
-------
    hfs_suspension_bias_matched_controls.csv
    fig_suspension_bias_control.png
"""
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hfs_config import (COLLAGEN_GROUPS, COLLAGEN_RE, HFS_PB_DIR, HFS_OUTPUT_DIR,
                        HFS_FIGURES_DIR, MIN_CPM_FOR_RATIO)

GREY, FOC, DEP = "#9aa5ad", "#0b5394", "#c1741a"

REGION_LABEL = {"acute kidney injury": "kidney, AKI",
                "chronic kidney disease": "kidney, CKD",
                "apex of heart": "heart, apex", "heart left ventricle": "heart, LV",
                "heart right ventricle": "heart, RV",
                "interventricular septum": "heart, septum",
                "left cardiac atrium": "heart, LA", "right cardiac atrium": "heart, RA",
                "donor": "heart, one donor (paired)"}


def _cpm(name):
    Z = np.load(HFS_PB_DIR / f"{name}.npz", allow_pickle=False)
    genes = np.array([str(x) for x in Z["genes"]])
    groups = [str(x) for x in Z["groups"]]
    S = Z["sum"].astype(np.float64)
    return S / S.sum(1, keepdims=True) * 1e6, genes, groups


def log_ratio_series(name, stratum):
    """Genome-wide log2(nucleus / cell) for one stratum, genes present in both."""
    cpm, genes, groups = _cpm(name)
    c = cpm[groups.index(f"cell | {stratum}")]
    n = cpm[groups.index(f"nucleus | {stratum}")]
    ok = (c > MIN_CPM_FOR_RATIO) & (n > MIN_CPM_FOR_RATIO)
    return pd.Series(np.log2(n[ok] / c[ok]), index=genes[ok])


def bias_table(name):
    """Per-stratum nucleus/cell fold-changes for ANTXR2, COL6 and their ratio."""
    cpm, genes, groups = _cpm(name)
    gi = {g: i for i, g in enumerate(genes)}
    col6 = [gi[g] for g in COLLAGEN_GROUPS["COL6"] if g in gi]
    all_col = [i for i, g in enumerate(genes) if re.fullmatch(COLLAGEN_RE, g)]

    out = []
    for stratum in sorted({g.split(" | ")[1] for g in groups}):
        try:
            c = cpm[groups.index(f"cell | {stratum}")]
            n = cpm[groups.index(f"nucleus | {stratum}")]
        except ValueError:
            continue   # stratum present in only one suspension type
        ac, an = c[gi["ANTXR2"]], n[gi["ANTXR2"]]
        l6c, l6n = c[col6].sum(), n[col6].sum()
        ok = (c > MIN_CPM_FOR_RATIO) & (n > MIN_CPM_FOR_RATIO)
        lr = np.log2(n[ok] / c[ok])
        s = pd.Series(lr, index=genes[ok])
        out.append(dict(stratum=stratum, n_genes=int(ok.sum()),
                        global_median_lr=float(np.median(lr)),
                        antxr2_fc=an / ac, col6_fc=l6n / l6c,
                        share_fc=(l6n / n[all_col].sum()) / (l6c / c[all_col].sum()),
                        load_fc=(l6n / an) / (l6c / ac),
                        antxr2_pctile=float((s < s["ANTXR2"]).mean() * 100)))
    return pd.DataFrame(out)


def build():
    K = bias_table("kidney_susp_pseudobulk").assign(tissue="kidney")
    H = bias_table("heart_susp_pseudobulk").assign(tissue="heart")
    HD = bias_table("heart_donor_susp_pseudobulk").assign(tissue="heart (1 donor, paired)")
    BI = pd.concat([K, H, HD], ignore_index=True)
    BI["matched_on"] = np.select(
        [BI.tissue == "kidney", BI.tissue == "heart"],
        ["donor set + disease + cell type + chemistry + study",
         "region + disease + cell type + chemistry + study"],
        default="same donor + cell type + chemistry + study")
    BI["load_understated_x"] = 1 / BI.load_fc
    return BI[["tissue", "stratum", "matched_on", "n_genes", "global_median_lr",
               "antxr2_fc", "antxr2_pctile", "col6_fc", "share_fc", "load_fc",
               "load_understated_x"]]


def figure(BI, tissue_table):
    """Four panels: genome-wide shift, ANTXR2 percentile, ratio decomposition,
    and the corrected interval for the two single-nucleus tendon rows."""
    try:                                     # project figure style if available
        from figure_style import apply_figure_style, panel_letter
        apply_figure_style(frame="open", sizes=(8, 7, 6))
    except ImportError:
        def panel_letter(ax, l):
            ax.text(-0.13, 1.06, l, transform=ax.transAxes, fontsize=9, fontweight="bold")

    SH = log_ratio_series("heart_susp_pseudobulk", "heart left ventricle")
    SK = log_ratio_series("kidney_susp_pseudobulk", "chronic kidney disease")

    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.34)
    ax = [fig.add_subplot(gs[i, j]) for i in (0, 1) for j in (0, 1)]

    # (a) the whole transcriptome moves
    for s, lab, c in [(SH, "heart, left ventricle", "#7aa6c2"), (SK, "kidney, CKD", "#b8b8b8")]:
        ax[0].hist(s.values, bins=70, range=(-8, 8), histtype="step", lw=1.2,
                   color=c, label=lab, density=True)
    ax[0].axvline(0, color=GREY, lw=0.7, ls=":")
    ymax = ax[0].get_ylim()[1]
    for g, c, frac in [("ANTXR2", FOC, 0.92), ("COL6A1", DEP, 0.62), ("COL1A1", DEP, 0.40)]:
        v = float(SH[g])
        ax[0].plot([v, v], [0, ymax * frac], color=c, lw=1.1)
        ax[0].plot([v], [ymax * frac], marker="v", ms=4, color=c)
        ax[0].text(v, ymax * frac + 0.006, f"$\\it{{{g}}}$", ha="center", va="bottom",
                   color=c, fontsize=6)
    ax[0].set_xlabel("nucleus / cell, log$_2$ fold-change")
    ax[0].set_ylabel("density of genes")
    ax[0].set_title("Every gene shifts; the receptor is\nnot an outlier", loc="left")
    ax[0].legend(frameon=False, loc="upper left", fontsize=6, handlelength=1.2)
    ax[0].set_xlim(-8, 8)
    ax[0].set_ylim(0, ymax * 1.22)

    # (b) where ANTXR2 sits in that distribution
    y = np.arange(len(BI))[::-1]
    ax[1].axvline(50, color=GREY, lw=0.8, ls="--")
    ax[1].scatter(BI.antxr2_pctile, y, s=26, color=FOC, zorder=3)
    ax[1].set_yticks(y)
    ax[1].set_yticklabels([REGION_LABEL.get(s, s) for s in BI.stratum], fontsize=6.5)
    ax[1].set_xlim(0, 100)
    ax[1].set_ylim(-0.8, len(BI) - 0.2)
    ax[1].set_xlabel("percentile of $\\it{ANTXR2}$ shift\namong all genes")
    ax[1].set_title("Upper-middle of the distribution,\nnot the tail", loc="left")
    ax[1].text(48, len(BI) - 0.45, "genome\nmedian", ha="right", va="top",
               color=GREY, fontsize=6)

    # (c) the ratio is biased from both ends
    xs = np.arange(3)
    for _, r in BI.iterrows():
        ax[2].plot(xs, np.log2([r.antxr2_fc, r.col6_fc, r.load_fc]),
                   color="#c8c8c8", lw=0.7, zorder=1)
    med = np.log2([BI.antxr2_fc.median(), BI.col6_fc.median(), BI.load_fc.median()])
    ax[2].plot(xs, med, color=FOC, lw=1.8, marker="o", ms=5, zorder=3)
    ax[2].axhline(0, color=GREY, lw=0.7, ls=":")
    ax[2].set_xticks(xs)
    ax[2].set_xticklabels(["$\\it{ANTXR2}$", "COL6\nsum", "COL6 per\n$\\it{ANTXR2}$"],
                          fontsize=6.5)
    ax[2].set_ylabel("nucleus / cell, log$_2$")
    ax[2].set_xlim(-0.35, 2.35)
    ax[2].set_title("Receptor up, collagen down —\nthe ratio is biased from both ends",
                    loc="left")
    ax[2].text(2.05, med[2] - 0.55, f"{1/BI.load_fc.median():.0f}x understated\n(median)",
               ha="right", va="top", color=FOC, fontsize=6)

    # (d) tendon, corrected
    wc = tissue_table[tissue_table.susp == "cell"]
    ref_core = wc[wc.involvement == "affected (core)"].col6_per_antxr2.median()
    ref_max = wc[wc.involvement == "not reported"].col6_per_antxr2.max()
    lo, hi = BI.load_understated_x.min(), BI.load_understated_x.max()
    tn = tissue_table[tissue_table.susp == "nucleus"]

    ax[3].axvspan(0, ref_max, color="#ececec", zorder=0)
    ax[3].axvline(ref_max, color=GREY, lw=0.9)
    ax[3].axvline(ref_core, color="#8c1d1d", lw=0.9, ls="--")
    for (_, r), yv in zip(tn.iterrows(), [1.0, 0.55]):
        ax[3].plot([r.col6_per_antxr2 * lo, r.col6_per_antxr2 * hi], [yv, yv],
                   color=FOC, lw=2.6, solid_capstyle="butt")
        ax[3].plot([r.col6_per_antxr2], [yv], marker="s", ms=5, mfc="white",
                   mec=GREY, mew=1.0)
        ax[3].text(r.col6_per_antxr2 * hi + 1.5, yv,
                   r.tissue.replace("tendon, ", "").capitalize() + " tendon",
                   va="center", fontsize=6.5)
    ax[3].set_ylim(0.25, 1.35)
    ax[3].set_yticks([])
    ax[3].set_xlim(0, 72)
    ax[3].set_xlabel("COL6 per $\\it{ANTXR2}$ (CPM ratio)")
    ax[3].set_title("Corrected, tendon cannot be placed", loc="left")
    ax[3].text(ref_max - 1.2, 0.33, "never-reported\nceiling", ha="right", va="bottom",
               color=GREY, fontsize=6)
    ax[3].text(ref_core + 1.2, 1.28, "core-affected\nmedian", ha="left", va="top",
               color="#8c1d1d", fontsize=6)
    ax[3].text(4, 0.70, "open square =\nuncorrected", fontsize=5.8, color=GREY,
               va="center")

    for a, l in zip(ax, "abcd"):
        panel_letter(a, l)
    out = HFS_FIGURES_DIR / "fig_suspension_bias_control.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    return out


def main():
    BI = build()
    BI.to_csv(HFS_OUTPUT_DIR / "hfs_suspension_bias_matched_controls.csv", index=False)
    T = pd.read_csv(HFS_OUTPUT_DIR / "hfs_col6_specificity_by_tissue.csv")
    out = figure(BI, T)
    print(f"{len(BI)} matched strata; ANTXR2 shift percentile "
          f"{BI.antxr2_pctile.min():.0f}-{BI.antxr2_pctile.max():.0f}; "
          f"load understated {BI.load_understated_x.min():.1f}-"
          f"{BI.load_understated_x.max():.1f}x\n-> {out}")
    return BI


if __name__ == "__main__":
    main()
