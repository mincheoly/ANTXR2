"""The gut is the inverse configuration, not a marginal anomaly.

Across tissues the clearance-arm reading is "high collagen VI per receptor =
high load". The intestine is a core-affected organ that sits at the *opposite*
corner of that plane: it has among the highest ANTXR2 of any tissue surveyed
and among the lowest collagen VI, giving it the lowest load-per-receptor of all
populations. Whatever the intestinal phenotype is, the collagen VI clearance
account does not explain it.

Panel b asks the follow-up question — where in the gut the receptor actually
sits — and answers stroma and glia, not epithelium, with every epithelial label
at the bottom of the ranking. That is in tension with the published mouse work
placing the Wnt role in epithelium cell-autonomously, and is the observation
that motivates the inflamed-vs-healthy epithelial arm: the receptor's Wnt
function is injury-conditional, so a healthy atlas cannot show it by
construction.

Outputs
-------
    hfs_gut_antxr2_by_compartment.csv   per-cell-type gut ranking with compartment
    hfs_antxr2_abundance_ranking.csv    whole-cell tissues ranked by receptor CPM
    fig_gut_anomaly.png
"""
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from hfs_config import (COLLAGEN_GROUPS, COLLAGEN_RE, HFS_PB_DIR, HFS_OUTPUT_DIR,
                        HFS_FIGURES_DIR, MIN_CELLS_PER_DONOR_GROUP)

# compartment assignment, first match wins; anything unmatched is immune / other
CLS = [("epithelial", r"enterocyte|epithelial cell|goblet|paneth|stem cell|tuft|"
                      r"enteroendocrine|transit amplif|BEST4|^M cell"),
       ("stromal / mural", r"stromal|myofibroblast|fibroblast|smooth muscle|pericyte|mesotheli"),
       ("endothelial", r"endothelial"),
       ("glial / neural", r"glial|neuron|neural")]

COL = {"stromal / mural": "#b2182b", "epithelial": "#2166ac", "endothelial": "#e08214",
       "glial / neural": "#5aae61", "immune / other": "#999999"}

MARK = {"affected (core)": ("o", 58, "#b2182b"), "occasional": ("^", 34, "#ef8a62"),
        "rare": ("s", 30, "#f4a582"), "not reported": ("D", 26, "#bbbbbb")}


def compartment(label):
    for name, pat in CLS:
        if re.search(pat, label, re.I):
            return name
    return "immune / other"


def load_celltype_cpm(name):
    """Per-cell-type CPM for ANTXR2, COL6 and all collagen from a pseudobulk npz."""
    Z = np.load(HFS_PB_DIR / f"{name}.npz", allow_pickle=False)
    genes = np.array([str(x) for x in Z["genes"]])
    key = "groups" if "groups" in Z else "celltypes"   # older builds used celltypes
    labels = [str(x) for x in Z[key]]
    S = Z["sum"].astype(np.float64)
    cpm = S / S.sum(1, keepdims=True) * 1e6
    gi = {g: i for i, g in enumerate(genes)}
    c6 = [gi[g] for g in COLLAGEN_GROUPS["COL6"] if g in gi]
    ca = [i for i, g in enumerate(genes) if re.fullmatch(COLLAGEN_RE, g)]
    return pd.DataFrame(dict(celltype=labels, ncell=Z["ncell"],
                             ANTXR2=cpm[:, gi["ANTXR2"]],
                             COL6=cpm[:, c6].sum(1),
                             COL_all=cpm[:, ca].sum(1))
                        ).assign(load=lambda d: d.COL6 / d.ANTXR2)


def figure(gut, wc):
    try:
        from figure_style import apply_figure_style, panel_letter
        apply_figure_style(frame="open", sizes=(8, 7, 6))
    except ImportError:
        def panel_letter(ax, l):
            ax.text(-0.13, 1.06, l, transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig, ax = plt.subplots(1, 2, figsize=(10.6, 5.0))

    # (a) receptor vs substrate across tissues
    for inv, (m, s, c) in MARK.items():
        d = wc[wc.involvement == inv]
        ax[0].scatter(d.ANTXR2, d.col6_cpm, marker=m, s=s, facecolor=c,
                      edgecolor="white", linewidth=0.5, label=inv, zorder=3)
    lab = {"Large_Intestine": (-6, -16, "right"), "Small_Intestine": (8, -14, "left"),
           "Skin": (-8, 10, "right"), "Gingiva*": (-8, 6, "right"),
           "synovium (JIA)": (8, 8, "left"), "Pancreas": (6, 10, "left"),
           "Trachea": (7, -12, "left")}
    for t, (dx, dy, ha) in lab.items():
        r = wc[wc.tissue == t]
        if len(r):
            r = r.iloc[0]
            ax[0].annotate(t.replace("_", " ").replace("*", ""), (r.ANTXR2, r.col6_cpm),
                           textcoords="offset points", xytext=(dx, dy), ha=ha,
                           fontsize=6.5, color="#333333")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("$\\it{ANTXR2}$ (CPM)")
    ax[0].set_ylabel("collagen VI (CPM)")
    ax[0].set_title("Gut sits opposite skin and gingiva:\nreceptor-rich, collagen VI-poor",
                    loc="left")
    ax[0].legend(frameon=False, fontsize=6, loc="lower left", handletextpad=0.3,
                 borderpad=0.2)
    panel_letter(ax[0], "a")

    # (b) where the receptor sits within the gut
    Gs = gut.sort_values("ANTXR2")
    yy = np.arange(len(Gs))
    ax[1].barh(yy, Gs.ANTXR2, color=[COL[c] for c in Gs.cls], height=0.72)
    ax[1].set_yticks(yy)
    ax[1].set_yticklabels(Gs.celltype, fontsize=5.4)
    ax[1].set_xlabel("$\\it{ANTXR2}$ (CPM), human intestine")
    ax[1].set_title("Within gut, the receptor is stromal —\nepithelium is the poorest "
                    "compartment", loc="left")
    ax[1].set_ylim(-0.8, len(Gs) - 0.2)
    panel_letter(ax[1], "b")

    fig.tight_layout()
    ax[1].legend(handles=[mpl.patches.Patch(facecolor=COL[n], label=n) for n, _ in CLS]
                 + [mpl.patches.Patch(facecolor=COL["immune / other"], label="immune / other")],
                 frameon=False, fontsize=6, loc="lower right", handletextpad=0.4,
                 borderpad=0.2, handlelength=1.1)

    out = HFS_FIGURES_DIR / "fig_gut_anomaly.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    return out


def main():
    gut = load_celltype_cpm("gut_pseudobulk")
    gut = gut[gut.ncell >= MIN_CELLS_PER_DONOR_GROUP].copy()
    gut["cls"] = gut.celltype.map(compartment)
    gut = gut.sort_values("ANTXR2", ascending=False).reset_index(drop=True)
    gut.insert(0, "rank", np.arange(1, len(gut) + 1))

    T = pd.read_csv(HFS_OUTPUT_DIR / "hfs_col6_specificity_by_tissue.csv")
    wc = T[T.susp == "cell"].sort_values("ANTXR2", ascending=False).copy()
    wc.insert(0, "antxr2_rank", np.arange(1, len(wc) + 1))

    gut.to_csv(HFS_OUTPUT_DIR / "hfs_gut_antxr2_by_compartment.csv", index=False)
    wc[["antxr2_rank", "tissue", "involvement", "susp", "assay", "n_cells", "donors",
        "ANTXR2", "col6_cpm", "col6_per_antxr2"]].to_csv(
        HFS_OUTPUT_DIR / "hfs_antxr2_abundance_ranking.csv", index=False)
    out = figure(gut, wc)

    epi = gut[gut.cls == "epithelial"]
    print(f"gut: {len(gut)} cell types >= {MIN_CELLS_PER_DONOR_GROUP} cells; "
          f"top = {gut.celltype.iloc[0]} ({gut.ANTXR2.iloc[0]:.0f} CPM); "
          f"epithelial median {epi.ANTXR2.median():.0f} CPM, "
          f"best epithelial rank {int(epi['rank'].min())}/{len(gut)}\n-> {out}")
    return gut, wc


if __name__ == "__main__":
    main()
