"""Figure 1 panels a-d from fig1_tables.py outputs (rebuilt 2026-11-22 on this
host from the atlas h5ads; replaces the placeholders in paper_figures_draft.py).

a  ANTXR2 across Tabula Sapiens cell types (10x 3' v3), donor-equal-weight;
   right: compartment level relative to stroma on three platforms
b  fibroblast ANTXR2 by tissue, in two comparability classes
c  receptor/substrate and paralog ratios + MRC2 on the same rows
d  ANTXR2 vs MRC2 detection across gut cell types (IBDverse, GCA)
e  gut differentiation ladder and within-population maturation (from the landscape outputs)

    -> /data/antxr2/fig1/panels/fig1.png
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import fig1_tables as T

warnings.filterwarnings("ignore")
TAB = T.TAB
OUT = T.ROOT / "panels"
OUT.mkdir(exist_ok=True)

INK, INK2, GRID, SURF, NEUTRAL = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#a3a29c"
COMP_COL = {"Stromal": "#2a78d6", "Immune": "#eb6834", "Epithelium": "#1baf7a", "Endothelium": "#eda100"}
plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "legend.frameon": False, "font.family": "DejaVu Sans",
})
LETTERS = []


def letter(ax, s):
    LETTERS.append((ax, s))


def place_letters(fig):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, s in LETTERS:
        tb = ax.get_tightbbox(r).transformed(inv)
        top = ax.title.get_window_extent(r).transformed(inv).y1 if ax.get_title(loc="left") else tb.y1
        fig.text(tb.x0 - 0.01, top + 0.004, s, fontsize=10, fontweight="bold", ha="left", va="bottom")


def grid(ax, axis="x"):
    ax.grid(axis=axis, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


LABEL_1A = {"neutrophil": "neutrophil", "fibroblast": "fibroblast",
            "enterocyte of epithelium proper of ileum": "ileal enterocyte",             "pericyte": "pericyte", "endothelial cell": "endothelial cell", "basal cell": "basal cell",
            "corneal epithelial cell": "corneal epithelial cell"}


def panel_a(ax, axc):
    t = pd.read_csv(TAB / "1a_celltypes.csv")
    t = t[t.cp10k_ANTXR2 > 0]
    x = np.arange(1, len(t) + 1)
    for comp, d in t.groupby(t.compartment.where(t.compartment.isin(COMP_COL), "Other")):
        c = COMP_COL.get(comp, NEUTRAL)
        xi = x[d.index]
        thick = ~d.thin.values
        ax.scatter(xi[thick], d.cp10k_ANTXR2[thick], s=11, color=c, lw=0, zorder=3,
                   label=comp if comp != "Other" else "Neural / germline")
        ax.scatter(xi[~thick], d.cp10k_ANTXR2[~thick], s=11, facecolor=SURF, edgecolor=c, lw=0.7, zorder=3)
    ax.set_yscale("log")
    ax.set_xlim(-2, len(t) + 3)
    lab = t[t.cell_type.isin(LABEL_1A) & ~t.thin]
    for i, r in lab.iterrows():
        last = r.cell_type == "corneal epithelial cell"
        ax.annotate(LABEL_1A[r.cell_type], (x[i], r.cp10k_ANTXR2), xytext=(-8 if last else 5, 13),
                    textcoords="offset points", fontsize=5.8, color=INK2, ha="right" if last else "left", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=NEUTRAL, shrinkA=0, shrinkB=2))
    ax.set_xlabel("Cell types ranked by ANTXR2 (Tabula Sapiens, 10x 3′ v3, 24 donors)")
    ax.set_ylabel("ANTXR2, donor-mean CP10k")
    ax.set_title(f"ANTXR2 is detected in {(t.cp10k_ANTXR2 > 0.05).sum()} of {len(t)} cell types", loc="left")
    grid(ax, "y")
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([], [], marker="o", ls="", mfc=SURF, mec=INK2, mew=0.7, ms=3.5))
    l.append(f"< {T.THIN_DONORS} donors")
    ax.legend(h, l, loc="lower left", ncol=2, handletextpad=0.2, columnspacing=0.8)

    # compartment / stroma on three platforms
    u = T.load()
    rows = []
    for name, d, key in [("TS 10x 3′ v3\n(intron-incl.)", u[(u.atlas == "tabula_sapiens") & (u.assay == "10x 3' v3")], "x0"),
                         ("TS Smart-seq2\n(full-length)", u[(u.atlas == "tabula_sapiens") & (u.assay == "Smart-seq2")], "x0"),
                         ("GCA adult gut\n(exon-only)", u[(u.atlas == "gca_healthy") & u.x1.isin(T.ADULT)], "x2")]:
        s, _ = T.summarise(d, [key])
        s = s.set_index(key).cp10k_ANTXR2
        s = s.rename({"Mesenchymal": "Stromal", "Myeloid": "Immune", "Epithelial": "Epithelium",
                      "Endothelial": "Endothelium"})
        ref = s["Stromal"]
        for comp in ["Immune", "Epithelium", "Endothelium"]:
            rows.append((name, comp, np.log2(s[comp] / ref)))
    r = pd.DataFrame(rows, columns=["platform", "comp", "lr"])
    plats = r.platform.unique()
    for j, comp in enumerate(["Immune", "Epithelium", "Endothelium"]):
        d = r[r.comp == comp].set_index("platform").reindex(plats)
        yy = np.arange(len(plats)) + (j - 1) * 0.22
        axc.hlines(yy, 0, d.lr, color=COMP_COL[comp], lw=1.6, zorder=3)
        axc.scatter(d.lr, yy, s=14, color=COMP_COL[comp], zorder=4, lw=0)
    axc.axvline(0, color=COMP_COL["Stromal"], lw=1.0)
    axc.set_yticks(range(len(plats)))
    axc.set_yticklabels(plats)
    axc.invert_yaxis()
    axc.set_xlabel("log2 (compartment / stroma)\nstroma = 0 (blue line)")
    axc.set_title("Immune rank depends\non platform", loc="left")
    grid(axc)


ROW_RENAME = {"IBDverse ileum": "Ileum (IBDverse)", "IBDverse rectum": "Rectum (IBDverse)",
              "Oral mucosa": "Oral mucosa / gingiva"}
METRICS = [("cp10k_ANTXR2", None, "ANTXR2 (CP10k)", "b"),
           ("gm_col6_antxr2", "lr_col6_antxr2", "COL6A1-3 / ANTXR2", "c"),
           ("gm_antxr1_antxr2", "lr_antxr1_antxr2", "ANTXR1 / ANTXR2", ""),
           ("cp10k_MRC2", None, "MRC2 (CP10k)", "")]


def panel_bc(fig, gs):
    t = pd.read_csv(TAB / "1bc_fibroblasts.csv")
    don = pd.read_csv(TAB / "1bc_fibroblasts_donors.csv")
    blocks = [("intron", "intron-inclusive"), ("exon", "exon-only")]
    sizes = [len(t[t.cls == c]) for c, _ in blocks]
    sub = gs.subgridspec(2, 4, height_ratios=sizes, hspace=0.18, wspace=0.12)
    axes = np.empty((2, 4), object)
    for bi, (cls, title) in enumerate(blocks):
        tb = t[t.cls == cls].sort_values(["atlas", "cp10k_ANTXR2"], ascending=[False, True]).reset_index(drop=True)
        tb = tb.iloc[::-1].reset_index(drop=True)
        for mi, (col, lr, xl, let) in enumerate(METRICS):
            ax = fig.add_subplot(sub[bi, mi], sharex=axes[0, mi] if bi else None)
            axes[bi, mi] = ax
            y = np.arange(len(tb))
            for i, r in tb.iterrows():
                d = don[(don.atlas == r.atlas) & (don.group == r.group)]
                v = 2 ** d[lr] if lr else d[col]
                ax.scatter(v, np.full(len(v), i) + np.random.default_rng(i).uniform(-0.18, 0.18, len(v)),
                           s=3, color=NEUTRAL, lw=0, alpha=0.6, zorder=2)
            ax.scatter(tb[col], y, s=20, color=INK, zorder=4, lw=0, marker="D" if True else "o")
            ax.set_xscale("log")
            ax.set_yticks(y)
            ax.set_ylim(len(tb) - 0.5, -0.5)
            if mi == 0:
                labels = [f"{ROW_RENAME.get(g, g)}  ({n})" for g, n in zip(tb.group, tb.n_donors)]
                ax.set_yticklabels(labels)
                for tl, h in zip(ax.get_yticklabels(), tb.hfs):
                    if h:
                        tl.set_fontweight("bold")
                yb = -0.25
                ax.text(0.057, yb, title, fontsize=6, color=INK2, style="italic", va="center", ha="left")
            else:
                ax.set_yticklabels([])
            # separate atlases inside a block
            brk = np.where(tb.atlas.values[1:] != tb.atlas.values[:-1])[0]
            for b in brk:
                ax.axhline(b + 0.5, color=INK2, lw=0.4, ls=(0, (2, 2)))
            grid(ax)
            if bi == 0:
                ax.set_title(xl)
                ax.tick_params(labelbottom=False)
                if let:
                    letter(ax, let)
    for mi in range(4):
        ax = axes[1, mi]
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    axes[0, 0].figure.text(0, 0, "")
    return axes


GUT_ORDER = {"ibdverse": ["Stem", "Enterocyte", "Colonocyte", "Secretory", "Myeloid", "Mesenchymal"],
             "gca_healthy": ["Epithelial_stem", "Transit_amplifying", "Absorptive", "Secretory", "Enteroendocrine",
                             "Macrophage", "Fibroblast"]}
GUT_LABEL = {"Epithelial_stem": "Stem", "Transit_amplifying": "Transit amplifying", "Mesenchymal": "Fibroblast",
             "Myeloid": "Myeloid", "Macrophage": "Macrophage"}


def panel_d(axs):
    t = pd.read_csv(TAB / "1d_gut_mrc2.csv")
    for ax, atlas, title in zip(axs, ["ibdverse", "gca_healthy"],
                                ["IBDverse, normal ileum + rectum\n(intron-incl.)",
                                 "Gut Cell Atlas, adult intestine\n(exon-only)"]):
        d = t[t.atlas == atlas].set_index("group").reindex(GUT_ORDER[atlas])
        y = np.arange(len(d))
        a, m = d.det_ANTXR2 * 100, d.det_MRC2 * 100
        zero = m < 0.06
        m = m.clip(lower=0.06)
        ax.hlines(y, np.minimum(a, m), np.maximum(a, m), color=GRID, lw=1.6, zorder=2)
        ax.scatter(a, y, s=22, color=INK, zorder=4, lw=0, label="ANTXR2")
        ax.scatter(m, y, s=22, facecolor=SURF, edgecolor=INK, lw=0.9, zorder=4, label="MRC2")
        for yi in y[zero.values]:
            ax.text(0.075, yi - 0.25, "0", fontsize=5.5, color=INK2, ha="left", va="bottom")
        ax.set_xscale("log")
        ax.set_xlim(0.05, 120)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{GUT_LABEL.get(g, g)}  ({n:.0f})" for g, n in zip(d.index, d.n_donors)])
        ax.invert_yaxis()
        n_epi = sum(g not in ("Mesenchymal", "Myeloid", "Macrophage", "Fibroblast") for g in d.index)
        ax.axhspan(-0.5, n_epi - 0.5, color="#f0efec", zorder=0)
        ax.text(1.02, (n_epi - 1) / 2, "epithelium", fontsize=6, color=INK2, rotation=270, va="center", ha="left",
                transform=matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData))
        ax.set_title(title, loc="left", fontsize=6.8, color=INK2)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        grid(ax)
        ax.set_xlabel("% cells with ≥1 UMI (donor mean)")
    axs[0].legend(loc="lower left", ncol=1, handletextpad=0.2, borderaxespad=0.2)


LANDSCAPE = Path("/data/antxr2/program_coexpr/output/landscape/differentiation")
STATES = ["normal", "cd_uninflamed", "cd_inflamed"]
STATE_LABEL = {"normal": "Normal", "cd_uninflamed": "CD uninflamed", "cd_inflamed": "CD inflamed"}
STATE_COL = {"normal": "#86b6ef", "cd_uninflamed": "#2a78d6", "cd_inflamed": "#104281"}
POP_SHORT = {
    "Enterocyte top villus": "Top-villus enterocyte", "Enterocyte_top_villus": "Top-villus enterocyte",
    "Enterocyte progenitor/precursor crypt OLFM4lo": "Crypt OLFM4lo",
    "Enterocyte_progenitor_precursor_crypt_OLFM4lo": "Crypt OLFM4lo",
    "Enterocyte progenitor crypt OLFM4hi": "Crypt OLFM4hi",
    "Enterocyte stem MKI67": "Stem MKI67", "Enterocyte_stem_MKI67": "Stem MKI67",
    "Enterocyte stem OLFM4lo LGR5hi": "Stem LGR5hi", "Enterocyte_stem_OLFM4lo_LGR5hi": "Stem LGR5hi",
}
OFF = {"normal": -0.22, "cd_uninflamed": 0, "cd_inflamed": 0.22}


def panel_e(ax, axw):
    """Gut absorptive ladder (IBDverse ileum): ANTXR2 excess over matched genes vs LGR5hi stem;
    right: within-population maturation (top vs bottom quartile) as a matched-gene percentile."""
    d = pd.read_csv(LANDSCAPE / "ladder.csv")
    d = d[d.branch == "absorptive"].copy()
    order = d.groupby("population").excess.mean().sort_values().index.tolist()
    for st in STATES:
        q = d[d.state == st].set_index("population").reindex(order)
        y = np.arange(len(order)) + OFF[st]
        ax.errorbar(q.excess, y, xerr=q.anchor_sd, fmt="none", ecolor=STATE_COL[st], elinewidth=0.8, zorder=2)
        ax.scatter(q.excess, y, s=22, color=STATE_COL[st], edgecolor=SURF, lw=0.6, zorder=3, label=STATE_LABEL[st])
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([POP_SHORT.get(p, p) for p in order])
    ax.axvline(0, color=INK2, lw=0.6, zorder=1)
    grid(ax)
    ax.set_xlabel("ANTXR2 log2 change vs LGR5hi stem,\nminus matched genes (± matched-gene SD)")
    ax.set_title("Down in dividing progenitors, up with differentiation", loc="left")
    ax.legend(loc="lower right", handletextpad=0.2, borderaxespad=0.2)

    w =pd.read_csv(LANDSCAPE / "within.csv")
    w = w[w.population.map(lambda p: p in POP_SHORT)]
    pops = w.groupby("population").excess.mean().sort_values().index.tolist()
    for st in STATES:
        q = w[w.state == st].set_index("population").reindex(pops)
        axw.scatter(q.antxr2_pct, np.arange(len(pops)) + OFF[st], s=22, color=STATE_COL[st], edgecolor=SURF,
                    lw=0.6, zorder=3)
    axw.set_yticks(range(len(pops)))
    axw.set_yticklabels([POP_SHORT[p] for p in pops])
    axw.axvline(50, color=INK2, lw=0.6, zorder=1)
    axw.set_xlim(0, 100)
    grid(axw)
    axw.set_xlabel("Within population: matched-gene\npercentile, top vs bottom maturation")
    axw.set_title("Within populations: holds,\nexcept inflamed top villus", loc="left")


FOOT_BC = ("Bold rows: tissues with HFS involvement fixed a priori (skin, gingiva, intestine).  Grey dots: donors; "
           "black: donor mean; (n) donors.\nRatios: per-donor pooled UMI ratios, geometric mean.  "
           "Values are not comparable between the two blocks.\n"
           "Intron-inclusive: Tabula Sapiens (Cell Ranger 7.0.1) and IBDverse (Cell Ranger 7.2.0), which count "
           "intronic reads by default.\nExon-only: Gut Cell Atlas (released STARsolo Gene counts) and the dermal "
           "fibroblast atlas (no long-gene excess vs GCA).")


def main():
    fig = plt.figure(figsize=(7.2, 12.2))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.75, 0.82, 0.78], hspace=0.55,
                          left=0.22, right=0.965, top=0.975, bottom=0.045)
    top = gs[0].subgridspec(1, 2, width_ratios=[2.6, 1], wspace=0.55)
    ax_a, ax_c = fig.add_subplot(top[0]), fig.add_subplot(top[1])
    panel_a(ax_a, ax_c)
    letter(ax_a, "a")
    bc = panel_bc(fig, gs[1])
    mid = gs[2].subgridspec(1, 2, wspace=0.75)
    axd = [fig.add_subplot(mid[0]), fig.add_subplot(mid[1])]
    panel_d(axd)
    letter(axd[0], "d")
    bot = gs[3].subgridspec(1, 2, width_ratios=[1.5, 1], wspace=0.75)
    axe = [fig.add_subplot(bot[0]), fig.add_subplot(bot[1])]
    panel_e(*axe)
    letter(axe[0], "e")
    fig.canvas.draw()
    y0 = bc[1, 0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted()).y0
    fig.text(0.04, y0 - 0.005, FOOT_BC, fontsize=5.6, color=INK2, va="top")
    place_letters(fig)
    fig.savefig(OUT / "fig1.png", dpi=220)
    print(f"-> {OUT / 'fig1.png'}")


if __name__ == "__main__":
    main()
