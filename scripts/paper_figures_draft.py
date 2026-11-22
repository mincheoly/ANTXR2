"""Draft paper figures (user, 2026-11-21: "show me drafts of the figure you want
here"). Layout follows planning_summaries/paper_story.md §3: three acts.

Figure 1 panels a-d need the cross-tissue cell-type means, which are no longer
on this host (/data/ANTXR2/hfs_atlas and celltype_expression are empty); they
are drawn as labelled placeholders. 1e and all of Figures 2-3 use real outputs.

Colors: reference palette of the dataviz skill (states = ordinal blue ramp from
step 250; categorical slots 1-4 in fixed order; diverging blue <-> red with a
gray midpoint). The palette validator could not be run (no JS runtime on host).

    -> OUT_DIR/paper_figures/fig{1,2,3}_draft.png

Usage: python paper_figures_draft.py
"""
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

O = Path("/data/antxr2/program_coexpr/output")
OUT = O / "paper_figures"
OUT.mkdir(exist_ok=True)

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
NEUTRAL = "#a3a29c"
STATES = ["normal", "cd_uninflamed", "cd_inflamed"]
STATE_LABEL = {"normal": "Normal", "cd_uninflamed": "CD uninflamed", "cd_inflamed": "CD inflamed"}
STATE_COL = {"normal": "#86b6ef", "cd_uninflamed": "#2a78d6", "cd_inflamed": "#104281"}
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # slots 1-4
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#f0efec", "#e34948"])

plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "legend.frameon": False, "font.family": "DejaVu Sans",
})

POP_SHORT = {
    "fibroblast_pooled": "Fibroblasts", "fibroblast (pooled)": "Fibroblasts",
    "Enterocyte_top_villus": "Top-villus enterocyte", "Enterocyte top villus": "Top-villus enterocyte",
    "Enterocyte_progenitor_precursor_crypt_OLFM4lo": "Crypt OLFM4lo",
    "Enterocyte progenitor/precursor crypt OLFM4lo": "Crypt OLFM4lo",
    "Enterocyte progenitor crypt OLFM4hi": "Crypt OLFM4hi",
    "Ileum_goblet_MUC2hi_BCAS1hi": "Goblet MUC2hi", "Ileum goblet MUC2hi BCAS1hi": "Goblet MUC2hi",
    "Ileum goblet MUC2hi SPINK4hi": "Goblet SPINK4hi", "Ileum goblet MUC2int": "Goblet MUC2int",
    "Enterocyte_stem_MKI67": "Stem MKI67", "Enterocyte stem MKI67": "Stem MKI67",
    "Enterocyte_stem_OLFM4lo_LGR5hi": "Stem LGR5hi", "Enterocyte stem OLFM4lo LGR5hi": "Stem LGR5hi",
}
IFN_POPS = ["fibroblast_pooled", "Enterocyte_top_villus", "Enterocyte_progenitor_precursor_crypt_OLFM4lo",
            "Ileum_goblet_MUC2hi_BCAS1hi", "Enterocyte_stem_MKI67", "Enterocyte_stem_OLFM4lo_LGR5hi"]


LETTERS = []


def letter(ax, s, **_):
    LETTERS.append((ax, s))


def place_letters(fig):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, s in LETTERS:
        tb = ax.get_tightbbox(r).transformed(inv)
        top = ax.title.get_window_extent(r).transformed(inv).y1 if ax.get_title(loc="left") else tb.y1
        fig.text(tb.x0 - 0.02, top + 0.004, s, fontsize=10, fontweight="bold", ha="left", va="bottom")
    LETTERS.clear()


def legend_below(ax, dy=-0.30, ncol=3, **kw):
    ax.legend(loc="upper left", bbox_to_anchor=(0, dy), ncol=ncol, columnspacing=0.8, handletextpad=0.3,
              borderaxespad=0, **kw)


def grid(ax, axis="x"):
    ax.grid(axis=axis, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def zero(ax, axis="x", v=0):
    (ax.axvline if axis == "x" else ax.axhline)(v, color=INK2, lw=0.6, ls="-", zorder=1)


def placeholder(ax, title, body):
    ax.set_axis_off()
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc="#f0efec", ec=NEUTRAL, lw=0.8, ls="--"))
    ax.text(0.5, 0.62, title, transform=ax.transAxes, ha="center", va="center", fontsize=7.5, fontweight="bold")
    ax.text(0.5, 0.36, body, transform=ax.transAxes, ha="center", va="center", fontsize=6.3, color=INK2,
            wrap=True)


def state_legend(ax, loc="lower right", below=None, **kw):
    for s in STATES:
        ax.scatter([], [], s=22, color=STATE_COL[s], label=STATE_LABEL[s])
    if below is not None:
        legend_below(ax, dy=below, **kw)
    else:
        ax.legend(loc=loc, handletextpad=0.2, borderaxespad=0.2, **kw)


# ---------------------------------------------------------------- Figure 1
def fig1():
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.15], hspace=0.65, wspace=0.75, top=0.92)
    ph = [
        ("ANTXR2 across 238 cell types", "Donor-equal-weight means, 3 collections;\nthin labels (n_donors) marked.\n"
         "Data to regenerate: compute_means.py\n(/data/ANTXR2 is empty on this host)"),
        ("Fibroblast expression across tissues", "Gut vs skin vs other stroma;\nexon-only vs intron-counted\natlases flagged"),
        ("Ratios, raw ingredients", "COL6 / ANTXR2, ANTXR1 / ANTXR2;\nno composite score"),
        ("MRC2 in gut epithelium", "At detection floor: the one clean\nabsence (partner-availability)"),
    ]
    axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 0:2]),
           fig.add_subplot(gs[1, 2])]
    for ax, (t, b), L in zip(axs, ph, "abcd"):
        placeholder(ax, t, b)
        letter(ax, L, dx=-0.02, dy=1.02)

    # 1e: differentiation ladder, excess over matched genes vs LGR5hi stem
    d = pd.read_csv(O / "landscape/differentiation/ladder.csv")
    d = d[d.branch == "absorptive"].copy()
    order = d.groupby("population").excess.mean().sort_values().index.tolist()
    ax = fig.add_subplot(gs[2, 0:2])
    off = {"normal": -0.22, "cd_uninflamed": 0, "cd_inflamed": 0.22}
    for s in STATES:
        q = d[d.state == s].set_index("population").reindex(order)
        y = np.arange(len(order)) + off[s]
        ax.errorbar(q.excess, y, xerr=q.anchor_sd, fmt="none", ecolor=STATE_COL[s], elinewidth=0.8, zorder=2)
        ax.scatter(q.excess, y, s=22, color=STATE_COL[s], edgecolor=SURF, lw=0.6, zorder=3)
    ax.set_yticks(range(len(order)), [POP_SHORT.get(p, p) for p in order])
    zero(ax)
    grid(ax)
    ax.set_xlabel("ANTXR2 log2 change vs LGR5hi stem, minus matched genes (± matched-gene SD)")
    ax.set_title("ANTXR2 falls in dividing progenitors and rises with differentiation", loc="left")
    state_legend(ax, loc="upper left")
    letter(ax, "e", dx=-0.3)

    ax = fig.add_subplot(gs[2, 2])
    w = pd.read_csv(O / "landscape/differentiation/within.csv")
    w = w[w.population.map(lambda p: p in POP_SHORT)]
    pops = w.groupby("population").excess.mean().sort_values().index.tolist()
    for s in STATES:
        q = w[w.state == s].set_index("population").reindex(pops)
        ax.scatter(q.antxr2_pct, np.arange(len(pops)) + off[s], s=22, color=STATE_COL[s], edgecolor=SURF,
                   lw=0.6, zorder=3)
    ax.set_yticks(range(len(pops)), [POP_SHORT[p] for p in pops])
    zero(ax, v=50)
    ax.set_xlim(0, 100)
    grid(ax)
    ax.set_xlabel("Within population: matched-gene\npercentile, top vs bottom maturation")
    ax.set_title("Same direction within populations", loc="left")
    fig.suptitle("Figure 1 — ANTXR2 is broadly expressed and fits the known mechanism, but varies by cell state",
                 x=0.02, ha="left", fontsize=8.5, fontweight="bold")
    place_letters(fig)
    fig.savefig(OUT / "fig1_draft.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def ifn_units():
    f = O / "landscape_ref/test_memento/Enterocyte_top_villus@normal__estimate__ifn__cells.csv"
    return pd.read_csv(f).unit.tolist()


def fig2():
    fig = plt.figure(figsize=(7.2, 10.4))
    gs = fig.add_gridspec(4, 6, height_ratios=[1, 1.55, 1, 1], hspace=0.95, wspace=2.2, top=0.94)
    ST = {"normal": "N", "cd_uninflamed": "U", "cd_inflamed": "I"}

    # 2a: growth axis survives removing MKI67+ cells
    ax = fig.add_subplot(gs[0, 0:3])
    c = pd.read_csv(O / "landscape/cycle_compare_cells.csv")
    c = c[c.lineage == "epithelial"]
    units = {"hallmark:E2F_TARGETS": "E2F targets", "hallmark:MYC_TARGETS_V1": "MYC targets V1",
             "hallmark:G2M_CHECKPOINT": "G2M checkpoint"}
    for i, u in enumerate(units):
        q = c[(c.unit == u)].dropna(subset=["main", "drop"])
        for _, r in q.iterrows():
            ax.plot([i - 0.18, i + 0.18], [r["main"], r["drop"]], color=STATE_COL[r["state"]], lw=0.7, alpha=0.8)
        cols_ = [STATE_COL[x] for x in q.state]
        ax.scatter(np.full(len(q), i - 0.18), q["main"], s=12, color=cols_, zorder=3)
        ax.scatter(np.full(len(q), i + 0.18), q["drop"], s=12, color=cols_, zorder=3, marker="D")
    ax.set_xticks(range(3), list(units.values()))
    ax.set_xlim(-0.5, 2.5)
    zero(ax, "y")
    grid(ax, "y")
    ax.set_ylabel("Coupling, z vs matched genes")
    ax.set_title("Growth coupling survives removing MKI67+ cells\n(left: all cells · right: MKI67+ removed)",
                 loc="left")
    state_legend(ax, below=-0.2)
    letter(ax, "a")

    # 2b: split-half reliability, three atlases
    ax = fig.add_subplot(gs[0, 3:6])
    rows = []
    for atlas, p in [("IBD", "landscape_ref"), ("GCA", "landscape_gca_ref"), ("Skin", "landscape_skin_ref")]:
        s = pd.read_csv(O / p / "donor_split/summary.csv")
        s["atlas"] = atlas
        rows.append(s)
    s = pd.concat(rows)
    s = s[s.population.str.contains("fibroblast|top villus|Enterocyte$|TA$", regex=True)
          & ~s.population.str.contains("ADAMDEC1")]
    pop = (s.population.str.replace("fibroblast", "fib.").str.replace("Enterocyte top villus", "top-villus")
           .str.replace("Enterocyte", "enterocyte"))
    s["lab"] = s.atlas + " · " + pop + np.where(s.atlas == "IBD", " · " + s.state.map(ST), "")
    s = s.sort_values("r_half_median")
    y = np.arange(len(s))
    col = {"IBD": CAT[0], "GCA": CAT[1], "Skin": CAT[2]}
    ax.hlines(y, s.r_half_q10, s.r_half_q90, color=[col[a] for a in s.atlas], lw=1.2)
    ax.scatter(s.r_half_median, y, s=20, color=[col[a] for a in s.atlas], edgecolor=SURF, lw=0.6, zorder=3)
    ax.set_yticks(y, s.lab, fontsize=5.6)
    zero(ax)
    grid(ax)
    ax.set_xlabel("Donor split-half agreement of the\nprogram profile (median, 10–90%)")
    ax.set_title("Reproducibility varies by atlas and cell type", loc="left")
    for a, lab in [("IBD", "IBDverse"), ("GCA", "Gut Cell Atlas"), ("Skin", "Skin")]:
        ax.scatter([], [], s=20, color=col[a], label=lab)
    legend_below(ax, dy=-0.44)
    letter(ax, "b")

    # 2c: landscape heatmap, anchor percentile - 50
    ax = fig.add_subplot(gs[1, 0:6])
    L = pd.read_csv(O / "landscape_ref/coupling.csv")
    keep_pop = ["Enterocyte stem OLFM4lo LGR5hi", "Enterocyte stem MKI67", "Enterocyte progenitor crypt OLFM4hi",
                "Enterocyte progenitor/precursor crypt OLFM4lo", "Enterocyte top villus",
                "Ileum goblet MUC2hi BCAS1hi", "fibroblast (pooled)"]
    L = L[L.population.isin(keep_pop) & (L.activity_sd >= 0.2)]
    L["col"] = L.population + "|" + L.state
    M = L.pivot_table(index="unit", columns="col", values="pct") - 50
    cols = [f"{p}|{s_}" for p in keep_pop for s_ in STATES]
    M = M.reindex(columns=cols)
    M = M[M.notna().mean(1) >= 0.7]
    ifn = [u for u in ifn_units() if u in M.index]
    ecm = [u for u in ["matrisome:NABA_COLLAGENS", "matrisome:NABA_CORE_MATRISOME"] if u in M.index]
    forced = ifn + ecm + [u for u in ["hallmark:E2F_TARGETS", "hallmark:MYC_TARGETS_V1"] if u in M.index]
    rest = M.drop(index=forced).abs().mean(1).sort_values(ascending=False).index[:14].tolist()
    H = M.loc[forced + rest]
    H = H.loc[H.mean(1).sort_values(ascending=False).index]
    im = ax.imshow(H.values, cmap=DIV, vmin=-50, vmax=50, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(H.shape[1]), [ST[c_.split("|")[1]] for c_ in H.columns], fontsize=5.8)
    for k, p in enumerate(keep_pop):
        ax.text(3 * k + 1, -0.8, POP_SHORT[p], ha="center", va="bottom", fontsize=6.2)
    labels = []
    for u in H.index:
        coll, name = u.split(":", 1)
        labels.append(f"{name.replace('_', ' ').title() if coll in ('hallmark', 'matrisome') else name}"
                      f"  ({coll})")
    ax.set_yticks(range(H.shape[0]), labels, fontsize=5.6)
    for t, u in zip(ax.get_yticklabels(), H.index):
        if u in ifn:
            t.set_fontweight("bold")
        if u in ecm:
            t.set_fontstyle("italic")
    for b in np.arange(3, H.shape[1], 3):
        ax.axvline(b - 0.5, color=SURF, lw=2)
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("ANTXR2 percentile among\n200 matched genes − 50", fontsize=6)
    cb.outline.set_visible(False)
    ax.set_title("Landscape excerpt (IBDverse ileum): 8 IFN programs in bold, ECM in italics; "
                 "N / U / I = normal / CD uninflamed / CD inflamed", loc="left", pad=14)
    letter(ax, "c")

    # 2d: collagen, pooled vs within author subtype (precision-weighted over states)
    ax = fig.add_subplot(gs[2, 0:3])
    genes = ["COL6A1", "COL6A2", "COL6A3", "COL1A1", "COL12A1", "COL4A5"]

    def pooled(dirn):
        fs = sorted(glob.glob(str(O / dirn / "fibroblast_pooled@*__genes__cells.csv")))
        d = pd.concat([pd.read_csv(f).assign(state=f) for f in fs])
        d = d[d.gene.isin(genes)].copy()
        d["w"] = 1 / d.excess_se ** 2
        g = d.groupby("gene").apply(lambda x: pd.Series({
            "excess": np.average(x.excess, weights=x.w), "se": 1 / np.sqrt(x.w.sum())}))
        return g.reindex(genes)

    a, b = pooled("landscape_ref/test_memento"), pooled("landscape_ref/test_memento/by_author_subtype")
    y = np.arange(len(genes))
    for g_, lab, dy, cc in [(a, "Mixed subtypes (donor × class)", -0.16, NEUTRAL),
                            (b, "Within author subtype", 0.16, CAT[0])]:
        ax.errorbar(g_.excess, y + dy, xerr=1.96 * g_.se, fmt="o", ms=4, color=cc, mec=SURF, mew=0.5,
                    elinewidth=0.8, label=lab)
    ax.set_yticks(y, genes)
    ax.invert_yaxis()
    zero(ax)
    grid(ax)
    ax.set_xlabel("Gut fibroblasts: excess correlation\nover matched genes (95% CI)")
    ax.set_title("Negative: collagen signal is subtype mix", loc="left")
    legend_below(ax, dy=-0.36, ncol=1)
    letter(ax, "d")

    # 2e: the bar is beatable: collagen-program percentile for ECM-role genes
    ax = fig.add_subplot(gs[2, 3:6])
    cs = pd.read_csv(O / "landscape_skin_ref/test_memento/fibroblast_A3@normal__controls.csv")
    cg = pd.read_csv(O / "landscape_ref/test_memento/fibroblast_ADAMDEC1_uninflamed@normal__controls.csv")
    u = "matrisome:NABA_COLLAGENS"
    s_ = cs[cs.unit == u].set_index("gene")
    g_ = cg[cg.unit == u].set_index("gene")
    order = s_.pct.sort_values().index.tolist()
    y = np.arange(len(order))
    for d_, lab, dy, cc in [(s_, "Skin A3 (70 donors)", -0.16, CAT[2]), (g_, "Gut ADAMDEC1 (19 donors)", 0.16, CAT[0])]:
        d_ = d_.reindex(order)
        sig = (d_.excess_q < 0.05).values
        ax.scatter(d_.pct, y + dy, s=22, color=cc, edgecolor=SURF, lw=0.5, label=lab, zorder=3)
        ax.scatter(d_.pct[sig], (y + dy)[sig], s=58, facecolor="none", edgecolor=INK, lw=0.7, zorder=4)
    ax.set_yticks(y, order)
    for t in ax.get_yticklabels():
        if t.get_text() == "ANTXR2":
            t.set_fontweight("bold")
    zero(ax, v=50)
    ax.set_xlim(0, 104)
    grid(ax)
    ax.scatter([], [], s=58, facecolor="none", edgecolor=INK, lw=0.7, label="excess q < 0.05")
    legend_below(ax, dy=-0.36, ncol=2)
    ax.set_xlabel("Collagen program: percentile\namong each gene's matched genes")
    ax.set_title("The bar is beatable: ECM-role genes clear it", loc="left")
    letter(ax, "e")

    # 2f: IFN coupling by population x state (8 programs; median and range)
    ax = fig.add_subplot(gs[3, 0:3])
    y0, yt, yl = 0, [], []
    for p in IFN_POPS:
        for k, s in enumerate(STATES):
            f = O / f"landscape_ref/test_memento/{p}@{s}__estimate__ifn__cells.csv"
            if not f.exists():
                continue
            d = pd.read_csv(f)
            yy = y0 + (k - 1) * 0.25
            ax.hlines(yy, d.excess.min(), d.excess.max(), color=STATE_COL[s], lw=1)
            ax.scatter(d.excess.median(), yy, s=22, color=STATE_COL[s], edgecolor=SURF, lw=0.6, zorder=3)
        yt.append(y0)
        yl.append(POP_SHORT[p])
        y0 += 1
    ax.set_yticks(yt, yl)
    ax.invert_yaxis()
    zero(ax)
    grid(ax)
    ax.set_xlabel("Coupling to 8 IFN programs, excess over\nmatched genes (dot: median; line: range)")
    ax.set_title("Positive: IFN coupling rises with inflammation", loc="left")
    state_legend(ax, below=-0.36)
    letter(ax, "f")

    # 2g: contamination filters
    ax = fig.add_subplot(gs[3, 3:6])
    ifn = ifn_units()
    filt = {"all": ("All cells", NEUTRAL), "no_immune_doublet": ("No immune doublets", CAT[0]),
            "low_ambient": ("Low ambient Ig", CAT[1])}
    short = {"fibroblast_pooled": "Fibroblasts", "Enterocyte_top_villus": "Top-villus",
             "Enterocyte_progenitor_precursor_crypt_OLFM4lo": "Crypt OLFM4lo"}
    y0, yt, yl = 0, [], []
    for p in IFN_POPS[:3]:
        for s in STATES:
            d = pd.read_csv(O / f"ifn_type_split/contam/landscape_ref/{p}@{s}.csv")
            d = d[d.unit.isin(ifn)]
            for k, (fk, (lab, cc)) in enumerate(filt.items()):
                v = d[d["filter"] == fk].excess.median()
                ax.scatter(v, y0 + (k - 1) * 0.22, s=18, color=cc, edgecolor=SURF, lw=0.5, zorder=3,
                           label=lab if y0 == 0 else None)
            yt.append(y0)
            yl.append(f"{short[p]} · {ST[s]}")
            y0 += 1
        y0 += 0.4
    ax.set_yticks(yt, yl, fontsize=5.8)
    ax.invert_yaxis()
    zero(ax)
    ax.set_xlim(left=-0.02)
    grid(ax)
    ax.set_xlabel("Median IFN-program excess correlation")
    ax.set_title("Not carried by immune-derived reads", loc="left")
    legend_below(ax, dy=-0.36, ncol=2)
    letter(ax, "g")

    fig.suptitle("Figure 2 — ANTXR2 joins programs by context: collagen a pertinent negative, interferon a "
                 "pertinent positive", x=0.02, ha="left", fontsize=8.5, fontweight="bold", y=0.985)
    place_letters(fig)
    fig.savefig(OUT / "fig2_draft.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def fig3():
    fig = plt.figure(figsize=(7.2, 8.4))
    gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 0.9], hspace=0.85, wspace=2.3, top=0.93)
    stims = {"alpha": ("IFN-α", CAT[0]), "beta": ("IFN-β", CAT[1]), "gamma": ("IFN-γ", CAT[2]),
             "lambda": ("IFN-λ", CAT[3])}

    # 3a: HTEC ciliated mean ratio vs control over time
    ax = fig.add_subplot(gs[0, 0:3])
    m = pd.read_csv(O / "htec_ifn/mean.csv")
    m = m[(m.ct == "ciliated") & (m.stim != "control")]
    g = m.groupby(["stim", "time"]).ratio.agg(["mean", "min", "max"]).reset_index()
    for st, (lab, cc) in stims.items():
        q = g[g.stim == st]
        ax.fill_between(q.time, np.log2(q["min"]), np.log2(q["max"]), color=cc, alpha=0.12, lw=0)
        ax.plot(q.time, np.log2(q["mean"]), color=cc, lw=2, marker="o", ms=4, mec=SURF, mew=0.6)
        nudge = {"alpha": 0.07, "lambda": -0.07}.get(st, 0)
        ax.text(q.time.iloc[-1] + 1.2, np.log2(q["mean"].iloc[-1]) + nudge, lab, color=INK, va="center",
                fontsize=6.3)
    zero(ax, "y")
    grid(ax, "y")
    ax.set_xlabel("Hours after stimulation")
    ax.set_ylabel("ANTXR2 per 10k, log2 vs control\n(line: donor mean; band: 2 donors)")
    ax.set_xticks([3, 6, 9, 24, 48])
    ax.set_xlim(0, 56)
    ax.set_title("Airway ciliated cells: induced by α / β", loc="left")
    letter(ax, "a", dx=-0.22)

    # 3b: HTEC coupling to IFN programs per stimulus
    ax = fig.add_subplot(gs[0, 3:6])
    cpl = pd.read_csv(O / "htec_ifn/coupling_ciliated.csv")
    ifn = ifn_units()
    cpl = cpl[cpl.unit.isin(ifn)]
    order = ["control", "alpha", "beta", "gamma", "lambda"]
    lab = {"control": "Unstim.", **{k: v[0] for k, v in stims.items()}}
    cc = {"control": NEUTRAL, **{k: v[1] for k, v in stims.items()}}
    for i, st in enumerate(order):
        q = cpl[cpl.stim == st]
        jit = np.linspace(-0.22, 0.22, len(q))
        sig = (q.excess_q < 0.05).values
        ax.scatter(i + jit[~sig], q.excess.values[~sig], s=16, facecolor="none", edgecolor=cc[st], lw=0.8)
        ax.scatter(i + jit[sig], q.excess.values[sig], s=16, color=cc[st], edgecolor=SURF, lw=0.4)
    ax.set_xticks(range(5), [lab[s] for s in order])
    zero(ax, "y")
    grid(ax, "y")
    ax.set_ylabel("Excess correlation over matched genes")
    ax.scatter([], [], s=16, color=INK2, label="q < 0.05")
    ax.scatter([], [], s=16, facecolor="none", edgecolor=INK2, label="n.s.")
    ax.legend(loc="upper right")
    ax.set_title("…and couples to IFN programs under α / β only", loc="left")
    letter(ax, "b", dx=-0.2)

    # 3c: type I vs gamma, between-sample means
    ax = fig.add_subplot(gs[1, 0:2])
    f = pd.read_csv(O / "ifn_type_split/donors/fits.csv")
    models = {"htec_ciliated": "HTEC ciliated", "Enterocyte_progenitor_precursor_crypt_OLFM4lo": "Crypt OLFM4lo",
              "Enterocyte_progenitor_precursor_crypt_OLFM4lo@cd_inflamed": "  · inflamed",
              "Enterocyte_stem_OLFM4lo_LGR5hi": "Stem LGR5hi",
              "Enterocyte_top_villus@cd_inflamed": "Top-villus · inflamed", "fibroblast_pooled": "Fibroblasts"}
    y = np.arange(len(models))
    for term, lab_, dy, c_ in [("IFN_I_biased", "Type I-biased set", -0.15, CAT[0]),
                               ("IFN_G_biased", "γ-biased set", 0.15, CAT[1])]:
        q = f[f.term == term].set_index("model").reindex(list(models))
        ax.errorbar(q.coef, y + dy, xerr=1.96 * q.se, fmt="o", ms=4, color=c_, mec=SURF, mew=0.5,
                    elinewidth=0.8, label=lab_)
    ax.set_yticks(y, list(models.values()))
    ax.invert_yaxis()
    zero(ax)
    grid(ax)
    ax.set_xlabel("log2 ANTXR2 per SD of set score\n(95% CI, donor-clustered)")
    ax.set_title("Means follow type I, not γ", loc="left")
    legend_below(ax, dy=-0.32, ncol=1)
    letter(ax, "c", dx=-0.75)

    # 3d: Shaw 2017, 10 species
    ax = fig.add_subplot(gs[1, 2:4])
    sh = pd.read_csv(O / "ifn_type_split/shaw2017/antxr_query.csv")
    sh = sh[sh["query"] == "ANTXR2"].copy()
    nice = {"homo_sapiens": "Human", "rattus_norvegicus": "Rat", "bos_taurus": "Cow", "ovis_aries": "Sheep",
            "sus_scrofa": "Pig", "equus_caballus": "Horse", "canis_familiaris": "Dog",
            "pteropus_vampyrus": "Flying fox", "myotis_lucifugus": "Little brown bat", "gallus_gallus": "Chicken"}
    sh["name"] = sh.species.map(nice)
    sh = sh.sort_values("log2fc")
    y = np.arange(len(sh))
    hum = (sh.species == "homo_sapiens").values
    colors = np.where(hum, CAT[1], NEUTRAL)
    ax.hlines(y, 0, sh.log2fc, color=colors, lw=1.4)
    sig = (sh.fdr < 0.05).values
    ax.scatter(sh.log2fc[sig], y[sig], s=26, color=colors[sig], edgecolor=SURF, lw=0.5, zorder=3)
    ax.scatter(sh.log2fc[~sig], y[~sig], s=26, facecolor=SURF, edgecolor=colors[~sig], lw=0.9, zorder=3)
    ax.set_yticks(y, sh.name)
    zero(ax)
    grid(ax)
    ax.set_xlabel("ANTXR2 log2FC, type I IFN\n(filled: FDR < 0.05)")
    ax.set_title("Shaw 2017: human only", loc="left")
    ax.text(0.97, 0.42, "present + expressed\nin all 10 species;\nmammal identity 86–91%;\ndN/dS ≈ 0.2",
            transform=ax.transAxes, ha="right", va="center", fontsize=5.6, color=INK2)
    letter(ax, "d", dx=-0.62)

    # 3e: Hagai paired bulk, 4 species
    ax = fig.add_subplot(gs[1, 4:6])
    b = pd.read_csv(O / "ifn_type_split/fibroblast/bulk_ifnb.csv")
    sp = ["human", "rhesus", "rat", "mouse"]
    for i, s_ in enumerate(sp):
        q = b[b.species == s_]
        c_ = CAT[1] if s_ == "human" else NEUTRAL
        ax.scatter(np.full(len(q), i) + np.linspace(-0.15, 0.15, len(q)), q.antxr2_lfc, s=18, color=c_,
                   edgecolor=SURF, lw=0.5, zorder=3)
        ax.hlines(q.antxr2_lfc.median(), i - 0.28, i + 0.28, color=INK, lw=1)
    ax.set_xticks(range(4), [s.capitalize() for s in sp])
    zero(ax, "y")
    grid(ax, "y")
    ax.set_ylabel("ANTXR2 log2FC, IFN-β 8 h\n(paired, per donor)")
    ax.set_title("Hagai 2018: human 6 / 6 up", loc="left")
    letter(ax, "e", dx=-0.42)

    # 3f: scale: ANTXR2 vs canonical ISGs in the same experiments
    ax = fig.add_subplot(gs[2, 0:3])
    mx = json.load(open(O / "ifn_type_split/shaw2017/q_MX1.json"))["orthoClusters"][0]["speciesToGenes"]
    mx1 = mx["homo_sapiens"][0][3]
    an = sh.loc[sh.species == "homo_sapiens", "log2fc"].iloc[0]
    hb = b[b.species == "human"]
    vals = [("Shaw · ANTXR2", an, CAT[1]), ("Shaw · MX1", mx1, NEUTRAL),
            ("Hagai · ANTXR2", hb.antxr2_lfc.median(), CAT[1]), ("Hagai · ISG score\n(MX1, IFIT1, RSAD2)",
                                                                  hb.isg_lfc.median(), NEUTRAL)]
    y = np.arange(len(vals))[::-1]
    for yy, (lab_, v, c_) in zip(y, vals):
        ax.barh(yy, v, color=c_, height=0.55)
        ax.text(v + 0.15, yy, f"{v:+.2f}  ({2 ** v:,.0f}×)" if v > 3 else f"{v:+.2f}  ({2 ** v:.1f}×)",
                va="center", fontsize=6.3, color=INK)
    ax.set_yticks(y, [v[0] for v in vals])
    ax.set_xlim(0, 12)
    grid(ax)
    ax.set_xlabel("log2 fold change, human fibroblasts, type I IFN")
    ax.set_title("“Weak” is a measurement", loc="left")
    letter(ax, "f", dx=-0.3)

    # 3g: promoter schematic + ChIP-Atlas occupancy
    ax = fig.add_subplot(gs[2, 3:6])
    ax.set_xlim(-1050, 560)
    ax.set_ylim(-0.4, 3.1)
    for yy, lab_ in [(2.2, "Human ANTXR2"), (1.2, "Mouse Antxr2")]:
        ax.hlines(yy, -1000, 500, color=INK2, lw=1.5)
        ax.text(-1000, yy + 0.28, lab_, fontsize=6.3, va="bottom")
        ax.annotate("", xy=(120, yy + 0.25), xytext=(0, yy), arrowprops=dict(arrowstyle="-|>", lw=0.8,
                    color=INK, connectionstyle="angle,angleA=90,angleB=0,rad=0"))
    ax.add_patch(Rectangle((72 - 15, 2.05), 30, 0.3, color=CAT[1]))
    ax.text(130, 2.45, "ISRE +72  AGTTTCTTTTTCT", ha="left", va="bottom", fontsize=5.5)
    ax.add_patch(Rectangle((-179 - 10, 2.08), 20, 0.24, fc=SURF, ec=CAT[1], lw=0.9))
    ax.text(-179, 2.45, "core −179", ha="center", fontsize=5.5, color=INK2)
    ax.text(-250, 0.92, "no ISRE match in −1 kb…+0.5 kb", fontsize=5.8, color=INK2, ha="center", va="top")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-1000, -500, 0, 500], ["−1 kb", "−500", "TSS", "+500"])
    ax.text(0.0, 0.1, "ChIP-Atlas STAT1 (percentile among genes): ANTXR2 58th · MX1, IFIT1 99th\n"
            "→ GAPDH-level occupancy, as expected for a ~2× responder",
            transform=ax.transAxes, fontsize=5.7, color=INK2, va="bottom")
    ax.set_title("Human-only ISRE; occupancy fits a weak ISG", loc="left")
    letter(ax, "g", dx=-0.08)

    fig.suptitle("Figure 3 — ANTXR2 is a weak, human- and type I-specific interferon-stimulated gene", x=0.02,
                 ha="left", fontsize=8.5, fontweight="bold", y=0.995)
    place_letters(fig)
    fig.savefig(OUT / "fig3_draft.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3):
        f()
        print(f"-> {OUT / (f.__name__ + '_draft.png')}")
