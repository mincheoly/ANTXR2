"""UMAPs: gut atlas -> gut endothelium -> gut lymphatic EC, plus the Tabula
Sapiens endothelium view that explains the annotation-granularity effect.

Figures (PNG + SVG) in <out>:
  A_gut_atlas.png          atlas UMAP (shipped X_umap): endothelium / LEC
                           highlighted on grey; ANTXR2 and PROX1
  B_ts_endothelium.png     Tabula Sapiens endothelium (shipped scVI compartment
                           UMAP): LEC-labelled vs generic "endothelial cell" vs
                           other EC labels; PROX1, CCL21, ANTXR2
  C_gut_endothelium.png    gut endothelium, UMAP recomputed on the atlas's
                           integrated X_pca restricted to endothelium
  D_gut_lec_subtypes.png   gut LEC, UMAP recomputed on X_pca restricted to LEC;
                           one panel per author subtype (small multiples) plus
                           stage / region / chemistry
  E_gut_lec_features.png   LEC feature plots
  F_lec_qc.png             per donor x subtype QC: genes, UMIs, % mito, doublet score

Embeddings use each atlas's own integrated space (gut: X_pca; TS: X_scvi /
shipped UMAP). A fresh PCA on the subset would mostly separate 3'/5'
chemistry and fetal/adult samples.

Expression = log1p(counts / library x 1e4). Zero cells drawn first in grey,
expressing cells on top in a single-hue ramp. Categorical colour is capped at
three hues per panel (scatter: every pair co-visible); more categories are
faceted.

Usage: python lymphatic_umaps.py <workspace dir> <out dir>
"""
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fib_cell_extract as F                       # noqa: E402
import stream_panel_extract as S                   # noqa: E402

PANEL = ["ANTXR2", "ANTXR1", "PROX1", "LYVE1", "CCL21", "TFF3", "MMRN1", "PDPN", "FLT4",
         "PECAM1", "CDH5", "PLVAP", "ADAMTS4", "FN1", "MFAP4", "PIEZO2", "CLDN11", "FOXC2",
         "GJA4", "ITGA9", "ACKR4", "STAB2", "MADCAM1", "ADGRG3", "MRC2", "CSF3", "CXCL1",
         "NFKBIZ", "TNFAIP3", "IGF1", "COL6A3", "BMP2", "MKI67", "TOP2A", "IL33", "ADM",
         "HEXB", "SNX2", "RELN", "PTPRC", "EPCAM", "DCN", "COL1A2"]
GUT_OBS = ["author_cell_type", "cell_type", "category", "Age_group", "tissue", "assay",
           "donor_id", "doublet_scores", "pct_counts_mt", "n_genes", "total_counts"]
TS_OBS = ["cell_type", "tissue_in_publication", "donor_id", "assay"]
LEC_LABEL = "endothelial cell of lymphatic vessel"

# reference palette (dataviz skill), light mode
SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GREY = "#d6d5d1"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
RAMP = LinearSegmentedColormap.from_list(
    "blue", ["#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"])
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK2, "font.size": 9,
                     "axes.titlesize": 10, "axes.titleweight": "normal"})


def load_panel(h5ad, cache, obs_cols):
    if not os.path.exists(cache):
        S.extract(h5ad, cache, panel=PANEL, obs_cols=obs_cols)
    z = np.load(cache, allow_pickle=False)
    expr = pd.DataFrame(np.log1p(z["counts"] / z["total"][:, None] * 1e4), columns=z["genes"].astype(str))
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    obs["lib"] = z["total"]
    return expr, obs


def obsm(h5ad, key):
    with h5py.File(h5ad, "r") as f:
        return f["obsm"][key][:]


def sub_umap(pcs, seed=0, n_neighbors=15):
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(np.zeros((pcs.shape[0], 1), np.float32))
    a.obsm["X_int"] = pcs.astype(np.float32)
    sc.pp.neighbors(a, use_rep="X_int", n_neighbors=n_neighbors, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    return a.obsm["X_umap"]


def clean(ax, title):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, color=INK, loc="left")


def pt(n):
    return float(np.clip(30000 / max(n, 1), 0.3, 12))


def categorical(ax, xy, masks, title):
    """masks: list of (label, bool mask, colour); unmasked cells are grey context."""
    any_m = np.zeros(len(xy), bool)
    for _, m, _ in masks:
        any_m |= m
    s = pt(len(xy))
    ax.scatter(*xy[~any_m].T, s=s, c=GREY, lw=0, rasterized=True)
    for lab, m, col in masks:
        ax.scatter(*xy[m].T, s=s * 1.6, c=col, lw=0, rasterized=True, label=f"{lab} ({m.sum():,})")
    # below the panel, never over the points
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 0), markerscale=4 / np.sqrt(max(s, 0.5)),
              fontsize=8, labelcolor=INK2, handletextpad=0.2, borderaxespad=0)
    clean(ax, title)


def feature(ax, xy, v, title):
    s = pt(len(xy))
    z = v <= 0
    ax.scatter(*xy[z].T, s=s, c=GREY, lw=0, rasterized=True)
    o = np.argsort(v[~z])
    sc_ = ax.scatter(*xy[~z][o].T, s=s * 1.6, c=v[~z][o], cmap=RAMP, lw=0, rasterized=True,
                     vmin=0, vmax=max(np.percentile(v[~z], 99), 1e-6) if (~z).any() else 1)
    cb = plt.colorbar(sc_, ax=ax, fraction=0.04, pad=0.01)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=INK2)
    clean(ax, f"{title}  ({(~z).mean() * 100:.0f}% >0)")


def label_clusters(ax, xy, labels, min_cells=50, size=8):
    """Direct text labels at each label's median position (identity by text,
    not colour); white halo keeps them legible over points."""
    from matplotlib import patheffects as pe
    span = np.ptp(xy[:, 1])
    placed = []
    for lab in sorted(pd.unique(labels), key=lambda l: -np.sum(labels == l)):
        m = labels == lab
        if m.sum() < min_cells:
            continue
        x, y = np.median(xy[m], axis=0)
        # nudge down past any label already placed too close (no overlapping text)
        while any(abs(x - px) < 0.22 * np.ptp(xy[:, 0]) and abs(y - py) < 0.045 * span for px, py in placed):
            y -= 0.05 * span
        placed.append((x, y))
        ax.text(x, y, f"{lab}", fontsize=size, color=INK, ha="center", va="center",
                path_effects=[pe.withStroke(linewidth=2.5, foreground=SURFACE)])


def save(fig, out, name):
    fig.savefig(os.path.join(out, name + ".png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(out, name + ".svg"), bbox_inches="tight")
    plt.close(fig)


def main(ws, out):
    os.makedirs(out, exist_ok=True)
    cache = os.path.join(out, "_cache")
    os.makedirs(cache, exist_ok=True)
    gut = os.path.join(ws, "data", "gut.h5ad")
    ts = os.path.join(ws, "data", "ts_endothelium.h5ad")
    ge, go = load_panel(gut, os.path.join(cache, "gut_panel.npz"), GUT_OBS)
    act = go.author_cell_type.astype(str)
    is_endo = (go.category == "Endothelial").values
    is_lec = act.str.startswith("LEC").values

    # A: whole gut atlas
    xy = obsm(gut, "X_umap")
    cat = go.category.astype(str).values
    fig, axs = plt.subplots(1, 4, figsize=(21, 5.6), layout="constrained")
    axs[0].scatter(*xy.T, s=pt(len(xy)), c=GREY, lw=0, rasterized=True)
    label_clusters(axs[0], xy, cat, size=10)
    clean(axs[0], f"Gut Cell Atlas (Elmentaite 2021), {len(xy):,} cells: atlas categories")
    categorical(axs[1], xy, [("blood-vessel endothelium", is_endo & ~is_lec, C1),
                             ("lymphatic endothelium", is_lec, C2)], "endothelium")
    feature(axs[2], xy, ge.ANTXR2.values, "ANTXR2")
    label_clusters(axs[2], xy, cat, size=8)
    feature(axs[3], xy, ge.PROX1.values, "PROX1")
    save(fig, out, "A_gut_atlas")

    # B: Tabula Sapiens endothelium, annotation granularity
    te, to = load_panel(ts, os.path.join(cache, "ts_endo_panel.npz"), TS_OBS)
    txy = obsm(ts, "X_umap_compartment_scvi_donorassay")
    ct = to.cell_type.astype(str).values
    fig, axs = plt.subplots(1, 4, figsize=(19, 5.4), layout="constrained")
    categorical(axs[0], txy, [("labelled 'endothelial cell' (generic)", ct == "endothelial cell", C1),
                              ("labelled lymphatic", ct == LEC_LABEL, C2)],
                f"Tabula Sapiens endothelium, {len(txy):,} cells")
    for ax, g in zip(axs[1:], ["PROX1", "CCL21", "ANTXR2"]):
        feature(ax, txy, te[g].values, g)
    save(fig, out, "B_ts_endothelium")
    gen = ct == "endothelial cell"
    lymph_like = (te.PROX1.values > 0) | (te.CCL21.values > 0)
    tab = (pd.DataFrame({"tissue": to.tissue_in_publication.astype(str).values, "generic": gen,
                         "labelled_lec": ct == LEC_LABEL, "prox1_or_ccl21": lymph_like})
           .groupby("tissue").apply(lambda d: pd.Series({
               "n_cells": len(d), "frac_generic_label": d.generic.mean(),
               "n_labelled_lec": int(d.labelled_lec.sum()),
               "n_generic_prox1_or_ccl21": int((d.generic & d.prox1_or_ccl21).sum())}))
           .sort_values("n_generic_prox1_or_ccl21", ascending=False))
    tab.to_csv(os.path.join(out, "B_ts_granularity_by_tissue.csv"))

    # C: gut endothelium, recomputed UMAP on integrated PCA
    pca = obsm(gut, "X_pca")
    e_idx = np.flatnonzero(is_endo)
    exy = sub_umap(pca[e_idx])
    esub = act.values[e_idx]
    # LEC subtypes are resolved in figure D; label their island once here
    elab = np.where(pd.Series(esub).astype(str).str.startswith("LEC").values, "lymphatic (LEC1-6)", esub)
    fig, axs = plt.subplots(1, 4, figsize=(21, 5.6), layout="constrained")
    axs[0].scatter(*exy.T, s=pt(len(exy)), c=GREY, lw=0, rasterized=True)
    label_clusters(axs[0], exy, elab, min_cells=40, size=8)
    clean(axs[0], f"Gut endothelium, {len(e_idx):,} cells: author subtypes")
    feature(axs[1], exy, ge.ANTXR2.values[e_idx], "ANTXR2")
    label_clusters(axs[1], exy, elab, min_cells=40, size=7)
    for ax, g in zip(axs[2:], ["PROX1", "PLVAP"]):
        feature(ax, exy, ge[g].values[e_idx], g)
    save(fig, out, "C_gut_endothelium")

    # G: ANTXR2 per endothelial subtype, one point per donor (pseudobulk)
    zc = np.load(os.path.join(cache, "gut_panel.npz"), allow_pickle=False)
    a2 = zc["counts"][:, list(zc["genes"].astype(str)).index("ANTXR2")][e_idx]
    lib = zc["total"][e_idx]
    stage_e = np.where(pd.Series(go.Age_group.values[e_idx]).str.contains("trim"), "fetal", "postnatal")
    d = pd.DataFrame({"subtype": esub, "stage": stage_e, "donor": go.donor_id.values[e_idx], "a2": a2,
                      "lib": lib, "pos": a2 > 0})
    g = d.groupby(["subtype", "stage", "donor"])
    pbk = pd.DataFrame({"n_cells": g.size(), "antxr2_cp10k": g.a2.sum() / g.lib.sum() * 1e4,
                        "antxr2_det": g.pos.mean()}).reset_index()
    pbk = pbk[pbk.n_cells >= 20]
    summ = (pbk.groupby(["subtype", "stage"]).agg(n_donors=("donor", "size"),
                                                   median_cp10k=("antxr2_cp10k", "median"),
                                                   median_det=("antxr2_det", "median"))
            .reset_index().sort_values("median_cp10k", ascending=False))
    pbk.to_csv(os.path.join(out, "G_endothelial_subtype_antxr2_by_donor.csv"), index=False)
    summ.to_csv(os.path.join(out, "G_endothelial_subtype_antxr2_summary.csv"), index=False)
    order = (pbk.groupby("subtype").antxr2_cp10k.median().sort_values(ascending=False).index.tolist())
    fig, ax = plt.subplots(figsize=(13, 4.8), layout="constrained")
    for i, sname in enumerate(order):
        for st, col in (("postnatal", C2), ("fetal", C1)):
            dd = pbk[(pbk.subtype == sname) & (pbk.stage == st)]
            if dd.empty:
                continue
            jit = np.random.default_rng(i).uniform(-0.18, 0.18, len(dd))
            ax.scatter(i + jit, dd.antxr2_cp10k, s=22, c=col, lw=0, label=st)
        ax.hlines(pbk[pbk.subtype == sname].antxr2_cp10k.median(), i - 0.3, i + 0.3, color=INK, lw=1.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=35, ha="right", color=INK2)
    ax.set_ylabel("ANTXR2 CP10K (donor pseudobulk)")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]:
        ax.spines[sp].set_color(GREY)
    ax.tick_params(colors=INK2)
    h, l = ax.get_legend_handles_labels()
    keep = dict(zip(l, h))
    ax.legend(keep.values(), keep.keys(), frameon=False, loc="upper right", labelcolor=INK2)
    ax.set_title("Gut endothelium: ANTXR2 per subtype, one point per donor (>= 20 cells); line = median",
                 loc="left", color=INK)
    save(fig, out, "G_endothelial_subtype_antxr2")
    print(summ.round(3).to_string(index=False))

    # D: gut LEC subtypes
    l_idx = np.flatnonzero(is_lec)
    lxy = sub_umap(pca[l_idx], n_neighbors=15)
    pd.DataFrame(lxy, columns=["umap1", "umap2"]).assign(
        subtype=act.values[l_idx], donor=go.donor_id.values[l_idx]).to_csv(
        os.path.join(out, "D_gut_lec_umap_coords.csv"), index=False)
    lsub = act.values[l_idx]
    subs = sorted(set(lsub))
    stage = np.where(pd.Series(go.Age_group.values[l_idx]).str.contains("trim"), "fetal", "postnatal")
    region = np.where(go.tissue.values[l_idx] == "mesenteric lymph node", "mesenteric lymph node", "gut wall")
    assay = go.assay.values[l_idx]
    fig, axs = plt.subplots(3, 3, figsize=(13, 14), layout="constrained")
    for ax, s in zip(axs.flat, subs):
        categorical(ax, lxy, [(s, lsub == s, C2)], s)
    categorical(axs.flat[6], lxy, [("fetal", stage == "fetal", C1), ("postnatal", stage == "postnatal", C2)],
                "developmental stage")
    categorical(axs.flat[7], lxy, [("mesenteric lymph node", region == "mesenteric lymph node", C1),
                                   ("gut wall", region == "gut wall", C2)], "region")
    categorical(axs.flat[8], lxy, [("10x 3' v2", assay == "10x 3' v2", C1), ("10x 5' v2", assay == "10x 5' v2", C2)],
                "chemistry")
    fig.suptitle(f"Gut lymphatic endothelium, {len(l_idx):,} cells (UMAP on the atlas's integrated PCA)",
                 x=0.01, ha="left", color=INK)
    save(fig, out, "D_gut_lec_subtypes")

    # E: LEC feature plots
    feats = ["ANTXR2", "PROX1", "LYVE1", "CCL21", "ADAMTS4", "FN1", "MFAP4", "COL6A3", "PIEZO2",
             "CSF3", "NFKBIZ", "IGF1", "CLDN11", "FOXC2", "GJA4", "ACKR4", "STAB2", "MADCAM1",
             "ADGRG3", "MRC2", "MKI67", "IL33", "ADM", "HEXB"]
    fig, axs = plt.subplots(4, 6, figsize=(24, 15), layout="constrained")
    for ax, g in zip(axs.flat, feats):
        feature(ax, lxy, ge[g].values[l_idx], g)
    fig.suptitle("Gut lymphatic endothelium: ANTXR2, LEC identity, LEC6 program, valve, other subtypes, partners",
                 x=0.01, ha="left", color=INK)
    save(fig, out, "E_gut_lec_features")

    # F: QC per donor x subtype (>= 10 cells), one measure per panel
    num = lambda c: pd.to_numeric(go[c].values[l_idx])   # panel cache stores obs as text
    q = pd.DataFrame({"subtype": lsub, "stage": stage, "donor": go.donor_id.values[l_idx],
                      "n_genes": num("n_genes"), "UMIs": num("total_counts"),
                      "pct_mito": num("pct_counts_mt"), "doublet_score": num("doublet_scores")})
    g = q.groupby(["subtype", "stage", "donor"])
    dq = g[["n_genes", "UMIs", "pct_mito", "doublet_score"]].median().assign(n=g.size()).reset_index()
    dq = dq[dq.n >= 10]
    dq.to_csv(os.path.join(out, "F_lec_qc_by_donor_subtype.csv"), index=False)
    order = [s for s in subs if s in set(dq.subtype)]
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.6), layout="constrained")
    for ax, m in zip(axs, ["n_genes", "UMIs", "pct_mito", "doublet_score"]):
        for i, s in enumerate(order):
            d = dq[dq.subtype == s]
            col = C2 if s.startswith("LEC6") else C1
            ax.scatter(np.full(len(d), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(d)),
                       d[m], s=22, c=col, lw=0)
            ax.hlines(d[m].median(), i - 0.3, i + 0.3, color=INK, lw=1.5)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([s.split(" ")[0] for s in order], color=INK2)
        ax.set_title(f"median {m} per donor", loc="left", color=INK)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color(GREY)
        ax.spines["bottom"].set_color(GREY)
        ax.tick_params(colors=INK2)
    fig.suptitle("Gut LEC quality metrics, one point per donor x subtype (>= 10 cells); LEC6 in orange",
                 x=0.01, ha="left", color=INK)
    save(fig, out, "F_lec_qc")
    print("figures written to", out)
    print(tab.head(12).to_string())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
