"""ECM-clearance-pathway heatmap: ANTXR1, ANTXR2, and the clearance machinery
(MRC2/CTSB/CTSK/MMP14/TIMP2/LAMP1) across a curated set of cell types, annotated
with GAPO/HFS disease-tissue involvement.

Two panels, stacked vertically: whole-body atlas (curated cell types) and skin
fibroblast atlas (all 13 subtypes, healthy/nonlesional). Independent log color
scales per panel (see cell_type_curation.SCALE_CAVEAT_NOTE). Row order within
each panel is hierarchical-clustering leaf order (average linkage, correlation
distance on the log-transformed gene profile) -- not category or a fixed
reading order; no dendrogram is drawn, only the resulting order is used. This
figure shows raw expression + independently-sourced disease annotations only --
see cell_type_curation.SCORE_CAVEAT_NOTE: no composite score is computed here,
ever.

Usage: python plot_ecm_clearance_heatmap.py
"""
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

from cell_type_curation import (
    COVERAGE_GAP_NOTE, NORMALIZATION_NOTE, SCALE_CAVEAT_NOTE, SCORE_CAVEAT_NOTE,
    SKIN_DISEASE_STATUS, SKIN_STATUS_COLORS, SKIN_SUBTYPES,
    WHOLE_BODY_CATEGORIES, WHOLE_BODY_CATEGORY_COLORS,
)
from config import FIGURES_DIR, GENE_PANEL, OUTPUT_DIR
from gene_panel_query import aggregate_skin, aggregate_whole_body, scan_row_groups_for_genes

WHOLE_BODY_PATH = os.path.join(OUTPUT_DIR, "combined_celltype_means.parquet")
SKIN_PATH = os.path.join(OUTPUT_DIR, "skin_fibroblast_celltype_means.parquet")

COVERAGE_GAP_KEYWORDS = re.compile(
    r"chondro|osteo|thyroid|adrenal|dental|tooth|synov|gonad", re.IGNORECASE
)


SCAN_CACHE_DIR = "/tmp/ecm_heatmap_scan_cache"


def build_whole_body_matrix():
    curated_cell_types = list(WHOLE_BODY_CATEGORIES.keys())
    df = scan_row_groups_for_genes(
        WHOLE_BODY_PATH, GENE_PANEL,
        cache_path=os.path.join(SCAN_CACHE_DIR, "whole_body.parquet"),
    )
    agg = aggregate_whole_body(df, cell_types=curated_cell_types)

    found = set(agg["cell_type"].unique())
    missing = set(curated_cell_types) - found
    if missing:
        print(f"  WARNING: curated whole-body cell types not found in scan: {sorted(missing)}")

    matrix = agg.pivot(index="cell_type", columns="gene", values="mean_expression")
    matrix = matrix.reindex(columns=GENE_PANEL)
    n_meta = agg.pivot(index="cell_type", columns="gene", values="n_donors")
    row_meta = pd.DataFrame({
        "n": n_meta.max(axis=1),
        "n_constant": (n_meta.nunique(axis=1) == 1),
        "category": matrix.index.map(WHOLE_BODY_CATEGORIES),
    }, index=matrix.index)
    return matrix, row_meta


def build_skin_matrix():
    df = scan_row_groups_for_genes(
        SKIN_PATH, GENE_PANEL,
        cache_path=os.path.join(SCAN_CACHE_DIR, "skin.parquet"),
    )
    agg = aggregate_skin(df, subtypes=SKIN_SUBTYPES)

    found = set(agg["fibroblast_subtype"].unique())
    missing = set(SKIN_SUBTYPES) - found
    if missing:
        print(f"  WARNING: curated skin subtypes not found in scan: {sorted(missing)}")

    matrix = agg.pivot(index="fibroblast_subtype", columns="gene", values="mean_expression")
    matrix = matrix.reindex(columns=GENE_PANEL)
    n_meta = agg.pivot(index="fibroblast_subtype", columns="gene", values="n_samples")
    row_meta = pd.DataFrame({
        "n": n_meta.max(axis=1),
        "n_constant": (n_meta.nunique(axis=1) == 1),
        "category": matrix.index.map(SKIN_DISEASE_STATUS),
    }, index=matrix.index)
    return matrix, row_meta


def log_transform(matrix):
    nonzero = matrix.values[matrix.values > 0]
    floor = nonzero.min() / 2 if nonzero.size else 1e-9
    log_matrix = np.log10(matrix.astype(float) + floor)
    return log_matrix, floor


def cluster_row_order(log_matrix):
    """Average-linkage, correlation-distance leaf order over rows (clusters by
    the *shape* of each row's gene profile, not overall magnitude -- a
    deliberate choice, not the only reasonable one; euclidean distance would
    instead group primarily by magnitude). No dendrogram is drawn; only the
    resulting leaf order is used to reindex the matrix.
    """
    from scipy.cluster.hierarchy import dendrogram

    if len(log_matrix) < 3:
        return list(log_matrix.index)
    dist = pdist(log_matrix.values, metric="correlation")
    if np.isnan(dist).any():
        # correlation distance is undefined for a constant row (zero variance
        # across genes) -- fall back to euclidean for robustness rather than
        # crashing on a curated row that happens to be flat.
        dist = pdist(log_matrix.values, metric="euclidean")
    Z = linkage(dist, method="average")
    order = dendrogram(Z, no_plot=True)["leaves"]
    return [log_matrix.index[i] for i in order]


def normalize_per_gene(log_matrix, method="minmax"):
    """Rescale each gene (column) independently, computed on the log-transformed
    values. Without this, a uniformly-high-expression gene like LAMP1 (a
    lysosomal housekeeping gene, broadly expressed everywhere) anchors the
    shared color scale and drowns out its own -- and every other gene's --
    relative variation across cell types; this is about comparing each gene's
    own pattern on equal footing, not about absolute expression level, which
    is why it's applied per gene (column here; conventionally "row" in the
    gene-as-row heatmap layout this project doesn't use) rather than per cell
    type. Only the returned matrix is used for cell *color* -- text
    annotations and the CSV outputs still carry the real log/raw values, never
    the normalized ones, so nothing about the underlying numbers is hidden.

    - "minmax": rescale each gene to [0, 1] using its own min/max across the
      shown cell types. Keeps the existing sequential single-hue colormap
      valid (0=that gene's lowest shown cell type, 1=its highest) -- default,
      since it doesn't require redesigning the color scheme.
    - "zscore": rescale each gene to mean 0, unit variance. This produces a
      genuinely diverging quantity (above/below that gene's own average) and
      would need a diverging colormap (two hues + neutral midpoint) to display
      correctly per the dataviz skill -- not implemented here to keep this a
      single, non-disruptive change; ask if you want this instead.
    """
    if method == "minmax":
        col_min = log_matrix.min(axis=0)
        col_max = log_matrix.max(axis=0)
        span = (col_max - col_min).replace(0, 1)  # constant column -> avoid /0, reads as 0
        return (log_matrix - col_min) / span
    elif method == "zscore":
        raise NotImplementedError(
            "zscore normalization needs a diverging colormap redesign -- see docstring"
        )
    raise ValueError(f"unknown normalize method {method!r}")


def draw_heatmap_panel(fig, gs_slice, log_matrix, display_matrix, row_meta, color_map, title, annotate_cells):
    inner = gs_slice.subgridspec(1, 3, width_ratios=[0.35, 10, 0.4], wspace=0.05)
    ax_strip = fig.add_subplot(inner[0, 0])
    ax_main = fig.add_subplot(inner[0, 1])
    ax_cbar = fig.add_subplot(inner[0, 2])

    n_rows, n_cols = log_matrix.shape
    strip_rgb = np.array([
        [int(color_map.get(row_meta.loc[r, "category"], ("#cccccc",))[0].lstrip("#")[i:i+2], 16) / 255
         for i in (0, 2, 4)]
        for r in log_matrix.index
    ]).reshape(n_rows, 1, 3)
    ax_strip.imshow(strip_rgb, aspect="auto")
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)

    # Color comes from display_matrix (per-gene-normalized, [0,1]); text
    # annotations and everything else read from log_matrix (the real values) --
    # normalization only ever changes what color a cell is painted, never what
    # number is reported for it.
    vmin, vmax = 0.0, 1.0
    im = ax_main.imshow(display_matrix.values, aspect="auto", cmap="Blues", vmin=vmin, vmax=vmax)
    ax_main.set_xticks(range(n_cols))
    labels = list(log_matrix.columns)
    ax_main.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    for tick, gene in zip(ax_main.get_xticklabels(), labels):
        if gene in ("ANTXR1", "ANTXR2"):
            tick.set_fontweight("bold")
    ax_main.set_yticks(range(n_rows))
    row_labels = [f"{r} (n={int(row_meta.loc[r, 'n'])})" for r in log_matrix.index]
    ax_main.set_yticklabels(row_labels, fontsize=8)
    ax_main.set_title(title, fontsize=11, loc="left", fontweight="600")

    if annotate_cells:
        for i in range(n_rows):
            for j in range(n_cols):
                real_val = log_matrix.values[i, j]
                rel = display_matrix.values[i, j]
                color = "white" if rel > 0.6 else "#333333"
                ax_main.text(j, i, f"{10**real_val:.1e}", ha="center", va="center",
                             fontsize=6, color=color)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("expression, scaled 0-1 per gene\n(this panel's cell types)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    for r in log_matrix.index:
        if not row_meta.loc[r, "n_constant"]:
            print(f"  NOTE: n_donors/n_samples not constant across genes for row {r!r}")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    out_dir = os.path.join(FIGURES_DIR, "ecm_clearance_heatmap")
    os.makedirs(out_dir, exist_ok=True)

    print("=== whole-body atlas ===")
    wb_matrix, wb_meta = build_whole_body_matrix()
    print("=== skin fibroblast atlas ===")
    skin_matrix, skin_meta = build_skin_matrix()

    print(f"\ngenes found: whole-body={sorted(wb_matrix.columns)}, skin={sorted(skin_matrix.columns)}")
    assert sorted(wb_matrix.columns.dropna()) == sorted(GENE_PANEL), "whole-body gene panel mismatch"
    assert sorted(skin_matrix.columns.dropna()) == sorted(GENE_PANEL), "skin gene panel mismatch"

    wb_log, wb_floor = log_transform(wb_matrix)
    skin_log, skin_floor = log_transform(skin_matrix)
    print(f"\nwhole-body: floor={wb_floor:.3e}, log range=({wb_log.values.min():.2f}, {wb_log.values.max():.2f})")
    print(f"skin: floor={skin_floor:.3e}, log range=({skin_log.values.min():.2f}, {skin_log.values.max():.2f})")

    mrc2_wb = wb_matrix["MRC2"]
    print(f"\nMRC2 whole-body curated range: {mrc2_wb.min():.3e} - {mrc2_wb.max():.3e}")
    mrc2_skin = skin_matrix["MRC2"]
    print(f"MRC2 skin range: {mrc2_skin.min():.3e} - {mrc2_skin.max():.3e}")

    ratio = np.log10(skin_matrix["ANTXR2"] + skin_floor) - np.log10(skin_matrix["ANTXR1"] + skin_floor)
    n_favor_2 = (ratio > 0).sum()
    print(f"\n[verification only, not written to file] skin ANTXR2>ANTXR1 in {n_favor_2}/{len(ratio)} subtypes")

    # cheap re-check against the curated list's own keys plus a keyword scan of
    # the categories used -- full 238-vocab keyword recheck done separately in
    # the exploration phase; here just confirm the coverage-gap note's keywords
    # don't appear in the curated list itself (sanity, not a full atlas
    # rescan, which would take ~1-2 min for no new information here).
    for ct in WHOLE_BODY_CATEGORIES:
        if COVERAGE_GAP_KEYWORDS.search(ct):
            print(f"  NOTE: curated cell type {ct!r} unexpectedly matches a coverage-gap keyword")

    wb_order = cluster_row_order(wb_log)
    skin_order = cluster_row_order(skin_log)
    wb_log = wb_log.loc[wb_order]
    wb_meta = wb_meta.loc[wb_order]
    skin_log = skin_log.loc[skin_order]
    skin_meta = skin_meta.loc[skin_order]

    # Per-gene (column) min-max normalization for cell *color* only -- see
    # normalize_per_gene's docstring. Requested because a uniformly-high gene
    # like LAMP1 was dominating the shared color scale and drowning out every
    # gene's own relative cross-cell-type pattern.
    wb_display = normalize_per_gene(wb_log)
    skin_display = normalize_per_gene(skin_log)
    print("\npre-normalization log10 range per gene, whole-body panel (shows the "
          "between-gene baseline-offset problem this normalization fixes):")
    print(pd.DataFrame({"min": wb_log.min(axis=0), "max": wb_log.max(axis=0)}).round(2).to_string())

    # Explicit spacer row (in the same height-ratio units as the row counts,
    # i.e. per-row-inch scaled) rather than a GridSpec `hspace` fraction --
    # hspace scales with each axes' own height, so a single fraction either
    # collided with the top panel's rotated x-tick labels (too small, given
    # the top panel is much taller than the bottom one) or left a huge gap
    # (too large). A fixed-size spacer avoids retuning this if row counts
    # change.
    ROW_IN_WB, ROW_IN_SKIN, SPACER_IN = 0.24, 0.30, 1.0
    fig_h = ROW_IN_WB * len(wb_log) + ROW_IN_SKIN * len(skin_log) + SPACER_IN + 3.2
    fig = plt.figure(figsize=(11, fig_h))
    fig.patch.set_facecolor("#fcfcfb")
    gs = fig.add_gridspec(
        3, 1,
        height_ratios=[ROW_IN_WB * len(wb_log), SPACER_IN, ROW_IN_SKIN * len(skin_log)],
        hspace=0, top=0.94, bottom=0.16, left=0.28, right=0.90,
    )

    draw_heatmap_panel(fig, gs[0], wb_log, wb_display, wb_meta, WHOLE_BODY_CATEGORY_COLORS,
                        "Whole-body atlas (curated cell types)", annotate_cells=False)
    draw_heatmap_panel(fig, gs[2], skin_log, skin_display, skin_meta, SKIN_STATUS_COLORS,
                        "Skin fibroblast atlas (healthy, nonlesional)", annotate_cells=True)

    fig.suptitle("ANTXR1 / ANTXR2 and the ECM-clearance pathway, by cell type",
                 fontsize=13, fontweight="700", x=0.28, ha="left")

    wb_handles = [Patch(facecolor=c[0], label=cat) for cat, c in WHOLE_BODY_CATEGORY_COLORS.items()]
    skin_handles = [Patch(facecolor=c[0], label=cat) for cat, c in SKIN_STATUS_COLORS.items()]
    fig.legend(handles=wb_handles, loc="upper left", bbox_to_anchor=(0.005, 0.98),
               fontsize=7, frameon=False, title="Whole-body category", title_fontsize=7.5)
    fig.legend(handles=skin_handles, loc="upper left", bbox_to_anchor=(0.005, 0.30),
               fontsize=7, frameon=False, title="Skin subtype status", title_fontsize=7.5)

    caption = "\n".join([SCORE_CAVEAT_NOTE, SCALE_CAVEAT_NOTE, NORMALIZATION_NOTE, COVERAGE_GAP_NOTE])
    fig.text(0.28, 0.01, caption, fontsize=6.5, color="#52514e", wrap=True,
              va="bottom", ha="left")

    png_path = os.path.join(out_dir, "heatmap.png")
    svg_path = os.path.join(out_dir, "heatmap.svg")
    fig.savefig(png_path, dpi=200, facecolor=fig.get_facecolor())
    fig.savefig(svg_path, facecolor=fig.get_facecolor())
    print(f"\nwrote {png_path}")
    print(f"wrote {svg_path}")

    def to_long(matrix, display, meta, id_col):
        m = matrix.copy()
        m.index.name = id_col
        long = m.reset_index().melt(id_vars=id_col, var_name="gene", value_name="mean_expression")
        long["log10_mean_expression"] = np.log10(long["mean_expression"] + (
            wb_floor if id_col == "cell_type" else skin_floor))
        d = display.copy()
        d.index.name = id_col
        norm_long = d.reset_index().melt(id_vars=id_col, var_name="gene", value_name="normalized_0_1")
        long = long.merge(norm_long, on=[id_col, "gene"])
        long = long.merge(meta.reset_index().rename(columns={"index": id_col}), on=id_col)
        return long

    wb_csv = to_long(wb_matrix, wb_display, wb_meta, "cell_type")
    skin_csv = to_long(skin_matrix, skin_display, skin_meta, "fibroblast_subtype")

    for path, df_, note in [
        (os.path.join(out_dir, "whole_body_curated.csv"), wb_csv, COVERAGE_GAP_NOTE),
        (os.path.join(out_dir, "skin_curated.csv"), skin_csv, COVERAGE_GAP_NOTE),
    ]:
        with open(path, "w") as f:
            f.write(f"# {SCORE_CAVEAT_NOTE}\n# {SCALE_CAVEAT_NOTE}\n# {NORMALIZATION_NOTE}\n# {note}\n")
            df_.to_csv(f, index=False)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
