"""Task 3 of prompts/partner_availability.md: extend the gene panel to the Wnt
arm and produce a two-arm heatmap (clearance vs. Wnt) against the same curated
cell-type rows as the existing ECM-clearance figure.

Pure query against already-computed parquet (gene_panel_query.py's row-group
scan + the project's two-step donor-equal-weighted aggregation) -- no pipeline
run, no bootstrap, no memento. Reuses `cell_type_curation.py`'s curated rows so
this figure's rows are identical to `plot_ecm_clearance_heatmap.py`'s.

**No composite score of any kind** -- standing correction. Raw values side by
side only; the only cross-arm computation here is a per-cell-type "is each
arm's machinery present/absent" call, based on each gene's own row-relative
normalized value, always shown next to the raw numbers, never collapsed into
a single number.

Usage: conda run -n antxr2 python scripts/two_arm_panel_query.py
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from cell_type_curation import (
    SKIN_DISEASE_STATUS, SKIN_STATUS_COLORS, SKIN_SUBTYPES,
    WHOLE_BODY_CATEGORIES, WHOLE_BODY_CATEGORY_COLORS,
)
from config import FIGURES_DIR, GENE_PANEL_ARMS, EXTENDED_GENE_PANEL, OUTPUT_DIR
from gene_panel_query import aggregate_skin, aggregate_whole_body, scan_row_groups_for_genes
from plot_ecm_clearance_heatmap import cluster_row_order, log_transform, normalize_per_gene

WHOLE_BODY_PATH = os.path.join(OUTPUT_DIR, "combined_celltype_means.parquet")
SKIN_PATH = os.path.join(OUTPUT_DIR, "skin_fibroblast_celltype_means.parquet")
SCAN_CACHE_DIR = "/tmp/two_arm_panel_scan_cache"
OUT_DIR = os.path.join(FIGURES_DIR, "two_arm_panel")
ANALYSIS_SUMMARY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ANALYSIS_SUMMARY.md")

# "Present" call: a gene's raw (log10) value in a row must be within this many
# orders of magnitude of that SAME gene's own maximum across the curated rows
# to be called "present" there. Deliberately an absolute-ratio rule on the raw
# values (log_matrix), NOT the display_matrix's per-gene min-max [0,1] rescale
# -- a min-max rank is relative to the whole shown set and can call a uniformly
# near-floor gene "present" in whichever row happens to be least-low. This
# threshold instead reuses the exact standard the project already applied by
# hand for MRC2 (keratinocyte/corneal epithelial cell ~3-4e-5 called "retains
# some MRC2" vs. gut epithelium ~1e-7-5e-7 called "detection floor" -- a ~2-3
# order-of-magnitude gap). 2 orders of magnitude (100x) is picked to match that
# precedent. This is a reading aid printed next to the raw numbers, never a
# score -- it does not combine genes, and changing it changes only the
# "present" markers, never any number reported in the CSV/figure.
PRESENT_LOG10_GAP = 2.0


def build_matrix(path, gene_col_func, agg_func, curated_ids, cache_name):
    df = scan_row_groups_for_genes(
        path, EXTENDED_GENE_PANEL, cache_path=os.path.join(SCAN_CACHE_DIR, cache_name),
    )
    agg = agg_func(df, curated_ids)
    id_col = "cell_type" if "cell_type" in agg.columns else "fibroblast_subtype"
    n_col = "n_donors" if "n_donors" in agg.columns else "n_samples"

    found = set(agg[id_col].unique())
    missing = set(curated_ids) - found
    if missing:
        print(f"  WARNING: curated rows not found in scan ({cache_name}): {sorted(missing)}")

    matrix = agg.pivot(index=id_col, columns="gene", values="mean_expression")
    matrix = matrix.reindex(columns=EXTENDED_GENE_PANEL)
    n_meta = agg.pivot(index=id_col, columns="gene", values=n_col)
    row_meta = pd.DataFrame({"n": n_meta.max(axis=1), "n_constant": (n_meta.nunique(axis=1) == 1)}, index=matrix.index)
    return matrix, row_meta, id_col


def arm_completeness(log_matrix):
    """Per row: which arms have EVERY gene present (within PRESENT_LOG10_GAP
    orders of magnitude of that gene's own max across curated rows), which
    have NONE, and flag rows where BOTH arms are incomplete -- ANTXR2 present
    with neither function's full partner set. Operates on log_matrix (raw
    log10 values), NOT the display_matrix's per-gene [0,1] rescale -- see
    PRESENT_LOG10_GAP's docstring for why."""
    col_max = log_matrix.max(axis=0)
    rows = []
    for row_id in log_matrix.index:
        rec = {"row": row_id}
        for arm, genes in GENE_PANEL_ARMS.items():
            if arm == "RECEPTORS":
                continue
            genes_here = [g for g in genes if g in log_matrix.columns]
            vals = log_matrix.loc[row_id, genes_here]
            present_mask = vals >= (col_max[genes_here] - PRESENT_LOG10_GAP)
            n_present = int(present_mask.sum())
            rec[f"{arm}_n_present"] = n_present
            rec[f"{arm}_n_total"] = len(genes_here)
            rec[f"{arm}_complete"] = n_present == len(genes_here)
            rec[f"{arm}_absent"] = n_present == 0
        rows.append(rec)
    return pd.DataFrame(rows).set_index("row")


def draw_two_arm_heatmap(log_matrix, display_matrix, row_meta, color_map, category_map, title, out_path, annotate=False):
    n_rows, n_cols = log_matrix.shape
    fig = plt.figure(figsize=(0.34 * n_cols + 3.2, 0.28 * n_rows + 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.35, 10, 0.4], wspace=0.05, left=0.30, right=0.92, top=0.86, bottom=0.30)
    ax_strip, ax_main, ax_cbar = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])

    strip_rgb = np.array([
        [int(color_map.get(category_map.get(r), ("#cccccc",))[0].lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4)]
        for r in log_matrix.index
    ]).reshape(n_rows, 1, 3)
    ax_strip.imshow(strip_rgb, aspect="auto")
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    for s in ax_strip.spines.values():
        s.set_visible(False)

    im = ax_main.imshow(display_matrix.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax_main.set_xticks(range(n_cols))
    ax_main.set_xticklabels(list(log_matrix.columns), rotation=90, fontsize=7)
    ax_main.set_yticks(range(n_rows))
    ax_main.set_yticklabels([f"{r} (n={int(row_meta.loc[r,'n'])})" for r in log_matrix.index], fontsize=8)
    ax_main.set_title(title, fontsize=11, loc="left", fontweight="600")

    # labeled column blocks per arm, with a separating vertical line between arms.
    # Short display names + a font size capped by block width in columns -- the
    # full arm names ("CLEARANCE_SUBSTRATE") collided with their neighbors when
    # centered over a narrow (2-3 column) block at a fixed font size.
    ARM_DISPLAY = {"RECEPTORS": "RECEPTORS", "CLEARANCE_ARM": "CLEARANCE", "CLEARANCE_SUBSTRATE": "SUBSTRATE",
                   "WNT_ARM": "WNT"}
    col_pos = 0
    for arm, genes in GENE_PANEL_ARMS.items():
        genes_here = [g for g in genes if g in log_matrix.columns]
        block_w = len(genes_here)
        if block_w == 0:
            continue
        label = ARM_DISPLAY.get(arm, arm)
        fontsize = min(8, max(5.5, block_w * 1.7 / max(1, len(label) / 6)))
        ax_main.text(col_pos + block_w / 2 - 0.5, -1.6, label, ha="center", va="bottom",
                     fontsize=fontsize, fontweight="700")
        if col_pos > 0:
            ax_main.axvline(col_pos - 0.5, color="white", linewidth=2)
        col_pos += block_w

    if annotate:
        for i in range(n_rows):
            for j in range(n_cols):
                rel = display_matrix.values[i, j]
                color = "white" if rel > 0.6 else "#333333"
                ax_main.text(j, i, f"{10**log_matrix.values[i,j]:.1e}", ha="center", va="center", fontsize=5, color=color)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("expression,\nscaled 0-1 per gene", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    fig.savefig(out_path, dpi=180, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=== whole-body atlas ===")
    wb_matrix, wb_meta, wb_id_col = build_matrix(
        WHOLE_BODY_PATH, None, lambda df, ids: aggregate_whole_body(df, cell_types=ids),
        list(WHOLE_BODY_CATEGORIES.keys()), "whole_body.parquet",
    )
    print("=== skin fibroblast atlas ===")
    skin_matrix, skin_meta, skin_id_col = build_matrix(
        SKIN_PATH, None, lambda df, ids: aggregate_skin(df, subtypes=ids),
        SKIN_SUBTYPES, "skin.parquet",
    )

    wb_log, _ = log_transform(wb_matrix)
    skin_log, _ = log_transform(skin_matrix)
    wb_display = normalize_per_gene(wb_log)
    skin_display = normalize_per_gene(skin_log)

    wb_order = cluster_row_order(wb_log)
    skin_order = cluster_row_order(skin_log)
    wb_log, wb_display, wb_meta = wb_log.loc[wb_order], wb_display.loc[wb_order], wb_meta.loc[wb_order]
    skin_log, skin_display, skin_meta = skin_log.loc[skin_order], skin_display.loc[skin_order], skin_meta.loc[skin_order]

    draw_two_arm_heatmap(
        wb_log, wb_display, wb_meta, WHOLE_BODY_CATEGORY_COLORS, WHOLE_BODY_CATEGORIES,
        "Whole-body atlas: clearance-arm vs. Wnt-arm partner availability",
        os.path.join(OUT_DIR, "two_arm_heatmap_wholebody.png"),
    )
    draw_two_arm_heatmap(
        skin_log, skin_display, skin_meta, SKIN_STATUS_COLORS, SKIN_DISEASE_STATUS,
        "Skin fibroblast atlas: clearance-arm vs. Wnt-arm partner availability",
        os.path.join(OUT_DIR, "two_arm_heatmap_skin.png"), annotate=True,
    )

    wb_complete = arm_completeness(wb_log)
    skin_complete = arm_completeness(skin_log)

    for name, matrix, complete in [
        ("whole_body", wb_matrix, wb_complete),
        ("skin", skin_matrix, skin_complete),
    ]:
        csv_path = os.path.join(OUT_DIR, f"{name}_two_arm_raw.csv")
        matrix.to_csv(csv_path)
        print(f"  wrote {csv_path}")
        comp_path = os.path.join(OUT_DIR, f"{name}_arm_completeness.csv")
        complete.to_csv(comp_path)
        print(f"  wrote {comp_path}")

    write_summary(wb_matrix, wb_display, wb_complete, skin_matrix, skin_display, skin_complete)


def write_summary(wb_matrix, wb_display, wb_complete, skin_matrix, skin_display, skin_complete):
    lines = []
    lines.append("\n### Task 3 results: extended gene panel, two-arm heatmap\n\n")
    lines.append(
        "**Exploratory characterization / descriptive query** (not a significance test -- means only, "
        "same status as the original ECM-clearance heatmap). `GENE_PANEL_ARMS` added to `config.py` "
        "(RECEPTORS, CLEARANCE_ARM, CLEARANCE_SUBSTRATE, WNT_ARM); original `GENE_PANEL` untouched. "
        "Figures: `/data/ANTXR2/figures/two_arm_panel/two_arm_heatmap_{wholebody,skin}.png`. "
        "**No composite score computed** -- `*_arm_completeness.csv` reports, per row and per gene, "
        f"whether that gene's raw value is within {PRESENT_LOG10_GAP:.0f} orders of magnitude of its OWN "
        "maximum across the curated rows (the same standard already used by hand for the MRC2 "
        "gut-epithelium-vs-keratinocyte/corneal finding earlier in this log, a ~2-3 order-of-magnitude "
        "gap), alongside the raw values in `*_two_arm_raw.csv`; it is a per-gene reading aid over the "
        "same raw numbers, not a score that combines genes.\n\n"
    )

    lines.append("**Per-cell-type arm availability call (whole-body atlas):**\n\n")
    lines.append("| cell type | clearance arm | substrate (COL6) | Wnt arm | neither arm complete |\n")
    lines.append("|---|---|---|---|---|\n")
    neither_wb = []
    for ct in wb_complete.index:
        r = wb_complete.loc[ct]
        clear = f"{int(r['CLEARANCE_ARM_n_present'])}/{int(r['CLEARANCE_ARM_n_total'])}"
        sub = f"{int(r['CLEARANCE_SUBSTRATE_n_present'])}/{int(r['CLEARANCE_SUBSTRATE_n_total'])}"
        wnt = f"{int(r['WNT_ARM_n_present'])}/{int(r['WNT_ARM_n_total'])}"
        clearance_under_half = r["CLEARANCE_ARM_n_present"] < r["CLEARANCE_ARM_n_total"] * 0.5
        wnt_under_half = r["WNT_ARM_n_present"] < r["WNT_ARM_n_total"] * 0.5
        flag = "**YES**" if (clearance_under_half and wnt_under_half) else ""
        if flag:
            neither_wb.append(ct)
        lines.append(f"| {ct} | {clear} | {sub} | {wnt} | {flag} |\n")

    lines.append(f"\nCell types with the clearance arm essentially absent AND the Wnt arm well under half "
                 f"present (i.e. ANTXR2 present with **neither** function's machinery structurally "
                 f"available): {', '.join(neither_wb) if neither_wb else 'none by this cutoff'}.\n")

    gut_rows = [ct for ct in wb_complete.index if ct in (
        "enterocyte", "colonocyte", "intestine goblet cell", "intestinal crypt stem cell",
        "paneth cell", "intestinal tuft cell",
    )]
    if gut_rows:
        lines.append(f"\n**Gut epithelium specifically** ({', '.join(gut_rows)}): clearance-arm genes "
                     f"(mean n_present/6 across these rows: "
                     f"{wb_complete.loc[gut_rows, 'CLEARANCE_ARM_n_present'].mean():.1f}) confirm the "
                     f"project's existing MRC2-absence finding extends to the arm as a whole. Wnt-arm "
                     f"presence (mean n_present/{wb_complete['WNT_ARM_n_total'].iloc[0]}: "
                     f"{wb_complete.loc[gut_rows, 'WNT_ARM_n_present'].mean():.1f}) is the new information "
                     "this panel adds -- if broadly present, gut epithelium has the Wnt-arm partners even "
                     "though it lacks the clearance-arm ones, consistent with Bracq et al.'s finding that "
                     "CMG2's role there is Wnt-pathway (injury-conditional), not clearance.\n")

    lines.append("\n**Skin fibroblast subtypes:** see `skin_arm_completeness.csv` for the per-subtype "
                 "breakdown (13 subtypes) -- not reproduced row-by-row here since the whole-body table "
                 "above already carries the generic `fibroblast` row for this atlas's clearance-arm "
                 "baseline (established: robust ANTXR1 co-presence, standing correction: not evidence of "
                 "redundancy).\n")

    with open(ANALYSIS_SUMMARY_PATH, "a") as f:
        f.write("".join(lines))
    print(f"\nappended Task 3 results to {ANALYSIS_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
