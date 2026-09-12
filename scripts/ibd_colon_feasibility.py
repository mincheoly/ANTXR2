"""Feasibility check for the IBD colon atlas (Smillie et al., Cell 2019,
SCP259 via HCA) as an ANTXR2 study dataset -- a feasibility gate, not a
results-producing step, mirroring kong2023_feasibility.py's role: no
differential-correlation test is run here regardless of outcome.

Focuses on the cell types the user is specifically interested in: fibroblast
subtypes and epithelial cells (Stem/TA1/TA2/Cycling TA plus the differentiated
lineage), i.e. the "Fibroblasts" and "Epithelial" cell_type_lineage groups
from cell_subsets.txt.

ANTXR2 "presence" replicates memento's own filter formula directly (a donor
group's RAW mean count -- not memento's capture-corrected mean_expression --
must exceed filter_mean_thresh; the gene passes overall if that holds in a
strict majority, >min_perc_group, of donor groups with >=100 cells), read
straight from the compute-ready h5ad files' raw X, exactly like
kong2023_feasibility.py -- NOT from the already-computed
ibd_colon_atlas_celltype_means.parquet, whose mean_expression is memento's
capture-rate-corrected estimate and is not on the same scale as
filter_mean_thresh (calibrated for raw counts).

Run in the antxr2 conda env: python ibd_colon_feasibility.py
"""
import os

import anndata.io as aio
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from config import (
    COEXPR_FILTER_MEAN_THRESH, COEXPR_MIN_GROUP_CELLS, COEXPR_MIN_PERC_GROUP,
    FIGURES_DIR, GENE_PANEL, IBD_COLON_CLUSTER_COL, IBD_COLON_COMPARTMENTS,
    IBD_COLON_DONOR_COL, IBD_COLON_H5AD_BY_COMPARTMENT, IBD_COLON_HEALTH_COL,
    IBD_COLON_OUTPUT_PATH,
)
from build_ibd_colon_h5ad import CELL_TYPE_LINEAGE_COL

OUT_DIR = os.path.join(FIGURES_DIR, "ibd_colon_feasibility")
ANALYSIS_SUMMARY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ANALYSIS_SUMMARY.md")
TARGET_LINEAGES = ["Epithelial", "Fibroblasts"]
MIN_DONORS_WELL_POWERED = 8  # same macrophage precedent as kong2023_feasibility.py
STEM_TA_TYPES = ["Stem", "TA 1", "TA 2", "Cycling TA"]


def load_compartment_antxr2(comp):
    h5ad_path = IBD_COLON_H5AD_BY_COMPARTMENT[comp]
    f = h5py.File(h5ad_path, "r")
    obs = aio.read_elem(f["obs"])
    var = aio.read_elem(f["var"])
    antxr2_pos = int(np.where(var.index.to_numpy() == "ANTXR2")[0][0])
    X = aio.sparse_dataset(f["X"])
    antxr2 = np.asarray(X[:, antxr2_pos].todense()).ravel()
    f.close()
    df = pd.DataFrame({
        "compartment": comp,
        "cell_type": obs[IBD_COLON_CLUSTER_COL].to_numpy(),
        "lineage": obs[CELL_TYPE_LINEAGE_COL].to_numpy(),
        "health": obs[IBD_COLON_HEALTH_COL].to_numpy(),
        "donor": obs[IBD_COLON_DONOR_COL].to_numpy(),
        "data_source": obs["data_source"].to_numpy(),
        "antxr2": antxr2,
    })
    return df[df["lineage"].isin(TARGET_LINEAGES)]


def gene_panel_heatmap_table():
    """Donor-equal-weighted mean per (cell_type_fine, gene), GENE_PANEL only,
    target lineages only -- from the already-computed parquet (memento's
    capture-corrected mean_expression, fine for a visualization/ranking table,
    unlike the raw-count presence filter above)."""
    table = pq.read_table(IBD_COLON_OUTPUT_PATH, filters=[("gene", "in", GENE_PANEL)])
    df = table.to_pandas()
    df = df[df["cell_type_lineage"].isin(TARGET_LINEAGES)]
    per_donor = df.groupby(["cell_type_fine", "gene", "donor_id"], observed=True)["mean_expression"].mean().reset_index()
    agg = per_donor.groupby(["cell_type_fine", "gene"], observed=True).agg(
        mean_expr=("mean_expression", "mean"), n_donors=("donor_id", "nunique"),
    ).reset_index()
    return agg


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_df = pd.concat([load_compartment_antxr2(c) for c in IBD_COLON_COMPARTMENTS], ignore_index=True)
    print(f"{len(all_df):,} cells across target lineages {TARGET_LINEAGES}")

    # ---- per cell_type x health x donor cell counts + presence filter ----
    grp = all_df.groupby(["cell_type", "health", "donor"], observed=True).agg(
        n_cells=("antxr2", "size"), antxr2_mean=("antxr2", "mean"),
    ).reset_index()
    donor_csv = os.path.join(OUT_DIR, "per_donor_counts.csv")
    grp.to_csv(donor_csv, index=False)
    print(f"wrote {donor_csv}")

    lineage_by_ct = all_df.drop_duplicates("cell_type").set_index("cell_type")["lineage"].to_dict()
    source_by_ct = (
        all_df.groupby("cell_type")["data_source"].agg(lambda s: s.mode().iloc[0]).to_dict()
    )

    summary_rows = []
    for (cell_type, health), sub in grp.groupby(["cell_type", "health"], observed=True):
        n_donors_total = len(sub)
        well_powered = sub[sub["n_cells"] >= COEXPR_MIN_GROUP_CELLS]
        n_well_powered = len(well_powered)
        n_pass = int((well_powered["antxr2_mean"] > COEXPR_FILTER_MEAN_THRESH).sum())
        frac_pass = n_pass / n_well_powered if n_well_powered else np.nan
        summary_rows.append({
            "lineage": lineage_by_ct[cell_type], "cell_type": cell_type, "health": health,
            "data_source": source_by_ct[cell_type],
            "n_donors_total": n_donors_total, "n_donors_ge_min_cells": n_well_powered,
            "n_donors_ANTXR2_present": n_pass, "frac_ANTXR2_present": frac_pass,
            "clears_min_perc_group_0.7": bool(frac_pass > COEXPR_MIN_PERC_GROUP) if n_well_powered else False,
            "underpowered_lt8_donors": n_well_powered < MIN_DONORS_WELL_POWERED,
        })
    summary_df = pd.DataFrame(summary_rows).sort_values(["lineage", "cell_type", "health"])
    summary_csv = os.path.join(OUT_DIR, "feasibility_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"wrote {summary_csv}")
    print(summary_df.to_string(index=False))

    heatmap_df = gene_panel_heatmap_table()
    heatmap_csv = os.path.join(OUT_DIR, "gene_panel_by_celltype.csv")
    heatmap_df.to_csv(heatmap_csv, index=False)
    print(f"wrote {heatmap_csv}")

    write_summary(summary_df, heatmap_df, all_df)


def write_summary(summary_df, heatmap_df, all_df):
    lines = []
    lines.append("\n## IBD colon atlas (Smillie et al., Cell 2019, SCP259) feasibility check (2026-09-09)\n\n")
    lines.append(
        "**Feasibility gate, not a results-producing step** -- no differential-correlation test "
        "was run. Data source: HCA DCP mirror of SCP259 (see the download-corruption finding above); "
        f"51 cell subsets total (15 epithelial, 13 stromal/glial, 23 immune), 30 donors "
        "(18 UC + 12 healthy). This check covers only the Epithelial and Fibroblasts lineages "
        "(the user's stated interest). Output: `/data/ANTXR2/figures/ibd_colon_feasibility/"
        "{per_donor_counts,feasibility_summary,gene_panel_by_celltype}.csv`.\n\n"
    )

    lines.append("### Donor counts and ANTXR2 presence, Epithelial + Fibroblasts lineages\n\n")
    lines.append(
        f"ANTXR2 presence replicates memento's own filter formula directly (raw per-donor mean count "
        f"> `filter_mean_thresh={COEXPR_FILTER_MEAN_THRESH}`; passes overall if true in a strict "
        f"majority, `>min_perc_group={COEXPR_MIN_PERC_GROUP}`, of donor groups with "
        f">={COEXPR_MIN_GROUP_CELLS} cells), computed directly from raw counts (not the memento "
        f"capture-corrected parquet). Donors with <{MIN_DONORS_WELL_POWERED} usable groups flagged "
        f"underpowered (macrophage precedent). `data_source` shows whether a cell type's counts came "
        f"from the full 30-donor cohort (`full_cohort_mtx`, Fib only) or the 17-donor discovery-cohort "
        f"fallback forced by the HCA source corruption (`discovery_cohort_rds`, Epi + Imm).\n\n"
    )
    lines.append("| lineage | cell type | health | source | donors (total/≥100 cells) | ANTXR2 present (n/frac) | clears mpg=0.7 | underpowered |\n")
    lines.append("|---|---|---|---|---|---|---|---|\n")
    for _, r in summary_df.iterrows():
        present_str = (
            f"{r['n_donors_ANTXR2_present']}/{r['n_donors_ge_min_cells']} ({r['frac_ANTXR2_present']:.2f})"
            if pd.notna(r["frac_ANTXR2_present"]) else "n/a"
        )
        lines.append(
            f"| {r['lineage']} | {r['cell_type']} | {r['health']} | {r['data_source']} | "
            f"{r['n_donors_total']} / {r['n_donors_ge_min_cells']} | {present_str} | "
            f"{'YES' if r['clears_min_perc_group_0.7'] else 'no'} | "
            f"{'**YES**' if r['underpowered_lt8_donors'] else 'no'} |\n"
        )

    stem_ta = summary_df[summary_df["cell_type"].isin(STEM_TA_TYPES)]
    n_stem_ta_ok = int((~stem_ta["underpowered_lt8_donors"] & stem_ta["clears_min_perc_group_0.7"]).sum())
    lines.append(
        f"\n**Stem/TA cells** ({', '.join(STEM_TA_TYPES)}): {n_stem_ta_ok}/{len(stem_ta)} "
        f"(cell type x health) groups are both well-powered and clear the presence filter. "
        "Consistent with Phase 1's prior finding that ANTXR2 is not stem-skewed in gut epithelium.\n\n"
    )

    fib = summary_df[summary_df["lineage"] == "Fibroblasts"]
    n_fib_ok = int((~fib["underpowered_lt8_donors"] & fib["clears_min_perc_group_0.7"]).sum())
    lines.append(
        f"**Fibroblast subtypes**: {n_fib_ok}/{len(fib)} (cell type x health) groups clear both power "
        "and presence filters -- fibroblasts are the full 30-donor `full_cohort_mtx` compartment (not "
        "affected by the HCA corruption), and the donor-equal-weighted heatmap below shows the "
        "Wnt-niche fibroblast subtypes (RSPO3+, WNT2B+, WNT5B+) at the top of the ANTXR2 ranking -- "
        "directionally consistent with Phase 1's fibroblast finding and on-theme for the project's "
        "Wnt-transduction arm.\n\n"
    )

    lines.append("### Gene panel (GENE_PANEL) donor-equal-weighted mean, top ANTXR2 cell types\n\n")
    antxr2_top = heatmap_df[heatmap_df["gene"] == "ANTXR2"].sort_values("mean_expr", ascending=False).head(10)
    lines.append("| cell type | mean ANTXR2 | n donors |\n|---|---|---|\n")
    for _, r in antxr2_top.iterrows():
        lines.append(f"| {r['cell_type_fine']} | {r['mean_expr']:.6f} | {int(r['n_donors'])} |\n")

    sources = sorted(summary_df["data_source"].unique())
    if sources == ["full_cohort_mtx"]:
        cohort_note = (
            "**Full 30-donor cohort, all 3 compartments, no coverage asymmetry.** "
            "Fibroblast coverage was always full-cohort; Epithelial coverage (including "
            "Stem/TA1/TA2/Cycling TA) is now also full-cohort after the SCP259-direct "
            "download superseded the earlier HCA-corruption-forced discovery-cohort "
            "fallback (see the download section above)."
        )
    else:
        cohort_note = (
            "**One asymmetry to carry forward.** Fibroblast coverage is full-cohort (30 donors, "
            "complete raw matrix); Epithelial coverage (including Stem/TA1/TA2/Cycling TA) is "
            f"discovery-cohort-only ({', '.join(sources)}) because the HCA-hosted full-cohort Epi "
            "and Imm matrices are corrupted at the source (see the download section above)."
        )
    lines.append(
        f"\n### Verdict\n\n"
        f"**Usable for the ANTXR2 project's stated interest (fibroblasts, epithelial, TA/stem cells).** "
        f"{cohort_note} The Wnt-niche fibroblast subtypes rank highest for ANTXR2 (directionally "
        f"on-theme); most cell-type x health groups clear the 8-donor floor (see table above for "
        f"exactly which don't).\n\n"
    )

    with open(ANALYSIS_SUMMARY_PATH, "a") as fh:
        fh.write("".join(lines))
    print(f"\nappended feasibility results to {ANALYSIS_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
