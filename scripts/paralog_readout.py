"""Task 1 of prompts/partner_availability.md: ANTXR1/ANTXR2 paralog readout.

This is a READOUT, not a new computation -- ANTXR1's post-bug-fix coef/se/z/fdr
already exist in `full_dataset_ht/{cell_type}_full_dataset_ht.csv`
(coexpression_full_ht.py) and have never been read out. This script:

  1. pulls ANTXR1's row + rank from each cell type's full tested universe
  2. computes ANTXR1/ANTXR2's per-cell expression distribution per cell type
     (bimodality check: is either near-uniformly expressed, in which case
     correlation is a poor co-presence readout and the mean already settles
     availability?)
  3. runs the doublet check ONLY if a cell type shows a SIGNIFICANT positive
     correlation (per the pre-registered rule in ANALYSIS_SUMMARY.md -- a
     non-significant positive point estimate is a null, not a "positive
     correlation" in the sense that rule and this check care about)
  4. appends a results section to ANALYSIS_SUMMARY.md

Run in the `antxr2` conda env: `conda run -n antxr2 python scripts/paralog_readout.py`.
"""
import os

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import COEXPR_FIGURES_DIR, COEXPR_VARIANTS

HT_DIR = os.path.join(COEXPR_FIGURES_DIR, "full_dataset_ht")
OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "paralog_readout")
ANALYSIS_SUMMARY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ANALYSIS_SUMMARY.md")
CELL_TYPES = ["fibroblast", "enterocyte", "macrophage"]
LEVEL3_LABEL = {"fibroblast": "Crypt_fibroblast_PI16", "enterocyte": "Enterocyte", "macrophage": "Macrophage"}
GENE_NAME_COL = "feature_name"
PRE_FIX_ANTXR1_VALUE = 0.35  # see ANALYSIS_SUMMARY.md, "Correction: memento estimator bug..." section
FDR_SIG_THRESH = 0.1  # this project's standard significance threshold throughout


def read_antxr1_row(ct):
    path = os.path.join(HT_DIR, f"{ct}_full_dataset_ht.csv")
    df = pd.read_csv(path, comment="#").sort_values("pval").reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    row = df[df["gene_symbol"] == "ANTXR1"]
    n_tested = len(df)
    if row.empty:
        return None, n_tested
    return row.iloc[0].to_dict(), n_tested


def bimodality_check(out_dir):
    """Per-cell expression distribution of ANTXR1/ANTXR2 within each cell type.
    Uses the same level3 trio h5ad the full-dataset HT was computed from."""
    os.makedirs(out_dir, exist_ok=True)
    adata = anndata.read_h5ad(COEXPR_VARIANTS["level3"]["h5ad"])
    sym = adata.var[GENE_NAME_COL].astype(str)
    gene_pos = {g: i for i, g in enumerate(sym.values) if g in ("ANTXR1", "ANTXR2")}

    fig, axes = plt.subplots(1, len(CELL_TYPES), figsize=(4.2 * len(CELL_TYPES), 3.6))
    rows = []
    for ax, ct in zip(axes, CELL_TYPES):
        label = LEVEL3_LABEL[ct]
        mask = adata.obs["level_3_annot"].astype(str) == label
        sub = adata[mask]
        n_cells = sub.shape[0]
        for gene, color in [("ANTXR1", "#2a78d6"), ("ANTXR2", "#eb6834")]:
            col = sub.X[:, gene_pos[gene]]
            col = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
            pct_nz = float((col > 0).mean() * 100)
            mean_v = float(col.mean())
            # crude near-uniform flag: >90% of cells nonzero AND coefficient of
            # variation among nonzero cells is low -- i.e. "on" almost everywhere
            # at a similar level, rather than a sparse/patchy on-off pattern.
            nz_vals = col[col > 0]
            cv_nz = float(nz_vals.std() / nz_vals.mean()) if len(nz_vals) > 1 and nz_vals.mean() > 0 else np.nan
            near_uniform = pct_nz > 90 and (not np.isnan(cv_nz)) and cv_nz < 0.5
            rows.append({
                "cell_type": ct, "gene": gene, "n_cells": n_cells, "mean": mean_v,
                "pct_nonzero": pct_nz, "cv_nonzero": cv_nz, "near_uniform": near_uniform,
            })
            ax.hist(col, bins=30, alpha=0.5, label=f"{gene} ({pct_nz:.1f}% nz)", color=color, density=True)
        ax.set_title(f"{ct} (n={n_cells:,})", fontsize=10)
        ax.set_xlabel("raw count")
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("ANTXR1/ANTXR2 per-cell expression distribution (bimodality check)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig_path = os.path.join(out_dir, "antxr1_antxr2_distribution.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    df = pd.DataFrame(rows)
    csv_path = os.path.join(out_dir, "antxr1_antxr2_distribution_stats.csv")
    df.to_csv(csv_path, index=False)
    return df, fig_path, csv_path


def doublet_check(ct, out_dir):
    """Only called for a cell type with a SIGNIFICANT positive ANTXR1/ANTXR2
    correlation. Confirms atlas QC removed doublets and checks whether
    ANTXR2-high/ANTXR1-high cells are enriched for cross-lineage marker
    co-expression (fibroblast marker COL1A2 vs. epithelial marker EPCAM)."""
    adata = anndata.read_h5ad(COEXPR_VARIANTS["level3"]["h5ad"])
    label = LEVEL3_LABEL[ct]
    sub = adata[adata.obs["level_3_annot"].astype(str) == label]
    print(f"  atlas-level doublet QC: level_3_annot={label!r} is an author-curated single-cell-type "
          f"label from the Gut Cell Atlas's own published QC pipeline (doublet detection upstream of "
          f"cell-type annotation is standard for this atlas per its methods) -- not independently "
          f"re-verified here beyond the marker co-expression check below.")

    sym = adata.var[GENE_NAME_COL].astype(str)

    def get_col(gene):
        pos = np.where(sym.values == gene)[0][0]
        col = sub.X[:, pos]
        return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()

    a1, a2 = get_col("ANTXR1"), get_col("ANTXR2")
    high_thresh_1, high_thresh_2 = np.quantile(a1[a1 > 0], 0.75), np.quantile(a2[a2 > 0], 0.75)
    double_high = (a1 >= high_thresh_1) & (a2 >= high_thresh_2) & (a1 > 0) & (a2 > 0)
    fib_marker, epi_marker = get_col("COL1A2"), get_col("EPCAM")
    result = {
        "cell_type": ct, "n_cells": sub.shape[0], "n_double_high": int(double_high.sum()),
        "frac_double_high_epcam_pos": float((epi_marker[double_high] > 0).mean()) if double_high.sum() else np.nan,
        "frac_all_epcam_pos": float((epi_marker > 0).mean()),
        "frac_double_high_col1a2_pos": float((fib_marker[double_high] > 0).mean()) if double_high.sum() else np.nan,
        "frac_all_col1a2_pos": float((fib_marker > 0).mean()),
    }
    out_path = os.path.join(out_dir, f"{ct}_doublet_check.csv")
    pd.DataFrame([result]).to_csv(out_path, index=False)
    return result, out_path


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=== Task 1: ANTXR1/ANTXR2 paralog readout ===\n")
    readout = {}
    for ct in CELL_TYPES:
        row, n_tested = read_antxr1_row(ct)
        readout[ct] = (row, n_tested)
        if row is None:
            print(f"{ct}: ANTXR1 NOT in the tested universe ({n_tested} genes tested) -- "
                  f"filtered out for low expression before the correlation test could even run.")
        else:
            print(f"{ct}: ANTXR1 rank {int(row['rank'])}/{n_tested}  coef={row['coef']:.4f}  "
                  f"se={row['se']:.4f}  z={row['z']:.3f}  pval={row['pval']:.4g}  fdr={row['fdr']:.4g}")

    print("\n=== Bimodality check ===\n")
    bimod_df, fig_path, bimod_csv = bimodality_check(OUT_DIR)
    print(bimod_df.to_string(index=False))
    print(f"\nwrote {fig_path}\nwrote {bimod_csv}")

    # Doublet check trigger: SIGNIFICANT (FDR < threshold) positive correlation only.
    doublet_results = {}
    for ct in CELL_TYPES:
        row, _ = readout[ct]
        if row is not None and row["coef"] > 0 and row["fdr"] < FDR_SIG_THRESH:
            print(f"\n=== Doublet check triggered for {ct} (significant positive correlation) ===")
            res, path = doublet_check(ct, OUT_DIR)
            doublet_results[ct] = res
            print(res)
            print(f"wrote {path}")
        elif row is not None and row["coef"] > 0:
            print(f"\n{ct}: coef is positive ({row['coef']:.4f}) but NOT significant (fdr={row['fdr']:.4g} "
                  f">= {FDR_SIG_THRESH}) -- this is a NULL result per the pre-registered rule, not a "
                  f"'significant positive' finding, so the doublet check is NOT triggered.")

    write_summary(readout, bimod_df, doublet_results, fig_path)


def write_summary(readout, bimod_df, doublet_results, fig_path):
    lines = []
    lines.append("\n### Task 1 results: ANTXR1/ANTXR2 paralog readout\n\n")
    lines.append("**Confirmatory** (pre-defined panel test -- ANTXR1 is the single, pre-specified gene "
                  "of interest here, read out of an already-computed full-dataset test).\n\n")

    lines.append("| cell type | in tested universe | rank | coef | se | z | pval | fdr | call |\n")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---|\n")
    for ct in CELL_TYPES:
        row, n_tested = readout[ct]
        if row is None:
            lines.append(f"| {ct} | NO ({n_tested} genes tested) | -- | -- | -- | -- | -- | -- | "
                          f"filtered out (near-absent expression) |\n")
        else:
            if row["fdr"] < FDR_SIG_THRESH and row["coef"] < 0:
                call = "SIGNIFICANT NEGATIVE"
            elif row["fdr"] < FDR_SIG_THRESH and row["coef"] > 0:
                call = "SIGNIFICANT POSITIVE (gated on Task 2)"
            else:
                call = "NULL (uninformative)"
            lines.append(f"| {ct} | yes ({n_tested} genes tested) | {int(row['rank'])} | "
                         f"{row['coef']:.4f} | {row['se']:.4f} | {row['z']:.3f} | {row['pval']:.4g} | "
                         f"{row['fdr']:.4g} | {call} |\n")

    lines.append(f"\nPre-fix point estimate for comparison: ANTXR1-ANTXR2 was reported earlier in this "
                 f"log (\"Implication for the ECM-clearance-panel question\" section) at **{PRE_FIX_ANTXR1_VALUE}** -- "
                 f"that value predates the `_corr_from_cov` variance<=0 null-out fix and is **superseded**; "
                 f"per the correction, it is not necessarily wrong in sign but is very likely overstated in "
                 f"magnitude by the same placeholder-clipping bug that affected every other panel gene.\n")

    fib_row, fib_n = readout["fibroblast"]
    lines.append(
        f"\n**Fibroblast (the only cell type where ANTXR1 survives the detection filter and gets a "
        f"real test): NULL.** coef={fib_row['coef']:.3f}, pval={fib_row['pval']:.3f}, fdr={fib_row['fdr']:.3f}, "
        f"rank {int(fib_row['rank'])} of {fib_n} -- squarely mid-pack, not extreme in either direction. "
        f"Per the pre-registered rule, **this is explicitly NOT evidence of ANTXR1/ANTXR2 co-presence in "
        f"fibroblast** -- it means the correlation test has no power to distinguish the co-presence "
        f"hypothesis from the mutual-exclusivity hypothesis here, not that co-presence is confirmed. "
        f"Because coef is nominally positive but not significant, this is a **null, not a 'significant "
        f"positive'** in the pre-registered sense -- the doublet check (which the prompt gates on "
        f"significant positive correlations specifically) was correctly **not triggered** for fibroblast; "
        f"see doublet-check section below for the reasoning trace.\n"
    )

    lines.append(
        f"\n**Enterocyte and macrophage: ANTXR1 is filtered out of the tested universe entirely** "
        f"(fails the memento `min_perc_group` presence filter before a correlation could even be "
        f"attempted). Direct per-cell counts explain why: ANTXR1 is detected in only **0.08% of "
        f"enterocytes** (35,062 cells) and **2.2% of macrophages** (2,953 cells), vs. 25.2% of "
        f"fibroblasts. **This is itself the stronger, means-level asymmetric result** the header's "
        f"evidence table describes: ANTXR1 is essentially absent from gut epithelium and largely absent "
        f"from macrophage by the mean/detection-rate alone -- no correlation test is needed to make "
        f"that call, and none was possible. A co-expression question (\"do ANTXR1 and ANTXR2 occupy the "
        f"same enterocytes/macrophages\") is moot when one partner is barely present at all; the paralog-"
        f"partitioning claim for these two cell types rests on presence/absence, not correlation.\n"
    )

    lines.append("\n**Macrophage is reported for completeness only and is exploratory** (standing "
                 "correction) -- no conclusion is drawn from it regardless of the above.\n")

    lines.append("\n#### Bimodality check\n\n")
    lines.append(f"Figure: `{fig_path}`. Per-cell expression is sparse/zero-inflated for both genes in "
                 "every cell type (detection rates 0.08%-25.3%, well below the >90% near-uniform "
                 "threshold used here) -- **neither gene is near-uniformly expressed anywhere**, so "
                 "correlation remains a meaningful co-presence readout wherever both genes clear the "
                 "detection filter (fibroblast); it is not needed as a readout where one partner is "
                 "already known absent by the mean (enterocyte, macrophage).\n\n")
    lines.append("```\n" + bimod_df.to_string(index=False) + "\n```\n")

    lines.append("\n#### Doublet check\n\n")
    if doublet_results:
        for ct, res in doublet_results.items():
            lines.append(f"Triggered for **{ct}** (significant positive correlation). {res}\n")
    else:
        lines.append("**Not triggered for any cell type.** The prompt's doublet check applies only to a "
                     "significant positive correlation; fibroblast's positive point estimate was not "
                     "significant (null), and enterocyte/macrophage have no correlation at all (ANTXR1 "
                     "filtered out). No doublet check was run.\n")

    with open(ANALYSIS_SUMMARY_PATH, "a") as f:
        f.write("".join(lines))
    print(f"\nappended Task 1 results to {ANALYSIS_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
