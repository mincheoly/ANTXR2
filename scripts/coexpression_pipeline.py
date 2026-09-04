"""ANTXR2 co-expression pipeline (Phase 2 test run): memento point-estimate
correlations on the macrophage/enterocyte/fibroblast trio.

Continues directly from `subset_coexpression_data.py` (see
planning_summaries/dataset_subsetting.md), which built the working set but
deliberately stopped before any correlation computation. This script:

  1. loads and validates the "level3" trio h5ad, drops any donor x cell_type
     group below COEXPR_MIN_GROUP_CELLS cells
  2. computes ANTXR2-vs-every-other-gene point-estimate correlations, per
     donor x cell_type group (memento's compute_2d_moments, no bootstrap)
  3. averages those across donors, per cell type (unweighted mean over
     non-NaN donor values -- same donor-equal-weighting convention as
     Phase 1's "binding aggregation convention" in ANALYSIS_SUMMARY.md)
  4. takes the top 50 genes by |donor-averaged correlation| per cell type,
     and their union (<=151 genes incl. ANTXR2)
  5. recomputes the full pairwise correlation matrix over that union panel,
     again per donor x cell_type group (memento.get_corr_matrix)
  6. averages those matrices across donors, per cell type

Point estimates only throughout -- no bootstrap, no p-values (memento's
ht_*/binary_test_* functions are never called). Writes a NEW h5ad
(elmentaite2021_trio_level3.h5ad, the source prepared by
subset_coexpression_data.py, is left untouched), one heatmap per cell type,
and appends a summary section to ANALYSIS_SUMMARY.md.

Gene identity gotcha: adata.var.index is Ensembl IDs, not symbols. ANTXR2's
Ensembl ID is resolved from var['feature_name'] at the start and used
throughout; symbols are only for display (top-gene tables, heatmap labels).

Run in the `antxr2` conda env: `conda run -n antxr2 python
scripts/coexpression_pipeline.py`.
"""
import itertools
import os
import sys
import time
import warnings
from contextlib import contextmanager
from importlib.metadata import version as pkg_version

import anndata
import matplotlib.pyplot as plt
import memento
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from config import (
    COEXPR_CELL_TYPES, COEXPR_DONOR_COL, COEXPR_FILTER_MEAN_THRESH,
    COEXPR_FIGURES_DIR, COEXPR_MIN_GROUP_CELLS, COEXPR_MIN_PERC_GROUP,
    COEXPR_OUTPUT_H5AD, COEXPR_TARGET_GENE, COEXPR_TOP_N_GENES, COEXPR_VARIANTS,
)

INPUT_H5AD = COEXPR_VARIANTS["level3"]["h5ad"]
GENE_NAME_COL = "feature_name"
Q_COLUMN = "capture_rate_pbmc"
ANALYSIS_SUMMARY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ANALYSIS_SUMMARY.md")

TIMINGS = {}


@contextmanager
def step_timer(name):
    t0 = time.time()
    print(f"\n=== {name} ===")
    yield
    dt = time.time() - t0
    TIMINGS[name] = dt
    print(f"  [{name}] done in {dt:.1f}s")


# ---------------------------------------------------------------------------
# Step 1: load + validate + filter small groups
# ---------------------------------------------------------------------------

def load_and_filter(h5ad_path):
    adata = anndata.read_h5ad(h5ad_path)
    print(f"  loaded {h5ad_path}: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    for col in ["cell_type", Q_COLUMN, COEXPR_DONOR_COL]:
        assert col in adata.obs.columns, f"obs missing required column {col!r}"

    sym_to_id = dict(zip(adata.var[GENE_NAME_COL].astype(str), adata.var.index))
    if COEXPR_TARGET_GENE not in sym_to_id:
        raise SystemExit(f"{COEXPR_TARGET_GENE} not found in var[{GENE_NAME_COL}]")
    target_id = sym_to_id[COEXPR_TARGET_GENE]
    print(f"  {COEXPR_TARGET_GENE} -> {target_id} (var.index is Ensembl ID, not symbol)")

    print("\n  cell counts per cell_type:")
    ct_counts = adata.obs["cell_type"].value_counts()
    print(ct_counts.to_string())

    # Raw-counts check: same heuristic as compute_means.py's resolve_raw_slot
    # (max > 20 and >99% of a sample look integer), cross-checked against this
    # file's own provenance record.
    X = adata.X
    sample = X.data[:5000] if X.nnz else np.array([])
    frac_int = float(np.mean(np.isclose(sample, np.round(sample)))) if sample.size else 1.0
    x_max = float(X.max())
    looks_raw = x_max > 20 and frac_int > 0.99
    prov_slot = adata.uns.get("provenance", {}).get("source_slot")
    print(f"\n  X: dtype={X.dtype}, max={x_max:.1f}, frac_integer(sample)={frac_int:.4f}, "
          f"looks_raw={looks_raw}, uns.provenance.source_slot={prov_slot!r}")
    if not looks_raw:
        raise SystemExit("adata.X does not look like raw counts -- aborting")

    print(f"\n  {Q_COLUMN} distribution per cell_type:")
    print(adata.obs.groupby("cell_type")[Q_COLUMN].describe().to_string())
    n_distinct = adata.obs[Q_COLUMN].nunique()
    print(f"  NOTE: {Q_COLUMN} takes only {n_distinct} distinct value(s) across all "
          f"cells -- it is a per-assay constant (10x 3' v2 vs 5' v2), not a "
          f"continuous per-cell estimate (see subset_coexpression_data.py for how "
          f"it was derived). memento still consumes it correctly via q_column, "
          f"taking a per-group mean internally.")
    if adata.obs[Q_COLUMN].isna().any():
        raise SystemExit(f"{Q_COLUMN} has NaNs")

    donor = adata.obs[COEXPR_DONOR_COL].astype(str)
    cell_type = adata.obs["cell_type"].astype(str)
    crosstab = pd.crosstab(donor, cell_type)
    print(f"\n  donor x cell_type crosstab ({crosstab.shape[0]} donors x "
          f"{crosstab.shape[1]} cell types):")
    print(crosstab.to_string())

    long = crosstab.stack().rename("n_cells").reset_index()
    long.columns = [COEXPR_DONOR_COL, "cell_type", "n_cells"]
    long = long[long["n_cells"] > 0].reset_index(drop=True)
    dropped = long[long["n_cells"] < COEXPR_MIN_GROUP_CELLS].reset_index(drop=True)
    kept = long[long["n_cells"] >= COEXPR_MIN_GROUP_CELLS].reset_index(drop=True)

    print(f"\n  {len(dropped)} of {len(long)} existing donor x cell_type groups fall "
          f"below the {COEXPR_MIN_GROUP_CELLS}-cell floor and are DROPPED:")
    if len(dropped):
        print(dropped.to_string(index=False))
    print(f"\n  {len(kept)} groups kept, by cell type:")
    print(kept.groupby("cell_type").size().to_string())

    keep_index = pd.MultiIndex.from_frame(kept[[COEXPR_DONOR_COL, "cell_type"]])
    pair_index = pd.MultiIndex.from_arrays([donor.values, cell_type.values])
    keep_mask = pair_index.isin(keep_index)
    adata_f = adata[keep_mask].copy()
    for col in ["cell_type", COEXPR_DONOR_COL]:
        if isinstance(adata_f.obs[col].dtype, pd.CategoricalDtype):
            adata_f.obs[col] = adata_f.obs[col].cat.remove_unused_categories()
    print(f"\n  {adata_f.shape[0]:,} cells remain after dropping small groups "
          f"(from {adata.shape[0]:,})")

    report = {
        "cell_counts": ct_counts,
        "capture_rate_describe": adata.obs.groupby("cell_type")[Q_COLUMN].describe(),
        "capture_rate_n_distinct": n_distinct,
        "crosstab": crosstab,
        "dropped_groups": dropped,
        "kept_groups": kept,
    }
    return adata_f, report, target_id, sym_to_id


# ---------------------------------------------------------------------------
# Step 2: ANTXR2 vs. every other filtered gene, per donor x cell_type group
# ---------------------------------------------------------------------------

def compute_antxr2_correlations(adata, target_id):
    adata = adata.copy()
    adata.X = sp.csr_matrix(adata.X)
    original_var_index = adata.var.index.copy()

    memento.setup_memento(
        adata, q_column=Q_COLUMN, filter_mean_thresh=COEXPR_FILTER_MEAN_THRESH,
        min_cell_count=COEXPR_MIN_GROUP_CELLS,
    )
    memento.create_groups(adata, label_columns=[COEXPR_DONOR_COL, "cell_type"])
    groups = adata.uns["memento"]["groups"]
    group_meta = memento.get_groups(adata)
    group_n_cells = {g: adata.uns["memento"]["group_cells"][g].shape[0] for g in groups}
    print(f"  {len(groups)} donor x cell_type groups after create_groups "
          f"(min_cell_count={COEXPR_MIN_GROUP_CELLS})")

    memento.compute_1d_moments(adata, min_perc_group=COEXPR_MIN_PERC_GROUP, filter_genes=True)

    # NB: adata.uns['memento']['gene_filter'] is never resliced after
    # compute_1d_moments's internal _inplace_subset_var call -- it stays keyed
    # on `original_var_index` (18,370 genes) even though adata.var has since
    # shrunk to the globally-filtered gene_list. Must index it with the ORIGINAL
    # index, not the post-filter one.
    gene_filter = adata.uns["memento"]["gene_filter"]
    per_group_gene_counts = {g: int(mask.sum()) for g, mask in gene_filter.items()}
    antxr2_orig_pos = original_var_index.get_loc(target_id)
    antxr2_survival = {g: bool(mask[antxr2_orig_pos]) for g, mask in gene_filter.items()}

    gene_list = adata.uns["memento"]["gene_list"]
    print(f"\n  global filtered gene_list: {len(gene_list)} / {len(original_var_index)} genes "
          f"(filter_mean_thresh={COEXPR_FILTER_MEAN_THRESH}, min_perc_group={COEXPR_MIN_PERC_GROUP})")
    print("  per-group genes passing filter_mean_thresh, and ANTXR2 survival:")
    for g in groups:
        surv = "OK" if antxr2_survival[g] else "low-expression"
        print(f"    {g:40s} n_cells={group_n_cells[g]:6d}  genes_passing={per_group_gene_counts[g]:6d}  ANTXR2={surv}")

    # NB: `antxr2_survival[g]` is False whenever ANTXR2 fails the per-group MEAN
    # threshold (obs_mean > filter_mean_thresh) -- this does NOT necessarily mean
    # its variance is zero. compute_2d_moments only NaNs a correlation when
    # var<=0 exactly (memento/estimator.py's _corr_from_cov); a gene that's just
    # barely/rarely detected in a group can still have var>0 (verified: e.g.
    # sg^F3^enterocyte has ANTXR2 mean=4.6e-6, var=8.3e-11 -- both nonzero) and so
    # still produces a NUMERIC correlation there, not NaN. That correlation is
    # real output, not an error, but it's derived from near-zero signal and is
    # correspondingly noisy -- flagged here as "low-expression" so it's visible,
    # not because it silently drops out downstream.
    low_expr_groups = [g for g, ok in antxr2_survival.items() if not ok]
    if low_expr_groups:
        print(f"\n  NOTE: ANTXR2 is below its OWN group's mean-expression threshold (but "
              f"survives the global >{COEXPR_MIN_PERC_GROUP*100:.0f}% threshold) in "
              f"{len(low_expr_groups)} group(s): {low_expr_groups} -- these groups still "
              f"produce a numeric correlation (memento only NaNs when variance is exactly "
              f"zero, which is rarer than failing the mean filter), but that correlation "
              f"is derived from near-undetected expression and should be treated as noisier "
              f"than groups where ANTXR2 clears its own filter.")

    if target_id not in gene_list:
        raise SystemExit("ANTXR2 filtered out of the GLOBAL gene_list -- cannot compute step 3")

    gene_pairs = [(target_id, g) for g in gene_list if g != target_id]
    print(f"\n  {len(gene_pairs)} gene pairs (ANTXR2 vs. every other filtered gene) x "
          f"{len(groups)} groups")
    memento.compute_2d_moments(adata, gene_pairs)
    moment_corr_df, _cell_counts = memento.get_2d_moments(adata)
    antxr2_gene_order = moment_corr_df["gene_2"].tolist()

    print("\n  NaN correlations per group (failed estimates, e.g. from sparsity):")
    by_donor = {ct: {} for ct in COEXPR_CELL_TYPES}
    for g in groups:
        donor = group_meta.loc[g, COEXPR_DONOR_COL]
        ct = group_meta.loc[g, "cell_type"]
        vals = moment_corr_df[g].values.astype(float)
        n_nan = int(np.isnan(vals).sum())
        print(f"    {g:40s} NaN: {n_nan}/{len(vals)}")
        by_donor.setdefault(ct, {})[donor] = {
            "n_cells": group_n_cells[g],
            "antxr2_vs_all": vals,
            "genes_passing_group_filter": per_group_gene_counts[g],
            "antxr2_survives_group_filter": antxr2_survival[g],
        }

    results = {
        "antxr2_gene_order": antxr2_gene_order,
        "by_donor": by_donor,
        "gene_list_size": len(gene_list),
        "gene_list_total": len(original_var_index),
    }
    return results


# ---------------------------------------------------------------------------
# Step 3: average ANTXR2-vs-all correlations across donors, per cell type
# ---------------------------------------------------------------------------

def average_across_donors_1d(results):
    antxr2_gene_order = results["antxr2_gene_order"]
    donor_avg = {}
    for ct in COEXPR_CELL_TYPES:
        donors = results["by_donor"].get(ct, {})
        if not donors:
            print(f"  {ct}: 0 surviving donor groups -- skipping")
            continue
        stacked = np.vstack([d["antxr2_vs_all"] for d in donors.values()])  # donors x genes
        mean_corr = np.nanmean(stacked, axis=0)
        n_donors = np.sum(~np.isnan(stacked), axis=0)
        print(f"  {ct}: {len(donors)} donor groups contribute; "
              f"genes with full donor coverage: {int((n_donors == len(donors)).sum())}/{len(antxr2_gene_order)}, "
              f"median n_donors/gene: {np.median(n_donors):.0f}")
        donor_avg[ct] = {
            "mean_corr": mean_corr,
            "n_donors": n_donors,
            "n_donor_groups_total": len(donors),
        }
    results["donor_averaged_1d"] = donor_avg
    return results


# ---------------------------------------------------------------------------
# Step 4: top-N genes per cell type + union panel
# ---------------------------------------------------------------------------

def select_top_genes(results, sym_map, target_id, top_n=COEXPR_TOP_N_GENES):
    antxr2_gene_order = np.array(results["antxr2_gene_order"])
    top_tables = {}
    for ct, d in results["donor_averaged_1d"].items():
        df = pd.DataFrame({
            "gene_id": antxr2_gene_order,
            "gene_symbol": [sym_map.get(g, g) for g in antxr2_gene_order],
            "mean_corr": d["mean_corr"],
            "n_donors": d["n_donors"],
        })
        df = df[df["gene_id"] != target_id].copy()
        df["abs_corr"] = df["mean_corr"].abs()
        df = df.dropna(subset=["mean_corr"])
        # deterministic tiebreak on abs_corr ties: gene_id ascending
        df = df.sort_values(["abs_corr", "gene_id"], ascending=[False, True]).reset_index(drop=True)
        top = df.head(top_n).copy()
        top.insert(0, "rank", np.arange(1, len(top) + 1))
        top_tables[ct] = top
        print(f"  {ct}: top {len(top)} genes selected "
              f"(|corr| range {top['abs_corr'].min():.3f} - {top['abs_corr'].max():.3f})")

    union_ids = sorted(set().union(*[set(t["gene_id"]) for t in top_tables.values()]) | {target_id})
    print(f"\n  union panel: {len(union_ids)} genes "
          f"({len(union_ids) - 1} from top-{top_n} lists + ANTXR2)")

    overlaps = {}
    for a, b in itertools.combinations(top_tables, 2):
        n = len(set(top_tables[a]["gene_id"]) & set(top_tables[b]["gene_id"]))
        overlaps[f"{a} vs {b}"] = n
        print(f"    overlap {a} vs {b}: {n} genes")

    results["top_genes"] = top_tables
    results["union_gene_ids"] = union_ids
    results["top_gene_overlap"] = overlaps
    return results


# ---------------------------------------------------------------------------
# Step 5: pairwise correlations among the union panel, per donor x cell_type group
# ---------------------------------------------------------------------------

def compute_pairwise_correlations(raw_adata, union_ids):
    sub = raw_adata[:, union_ids].copy()
    sub.X = sp.csr_matrix(sub.X)
    assert list(sub.var.index) == union_ids

    memento.setup_memento(
        sub, q_column=Q_COLUMN, filter_mean_thresh=COEXPR_FILTER_MEAN_THRESH,
        min_cell_count=COEXPR_MIN_GROUP_CELLS,
    )
    memento.create_groups(sub, label_columns=[COEXPR_DONOR_COL, "cell_type"])
    groups = sub.uns["memento"]["groups"]
    group_meta = memento.get_groups(sub)
    group_n_cells = {g: sub.uns["memento"]["group_cells"][g].shape[0] for g in groups}
    print(f"  {len(groups)} donor x cell_type groups (should match step 2's groups exactly)")

    # filter_genes=False deliberately: none of the curated union genes should be
    # dropped even if a specific group's own mean filter would otherwise exclude
    # them -- a gene failing in one group still gets a NaN-masked entry there via
    # memento's own var<=0 -> NaN handling in _corr_from_cov, not silent removal.
    memento.compute_1d_moments(sub, min_perc_group=COEXPR_MIN_PERC_GROUP, filter_genes=False)
    assert list(sub.var.index) == union_ids, "filter_genes=False should not reorder/subset genes"

    by_donor = {ct: {} for ct in COEXPR_CELL_TYPES}
    for g in groups:
        donor = group_meta.loc[g, COEXPR_DONOR_COL]
        ct = group_meta.loc[g, "cell_type"]
        mat = memento.get_corr_matrix(sub, g)
        n_nan = int(np.isnan(mat).sum())
        print(f"    {g:40s} n_cells={group_n_cells[g]:6d}  NaN entries: {n_nan}/{mat.size}")
        by_donor.setdefault(ct, {})[donor] = {"pairwise_matrix": mat}

    return {"union_gene_order": union_ids, "by_donor": by_donor}


# ---------------------------------------------------------------------------
# Step 6: average pairwise matrices across donors, per cell type
# ---------------------------------------------------------------------------

def average_pairwise_across_donors(pairwise_results):
    donor_avg = {}
    for ct in COEXPR_CELL_TYPES:
        donors = pairwise_results["by_donor"].get(ct, {})
        if not donors:
            print(f"  {ct}: 0 surviving donor groups -- skipping")
            continue
        stacked = np.stack([d["pairwise_matrix"] for d in donors.values()], axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN slices are expected/reported below
            mean_mat = np.nanmean(stacked, axis=0)
        n_donors_mat = np.sum(~np.isnan(stacked), axis=0)
        all_nan_frac = float(np.isnan(mean_mat).mean())
        print(f"  {ct}: {len(donors)} donor groups averaged; "
              f"{all_nan_frac*100:.1f}% of matrix entries are NaN in every donor group")
        donor_avg[ct] = {
            "pairwise_matrix_mean": mean_mat,
            "pairwise_matrix_n_donors": n_donors_mat,
            "n_donor_groups_total": len(donors),
        }
    return donor_avg


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def cluster_order_from_corr(mat):
    m = mat.copy()
    np.fill_diagonal(m, 1.0)
    m = np.nan_to_num(m, nan=0.0)  # missing pairs treated as uncorrelated for clustering only
    dist = 1 - m
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    dist[dist < 0] = 0
    condensed = squareform(dist, checks=False)
    if mat.shape[0] < 3:
        return list(range(mat.shape[0]))
    Z = linkage(condensed, method="average")
    return dendrogram(Z, no_plot=True)["leaves"]


def plot_heatmaps(donor_avg_2d, union_ids, sym_map, target_id, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    labels_full = [sym_map.get(g, g) for g in union_ids]
    for ct, d in donor_avg_2d.items():
        mat = d["pairwise_matrix_mean"]
        order = cluster_order_from_corr(mat)
        mat_o = mat[np.ix_(order, order)]
        labels_o = [labels_full[i] for i in order]

        n = len(labels_o)
        fig_size = max(6, n * 0.14)
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        im = ax.imshow(mat_o, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        fontsize = 6 if n > 80 else 7
        ax.set_xticklabels(labels_o, rotation=90, fontsize=fontsize)
        ax.set_yticklabels(labels_o, fontsize=fontsize)
        for tick_x, tick_y, lab in zip(ax.get_xticklabels(), ax.get_yticklabels(), labels_o):
            if lab == sym_map.get(target_id, target_id):
                tick_x.set_fontweight("bold")
                tick_y.set_fontweight("bold")
        ax.set_title(
            f"{ct}: donor-averaged pairwise correlation, union top-gene panel "
            f"(n={d['n_donor_groups_total']} donor groups)",
            fontsize=10,
        )
        cbar = fig.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label("memento point-estimate correlation", fontsize=8)
        fig.tight_layout()
        path = os.path.join(out_dir, f"{ct}_pairwise_correlation.png")
        fig.savefig(path, dpi=180)
        plt.close(fig)
        print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def write_summary(report, results, donor_avg_2d, target_id, sym_map, timings, versions):
    n_donors_total = report["crosstab"].shape[0]
    lines = []
    lines.append("\n## Phase 2: co-expression analysis (test run)\n")
    lines.append(
        f"Ran `scripts/coexpression_pipeline.py` on the `level3` trio working set "
        f"({INPUT_H5AD}) -- the test run planned in `planning_summaries/coexpression.md`. "
        f"Point estimates only (memento's `compute_1d_moments`/`compute_2d_moments`/"
        f"`get_corr_matrix`), no bootstrap/hypothesis testing. Output: "
        f"`{COEXPR_OUTPUT_H5AD}` (a NEW file -- the source h5ad is untouched), figures "
        f"in `{COEXPR_FIGURES_DIR}/`.\n"
    )

    lines.append("### Parameters\n")
    lines.append(f"- Donor column: `obs['{COEXPR_DONOR_COL}']` (not `donor_id` -- see donor-identity note in "
                  f"the Phase 2 prep section above)\n")
    lines.append(f"- Capture rate: `obs['{Q_COLUMN}']`, passed via memento's `q_column` "
                  f"(NOTE: only {report['capture_rate_n_distinct']} distinct values across all cells -- "
                  f"a per-assay constant, not continuous per-cell; expected, see prep notes)\n")
    lines.append(f"- Minimum group size: {COEXPR_MIN_GROUP_CELLS} cells per donor x cell_type group\n")
    lines.append(f"- Gene filter (memento defaults, recorded explicitly): "
                  f"`filter_mean_thresh={COEXPR_FILTER_MEAN_THRESH}`, `min_perc_group={COEXPR_MIN_PERC_GROUP}`\n")
    lines.append(f"- Ranking: top {COEXPR_TOP_N_GENES} genes per cell type by |donor-averaged correlation|\n")
    lines.append(f"- Averaging across donors: plain unweighted mean over non-NaN donor values "
                  f"(same donor-equal-weighting convention as Phase 1's binding aggregation)\n")

    lines.append("\n### Donor x cell_type groups\n")
    lines.append(f"{n_donors_total} donors total. {len(report['dropped_groups'])} of "
                  f"{len(report['dropped_groups']) + len(report['kept_groups'])} existing groups fell below "
                  f"the {COEXPR_MIN_GROUP_CELLS}-cell floor and were dropped:\n\n")
    lines.append("```\n" + report["dropped_groups"].to_string(index=False) + "\n```\n")
    kept_by_ct = report["kept_groups"].groupby("cell_type").size()
    lines.append("\nGroups kept, by cell type:\n\n```\n" + kept_by_ct.to_string() + "\n```\n")
    lines.append(
        "\n**Macrophage is severely thinned by the 100-cell floor** -- only "
        f"{kept_by_ct.get('macrophage', 0)} of {report['crosstab'].shape[0]} donors survive for "
        f"macrophage (vs. {kept_by_ct.get('fibroblast', 0)} for fibroblast, "
        f"{kept_by_ct.get('enterocyte', 0)} for enterocyte). Macrophage's donor-averaged "
        "correlations rest on a much smaller donor sample than the other two cell types -- "
        "treat any macrophage-specific finding here as exploratory, not confirmatory.\n"
    )

    lines.append("\n### Gene filtering\n")
    lines.append(f"Global filtered gene list: {results['gene_list_size']} / "
                  f"{results['gene_list_total']} genes pass "
                  f"(`filter_mean_thresh > {COEXPR_FILTER_MEAN_THRESH}` in "
                  f">{COEXPR_MIN_PERC_GROUP*100:.0f}% of donor x cell_type groups). "
                  f"ANTXR2 survives the global filter (required for step 2 to run at all); "
                  "per-group survival against ANTXR2's OWN group's mean-expression threshold "
                  "is logged at run time. **Caveat, not a NaN case**: a group where ANTXR2 "
                  "fails its own local mean filter still produces a numeric correlation there "
                  "(memento only emits NaN when variance is exactly zero, which is stricter "
                  "than failing the mean filter -- verified directly, e.g. "
                  "sg^F3^enterocyte has ANTXR2 mean=4.6e-6, var=8.3e-11, both nonzero). That "
                  "correlation is real output, not dropped, but is derived from near-"
                  "undetected expression and is noisier than groups where ANTXR2 clears its "
                  "own filter -- most relevant to the macrophage top-gene list below, whose "
                  "correlation magnitudes (up to 1.000) likely reflect this combined with the "
                  "small donor count.\n")

    lines.append("\n### Top ANTXR2-correlated genes per cell type (donor-averaged)\n")
    for ct, top in results["top_genes"].items():
        lines.append(f"\n**{ct}** (top 10 of {len(top)} shown, by |mean correlation|):\n\n")
        show = top.head(10)[["rank", "gene_symbol", "mean_corr", "n_donors"]]
        lines.append("```\n" + show.to_string(index=False) + "\n```\n")
    lines.append(f"\nUnion panel: {len(results['union_gene_ids'])} genes "
                  f"({len(results['union_gene_ids']) - 1} unique top-{COEXPR_TOP_N_GENES} genes + ANTXR2).\n")
    lines.append("\nPairwise overlap between cell types' top-gene lists:\n\n```\n"
                  + "\n".join(f"{k}: {v}" for k, v in results["top_gene_overlap"].items()) + "\n```\n")

    lines.append("\n### Package versions\n")
    lines.append("```\n" + "\n".join(f"{k}: {v}" for k, v in versions.items()) + "\n```\n")

    lines.append("\n### Runtime (test-run trio: 56,882 cells x 18,370 genes)\n")
    lines.append("```\n" + "\n".join(f"{k}: {v:.1f}s" for k, v in timings.items())
                  + f"\ntotal: {sum(timings.values()):.1f}s\n```\n")
    lines.append(
        "\nStep 2 (ANTXR2 vs. all filtered genes) and step 5 (pairwise over the union "
        "panel) are the two that will scale with dataset size when this moves beyond the "
        "trio test run -- step 2 scales roughly with n_filtered_genes x n_groups, step 5 "
        f"with n_union_genes^2 x n_groups (here n_union_genes is fixed at "
        f"{len(results['union_gene_ids'])} regardless of dataset size, so step 5's cost is "
        "mostly driven by n_groups, i.e. how many donor x cell_type combinations exist).\n"
    )

    with open(ANALYSIS_SUMMARY_PATH, "a") as f:
        f.write("".join(lines))
    print(f"\n  appended Phase 2 summary to {ANALYSIS_SUMMARY_PATH}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print(f"=== ANTXR2 co-expression pipeline ===")
    print(f"input: {INPUT_H5AD}")
    print(f"output: {COEXPR_OUTPUT_H5AD}")

    versions = {
        "anndata": pkg_version("anndata"),
        "scanpy": pkg_version("scanpy"),
        "memento-de": pkg_version("memento-de"),
    }
    print(f"package versions: {versions}")

    with step_timer("1_load_and_filter"):
        raw_adata, report, target_id, _sym_to_id = load_and_filter(INPUT_H5AD)
        sym_map = dict(zip(raw_adata.var.index, raw_adata.var[GENE_NAME_COL].astype(str)))

    with step_timer("2_compute_antxr2_correlations"):
        results = compute_antxr2_correlations(raw_adata, target_id)

    with step_timer("3_average_across_donors_1d"):
        results = average_across_donors_1d(results)

    with step_timer("4_select_top_genes"):
        results = select_top_genes(results, sym_map, target_id)

    with step_timer("5_compute_pairwise_correlations"):
        pairwise_results = compute_pairwise_correlations(raw_adata, results["union_gene_ids"])

    with step_timer("6_average_pairwise_across_donors"):
        donor_avg_2d = average_pairwise_across_donors(pairwise_results)

    with step_timer("7_write_output_and_figures"):
        # Build the uns['memento_correlations'] structure per the spec's schema:
        #   by_donor[cell_type][donor] -> step 3 antxr2_vs_all + step 5 pairwise_matrix
        #   donor_averaged[cell_type]  -> step 3b averaged antxr2_vs_all (+n_donors)
        #                                  + step 5b averaged pairwise_matrix (+n_donors)
        by_donor = {}
        for ct in COEXPR_CELL_TYPES:
            by_donor[ct] = {}
            d1 = results["by_donor"].get(ct, {})
            d2 = pairwise_results["by_donor"].get(ct, {})
            for donor in set(d1) | set(d2):
                by_donor[ct][donor] = {
                    **({"n_cells": d1[donor]["n_cells"],
                        "antxr2_vs_all": d1[donor]["antxr2_vs_all"]} if donor in d1 else {}),
                    **({"pairwise_matrix": d2[donor]["pairwise_matrix"]} if donor in d2 else {}),
                }

        donor_averaged = {}
        for ct in COEXPR_CELL_TYPES:
            donor_averaged[ct] = {}
            if ct in results["donor_averaged_1d"]:
                donor_averaged[ct]["antxr2_vs_all_mean"] = results["donor_averaged_1d"][ct]["mean_corr"]
                donor_averaged[ct]["antxr2_vs_all_n_donors"] = results["donor_averaged_1d"][ct]["n_donors"]
                donor_averaged[ct]["top_genes"] = results["top_genes"][ct]
            if ct in donor_avg_2d:
                donor_averaged[ct]["pairwise_matrix_mean"] = donor_avg_2d[ct]["pairwise_matrix_mean"]
                donor_averaged[ct]["pairwise_matrix_n_donors"] = donor_avg_2d[ct]["pairwise_matrix_n_donors"]

        raw_adata.uns["memento_correlations"] = {
            "params": {
                "donor_col": COEXPR_DONOR_COL,
                "q_column": Q_COLUMN,
                "target_gene": COEXPR_TARGET_GENE,
                "target_gene_id": target_id,
                "min_group_cells": COEXPR_MIN_GROUP_CELLS,
                "filter_mean_thresh": COEXPR_FILTER_MEAN_THRESH,
                "min_perc_group": COEXPR_MIN_PERC_GROUP,
                "top_n_genes": COEXPR_TOP_N_GENES,
                "estimator_type": "hyper_relative",
                "point_estimates_only": True,
                "package_versions": versions,
                "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
            "gene_symbol": sym_map,
            "antxr2_gene_order": results["antxr2_gene_order"],
            "union_gene_order": results["union_gene_ids"],
            "by_donor": by_donor,
            "donor_averaged": donor_averaged,
            "dropped_groups": report["dropped_groups"],
            "kept_groups": report["kept_groups"],
            "top_gene_overlap": results["top_gene_overlap"],
        }

        os.makedirs(os.path.dirname(COEXPR_OUTPUT_H5AD), exist_ok=True)
        print(f"  writing {COEXPR_OUTPUT_H5AD} ...")
        raw_adata.write_h5ad(COEXPR_OUTPUT_H5AD, compression="gzip", compression_opts=4)
        size = os.path.getsize(COEXPR_OUTPUT_H5AD)
        print(f"  wrote {size/1e9:.2f} GB")

        plot_heatmaps(donor_avg_2d, results["union_gene_ids"], sym_map, target_id, COEXPR_FIGURES_DIR)

        for ct, top in results["top_genes"].items():
            csv_path = os.path.join(COEXPR_FIGURES_DIR, f"{ct}_top_genes.csv")
            top.to_csv(csv_path, index=False)
            print(f"  wrote {csv_path}")

    write_summary(report, results, donor_avg_2d, target_id, sym_map, TIMINGS, versions)

    print(f"\n=== done. total runtime {sum(TIMINGS.values()):.1f}s ===")


if __name__ == "__main__":
    main()
