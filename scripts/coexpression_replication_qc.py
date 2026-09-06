"""Donor-replication sanity check for the ANTXR2 co-expression correlations
(Phase 2 QC), per user request in this session.

Wording precision (per user instruction): each donor x cell_type group in
`coexpression_pipeline.py`'s output already holds a memento point-estimate of
the SINGLE-CELL correlation between ANTXR2 and every other filtered gene, one
vector per donor. What this script computes is the correlation *between
donors' single-cell correlation estimates* -- i.e. do two independently
estimated donors' ~3,746-long "ANTXR2 vs. gene" correlation vectors agree with
each other. That second-order agreement is what would make the project's
donor-averaged tables (top-gene lists, GSEA) a meaningful summary rather than
an average of independent noise.

Scope decision from this session's investigation (see ANALYSIS_SUMMARY.md):
chemistry (10x 5' v2 vs 3' v2) is completely confounded with donor identity in
this working set -- no donor was profiled with both chemistries, and
macrophage is 100% 5' v2. A true same-donor technical replicate does not
exist here. Per user decision, this script reports BIOLOGICAL (donor-vs-donor)
replication only; chemistry composition is recorded as a caveat column, not a
formal split. Macrophage is included but flagged exploratory throughout (5
donors, single chemistry -- same caveat as everywhere else in this project).

No new memento computation -- pure analysis of the already-verified per-donor
correlation vectors in
`.uns['memento_correlations']['by_donor'][cell_type][donor]['antxr2_vs_all']`
written by coexpression_pipeline.py.

Run in the `antxr2` conda env: `conda run -n antxr2 python
scripts/coexpression_replication_qc.py`.
"""
import os

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from config import COEXPR_CELL_TYPES, COEXPR_DONOR_COL, COEXPR_FIGURES_DIR, COEXPR_OUTPUT_H5AD
from coexpression_pipeline import CELLTYPE_COLOR, cluster_order_from_corr

OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "replication_qc")
N_PERMUTATIONS = 200
TOP_N_CLEAN = 8
RNG_SEED = 0
EXPLORATORY_CELL_TYPES = {"macrophage"}

CSV_CAVEAT_HEADER = """\
# Donor-replication QC for ANTXR2's memento single-cell correlation estimates.
# Each value is the Pearson correlation BETWEEN two donors' independently
# estimated ~3,746-gene "ANTXR2 vs. gene" single-cell correlation vectors --
# i.e. a correlation of correlation estimates, not a single-cell correlation
# itself.
# CAVEAT 1: chemistry (10x 5' v2 vs 3' v2) is fully confounded with donor
# identity in this dataset -- no donor was profiled with both chemistries, and
# macrophage is 100% 5' v2. This file measures BIOLOGICAL (donor-vs-donor)
# replication only; it cannot and does not test a technical-replicate effect.
# CAVEAT 2: macrophage is exploratory-only (5 donors, single chemistry, thin
# power) -- same caveat as everywhere else this project reports macrophage.
# CAVEAT 3: n_genes_overlap is the number of genes with a non-NaN estimate in
# BOTH donors of a pair (pairwise-complete Pearson); it varies per pair and is
# recorded per row rather than assumed constant.
# CAVEAT 4: donors whose antxr2_vs_all vector was 100% NaN (ANTXR2 itself has
# zero variance in that donor's group -- too few cells) are dropped entirely
# before this table is built, since they cannot contribute a comparison; see
# replication_summary_stats.csv's n_donors_dropped_all_nan for counts per
# cell type and the script's console log for which donors.
"""


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_by_donor(h5ad_path):
    adata = anndata.read_h5ad(h5ad_path, backed="r")
    mc = adata.uns["memento_correlations"]
    obs = adata.obs
    return mc["by_donor"], mc["antxr2_gene_order"], obs


def build_gene_by_donor_matrix(by_donor_ct, antxr2_gene_order):
    donors = sorted(d for d, v in by_donor_ct.items() if "antxr2_vs_all" in v)
    df = pd.DataFrame({d: by_donor_ct[d]["antxr2_vs_all"] for d in donors}, index=antxr2_gene_order)
    n_cells = pd.Series({d: by_donor_ct[d]["n_cells"] for d in donors})
    return df, n_cells


def drop_uninformative_donors(df, n_cells, ct):
    """A donor whose antxr2_vs_all vector is entirely NaN cannot contribute any
    pairwise comparison -- this happens when ANTXR2 itself has var<=0 in that
    donor's group, which nulls EVERY entry (ANTXR2 is one of the two genes in
    every pair) via the estimator-bug fix in coexpression_pipeline.py. Left in,
    these donors inflate n_donors and show up as all-blank heatmap rows/pairplot
    panels; dropped here and reported explicitly rather than silently."""
    all_nan = df.columns[df.isna().all(axis=0)]
    if len(all_nan):
        dropped_n_cells = n_cells.reindex(all_nan)
        print(f"  DROPPED {len(all_nan)} donor(s) with 100% NaN antxr2_vs_all "
              f"(ANTXR2 has zero variance in their group -- too few cells for any "
              f"gene-level signal): {dict(dropped_n_cells)}")
        df = df.drop(columns=all_nan)
        n_cells = n_cells.drop(index=all_nan)
    return df, n_cells, list(all_nan)


# ---------------------------------------------------------------------------
# Pairwise correlation of per-donor correlation vectors
# ---------------------------------------------------------------------------

def pairwise_pearson(x, y):
    mask = ~np.isnan(x) & ~np.isnan(y)
    n = int(mask.sum())
    if n < 2:
        return np.nan, n
    r = float(np.corrcoef(x[mask], y[mask])[0, 1])
    return r, n


def compute_pairwise(df):
    """Exact-match source of truth: both the CSV/heatmap/summary stats AND the
    pairplot panel annotations call this, so they can never silently disagree."""
    donors = df.columns.tolist()
    n = len(donors)
    corr_mat = pd.DataFrame(np.eye(n), index=donors, columns=donors)
    overlap_mat = pd.DataFrame(np.nan, index=donors, columns=donors)
    rows = []
    for i, a in enumerate(donors):
        overlap_mat.loc[a, a] = int(df[a].notna().sum())
        for b in donors[i + 1:]:
            r, n_overlap = pairwise_pearson(df[a].values, df[b].values)
            corr_mat.loc[a, b] = corr_mat.loc[b, a] = r
            overlap_mat.loc[a, b] = overlap_mat.loc[b, a] = n_overlap
            rows.append({"donor_a": a, "donor_b": b, "n_genes_overlap": n_overlap, "pearson_r": r})
    return corr_mat, overlap_mat, pd.DataFrame(rows)


def permutation_null(df, n_perm=N_PERMUTATIONS, seed=RNG_SEED):
    """Independently shuffle each donor's gene axis (preserves that donor's own
    value distribution and NaN count exactly, destroys true gene-to-gene
    correspondence across donors), recompute the full pairwise-correlation
    matrix, pool off-diagonal values across permutations."""
    # NB: shuffling each donor's column independently moves each donor's OWN
    # NaN positions around at random, so two donors' non-NaN overlap for a
    # given pair becomes a random draw each permutation (unlike the real data,
    # where NaN positions are gene-identity-aligned and overlap is fixed). A
    # small fraction of permuted pairs land on very low overlap by chance and
    # legitimately produce NaN (matching pairwise_pearson's own n<2 -> NaN
    # rule) -- filtered out below with nan-aware reduction, not an error.
    rng = np.random.default_rng(seed)
    arr = df.values
    n_genes, n_donors = arr.shape
    iu = np.triu_indices(n_donors, k=1)
    null_vals = np.empty(n_perm * len(iu[0]))
    for p in range(n_perm):
        shuffled = np.empty_like(arr)
        for j in range(n_donors):
            shuffled[:, j] = rng.permutation(arr[:, j])
        shuffled_df = pd.DataFrame(shuffled, columns=df.columns)
        r_mat = shuffled_df.corr(method="pearson").values
        null_vals[p * len(iu[0]):(p + 1) * len(iu[0])] = r_mat[iu]
    n_nan = int(np.isnan(null_vals).sum())
    if n_nan:
        print(f"    (permutation null: {n_nan}/{null_vals.size} pooled values were NaN "
              f"from low-overlap draws -- dropped before summarizing)")
    return null_vals[~np.isnan(null_vals)]


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------

def plot_heatmap(corr_mat, n_cells, ct, out_dir, exploratory=False):
    order_idx = cluster_order_from_corr(corr_mat.values)
    donors = corr_mat.index.to_numpy()[order_idx]
    mat_o = corr_mat.values[np.ix_(order_idx, order_idx)]
    labels = [f"{d} (n={n_cells[d]:,})" for d in donors]

    n = len(donors)
    fig_w = max(6.5, n * 0.35)
    fig, ax = plt.subplots(figsize=(fig_w, fig_w))
    im = ax.imshow(mat_o, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    fontsize = 6 if n > 20 else 8
    ax.set_xticklabels(labels, rotation=90, fontsize=fontsize)
    ax.set_yticklabels(labels, fontsize=fontsize)
    title = f"{ct}: donor-pairwise correlation of ANTXR2-vs-gene estimates (n={n} donors)"
    if exploratory:
        title += "\nEXPLORATORY: single chemistry, thin donor count -- see caveats"
    fig.suptitle(title, fontsize=10, wrap=True)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label("Pearson r (correlation of correlation estimates)", fontsize=8)
    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    path = os.path.join(out_dir, f"{ct}_donor_pairwise_heatmap.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# Combined real-vs-null distribution figure
# ---------------------------------------------------------------------------

def plot_distributions(pairs_by_ct, null_by_ct, out_dir):
    cts = [ct for ct in COEXPR_CELL_TYPES if ct in pairs_by_ct]
    fig, axes = plt.subplots(1, len(cts), figsize=(5.5 * len(cts), 4.5), squeeze=False)
    axes = axes[0]
    for ax, ct in zip(axes, cts):
        real = pairs_by_ct[ct]["pearson_r"].dropna().values
        null = null_by_ct[ct]
        color = CELLTYPE_COLOR.get(ct, "#555555")
        bins = np.linspace(-1, 1, 41)
        ax.hist(null, bins=bins, density=True, color="#999999", alpha=0.45, label="permutation null")
        ax.hist(real, bins=bins, density=True, color=color, alpha=0.6, label="observed donor pairs")
        ax.axvline(np.median(real), color=color, linestyle="-", linewidth=2)
        ax.axvline(np.median(null), color="#444444", linestyle="--", linewidth=1.5)
        ax.axvline(0, color="black", linewidth=0.8)
        title = f"{ct}\nn_pairs={len(real)}, median r={np.median(real):.3f} (null median={np.median(null):.3f})"
        if ct in EXPLORATORY_CELL_TYPES:
            title += "\nEXPLORATORY (n=5 donors, single chemistry)"
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Pearson r between two donors' ANTXR2-vs-gene vectors")
        ax.set_xlim(-1, 1)
        if ax is axes[0]:
            ax.set_ylabel("density")
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Biological-replicate check: do independent donors' single-cell ANTXR2 "
                 "correlation estimates agree? (observed vs. gene-shuffled null)", fontsize=11)
    fig.tight_layout()
    path = os.path.join(out_dir, "donor_replication_distributions.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# Pairplot (R pairs()/ggpairs-style grid)
# ---------------------------------------------------------------------------

def plot_pairgrid(df, donor_order, out_path, title, mode, color):
    n = len(donor_order)
    panel_size = 0.6 if mode == "hex" else 1.15
    fig, axes = plt.subplots(n, n, figsize=(max(4, n * panel_size), max(4, n * panel_size)))
    if n == 1:
        axes = np.array([[axes]])
    for i, di in enumerate(donor_order):
        for j, dj in enumerate(donor_order):
            ax = axes[i, j]
            if j > i:
                ax.axis("off")
                continue
            if i == j:
                vals = df[di].dropna().values
                ax.hist(vals, bins=25, color="#777777", alpha=0.8)
                ax.set_yticks([])
                ax.set_xlim(-1, 1)
            else:
                x = df[dj].values
                y = df[di].values
                mask = ~np.isnan(x) & ~np.isnan(y)
                r, _ = pairwise_pearson(x, y)
                if mode == "hex":
                    ax.hexbin(x[mask], y[mask], gridsize=22, cmap="Blues", extent=(-1, 1, -1, 1), mincnt=1)
                else:
                    ax.scatter(x[mask], y[mask], s=4, alpha=0.2, color=color, linewidths=0)
                ax.plot([-1, 1], [-1, 1], color="black", linewidth=0.6, linestyle="--", alpha=0.6)
                ax.set_xlim(-1, 1)
                ax.set_ylim(-1, 1)
                ax.text(0.04, 0.86, f"r={r:.2f}", transform=ax.transAxes, fontsize=6)
            fontsize = 5 if n > 15 else 7
            if i == n - 1:
                ax.set_xlabel(dj, fontsize=fontsize, rotation=90)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(di, fontsize=fontsize, rotation=0, ha="right", va="center")
            else:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=4)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    dpi = 130 if mode == "hex" else 150
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"  wrote {out_path}")


def cluster_order_labels(corr_mat):
    order_idx = cluster_order_from_corr(corr_mat.values)
    return corr_mat.index.to_numpy()[order_idx].tolist()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"loading {COEXPR_OUTPUT_H5AD} ...")
    by_donor, antxr2_gene_order, obs = load_by_donor(COEXPR_OUTPUT_H5AD)

    donor_assay = obs[[COEXPR_DONOR_COL, "assay"]].drop_duplicates().set_index(COEXPR_DONOR_COL)["assay"]

    summary_rows = []
    all_pairs = []
    pairs_by_ct, null_by_ct = {}, {}
    per_ct_matrix = {}

    for ct in COEXPR_CELL_TYPES:
        if ct not in by_donor or not by_donor[ct]:
            print(f"\n{ct}: no data -- skipping")
            continue
        print(f"\n=== {ct} ===")
        df, n_cells = build_gene_by_donor_matrix(by_donor[ct], antxr2_gene_order)
        df, n_cells, dropped_donors = drop_uninformative_donors(df, n_cells, ct)
        n_donors = df.shape[1]
        print(f"  {n_donors} usable donors ({len(dropped_donors)} dropped), {df.shape[0]} genes")
        per_ct_matrix[ct] = (df, n_cells)

        corr_mat, overlap_mat, pairs_df = compute_pairwise(df)
        pairs_df.insert(0, "cell_type", ct)
        pairs_by_ct[ct] = pairs_df
        all_pairs.append(pairs_df)

        null_vals = permutation_null(df)
        null_by_ct[ct] = null_vals

        exploratory = ct in EXPLORATORY_CELL_TYPES
        plot_heatmap(corr_mat, n_cells, ct, OUT_DIR, exploratory=exploratory)

        chem_counts = donor_assay.reindex(df.columns).value_counts().to_dict()
        real_vals = pairs_df["pearson_r"].dropna().values
        null_p975 = float(np.percentile(null_vals, 97.5))
        frac_above_null = float((real_vals > null_p975).mean()) if len(real_vals) else np.nan
        # One-sided: are real donor-pair correlations stochastically greater than
        # the permutation null, not just "is the median outside an eyeballed band"?
        mwu_p = float(mannwhitneyu(real_vals, null_vals, alternative="greater").pvalue) if len(real_vals) else np.nan
        summary_rows.append({
            "cell_type": ct,
            "n_donors": n_donors,
            "n_donors_dropped_all_nan": len(dropped_donors),
            "n_pairs": len(pairs_df),
            "median_r_real": float(np.median(real_vals)) if len(real_vals) else np.nan,
            "mean_r_real": float(np.mean(real_vals)) if len(real_vals) else np.nan,
            "median_r_null": float(np.median(null_vals)),
            "null_p2.5": float(np.percentile(null_vals, 2.5)),
            "null_p97.5": null_p975,
            "frac_real_pairs_above_null_p97.5": frac_above_null,
            "mannwhitney_p_real_gt_null": mwu_p,
            "exploratory_flag": exploratory,
            "donor_chemistry_composition": "; ".join(f"{k}:{v}" for k, v in sorted(chem_counts.items())),
        })
        print(f"  real median r={np.median(real_vals):.3f}, mean r={np.mean(real_vals):.3f} "
              f"vs. null median r={np.median(null_vals):.3f} "
              f"(null 95% band [{np.percentile(null_vals, 2.5):.3f}, {null_p975:.3f}])")
        print(f"  {frac_above_null*100:.1f}% of real donor pairs exceed the null's 97.5th "
              f"percentile; Mann-Whitney U (real > null) p={mwu_p:.2e}")
        print(f"  donor chemistry composition: {chem_counts}")

        # pairplots
        full_order = cluster_order_labels(corr_mat)
        color = CELLTYPE_COLOR.get(ct, "#555555")
        plot_pairgrid(
            df, full_order,
            os.path.join(OUT_DIR, f"{ct}_donor_pairplot_full.png"),
            f"{ct}: all {n_donors} donors, ANTXR2-vs-gene correlation estimates "
            f"(hexbin density, lower triangle)" + ("  [EXPLORATORY]" if exploratory else ""),
            mode="hex", color=color,
        )

        top_donors = n_cells.sort_values(ascending=False).head(TOP_N_CLEAN).index.tolist()
        top_corr, _, _ = compute_pairwise(df[top_donors])
        top_order = cluster_order_labels(top_corr)
        plot_pairgrid(
            df, top_order,
            os.path.join(OUT_DIR, f"{ct}_donor_pairplot_top8.png"),
            f"{ct}: top {len(top_donors)} donors by n_cells (scatter, lower triangle)"
            + ("  [EXPLORATORY]" if exploratory else ""),
            mode="scatter", color=color,
        )

    plot_distributions(pairs_by_ct, null_by_ct, OUT_DIR)

    pairs_all = pd.concat(all_pairs, ignore_index=True)
    csv_path = os.path.join(OUT_DIR, "donor_pairwise_correlations.csv")
    with open(csv_path, "w") as f:
        f.write(CSV_CAVEAT_HEADER)
        pairs_all.to_csv(f, index=False)
    print(f"\nwrote {csv_path}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUT_DIR, "replication_summary_stats.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"wrote {summary_path}")

    print("\n=== summary ===")
    print(summary_df.to_string(index=False))
    if "macrophage" in summary_df["cell_type"].values:
        print("\nNOTE: macrophage numbers above are exploratory only -- 5 donors, "
              "100% single chemistry (10x 5' v2), thin power. Do not treat a weak "
              "macrophage result as invalidating fibroblast/enterocyte, and do not "
              "treat a strong one as confirmatory on its own.")


if __name__ == "__main__":
    main()
