"""Discovery/replication check for ANTXR2 co-expression correlations, using a
proper one-sample significance test rather than a magnitude-based gene
selection hack -- per user direction this session.

Background (see ANALYSIS_SUMMARY.md "Donor-replication sanity check" section
for the full arc): two earlier attempts at a gene-selection strategy for this
QC were tried and found wanting -- (1) selecting top genes by the CROSS-DONOR
MEAN correlation is circular, since the same donors being tested for
agreement were also used to pick which genes to test; (2) selecting each
donor's own top-N genes and taking the union removes circularity but is
noise-starved -- a hard top-N cut on one noisy donor's point estimate is
dominated by that donor's own sampling noise, understating real signal
(overlap between donors' own top-50 lists was only ~11-15%).

The fix: split donors into two independent halves. In the DISCOVERY half, run
memento's own one-sample bootstrap hypothesis test -- `ht_2d_moments` with a
constant (all-ones) treatment column, which collapses its usual two-group
regression machinery into exactly a one-sample test of "is this correlation
significantly different from zero", using memento's own per-cell bootstrap
for the standard error rather than a hard magnitude cutoff (see
`memento/hypothesis_test.py::_regress_2d`: `if (treatment==1).mean()==1:
corr_coef = np.average(boot_corr, ..., weights=Nc_list)`). Genes passing an
FDR threshold there are then evaluated for replication ONLY on the untouched
REPLICATION half: sign concordance (with a binomial test), an effect-size
replication scatter (discovery bootstrap coefficient vs. replication-half
mean correlation), and the cross-half donor-pairwise correlation-of-
correlations restricted to these genes. No donor is ever used for both gene
selection and its own replication test.

Point estimates for the replication half are read from the already-computed,
already-verified `.uns['memento_correlations']['by_donor']` in
elmentaite2021_trio_level3_coexpr.h5ad (the same source
coexpression_replication_qc.py uses) -- only the discovery half needs fresh
memento computation here, via the INPUT (pre-correlation) working set h5ad.

Run in the `antxr2` conda env: `conda run -n antxr2 python
scripts/coexpression_discovery_replication.py [--cell-types fibroblast,enterocyte,macrophage] [--num-boot 5000]`.
"""
import argparse
import os
import warnings

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import memento
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from config import (
    COEXPR_FIGURES_DIR, COEXPR_MIN_GROUP_CELLS, COEXPR_OUTPUT_H5AD,
    COEXPR_TARGET_GENE, COEXPR_VARIANTS,
)
from coexpression_pipeline import CELLTYPE_COLOR
from coexpression_replication_qc import (
    build_gene_by_donor_matrix, drop_uninformative_donors, load_by_donor, pairwise_pearson,
)

INPUT_H5AD = COEXPR_VARIANTS["level3"]["h5ad"]
GENE_NAME_COL = "feature_name"
Q_COLUMN = "capture_rate_pbmc"
DONOR_COL = "donor"
LEVEL3_LABEL = {"fibroblast": "Crypt_fibroblast_PI16", "enterocyte": "Enterocyte", "macrophage": "Macrophage"}
FILTER_MEAN_THRESH = 0.07
DEFAULT_MIN_PERC_GROUP = 0.7
SHRINKAGE, TRIM_PERCENT = 0.0, 0.5  # this pipeline's adopted defaults, see ANALYSIS_SUMMARY.md
CLIP_THRESH = 0.999  # matches coexpression_replication_qc.py's clip-artifact exclusion
OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "replication_qc")

CSV_CAVEAT_HEADER = """\
# ANTXR2-vs-gene correlations found significant (FDR<{fdr}) by memento's one-sample
# bootstrap test (ht_2d_moments with a constant treatment column) on a DISCOVERY half of
# donors, evaluated for replication on a completely separate REPLICATION half.
# coef/se/pval/fdr are from the discovery half only. replication_mean_corr is the plain
# (clip-excluded, |r|>=0.999 treated as an estimator artifact) mean of the already-computed
# per-donor point estimates across replication-half donors -- no replication donor was used
# in gene selection, so this is a non-circular replication check.
# See ANALYSIS_SUMMARY.md's "Donor-replication sanity check" section for the full method and
# why this replaced two earlier, flawed selection strategies (circular mean-based top-N, and
# noise-starved per-donor-top-N union).
"""


def split_donors(sub_full, donor_col):
    n_cells_by_donor = sub_full.obs[donor_col].value_counts()
    usable = n_cells_by_donor[n_cells_by_donor >= COEXPR_MIN_GROUP_CELLS].index.tolist()
    donors_sorted = n_cells_by_donor.reindex(usable).sort_values(ascending=False).index.tolist()
    return donors_sorted[0::2], donors_sorted[1::2]  # round-robin by cell count -> balanced power


def run_discovery_ht(full, ct, disc, num_boot, num_cpus, random_state, target_symbol=None):
    """target_symbol defaults to COEXPR_TARGET_GENE (ANTXR2) -- overridable so
    the identical pipeline can be run for an arbitrary anchor gene (see
    anchor_gene_control.py, prompts/partner_availability.md Task 2)."""
    target_symbol = target_symbol or COEXPR_TARGET_GENE
    label = LEVEL3_LABEL[ct]
    adata = full[(full.obs["level_3_annot"].astype(str) == label) & (full.obs[DONOR_COL].isin(disc))].copy()
    adata.obs[DONOR_COL] = adata.obs[DONOR_COL].astype(str)
    sym_to_id = dict(zip(adata.var[GENE_NAME_COL].astype(str), adata.var.index))
    target_id = sym_to_id[target_symbol]

    adata.X = sp.csr_matrix(adata.X)
    memento.setup_memento(adata, q_column=Q_COLUMN, filter_mean_thresh=FILTER_MEAN_THRESH,
                           min_cell_count=COEXPR_MIN_GROUP_CELLS, shrinkage=SHRINKAGE, trim_percent=TRIM_PERCENT)
    memento.create_groups(adata, label_columns=[DONOR_COL])
    groups = adata.uns["memento"]["groups"]
    adata_pre_filter = adata.copy()

    # ANTXR2 is borderline-expressed in some cell types (documented elsewhere in this
    # project); if it fails the global min_perc_group presence filter with only half the
    # donors, lower the threshold and retry rather than aborting -- same posture as
    # coexpression_pipeline.py's own documented deviations.
    mpg = DEFAULT_MIN_PERC_GROUP
    while True:
        adata_f = adata_pre_filter.copy()
        memento.compute_1d_moments(adata_f, min_perc_group=mpg, filter_genes=True)
        gene_list = adata_f.uns["memento"]["gene_list"]
        present = target_id in gene_list
        print(f"  {len(groups)} discovery donor groups, {len(gene_list)} genes pass global filter "
              f"(min_perc_group={mpg:.2f}), {target_symbol} present: {present}")
        if present or mpg <= 0.15:
            break
        mpg -= 0.1
    if target_id not in gene_list:
        raise SystemExit(f"{ct}: {target_symbol} filtered out of discovery-half gene list even at min_perc_group=0.15")

    test_genes = [g for g in gene_list if g != target_id]
    gene_pairs = [(target_id, g) for g in test_genes]
    memento.compute_2d_moments(adata_f, gene_pairs)

    # constant treatment=1 for every group collapses ht_2d_moments's two-group regression
    # into a one-sample bootstrap test of "correlation != 0" -- see module docstring.
    treatment = pd.DataFrame({"one": np.ones(len(groups))})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        memento.ht_2d_moments(adata_f, treatment=treatment, num_boot=num_boot, num_cpus=num_cpus,
                               verbose=0, random_state=random_state)
    ht = adata_f.uns["memento"]["2d_ht"]
    ht_df = pd.DataFrame({
        "gene_id": test_genes, "coef": ht["corr_coef"], "se": ht["corr_se"], "pval": ht["corr_asl"],
    }).dropna(subset=["pval"])
    ht_df["fdr"] = multipletests(ht_df["pval"].values, method="fdr_bh")[1]
    ht_df = ht_df.sort_values("pval").reset_index(drop=True)
    return ht_df, mpg


def evaluate_replication(ct, disc, repl, ht_df, by_donor, gene_order, fdr_thresh):
    df, n_cells = build_gene_by_donor_matrix(by_donor[ct], gene_order)
    df, n_cells, dropped = drop_uninformative_donors(df, n_cells, ct)
    disc_use = [d for d in disc if d in df.columns]
    repl_use = [d for d in repl if d in df.columns]

    arr_repl = df[repl_use].values
    arr_repl_clean = np.where(np.abs(arr_repl) >= CLIP_THRESH, np.nan, arr_repl)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean_repl_all = pd.Series(np.nanmean(arr_repl_clean, axis=1), index=df.index)

    sig = ht_df[ht_df["fdr"] < fdr_thresh].set_index("gene_id")
    common = sig.index.intersection(df.index)
    eval_df = pd.DataFrame({
        "coef_disc": sig.loc[common, "coef"], "mean_repl": mean_repl_all.loc[common],
    }).dropna()

    sub = df.loc[common]
    rows = []
    for da in disc_use:
        for db in repl_use:
            r, n = pairwise_pearson(sub[da].values, sub[db].values)
            if not np.isnan(r):
                rows.append({"donor_a": da, "donor_b": db, "n_overlap": n, "r": r})
    cross_df = pd.DataFrame(rows)

    return dict(disc_use=disc_use, repl_use=repl_use, sig=sig, common=common,
                eval_df=eval_df, cross_df=cross_df)


TMP_PKL = lambda ct: os.path.join(OUT_DIR, f"_tmp_discovery_repl_{ct}.pkl")


def run_one(ct, args, full, by_donor, gene_order):
    """Compute discovery-half hypothesis test + replication-half evaluation for
    ONE cell type and cache the result to disk. Split out from an all-cell-types
    main() so each cell type can be run as its own process invocation -- the
    full run (3 cell types x num_boot=5000) is a multi-minute job, and this
    project's tooling was unreliable running it as a single long-lived
    background process, but reliable running it cell-type-at-a-time in the
    foreground."""
    import pickle

    print(f"\n{'='*70}\n{ct}\n{'='*70}")
    label = LEVEL3_LABEL[ct]
    sub_full = full[full.obs["level_3_annot"].astype(str) == label]
    disc, repl = split_donors(sub_full, DONOR_COL)
    print(f"discovery n={len(disc)} donors, replication n={len(repl)} donors")

    ht_df, mpg_used = run_discovery_ht(full, ct, disc, args.num_boot, args.num_cpus, args.random_state)
    n_sig = int((ht_df["fdr"] < args.fdr_thresh).sum())
    print(f"{len(ht_df)} genes with a valid one-sample test (min_perc_group={mpg_used:.2f}); "
          f"{n_sig} pass FDR<{args.fdr_thresh}")

    result = evaluate_replication(ct, disc, repl, ht_df, by_donor, gene_order, args.fdr_thresh)
    eval_df, cross_df = result["eval_df"], result["cross_df"]
    n_eval = len(eval_df)

    sym_map = dict(zip(full.var.index, full.var[GENE_NAME_COL].astype(str)))
    sig_out = ht_df[ht_df["fdr"] < args.fdr_thresh].copy()
    sig_out.insert(0, "cell_type", ct)
    sig_out.insert(2, "gene_symbol", sig_out["gene_id"].map(sym_map))
    sig_out = sig_out.merge(
        eval_df.rename(columns={"mean_repl": "replication_mean_corr"})[["replication_mean_corr"]],
        left_on="gene_id", right_index=True, how="left")
    sig_out["in_replication_eval"] = sig_out["replication_mean_corr"].notna()

    row = {
        "cell_type": ct, "n_discovery_donors": len(disc), "n_replication_donors": len(repl),
        "min_perc_group_used": mpg_used, "n_genes_tested": len(ht_df), "n_significant_fdr": n_sig,
        "n_evaluable_in_replication": n_eval,
    }
    if n_eval >= 2:
        same_sign = int((np.sign(eval_df["coef_disc"]) == np.sign(eval_df["mean_repl"])).sum())
        bt = binomtest(same_sign, n_eval, 0.5, alternative="greater")
        r_repl = float(np.corrcoef(eval_df["coef_disc"], eval_df["mean_repl"])[0, 1])
        row.update({
            "sign_concordance_frac": same_sign / n_eval, "sign_concordance_n": same_sign,
            "sign_concordance_binom_p": bt.pvalue, "effect_size_replication_r": r_repl,
            "cross_half_median_r": cross_df["r"].median() if len(cross_df) else np.nan,
            "cross_half_n_pairs": len(cross_df),
        })
        print(f"sign concordance: {same_sign}/{n_eval} ({same_sign/n_eval*100:.1f}%), "
              f"binomial p={bt.pvalue:.2e}; effect-size replication r={r_repl:.3f}; "
              f"cross-half median r={row['cross_half_median_r']:.3f} (n={len(cross_df)} pairs)")
    else:
        print(f"only {n_eval} evaluable gene(s) -- too few for a meaningful replication check")

    cache = dict(ct=ct, disc=disc, repl=repl, disc_use=result["disc_use"], repl_use=result["repl_use"],
                 eval_df=eval_df, row=row, sig_out=sig_out, fdr_thresh=args.fdr_thresh)
    with open(TMP_PKL(ct), "wb") as f:
        pickle.dump(cache, f)
    print(f"cached intermediate result to {TMP_PKL(ct)}")


def assemble(cell_types, fdr_thresh_fallback):
    """Build the combined figure + summary CSVs from each cell type's cached
    intermediate result (written by run_one). Requires run_one to have been
    run for every cell type in `cell_types` first."""
    import pickle

    caches = {}
    for ct in cell_types:
        path = TMP_PKL(ct)
        if not os.path.exists(path):
            raise SystemExit(f"missing intermediate result for {ct!r} ({path}) -- run "
                              f"`--cell-types {ct} --assemble-only=false` first")
        with open(path, "rb") as f:
            caches[ct] = pickle.load(f)

    fig, axes = plt.subplots(1, len(cell_types), figsize=(5.3 * len(cell_types), 5.5), squeeze=False)
    axes = axes[0]
    summary_rows, all_sig_rows = [], []
    fdr_thresh = fdr_thresh_fallback

    for ax, ct in zip(axes, cell_types):
        c = caches[ct]
        fdr_thresh = c["fdr_thresh"]
        eval_df, row = c["eval_df"], c["row"]
        all_sig_rows.append(c["sig_out"])
        summary_rows.append(row)
        n_eval = row["n_evaluable_in_replication"]

        if n_eval >= 2:
            same_sign, r_repl = row["sign_concordance_n"], row["effect_size_replication_r"]
            color = CELLTYPE_COLOR[ct]
            lim = max(0.5, eval_df.abs().values.max() * 1.1)
            ax.scatter(eval_df["coef_disc"], eval_df["mean_repl"], s=28, alpha=0.65, color=color)
            ax.plot([-lim, lim], [-lim, lim], "k--", linewidth=1, alpha=0.5, label="y = x")
            ax.axhline(0, color="#ccc", linewidth=0.8, zorder=0)
            ax.axvline(0, color="#ccc", linewidth=0.8, zorder=0)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xlabel(f"discovery half (n={len(c['disc_use'])} donors): bootstrap coef")
            ax.set_ylabel(f"replication half (n={len(c['repl_use'])} donors): mean corr")
            title = (f"{ct}: {n_eval} FDR<{fdr_thresh}-significant genes\n"
                     f"replication r={r_repl:.3f}, sign concordance {same_sign}/{n_eval} "
                     f"({same_sign/n_eval*100:.0f}%)")
            if len(c["disc"]) < 6 or len(c["repl"]) < 6:
                title += "\nEXPLORATORY (too few donors for a robust split)"
            ax.set_title(title, fontsize=10)
            ax.legend(frameon=False, fontsize=8, loc="upper left")
        else:
            ax.text(0.5, 0.5, f"{ct}\nonly {n_eval} evaluable gene(s)", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle("Discovery (memento one-sample bootstrap test, FDR-selected) -> replication effect-size check\n"
                 "(non-circular: genes nominated on one donor half, evaluated only on the other)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig_path = os.path.join(OUT_DIR, "discovery_replication_ht.png")
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    print(f"wrote {fig_path}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUT_DIR, "discovery_replication_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"wrote {summary_path}")

    sig_all = pd.concat(all_sig_rows, ignore_index=True)
    sig_path = os.path.join(OUT_DIR, "discovery_significant_genes.csv")
    with open(sig_path, "w") as f:
        f.write(CSV_CAVEAT_HEADER.format(fdr=fdr_thresh))
        sig_all.to_csv(f, index=False)
    print(f"wrote {sig_path}")

    print("\n=== summary ===")
    print(summary_df.to_string(index=False))

    for ct in cell_types:
        os.remove(TMP_PKL(ct))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-types", default="fibroblast,enterocyte,macrophage",
                     help="comma-separated; each can be run as its own invocation (see run_one)")
    ap.add_argument("--num-boot", type=int, default=5000)
    ap.add_argument("--num-cpus", type=int, default=12)
    ap.add_argument("--fdr-thresh", type=float, default=0.1)
    ap.add_argument("--random-state", type=int, default=42)
    ap.add_argument("--assemble-only", action="store_true",
                     help="skip computation; build final outputs from cached per-cell-type results")
    ap.add_argument("--skip-assemble", action="store_true",
                     help="compute and cache this invocation's --cell-types only; don't build the "
                          "final combined figure/CSVs yet (use when calling once per cell type -- "
                          "run a final bare/--assemble-only call across ALL cell types afterwards)")
    args = ap.parse_args()
    cell_types = args.cell_types.split(",")

    os.makedirs(OUT_DIR, exist_ok=True)

    if not args.assemble_only:
        print(f"loading {INPUT_H5AD} ...")
        full = anndata.read_h5ad(INPUT_H5AD)
        by_donor, gene_order, obs = load_by_donor(COEXPR_OUTPUT_H5AD)
        for ct in cell_types:
            run_one(ct, args, full, by_donor, gene_order)

    if args.skip_assemble:
        return
    assemble(cell_types, args.fdr_thresh)


if __name__ == "__main__":
    main()
