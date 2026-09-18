"""ANTXR2 x Wnt-arm co-expression per gut cell type, healthy vs inflamed.

Motivation
----------
The cross-tissue arm (HFS_ATLAS_NOTES.md) found that in healthy human intestine
ANTXR2 is stromal and glial, with every epithelial label at the bottom of the
per-cell-type ranking. That is in tension with the mouse work placing the Wnt
role in epithelium cell-autonomously. The reconciling hypothesis is that the
Wnt function is **injury-conditional** — in which case a healthy atlas cannot
show it by construction, and the test has to be healthy vs inflamed in the same
tissue.

This script computes, for each (state, donor, cell_type) group, the memento
point-estimate correlation between ANTXR2 and each gene of
``GENE_PANEL_ARMS["WNT_ARM"]``, averages across donors with equal weight (the
project's binding aggregation convention), and contrasts states.

Status: EXPLORATORY
-------------------
Point estimates only, no bootstrap hypothesis test. Per the project's standing
rule, a **positive** correlation is not a claim until it clears the anchor-gene
control, because genes co-vary for reasons that have nothing to do with the
pathway. This script deliberately does **not** reimplement that control: the
confirmatory route is ``coexpression_discovery_replication.py`` +
``anchor_gene_control.py``, whose outcome differs by cell type and must not be
applied uniformly. Treat everything here as a screen that decides which
cell_type x gene pairs are worth the bootstrap.

Anti-correlations are the stronger result class here, as elsewhere in this
project: shared-activation confounding can manufacture co-expression but not
mutual exclusivity.

Usage
-----
    python scripts/wnt_state_coexpression.py --h5ad <colon.h5ad> \
        --state-col Health --donor-col Subject --celltype-col Cluster

The h5ad must carry raw counts, a per-cell capture-rate column (see
``capture_rate.py``), and a state column whose values map onto
{healthy, non-inflamed, inflamed} via ``--state-map``.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from config import GENE_PANEL_ARMS

TARGET = "ANTXR2"
WNT = GENE_PANEL_ARMS["WNT_ARM"]

DEFAULT_STATE_MAP = {"Healthy": "healthy", "Non-inflamed": "non_inflamed",
                     "Inflamed": "inflamed"}


# --------------------------------------------------------------- memento
def patch_memento():
    """Repair estimator._corr_from_cov (see ANALYSIS_SUMMARY.md).

    The shipped routine initialises the correlation array to 5.0, overwrites
    only entries with positive variance, then clips to [-1, 1] unconditionally
    — so any untouched entry silently becomes a perfect correlation of 1.0
    rather than NaN. Verified still present in the installed build.
    """
    from memento import estimator

    def _corr_from_cov(cov, var_1, var_2, boot=False):
        out = np.full(cov.shape, np.nan)
        ok = (var_1 > 0) & (var_2 > 0)
        out[ok] = cov[ok] / np.sqrt(var_1[ok] * var_2[ok])
        return np.clip(out, -1.0, 1.0, out=out, where=~np.isnan(out))

    estimator._corr_from_cov = _corr_from_cov
    return _corr_from_cov


def group_correlations(adata, target, partners, donor_col, celltype_col,
                       q_column, min_group_cells=50, min_perc_group=0.1,
                       filter_mean_thresh=0.07, shrinkage=0.5, trim_percent=0.1):
    """memento point-estimate corr(target, partner) per donor x cell_type."""
    import memento

    adata = adata.copy()
    adata.X = sp.csr_matrix(adata.X)
    original_var_index = adata.var.index.copy()

    memento.setup_memento(adata, q_column=q_column,
                          filter_mean_thresh=filter_mean_thresh,
                          min_cell_count=min_group_cells,
                          shrinkage=shrinkage, trim_percent=trim_percent)
    memento.create_groups(adata, label_columns=[donor_col, celltype_col])
    memento.compute_1d_moments(adata, min_perc_group=min_perc_group,
                               filter_genes=True)

    # NB compute_2d_moments maps pairs through `adata.var.index`, so pairs must
    # be NAMES from the post-filter var index -- not positional indices.
    gene_list = list(adata.var.index)
    if target not in gene_list:
        raise SystemExit(f"{target} did not survive the expression filter")
    keep = [g for g in partners if g in gene_list]
    dropped = sorted(set(partners) - set(keep))

    memento.compute_2d_moments(adata, [(target, g) for g in keep])

    # NB do NOT use get_2d_moments(groupby=...): it averages the per-group
    # correlations WEIGHTED BY CELL COUNT, which violates this project's
    # binding donor-equal-weight convention. Take the raw per-group columns
    # and average them unweighted below.
    moments, cell_counts = memento.get_2d_moments(adata, groupby=None)

    # group columns are 'sg^<donor>^<cell_type>' (create_groups delimiter '^')
    long = []
    for col in moments.columns:
        if not col.startswith("sg^"):
            continue
        parts = col.split("^")[1:]
        if len(parts) != 2:
            raise SystemExit(f"unexpected group key {col!r}: expected "
                             f"sg^<donor>^<cell_type>")
        donor, ct = parts
        for g, v in zip(moments["gene_2"], moments[col]):
            long.append(dict(donor=donor, cell_type=ct, partner=g,
                             corr=float(v), n_cells=cell_counts.get(col)))

    # gene_filter stays keyed on the ORIGINAL var index even after
    # compute_1d_moments subsets adata.var -- index it with original_var_index
    # if you need per-group survival flags.
    _ = original_var_index
    return pd.DataFrame(long), keep, dropped


# --------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--state-col", required=True)
    ap.add_argument("--donor-col", required=True)
    ap.add_argument("--celltype-col", required=True)
    ap.add_argument("--q-column", default="q")
    ap.add_argument("--state-map", default=None,
                    help="JSON mapping raw state values onto "
                         "healthy/non_inflamed/inflamed")
    ap.add_argument("--min-group-cells", type=int, default=50)
    ap.add_argument("--out-dir", default="output/wnt_state")
    args = ap.parse_args()

    import anndata as ad

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    smap = json.loads(args.state_map) if args.state_map else DEFAULT_STATE_MAP

    patch_memento()
    adata = ad.read_h5ad(args.h5ad)
    adata.obs["_state"] = adata.obs[args.state_col].map(smap)
    unmapped = adata.obs["_state"].isna().sum()
    if unmapped:
        print(f"WARNING {unmapped} cells have a state value outside --state-map "
              f"and are dropped: "
              f"{sorted(set(adata.obs[args.state_col][adata.obs['_state'].isna()]))}")
        adata = adata[adata.obs["_state"].notna()].copy()

    rows = []
    for state in ["healthy", "non_inflamed", "inflamed"]:
        sub = adata[adata.obs["_state"] == state]
        if sub.n_obs == 0:
            print(f"state {state}: absent from this dataset, skipping")
            continue
        print(f"\n=== {state}: {sub.n_obs} cells, "
              f"{sub.obs[args.donor_col].nunique()} donors ===", flush=True)

        df, kept, dropped = group_correlations(
            sub, TARGET, WNT, args.donor_col, args.celltype_col, args.q_column,
            min_group_cells=args.min_group_cells)
        if dropped:
            print(f"  Wnt genes below expression filter: {dropped}")
        print(f"  {len(kept)} Wnt partners x "
              f"{df[['donor', 'cell_type']].drop_duplicates().shape[0]} groups")
        df.insert(0, "state", state)
        rows.append(df)

    if not rows:
        raise SystemExit("no state produced any group -- check --state-map")
    long = pd.concat(rows, ignore_index=True)
    long.to_csv(out / "wnt_corr_by_donor.csv", index=False)

    # donor-equal-weight: each donor contributes one value per cell type
    avg = (long.groupby(["state", "cell_type", "partner"])["corr"]
           .agg(mean_corr="mean", sd="std", n_donors="count").reset_index())
    avg.to_csv(out / "wnt_corr_donor_averaged.csv", index=False)

    # state contrast, only where both arms have >= 2 donors
    wide = avg.pivot_table(index=["cell_type", "partner"], columns="state",
                           values="mean_corr")
    nd = avg.pivot_table(index=["cell_type", "partner"], columns="state",
                         values="n_donors")
    if {"healthy", "inflamed"}.issubset(wide.columns):
        delta = (wide["inflamed"] - wide["healthy"]).rename("delta_infl_minus_healthy")
        nd = nd.reindex(wide.index).fillna(0)
        # A Wnt gene can pass memento's expression filter in one state and not
        # the other. Contrasting those is comparing a correlation against a
        # missing one, so require BOTH arms present with >=2 donors; genes
        # dropped asymmetrically are reported separately, not silently.
        both = wide["healthy"].notna() & wide["inflamed"].notna()
        ok = both & (nd["healthy"] >= 2) & (nd["inflamed"] >= 2)
        contrast = pd.concat([wide, delta, nd.add_prefix("n_donors_")],
                             axis=1).loc[ok].reset_index()
        contrast.sort_values("delta_infl_minus_healthy").to_csv(
            out / "wnt_state_contrast.csv", index=False)
        asym = pd.concat([wide, nd.add_prefix("n_donors_")], axis=1).loc[~both].reset_index()
        asym.to_csv(out / "wnt_state_asymmetric.csv", index=False)
        print(f"\n{len(contrast)} cell_type x Wnt-gene contrasts with >=2 donors "
              f"in both arms -> {out/'wnt_state_contrast.csv'}")
        print(f"{len(asym)} pairs present in only one state (filter-asymmetric) "
              f"-> {out/'wnt_state_asymmetric.csv'}")
    else:
        print("\nOnly one state present — no contrast written. "
              "The healthy arm alone cannot address the injury-conditional "
              "hypothesis; it is a baseline, not a result.")

    print(f"{len(long)} donor-level correlations, "
          f"{avg.cell_type.nunique()} cell types -> {out}")


if __name__ == "__main__":
    main()
