"""One-sample memento test of corr(Wnt pathway gene, Wnt output gene) per cell type.

memento's own one-sample path: when the treatment column is all ones,
_regress_2d returns np.average(boot_corr, axis=0, weights=Nc_list) -- the
CELL-COUNT-WEIGHTED mean across donor groups -- and _compute_asl turns the
bootstrap distribution into a two-sided p-value. This is deliberately NOT the
donor-equal-weight average used earlier in the project.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, memento


def pathway_output_onesample(adata_ct, pathway, output, donor_col="SubjState",
                             q_col="q", num_boot=5000, num_cpus=10, seed=5,
                             min_group_cells=50, min_perc_group=0.1,
                             filter_mean_thresh=0.07, shrinkage=0.5,
                             trim_percent=0.1):
    """Return (ht_result, per_group_long, n_groups) for one cell type."""
    ad = adata_ct.copy()
    ad.X = sp.csr_matrix(ad.X)

    memento.setup_memento(ad, q_column=q_col,
                          filter_mean_thresh=filter_mean_thresh,
                          min_cell_count=min_group_cells,
                          shrinkage=shrinkage, trim_percent=trim_percent)
    memento.create_groups(ad, label_columns=[donor_col])
    memento.compute_1d_moments(ad, min_perc_group=min_perc_group,
                               filter_genes=True)

    # pairs must be NAMES from the POST-filter var index
    gl = list(ad.var.index)
    P = [g for g in pathway if g in gl]
    O = [g for g in output if g in gl]
    if not P or not O:
        return None, None, 0
    pairs = [(a, b) for a in P for b in O]
    memento.compute_2d_moments(ad, pairs)

    groups = list(ad.uns["memento"]["groups"])
    # all-ones treatment -> memento's one-sample, cell-count-weighted branch
    treatment = pd.DataFrame(np.ones((len(groups), 1)), columns=["intercept"])
    memento.ht_2d_moments(ad, treatment=treatment, num_boot=num_boot,
                          num_cpus=num_cpus, verbose=0, random_state=seed)
    ht = memento.get_2d_ht_result(ad).copy()

    # per-donor point estimates alongside, for the unweighted comparison
    moments, cell_counts = memento.get_2d_moments(ad, groupby=None)
    rows = []
    for col in moments.columns:
        if not col.startswith("sg^"):
            continue
        donor = col.split("^", 1)[1]
        for g1, g2, v in zip(moments["gene_1"], moments["gene_2"], moments[col]):
            rows.append(dict(donor=donor, gene_1=g1, gene_2=g2,
                             corr=float(v), n_cells=cell_counts.get(col)))
    return ht, pd.DataFrame(rows), len(groups)
