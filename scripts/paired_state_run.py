"""Paired within-patient inflamed-vs-non-inflamed contrast for the Wnt
pathway x output panel, with ANTXR2 and detection-matched controls inserted.

Design notes
------------
* UC patients only. Healthy donors are excluded by construction: the contrast
  is within-patient, so a donor contributing only one state carries no
  information once patient dummies are in the covariate matrix.
* Group = patient x state (SubjState). A patient enters a cell type only if it
  has >= MIN_CELLS in BOTH states there.
* treatment = inflamed indicator; covariate = patient dummies (drop_first).
  memento's _regress_2d residualizes boot_corr AND treatment on covariate with
  cell-count weights, so this is a genuine within-patient contrast rather than
  a cross-sectional one with patient as a nuisance term.
* MIN_CELLS = 100, not the project's earlier 250. The 250 floor was derived at
  q=0.10; at q=0.15 the pegged/undefined rates measured on this study's own
  per-donor estimates are 1.8%/4.1% at 100 cells, both inside the project's
  5% tolerance. 250 would leave 6/3/1 paired patients and is not runnable.
* Controls are rematched inside the UC-only population: ANTXR2 detection is
  higher in UC than pooled across states, so pooled-matched controls would be
  mismatched here.
"""
import os, time, numpy as np, pandas as pd, scipy.sparse as sp, memento

OUT = "paired_state"
MIN_CELLS = 100


def paired_state_contrast(adata_ct, gene1_list, output, donor_col="SubjState",
                          subject_col="Subject", state_col="Health",
                          inflamed_label="Inflamed", q_col="q",
                          num_boot=5000, num_cpus=10, seed=5,
                          min_group_cells=MIN_CELLS, min_perc_group=0.1,
                          filter_mean_thresh=0.07, shrinkage=0.5,
                          trim_percent=0.1):
    ad = adata_ct.copy()
    ad.X = sp.csr_matrix(ad.X)

    memento.setup_memento(ad, q_column=q_col,
                          filter_mean_thresh=filter_mean_thresh,
                          min_cell_count=min_group_cells,
                          shrinkage=shrinkage, trim_percent=trim_percent)
    memento.create_groups(ad, label_columns=[donor_col])
    memento.compute_1d_moments(ad, min_perc_group=min_perc_group,
                               filter_genes=True)

    gl = list(ad.var.index)
    G1 = [g for g in gene1_list if g in gl]
    O = [g for g in output if g in gl]
    if not G1 or not O:
        return None, None, 0, []
    pairs = [(a, b) for a in G1 for b in O]
    memento.compute_2d_moments(ad, pairs)

    groups = list(ad.uns["memento"]["groups"])
    # group name is "sg^<SubjState>"; SubjState is "<Subject>@<Health>"
    meta = []
    for gname in groups:
        val = gname.split("^", 1)[1]
        subj, state = val.rsplit("@", 1)
        meta.append((gname, subj, state))
    M = pd.DataFrame(meta, columns=["group", "subject", "state"])

    treatment = pd.DataFrame(
        {"inflamed": (M.state == inflamed_label).astype(float).values})
    covariate = pd.get_dummies(M.subject, drop_first=True).astype(float)
    covariate.index = treatment.index

    memento.ht_2d_moments(ad, treatment=treatment, covariate=covariate,
                          num_boot=num_boot, num_cpus=num_cpus, verbose=0,
                          random_state=seed)
    ht = memento.get_2d_ht_result(ad).copy()

    moments, cell_counts = memento.get_2d_moments(ad, groupby=None)
    rows = []
    for col in moments.columns:
        if not col.startswith("sg^"):
            continue
        val = col.split("^", 1)[1]
        subj, state = val.rsplit("@", 1)
        for g1, g2, v in zip(moments["gene_1"], moments["gene_2"], moments[col]):
            rows.append(dict(subject=subj, state=state, gene_1=g1, gene_2=g2,
                             corr=float(v), n_cells=cell_counts.get(col)))
    return ht, pd.DataFrame(rows), len(groups), list(M.subject.unique())


def run(E, uc_mask, cl, run_cts, pathway, output, ctrl_df, min_cells=MIN_CELLS):
    os.makedirs(OUT, exist_ok=True)
    for c in run_cts:
        t0 = time.time()
        sel = uc_mask & (cl == c)
        sub = E[sel].copy()
        # keep only patients with >= min_cells in BOTH states
        n = sub.obs.groupby(["Subject", "Health"], observed=True).size().unstack(fill_value=0)
        keep = n[(n.get("Inflamed", 0) >= min_cells) &
                 (n.get("Non-inflamed", 0) >= min_cells)].index
        sub = sub[sub.obs.Subject.isin(keep)].copy()
        ctrls = ctrl_df.query("cell_type == @c").control.tolist()
        g1 = list(pathway) + ["ANTXR2"] + ctrls
        ht, pg, ng, subs = paired_state_contrast(sub, g1, output,
                                                 min_group_cells=min_cells)
        if ht is None:
            print(f"{c}: no testable pairs", flush=True)
            continue
        ht["cell_type"] = c
        ht["klass"] = np.where(ht.gene_1 == "ANTXR2", "antxr2",
                               np.where(ht.gene_1.isin(ctrls), "control", "pathway"))
        pg["cell_type"] = c
        tag = c.replace(" ", "_").replace("+", "p")
        ht.to_csv(f"{OUT}/{tag}_ht.csv", index=False)
        pg.to_csv(f"{OUT}/{tag}_perdonor.csv", index=False)
        print(f"{c}: {len(keep)} paired patients, {ng} groups, {len(ht)} pairs, "
              f"{time.time()-t0:.0f}s", flush=True)
    return sorted(os.listdir(OUT))
