"""Replicate the LEC6 / ANTXR2 observation in Kong et al. 2023 (CELLxGENE Census).

Kong et al. 2023 (Immunity; Crohn's disease, terminal ileum and colon; 10x 3'),
stromal datasets only, which contain the endothelium. Independent of the
Elmentaite Gut Cell Atlas where LEC6 was defined (different cohort, lab and
chemistry). Census carries only the standardized cell_type.

LABEL ERROR, verified by markers (2026-10-05): in these two stromal datasets
Kong's lymphatic endothelial cells are labelled "lymphocyte" (4,248 cells).
They are not lymphocytes: PTPRC 2.1%, CD3E 0.9%, MS4A1 0.4% detected (the
same background as "endothelial cell"), while CCL21 85%, TFF3 86-96%, LYVE1
51-83%, MMRN1 52-69%, PROX1 39%, RELN 28%, FLT4 22%. The "endothelial cell"
label contains essentially no lymphatic cells (PROX1 <1%, CCL21 2%). So here
LEC = stromal "lymphocyte", BEC = "endothelial cell". No subtype labels are
used; the LEC6 state is probed with Elmentaite's LEC6 gene program.

Steps
  1. Census pull: "endothelial cell" + "lymphocyte" (= LEC, see above) in the
     two stromal datasets, all genes (cached h5ad, regenerable); fibroblasts
     for ANTXR2 / library only.
  2. LEC = "lymphocyte" label; marker table written for the record.
  3. Means-level, donor-paired (>= MIN_CELLS per pool):
     a. ANTXR2 pseudobulk LEC vs BEC, per segment, per disease.
     b. Within LEC, cells ranked by a program score (mean z of log CP10K over
        the program genes, within donor); ANTXR2 pseudobulk in the top vs
        bottom quartile, per donor. Programs: LEC6 (top Elmentaite
        LEC6-vs-other genes up in all 3 paired donors, ANTXR2 excluded), and
        two alternative-state controls: valve (FOXC2, GJA4, CLDN11, ITGA9) and
        capillary (CCL21, LYVE1, TFF3, PDPN). No lymph-node control: Kong
        sampled gut wall only.
  4. Leiden clusters within LEC with LEC6-gene detection, ANTXR2 and median
     UMI per cluster (kong_lec_clusters.csv); LEC UMAP coordinates and feature
     values for plotting.

Usage: python kong_lymphatic.py <workspace dir> <out dir>
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CENSUS = "2025-11-08"
DATASETS = {"6cf3634d-e911-44ad-bf52-c747a9af3c01": "colon",
            "0f4865d5-8000-4f68-8ac7-f5efea9e5e70": "TI"}
OBS = ["dataset_id", "donor_id", "disease", "tissue", "assay", "cell_type", "raw_sum", "nnz"]
MIN_CELLS = 20
LEC_GENES = ["PROX1", "CCL21", "TFF3", "MMRN1", "RELN", "PDPN", "LYVE1"]
BEC_GENES = ["PLVAP", "ACKR1"]
PROGRAMS = {
    "valve": ["FOXC2", "GJA4", "CLDN11", "ITGA9"],
    "capillary": ["CCL21", "LYVE1", "PDPN", "TFF3"],
}


def pull(cache):
    import cellxgene_census
    if os.path.exists(cache):
        import anndata as ad
        return ad.read_h5ad(cache)
    ids = list(DATASETS)
    with cellxgene_census.open_soma(census_version=CENSUS) as c:
        a = cellxgene_census.get_anndata(
            c, organism="Homo sapiens", obs_value_filter=f"dataset_id in {ids} and cell_type in ['endothelial cell', 'lymphocyte']",
            column_names={"obs": OBS, "var": ["feature_name"]})
    a.var_names = a.var.feature_name.astype(str).values
    a.var_names_make_unique()
    a.obs["segment"] = a.obs.dataset_id.map(DATASETS).astype(str)
    a.write_h5ad(cache, compression="gzip")
    return a


def pull_fibroblasts(cache):
    import cellxgene_census
    if os.path.exists(cache):
        return pd.read_parquet(cache)
    ids = list(DATASETS)
    with cellxgene_census.open_soma(census_version=CENSUS) as c:
        a = cellxgene_census.get_anndata(
            c, organism="Homo sapiens", obs_value_filter=f"dataset_id in {ids} and cell_type == 'fibroblast'",
            var_value_filter="feature_name in ['ANTXR2']", column_names={"obs": OBS, "var": ["feature_name"]})
    d = a.obs[["dataset_id", "donor_id", "disease", "raw_sum"]].copy()
    d["antxr2"] = np.asarray(a.X.sum(1)).ravel()
    d["segment"] = d.dataset_id.map(DATASETS).astype(str)
    d.to_parquet(cache)
    return d


def pseudobulk(counts, lib):
    return counts.sum() / lib.sum() * 1e4


def main(ws, out):
    import scanpy as sc
    os.makedirs(out, exist_ok=True)
    cache = os.path.join(out, "_cache")
    os.makedirs(cache, exist_ok=True)
    a = pull(os.path.join(cache, "kong_endothelium.h5ad"))
    fib = pull_fibroblasts(os.path.join(cache, "kong_fibroblast_antxr2.parquet"))
    lec6 = pd.read_csv(os.path.join(ws, "lymphatic", "LEC6_vs_otherLEC_paired.csv"))
    prog6 = [g for g in lec6[(lec6.n_up == lec6.n_units) & (lec6.gene != "ANTXR2")].head(30).gene
             if g in set(a.var_names)]
    programs = {"LEC6": prog6, **PROGRAMS}

    lib = np.asarray(a.X.sum(1)).ravel()
    a.obs["lib"] = lib
    a2 = np.asarray(a[:, "ANTXR2"].X.todense()).ravel()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    for name, genes in [("lec_score", LEC_GENES), ("bec_score", BEC_GENES)]:
        a.obs[name] = np.asarray(a[:, [g for g in genes if g in a.var_names]].X.mean(1)).ravel()

    # LEC = "lymphocyte" label (see module docstring); marker table for the record
    a.obs["is_lec"] = (a.obs.cell_type.astype(str) == "lymphocyte").values
    mk = [g for g in LEC_GENES + BEC_GENES + ["PTPRC", "CD3E", "MS4A1"] if g in a.var_names]
    det = pd.DataFrame((a.layers["counts"][:, [a.var_names.get_loc(g) for g in mk]] > 0).toarray(), columns=mk)
    tab = det.groupby(np.where(a.obs.is_lec, "LEC ('lymphocyte' label)", "BEC ('endothelial cell')")).mean()
    tab["n"] = a.obs.is_lec.map({True: "LEC ('lymphocyte' label)", False: "BEC ('endothelial cell')"}).value_counts()
    tab.to_csv(os.path.join(out, "kong_lec_label_markers.csv"))
    print(tab.round(3).to_string())

    # (a) LEC vs BEC, donor-paired, per segment x disease
    o = a.obs.assign(a2=a2)
    rows = []
    for (seg, dis, don), d in o.groupby(["segment", "disease", "donor_id"], observed=True):
        l, b = d[d.is_lec], d[~d.is_lec]
        if len(l) < MIN_CELLS or len(b) < MIN_CELLS:
            continue
        rows.append(dict(segment=seg, disease=dis, donor=don, n_lec=len(l), n_bec=len(b),
                         lec_cp10k=pseudobulk(l.a2, l.lib), bec_cp10k=pseudobulk(b.a2, b.lib),
                         lec_det=(l.a2 > 0).mean(), bec_det=(b.a2 > 0).mean()))
    lb = pd.DataFrame(rows)
    if lb.empty:
        raise SystemExit("no donor has >= MIN_CELLS LEC and BEC")
    fb = fib.groupby(["segment", "disease", "donor_id"], observed=True).apply(
        lambda d: pd.Series({"fib_cp10k": pseudobulk(d.antxr2, d.raw_sum), "n_fib": len(d)})).reset_index()
    lb = lb.merge(fb.rename(columns={"donor_id": "donor"}), on=["segment", "disease", "donor"], how="left")
    lb["log2_lec_bec"] = np.log2((lb.lec_cp10k + 0.05) / (lb.bec_cp10k + 0.05))
    lb.to_csv(os.path.join(out, "kong_lec_vs_bec_by_donor.csv"), index=False)

    # (b) program-quartile test within LEC, per donor
    L = a[a.obs.is_lec].copy()
    la2 = a2[a.obs.is_lec.values]
    llib = L.obs.lib.values
    qrows = []
    for pname, genes in programs.items():
        genes = [g for g in genes if g in L.var_names]
        X = np.asarray(L[:, genes].X.todense())
        for don in L.obs.donor_id.unique():
            m = (L.obs.donor_id == don).values
            if m.sum() < 4 * MIN_CELLS:
                continue
            Z = (X[m] - X[m].mean(0)) / (X[m].std(0) + 1e-9)
            sc_ = Z.mean(1)
            q1, q3 = np.quantile(sc_, [0.25, 0.75])
            top, bot = sc_ >= q3, sc_ <= q1
            qrows.append(dict(program=pname, donor=don, segment=L.obs.segment.values[m][0],
                              disease=L.obs.disease.values[m][0], n=int(m.sum()),
                              top_cp10k=pseudobulk(la2[m][top], llib[m][top]),
                              bot_cp10k=pseudobulk(la2[m][bot], llib[m][bot]),
                              top_det=(la2[m][top] > 0).mean(), bot_det=(la2[m][bot] > 0).mean()))
    qt = pd.DataFrame(qrows)
    qt["log2_top_bot"] = np.log2((qt.top_cp10k + 0.05) / (qt.bot_cp10k + 0.05))
    qt.to_csv(os.path.join(out, "kong_lec_program_quartiles.csv"), index=False)

    # LEC UMAP (fresh embedding within LEC) for figures
    sc.pp.highly_variable_genes(L, n_top_genes=2000, batch_key="segment")
    sc.pp.pca(L, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(L, n_neighbors=15, random_state=0)
    sc.tl.umap(L, random_state=0)
    # clusters within LEC: is there an LEC6-like state, or only depth structure?
    sc.tl.leiden(L, resolution=0.6, random_state=0, flavor="igraph", n_iterations=2)
    Lc = L.layers["counts"]
    gi = {g: i for i, g in enumerate(L.var_names)}
    crow = []
    for c in L.obs.leiden.cat.categories:
        m = (L.obs.leiden == c).values
        d = L.obs[m]
        det = lambda g: float((Lc[m][:, gi[g]].toarray() > 0).mean()) if g in gi else np.nan
        crow.append(dict(cluster=c, n=int(m.sum()), donors=d.donor_id.nunique(),
                         frac_colon=(d.segment == "colon").mean(), frac_crohn=(d.disease == "Crohn disease").mean(),
                         median_UMI=float(np.median(llib[m])),
                         LEC6_score=float(np.asarray(L[m][:, [g for g in prog6]].X.mean(1)).mean()),
                         ANTXR2_cp10k=pseudobulk(la2[m], llib[m]), ANTXR2_det=float((la2[m] > 0).mean()),
                         **{g: det(g) for g in ["ADAMTS4", "CSF3", "COL6A3", "MFAP4", "FN1", "IGF1",
                                                "LYVE1", "CCL21", "FOXC2", "CLDN11", "PROX1"]}))
    pd.DataFrame(crow).sort_values("median_UMI").to_csv(os.path.join(out, "kong_lec_clusters.csv"), index=False)
    for pname, genes in programs.items():
        L.obs[f"{pname}_score"] = np.asarray(L[:, [g for g in genes if g in L.var_names]].X.mean(1)).ravel()
    L.obs.assign(umap1=L.obsm["X_umap"][:, 0], umap2=L.obsm["X_umap"][:, 1], antxr2_counts=la2).to_csv(
        os.path.join(out, "kong_lec_cells.csv"))
    feats = ["ANTXR2", "PROX1", "LYVE1", "CCL21", "ADAMTS4", "FN1", "MFAP4", "COL6A3", "CSF3",
             "PIEZO2", "IGF1", "CLDN11", "FOXC2", "ACKR4", "MRC2"]
    pd.DataFrame(np.asarray(L[:, [f for f in feats if f in L.var_names]].X.todense()),
                 columns=[f for f in feats if f in L.var_names], index=L.obs_names).to_csv(
        os.path.join(out, "kong_lec_features_log1p.csv"))
    print(f"LEC cells: {L.n_obs} ({L.obs.segment.value_counts().to_dict()})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
