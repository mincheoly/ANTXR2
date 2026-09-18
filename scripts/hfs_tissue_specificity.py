"""Collagen VI burden per unit ANTXR2 across fibroblast populations.

The question this arm asks is not "where is ANTXR2 expressed" but "where is
there the most collagen VI per receptor to clear". Two readouts are produced
per tissue, both from library-size-normalised pseudobulk CPM:

``col6_frac_of_collagen``  COL6 as a share of the tissue's total collagen
                           output — a composition measure, robust to how much
                           collagen the tissue makes in absolute terms.
``col6_per_antxr2``        COL6 CPM divided by ANTXR2 CPM — the substrate:
                           receptor ratio. This is the load-per-receptor
                           metric referred to elsewhere in the project.

Suspension type is carried as a column, never averaged over. Single-nucleus
populations are written out but must not be compared to whole-cell ones on
``col6_per_antxr2`` without the correction from ``hfs_suspension_bias.py`` —
see HFS_ATLAS_NOTES.md.

Outputs
-------
    hfs_col6_specificity_by_tissue.csv   one row per tissue (pooled fibroblasts)
    hfs_col6_specificity_by_donor.csv    one row per tissue x label x donor
"""
import re

import numpy as np
import pandas as pd

from hfs_config import (COLLAGEN_GROUPS, COLLAGEN_RE, FIBROBLAST_LABELS, ORAL_SITES,
                        TISSUE_INVOLVEMENT, HFS_PB_DIR, HFS_OUTPUT_DIR,
                        MIN_CELLS_PER_DONOR_GROUP, MIN_CELLS_SMALL_ATLAS,
                        MIN_CELLS_POOLED)

GRP = COLLAGEN_GROUPS


def _load(name):
    Z = np.load(HFS_PB_DIR / f"{name}.npz", allow_pickle=False)
    genes = np.array([str(x) for x in Z["genes"]])
    groups = np.array([str(x) for x in Z["groups"]])
    return Z["sum"].astype(np.float64), genes, groups, Z["ncell"]


def _collagen_index(genes):
    gi = {g: i for i, g in enumerate(genes)}
    all_col = [i for i, g in enumerate(genes) if re.fullmatch(COLLAGEN_RE, g)]
    col6 = [gi[g] for g in GRP["COL6"] if g in gi]
    return gi, all_col, col6


def _metrics(summed, genes, gi, all_col, col6):
    """CPM-space collagen metrics for one pooled group's summed counts."""
    tot = summed.sum()
    cpm = summed / tot * 1e6
    col_all = cpm[all_col].sum()
    out = dict(col6_cpm=cpm[col6].sum(), col_all_cpm=col_all,
               ANTXR2=cpm[gi["ANTXR2"]] if "ANTXR2" in gi else np.nan)
    out["col6_frac_of_collagen"] = out["col6_cpm"] / col_all
    out["col6_per_antxr2"] = out["col6_cpm"] / out["ANTXR2"]
    for name, members in GRP.items():
        idx = [gi[g] for g in members if g in gi]
        out[f"{name}_pct_of_collagen"] = cpm[idx].sum() / col_all * 100
    out["other_pct_of_collagen"] = 100 - sum(out[f"{k}_pct_of_collagen"] for k in GRP)
    return out


# ------------------------------------------------------- Tabula Sapiens
def tabula_sapiens():
    S, genes, groups, N = _load("ts_donor_pseudobulk")
    gi, all_col, col6 = _collagen_index(genes)
    tot = S.sum(1)

    donor = pd.DataFrame({
        "tissue": [g.split(" | ")[0] for g in groups],
        "cell_type": [g.split(" | ")[1] for g in groups],
        "donor": [g.split(" | ")[2] for g in groups],
        "n_cells": N,
        "col6_cpm": S[:, col6].sum(1) / tot * 1e6,
        "col_all_cpm": S[:, all_col].sum(1) / tot * 1e6,
        "ANTXR2": S[:, gi["ANTXR2"]] / tot * 1e6})
    donor["col6_frac_of_collagen"] = donor.col6_cpm / donor.col_all_cpm
    donor["col6_per_antxr2"] = donor.col6_cpm / donor.ANTXR2
    donor = donor[donor.n_cells >= MIN_CELLS_PER_DONOR_GROUP].copy()
    donor["involvement"] = donor.tissue.map(TISSUE_INVOLVEMENT).fillna("not reported")

    tis = np.array([g.split(" | ")[0] for g in groups])
    is_fib = np.array([g.split(" | ")[1] in FIBROBLAST_LABELS for g in groups])
    rows = []
    for t in sorted(set(tis[is_fib])):
        j = np.where(is_fib & (tis == t))[0]
        if N[j].sum() < MIN_CELLS_POOLED:
            continue
        r = dict(tissue=t, n_cells=int(N[j].sum()),
                 labels="+".join(sorted({groups[k].split(" | ")[1] for k in j})),
                 donors=len({groups[k].split(" | ")[2] for k in j}))
        r.update(_metrics(S[j].sum(0), genes, gi, all_col, col6))
        rows.append(r)
    pooled = pd.DataFrame(rows)
    pooled["involvement"] = pooled.tissue.map(TISSUE_INVOLVEMENT).fillna("not reported")
    pooled["source"] = "Tabula Sapiens"
    pooled["susp"] = "cell"
    pooled["assay"] = "10x 3' v3"
    pooled["tissue_label"] = pooled.tissue.str.replace("_", " ", regex=False)
    return pooled, donor


# ------------------------------------------------------- oral atlas
def oral():
    S, genes, groups, N = _load("oral_v3_donor_pseudobulk")
    gi, all_col, col6 = _collagen_index(genes)
    display = {"gingiva": "Gingiva*", "buccal mucosa": "Buccal mucosa*",
               "hard palate": "Hard palate*"}
    inv = {"Gingiva*": "affected (core)", "Buccal mucosa*": "occasional",
           "Hard palate*": "occasional"}
    site = np.array([g.split(" | ")[0] for g in groups])
    rows = []
    for s in ORAL_SITES:
        j = [k for k in np.where(site == s)[0] if N[k] >= MIN_CELLS_PER_DONOR_GROUP]
        if not j:
            continue
        r = dict(tissue=display[s], n_cells=int(N[j].sum()), donors=len(j),
                 labels="fibroblast")
        r.update(_metrics(S[j].sum(0), genes, gi, all_col, col6))
        rows.append(r)
    df = pd.DataFrame(rows)
    df["involvement"] = df.tissue.map(inv)
    df["source"] = "Human Oral & Craniofacial Cell Atlas"
    df["susp"] = "cell"
    df["assay"] = "10x 3' v3"
    df["tissue_label"] = df.tissue
    return df


# ------------------------------------------------------- joint / tendon
def joint_tissues():
    """Synovium and the two tendon atlases.

    The tendon atlases are single-nucleus; that is recorded here and acted on in
    hfs_suspension_bias.py rather than silently pooled with whole-cell tissues.
    """
    meta = {
        "synovium_pseudobulk": dict(source="JIA synovium atlas", susp="cell",
                                    assay="10x 5' v2", display="Synovium (JIA)\u2020",
                                    keep=None),
        "tendon_ach_pseudobulk": dict(source="Achilles tendon atlas", susp="nucleus",
                                      assay="10x 3' v3", display="Tendon, Achilles\u2021",
                                      keep=None),
        "tendon_quad_pseudobulk": dict(source="Quadriceps tendon atlas", susp="nucleus",
                                       assay="10x 3' v3", display="Tendon, quadriceps\u2021",
                                       keep=lambda g: "fibroblast" in g),
    }
    rows = []
    for name, m in meta.items():
        S, genes, groups, N = _load(name)
        gi, all_col, col6 = _collagen_index(genes)
        tissues = np.array([g.split(" | ")[0] for g in groups])
        for t in sorted(set(tissues)):
            j = [k for k, g in enumerate(groups)
                 if g.split(" | ")[0] == t and N[k] >= MIN_CELLS_SMALL_ATLAS
                 and (m["keep"] is None or m["keep"](g))]
            if not j or S[j].sum() <= 0:
                continue
            r = dict(tissue=t, n_cells=int(N[j].sum()),
                     donors=len({groups[k].split(" | ")[-1] for k in j}),
                     labels="fibroblast", source=m["source"], susp=m["susp"],
                     assay=m["assay"], tissue_label=m["display"],
                     involvement="affected (core)")
            r.update(_metrics(S[j].sum(0), genes, gi, all_col, col6))
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    ts_pooled, ts_donor = tabula_sapiens()
    parts = [ts_pooled, oral(), joint_tissues()]
    cols = sorted(set().union(*[set(p.columns) for p in parts]))
    allt = (pd.concat([p.reindex(columns=cols) for p in parts], ignore_index=True)
            .sort_values("col6_frac_of_collagen", ascending=False)
            .reset_index(drop=True))
    allt["nucleus"] = allt.susp.eq("nucleus")

    lead = ["tissue", "tissue_label", "involvement", "source", "susp", "assay",
            "n_cells", "donors", "labels", "ANTXR2", "col6_cpm", "col_all_cpm",
            "col6_frac_of_collagen", "col6_per_antxr2"]
    allt = allt[lead + [c for c in allt.columns if c not in lead]]

    allt.to_csv(HFS_OUTPUT_DIR / "hfs_col6_specificity_by_tissue.csv", index=False)
    ts_donor.to_csv(HFS_OUTPUT_DIR / "hfs_col6_specificity_by_donor.csv", index=False)
    print(f"{len(allt)} tissue rows ({allt.nucleus.sum()} single-nucleus), "
          f"{len(ts_donor)} donor rows -> {HFS_OUTPUT_DIR}")
    return allt, ts_donor


if __name__ == "__main__":
    main()
