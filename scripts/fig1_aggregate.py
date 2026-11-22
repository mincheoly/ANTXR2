"""Figure 1 (Act 1): donor-level aggregation of the per-cell panel extracts
(fig1_extract.py), then donor-equal-weight summaries.

Unit of aggregation is atlas | donor | tissue | cell_type | annotation | assay |
suspension | disease. Per unit and gene: pooled CP10k (sum counts / sum library
size x 1e4), mean of per-cell CP10k (the memento-style mean up to the constant
1/q), and detection (fraction of cells with >=1 UMI). Summaries then average
units' values within donor first, then across donors (donor-equal weight).

    -> /data/antxr2/fig1/donor_units.parquet
Usage: python fig1_aggregate.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/data/antxr2/fig1")
CELLS = ROOT / "cells"

# atlas -> obs column holding the author's fine annotation (kept beside the ontology label)
ANNOT = {
    "tabula_sapiens": "free_annotation",
    "cross_tissue_immune": "Manually_curated_celltype",
    "gca_healthy": "level_3_annot",
    "gca_extended": "level_3_annot",
    "skin_fibroblast": "cluster_robust",
    "ibdverse": "author_cell_type",
}
# atlas -> obs column for tissue at the resolution we want
TISSUE = {"gca_healthy": "organ_groups", "gca_extended": "organ_groups"}
# atlas -> extra per-unit grouping columns (disease state beyond `disease`)
EXTRA = {"ibdverse": ["author_inflammation_status", "author_cell_type_category"],
         "skin_fibroblast": ["Author", "development_stage"],
         "gca_extended": ["study", "age_unified", "level_1_annot", "level_2_annot"],
         "gca_healthy": ["study", "age_unified", "level_1_annot", "level_2_annot"],
         "tabula_sapiens": ["compartment", "tissue_in_publication"]}
KEYS = ["atlas", "donor", "tissue", "cell_type", "annot", "assay", "suspension", "disease", "extra"]


def units(atlas):
    df = pd.read_parquet(CELLS / f"{atlas}.parquet")
    meta = json.loads((CELLS / f"{atlas}.meta.json").read_text())
    genes = meta["genes"]
    # No is_primary_data filter: the flag de-duplicates cells across datasets of
    # a collection, and GCA "Healthy reference" is entirely is_primary_data=False
    # (its primary copies sit in "Extended - All genes"). Each atlas here is
    # aggregated on its own, so within-file duplication does not arise; overlap
    # between gca_healthy and gca_extended is handled by never pooling the two.
    df = df[df.total > 0]
    donor = df["donorID_unified"] if "donorID_unified" in df else df["donor_id"]
    k = pd.DataFrame({
        "atlas": atlas,
        "donor": donor.astype(str).values,
        "tissue": df[TISSUE.get(atlas, "tissue")].astype(str).values,
        "cell_type": df["cell_type"].astype(str).values,
        "annot": df[ANNOT[atlas]].astype(str).values if ANNOT[atlas] in df else "",
        "assay": df["assay"].astype(str).values,
        "suspension": df["suspension_type"].astype(str).values,
        "disease": df["disease"].astype(str).values,
        "extra": (df[EXTRA[atlas]].astype(str).agg("|".join, axis=1).values if atlas in EXTRA else ""),
    })
    C = df[genes].to_numpy(np.float64)
    tot = df["total"].to_numpy(np.float64)
    N = C / tot[:, None] * 1e4
    frame = pd.concat([k.reset_index(drop=True),
                       pd.DataFrame(C, columns=[f"sum_{g}" for g in genes]),
                       pd.DataFrame(N, columns=[f"cp10k_{g}" for g in genes]),
                       pd.DataFrame((C > 0).astype(np.float32), columns=[f"det_{g}" for g in genes])], axis=1)
    frame["total"] = tot
    frame["n_cells"] = 1
    agg = {**{f"sum_{g}": "sum" for g in genes}, **{f"cp10k_{g}": "mean" for g in genes},
           **{f"det_{g}": "mean" for g in genes}, "total": "sum", "n_cells": "sum"}
    u = frame.groupby(KEYS, observed=True, sort=False).agg(agg).reset_index()
    for g in genes:
        u[f"pooled_{g}"] = u[f"sum_{g}"] / u.total * 1e4
    u = u.rename(columns={f"cp10k_{g}": f"mean_{g}" for g in genes})
    u["umi_per_cell"] = u.total / u.n_cells
    print(f"{atlas}: {len(df)} cells -> {len(u)} units, {u.donor.nunique()} donors", flush=True)
    return u


if __name__ == "__main__":
    atl = sorted(p.stem for p in CELLS.glob("*.parquet"))
    out = pd.concat([units(a) for a in atl], ignore_index=True)
    out.to_parquet(ROOT / "donor_units.parquet", index=False)
    print(f"-> {ROOT / 'donor_units.parquet'} ({len(out)} units)")
