"""Build every pseudobulk matrix the HFS cross-tissue arm consumes.

Each entry below records the cell filter and the grouping key for one atlas.
The grouping key is load-bearing: downstream scripts split it on " | " and rely
on field order, so changing a key here silently changes what the metrics mean.

    ts_donor        Tabula Sapiens fibroblasts   tissue | cell_type | donor
    gut             intestine, all cell types    cell_type
    skin / oral     all cell types               cell_type
    oral_v3_donor   oral fibroblasts, 3' v3      site | donor
    synovium        JIA synovial fibroblasts     "synovium (JIA)" | donor
    tendon_ach      Achilles tendon fibroblasts  "tendon, Achilles" | donor
    tendon_quad     quadriceps stromal cells     "tendon, quadriceps" | cell_type | donor
    kidney_susp     kidney interstitial fibro.   suspension_type | disease
    heart_susp      cardiac fibroblasts          suspension_type | region

The two ``*_susp`` matrices exist only for the suspension-type control: they are
the datasets that contain both single-cell and single-nucleus arms of the same
tissue, so a nucleus/cell bias can be measured rather than assumed.

Usage
-----
    python scripts/hfs_build_pseudobulk.py            # everything missing
    python scripts/hfs_build_pseudobulk.py gut skin   # named targets, forced
"""
import sys
import time

import h5py
import numpy as np
import pandas as pd

from hfs_config import (PANEL, FIBROBLAST_LABELS, ORAL_SITES, HFS_RAW_DIR, HFS_PB_DIR,
                        PSEUDOBULK_TARGET_NNZ)
from hfs_pseudobulk import read_obs_column, extract, pseudobulk

NNZ = PSEUDOBULK_TARGET_NNZ


def _obs(path, cols):
    with h5py.File(path, "r") as f:
        g = f["obs"]
        out = {}
        for c in cols:
            if c in g:
                v = read_obs_column(g, c)
                if v is not None:
                    out[c] = v
    return pd.DataFrame(out)


def _raw(key):
    return HFS_RAW_DIR / f"{key}.h5ad"


def _out(name):
    return HFS_PB_DIR / f"{name}.npz"


# --------------------------------------------------------------- builders
def build_ts_donor():
    """Tabula Sapiens stromal: 3' v3 fibroblasts, per tissue x label x donor."""
    p = _raw("ts_stromal")
    ob = _obs(p, ["assay", "tissue_in_publication", "cell_type", "donor_id"])
    mask = ((ob.assay == "10x 3' v3") & ob.cell_type.isin(FIBROBLAST_LABELS)).values
    key = (ob.tissue_in_publication.astype(str) + " | " + ob.cell_type.astype(str)
           + " | " + ob.donor_id.astype(str)).values
    return pseudobulk(p, mask, key, _out("ts_donor_pseudobulk"), NNZ)


def build_gut():
    """Intestine: adult, healthy, primary cells only; grouped by cell type.

    Every cell type is kept (not just fibroblasts) because the gut result turns
    on where ANTXR2 sits *within* the tissue, epithelium included.
    """
    p = _raw("gut")
    extract(p, _out("gut_panel"), PANEL, target_nnz=NNZ)
    ob = _obs(p, ["cell_type", "disease", "is_primary_data", "Age_group"])
    mask = ((ob.Age_group == "Adult") & (ob.disease == "normal")
            & (ob.is_primary_data.astype(str) == "True")).values
    return pseudobulk(p, mask, ob.cell_type.values.astype(str),
                      _out("gut_pseudobulk"), NNZ, count_nnz=True)


def build_simple(key, out_name):
    """Skin / oral: all cells, grouped by cell type (whole-tissue reference)."""
    p = _raw(key)
    extract(p, _out(f"{key}_panel"), PANEL, target_nnz=NNZ)
    ob = _obs(p, ["cell_type"])
    mask = np.ones(len(ob), bool)
    return pseudobulk(p, mask, ob.cell_type.values.astype(str), _out(out_name), NNZ,
                      count_nnz=True)


def build_oral_v3_donor():
    """Oral atlas: healthy 3' v3 fibroblasts at three sites, per donor.

    Chemistry is pinned to 3' v3 so the oral sites are comparable to Tabula
    Sapiens rather than confounded with capture chemistry.
    """
    p = _raw("oral")
    ob = _obs(p, ["assay", "tissue", "cell_type", "donor_id", "sampled_site_condition"])
    mask = ((ob.assay == "10x 3' v3") & ob.tissue.isin(ORAL_SITES)
            & (ob.cell_type == "fibroblast")
            & (ob.sampled_site_condition == "healthy")).values
    key = (ob.tissue.astype(str) + " | " + ob.donor_id.astype(str)).values
    return pseudobulk(p, mask, key, _out("oral_v3_donor_pseudobulk"), NNZ)


def build_synovium():
    p = _raw("synovium")
    ob = _obs(p, ["cell_type", "tissue", "assay", "donor_id"])
    mask = ((ob.cell_type == "fibroblast") & (ob.tissue == "layer of synovial tissue")).values
    key = ("synovium (JIA) | " + ob.donor_id.astype(str)).values
    return pseudobulk(p, mask, key, _out("synovium_pseudobulk"), NNZ)


def build_tendon_ach():
    p = _raw("tendon_ach")
    ob = _obs(p, ["cell_type", "donor_id"])
    mask = (ob.cell_type == "fibroblast").values
    key = ("tendon, Achilles | " + ob.donor_id.astype(str)).values
    return pseudobulk(p, mask, key, _out("tendon_ach_pseudobulk"), NNZ)


def build_tendon_quad():
    p = _raw("tendon_quad")
    ob = _obs(p, ["cell_type", "donor_id", "disease"])
    mask = (ob.cell_type.isin(["fibroblast", "osteoblast", "stromal cell"])
            & (ob.disease == "normal")).values
    key = ("tendon, quadriceps | " + ob.cell_type.astype(str) + " | "
           + ob.donor_id.astype(str)).values
    return pseudobulk(p, mask, key, _out("tendon_quad_pseudobulk"), NNZ)


def build_kidney_susp():
    """Kidney interstitial fibroblasts, restricted to donors sampled both ways.

    Restricting to shared donors is what makes this a control rather than
    another confounded comparison: suspension type varies, donor does not.
    """
    p = _raw("kidney")
    ob = _obs(p, ["cell_type", "suspension_type", "assay", "donor_id", "disease"])
    fib = ob.cell_type == "kidney interstitial fibroblast"
    shared = set.intersection(*[set(ob.donor_id[fib & (ob.suspension_type == s)])
                                for s in ("cell", "nucleus")])
    mask = (fib & ob.donor_id.isin(shared)).values
    key = (ob.suspension_type.astype(str) + " | " + ob.disease.astype(str)).values
    print(f"  kidney: {len(shared)} donors sampled in both suspension types")
    return pseudobulk(p, mask, key, _out("kidney_susp_pseudobulk"), NNZ)


def build_heart_susp():
    """Cardiac fibroblasts, 3' v3, healthy; suspension x anatomical region."""
    p = _raw("heart")
    ob = _obs(p, ["cell_type", "suspension_type", "assay", "donor_id", "disease", "tissue"])
    mask = ((ob.cell_type == "fibroblast") & (ob.assay == "10x 3' v3")
            & (ob.disease == "normal")).values
    key = (ob.suspension_type.astype(str) + " | " + ob.tissue.astype(str)).values
    return pseudobulk(p, mask, key, _out("heart_susp_pseudobulk"), NNZ)


def build_heart_donor_susp():
    """Within-individual check: the one donor contributing cells and nuclei."""
    p = _raw("heart")
    ob = _obs(p, ["cell_type", "suspension_type", "assay", "donor_id", "disease"])
    fib = ((ob.cell_type == "fibroblast") & (ob.assay == "10x 3' v3")
           & (ob.disease == "normal"))
    shared = set.intersection(*[set(ob.donor_id[fib & (ob.suspension_type == s)])
                                for s in ("cell", "nucleus")])
    if not shared:
        print("  heart: no donor sampled both ways, skipping within-donor check")
        return {}
    d = sorted(shared)[0]
    mask = (fib & (ob.donor_id == d)).values
    key = (ob.suspension_type.astype(str) + " | donor").values
    print(f"  heart: within-donor check on {d}")
    return pseudobulk(p, mask, key, _out("heart_donor_susp_pseudobulk"), NNZ)


BUILDERS = {
    "ts_donor": build_ts_donor,
    "gut": build_gut,
    "skin": lambda: build_simple("skin", "skin_pseudobulk"),
    "oral": lambda: build_simple("oral", "oral_pseudobulk"),
    "oral_v3_donor": build_oral_v3_donor,
    "synovium": build_synovium,
    "tendon_ach": build_tendon_ach,
    "tendon_quad": build_tendon_quad,
    "kidney_susp": build_kidney_susp,
    "heart_susp": build_heart_susp,
    "heart_donor_susp": build_heart_donor_susp,
}

OUTPUTS = {
    "ts_donor": "ts_donor_pseudobulk", "gut": "gut_pseudobulk",
    "skin": "skin_pseudobulk", "oral": "oral_pseudobulk",
    "oral_v3_donor": "oral_v3_donor_pseudobulk", "synovium": "synovium_pseudobulk",
    "tendon_ach": "tendon_ach_pseudobulk", "tendon_quad": "tendon_quad_pseudobulk",
    "kidney_susp": "kidney_susp_pseudobulk", "heart_susp": "heart_susp_pseudobulk",
    "heart_donor_susp": "heart_donor_susp_pseudobulk",
}


def main(targets=None):
    force = bool(targets)
    targets = targets or list(BUILDERS)
    for t in targets:
        if t not in BUILDERS:
            raise SystemExit(f"unknown target {t!r}\nknown: {list(BUILDERS)}")
        if not force and _out(OUTPUTS[t]).exists():
            print(f"skip  {t}", flush=True)
            continue
        t0 = time.time()
        print(f"build {t}", flush=True)
        n = BUILDERS[t]()
        kept = sum(1 for v in n.values() if v >= 25) if isinstance(n, dict) else "?"
        print(f"  done {t}: {len(n) if hasattr(n, '__len__') else '?'} groups "
              f"({kept} with >=25 cells) in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or None)
