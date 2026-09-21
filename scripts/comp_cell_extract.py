"""Stream epithelial / endothelial / immune cells x ALL genes out of the fig5
atlases, one CSR .npz per (tissue, compartment).

Cell masks reuse the fig5 population definitions exactly as the fibroblast
extraction did -- only the cell_type set changes:

    oral         10x 3' v3, sampled_site_condition=healthy, the three sites
    synovium     tissue = layer of synovial tissue (blood and synovial fluid
                 excluded, so immune cells are tissue-resident)
    tendon_ach   all cells (the file is entirely disease=normal)
    tendon_quad  disease = normal
    TS comp.     10x 3' v3, grouped by tissue_in_publication

Every cell of the compartment is kept and its cell_type carried in obs, so the
grouping decision (which subtypes become memento groups) is made at run time by
antxr2_coexpr.py rather than baked into the cache.
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fib_cell_extract as F                       # noqa: E402
from compartment_map import LOOKUP                 # noqa: E402

OUT = r"C:\Data\ANTXR2_workspace\comp_cells"
COMPARTMENTS = ["epithelial", "endothelial", "immune"]
MIN_CELLS = 200                                    # per (tissue, compartment)
CAP = 25_000                                       # cells per (tissue, compartment)
SEED = 5

# Why a cap: TS immune tissues run to >100k cells, whose CSR blocks are ~10x the
# largest fibroblast population (ts_Fat, 25.7k cells / 102M nnz) and would both
# dominate RAM and skew power across compartments. Subsampling is stratified by
# donor x cell_type so composition is preserved, and cells enter memento's
# moment estimates sublinearly, so the precision cost is modest.


def cap_mask(m, donor, ct, cap=CAP, seed=SEED):
    """Stratified down-sample of a population mask to `cap` cells."""
    idx = np.where(m)[0]
    if len(idx) <= cap:
        return m
    rng = np.random.default_rng(seed)
    key = pd.Series([f"{d}|{c}" for d, c in zip(donor[idx], ct[idx])])
    frac = cap / len(idx)
    keep = []
    for _, grp in key.groupby(key, sort=True):
        g = idx[grp.index.values]
        n = len(g)
        take = min(n, max(int(round(n * frac)), min(n, 100)))
        keep.append(rng.choice(g, size=take, replace=False))
    out = np.zeros_like(m)
    out[np.concatenate(keep)] = True
    return out


def _obs(path, cols):
    with h5py.File(path, "r") as f:
        have = [c for c in cols if c in f["obs"]]
        return pd.DataFrame({c: F.rc(f["obs"], c) for c in have})


def source_masks(tag):
    """-> (path, base_mask, tissue_label_per_cell, obs DataFrame)"""
    path = os.path.join(F.DATA, {"oral": "oral.h5ad", "synovium": "synovium.h5ad",
                                 "tendon_ach": "tendon_ach.h5ad",
                                 "tendon_quad": "tendon_quad.h5ad",
                                 "ts_epithelium": "ts_epithelium.h5ad",
                                 "ts_endothelium": "ts_endothelium.h5ad",
                                 "ts_immune": "ts_immune.h5ad"}[tag])
    if tag == "oral":
        ob = _obs(path, ["assay", "tissue", "cell_type", "donor_id", "sampled_site_condition"])
        m = ((ob.assay == "10x 3' v3") & ob.tissue.isin(F.ORAL_SITES) &
             (ob.sampled_site_condition == "healthy")).values
        lab = np.array(["oral_" + t.replace(" ", "_") for t in ob.tissue])
    elif tag == "synovium":
        ob = _obs(path, ["cell_type", "donor_id", "tissue"])
        m = (ob.tissue == "layer of synovial tissue").values
        lab = np.array(["syn_synovium_JIA"] * len(ob))
    elif tag == "tendon_ach":
        ob = _obs(path, ["cell_type", "donor_id", "disease"])
        m = np.ones(len(ob), bool)
        lab = np.array(["tach_tendon_Achilles"] * len(ob))
    elif tag == "tendon_quad":
        ob = _obs(path, ["cell_type", "donor_id", "disease"])
        m = (ob.disease == "normal").values
        lab = np.array(["tquad_tendon_quadriceps"] * len(ob))
    else:
        ob = _obs(path, ["assay", "tissue_in_publication", "cell_type", "donor_id"])
        m = (ob.assay == "10x 3' v3").values
        lab = np.array(["ts_" + t.replace(" ", "_") for t in ob.tissue_in_publication])
    return path, m, lab, ob


def compartment_of_labels(labels):
    """Unmapped labels are reported, not silently dropped."""
    out = np.array([LOOKUP.get(x, "UNMAPPED") for x in labels])
    return out


def run(tags):
    os.makedirs(OUT, exist_ok=True)
    report = []
    for tag in tags:
        path, base, lab, ob = source_masks(tag)
        cts = ob.cell_type.values
        if tag.startswith("ts_"):
            comp = np.array([tag.split("_", 1)[1].replace("epithelium", "epithelial")
                             .replace("endothelium", "endothelial")] * len(ob))
        else:
            comp = compartment_of_labels(cts)
            unmapped = sorted(set(cts[comp == "UNMAPPED"]))
            if unmapped:
                raise KeyError(f"{tag}: unmapped cell_type {unmapped}")
        for c in COMPARTMENTS:
            m = base & (comp == c)
            if m.sum() < MIN_CELLS:
                continue
            keep = pd.Series(lab[m]).value_counts()
            keep = set(keep[keep >= MIN_CELLS].index)
            if not keep:
                continue
            m = m & np.isin(lab, list(keep))
            for t in sorted(keep):                 # cap each tissue separately
                mt = m & (lab == t)
                capped = cap_mask(mt, ob.donor_id.values, cts)
                if capped.sum() < mt.sum():
                    print(f"  cap {c}/{t}: {mt.sum()} -> {capped.sum()} cells", flush=True)
                m = m & ~(mt & ~capped)
            written = F.extract(path, m, lab, {"donor": ob.donor_id.values,
                                               "cell_type": cts}, f"{c}", out_dir=OUT)
            for pop, nc, nnz, fn in written:
                report.append(dict(source=tag, compartment=c, population=pop,
                                   n_cells=nc, nnz=nnz, file=fn))
                print(f"{tag:14s} {c:12s} {pop:28s} cells={nc:7d} nnz={nnz:12d}", flush=True)
    pd.DataFrame(report).to_csv(os.path.join(OUT, "extract_report.csv"), index=False)
    return report


if __name__ == "__main__":
    run(sys.argv[1:] or ["oral", "synovium", "tendon_ach", "tendon_quad",
                         "ts_endothelium", "ts_epithelium", "ts_immune"])
