"""Figure 1 (Act 1) inputs: per-cell counts for a small gene panel plus library
size, streamed out of each atlas h5ad without loading it (block reads of the raw
CSR, bounded by nonzero count). Replaces the genome-wide memento cell-type means
of compute_means.py, whose outputs (/data/ANTXR2/celltype_expression) are gone:
Figure 1 needs only a handful of genes, and per-cell values keep the donor and
annotation columns so any grouping can be redone downstream (fig1_aggregate.py).

    -> /data/antxr2/fig1/cells/<atlas>.parquet   one row per cell:
       categorical obs columns, total (raw UMI sum over all genes), n_genes,
       one count column per panel gene present in the file
       /data/antxr2/fig1/cells/<atlas>.meta.json  file identity, matrix slot, genes missing

Usage: python fig1_extract.py [atlas ...]     (default: every atlas whose file exists)
"""
import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

RAW = Path("/data/antxr2/raw")
OUT = Path("/data/antxr2/fig1/cells")

# atlas -> (path, expected CELLxGENE dataset_version_id or None)
ATLASES = {
    "tabula_sapiens": (RAW / "tabula_sapiens/1c88f927-bcbb-4bb1-9881-281842945a2d.h5ad",
                       "1c88f927-bcbb-4bb1-9881-281842945a2d"),
    "cross_tissue_immune": (RAW / "cross_tissue_immune/09a381ce-feaf-48c8-a928-58f6768f4265.h5ad",
                            "09a381ce-feaf-48c8-a928-58f6768f4265"),
    "gca_healthy": (RAW / "gut_cell_atlas/healthy_reference.h5ad", "b801706d-6345-4af2-b3fe-8e30e0e6150b"),
    "gca_extended": (RAW / "gut_cell_atlas/extended_plus.h5ad", None),
    "skin_fibroblast": (RAW / "skin/dermal_fibroblast_atlas.h5ad", None),
    "ibdverse": (RAW / "ibdverse/ibdverse.h5ad", None),
}

PANEL = (
    # receptors and the HFS clearance arm with its substrate
    ["ANTXR1", "ANTXR2", "MRC2", "CTSB", "CTSK", "MMP14", "TIMP2", "LAMP1",
     "COL6A1", "COL6A2", "COL6A3"]
    # fibrillar collagens (the Act 2 pertinent negative)
    + ["COL1A1", "COL1A2", "COL3A1"]
    # Wnt receptor arm
    + ["LRP5", "LRP6"] + [f"FZD{i}" for i in range(1, 11)] + ["CTNNB1", "TCF7L2", "LGR5", "RNF43", "ZNRF3", "AXIN2"]
    # compartment / state markers
    + ["EPCAM", "PTPRC", "PECAM1", "PDGFRA", "PDGFRB", "DCN", "LUM", "VIM", "ACTA2", "MKI67", "OLFM4",
       "FABP6", "APOA4", "KRT20", "CEACAM5", "MUC2", "KRT14", "KRT5"]
    # type I interferon response (Act 3 bridge)
    + ["MX1", "ISG15", "IFIT1", "IFIT3", "IFI6", "IFI44L", "STAT1", "IRF7", "IFNAR1", "IFNAR2"]
    + ["ACTB", "GAPDH"]
)
TARGET_NNZ = 8_000_000
MAX_CATS = 20_000          # obs columns with more categories than this are per-cell ids, dropped


def read_col(g, key):
    o = g[key]
    if isinstance(o, h5py.Group):
        if not ("categories" in o and "codes" in o):
            return None
        cats = o["categories"][:]
        if len(cats) > MAX_CATS:
            return None
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        codes = o["codes"][:]
        return pd.Categorical.from_codes(codes, cats)  # code -1 -> NaN
    if o.ndim != 1:
        return None
    v = o[:]
    if v.dtype.kind in "SO":
        v = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])
        u = pd.unique(v)
        return pd.Categorical(v) if len(u) <= MAX_CATS else None
    return v


def file_version_id(f):
    if "uns" not in f or "citation" not in f["uns"]:
        return None
    v = f["uns"]["citation"][()]
    ids = re.findall(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
                     v.decode() if isinstance(v, bytes) else str(v))
    return ids[0] if ids else None


def gene_names(f, root):
    vg = f[f"{root}/var"] if root else f["var"]
    for k in ("feature_name", "gene_symbols", "gene_name"):
        if k in vg:
            v = vg[k]
            if isinstance(v, h5py.Group):
                cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in v["categories"][:]])
                return cats[v["codes"][:]]
            return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v[:]])
    idx = vg.attrs["_index"]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in vg[idx][:]])


def extract(name):
    path, expect = ATLASES[name]
    t0 = time.time()
    f = h5py.File(path, "r")
    vid = file_version_id(f)
    if expect and vid != expect:
        raise ValueError(f"{name}: file declares {vid}, expected {expect}")
    root = "raw" if ("raw" in f and "X" in f["raw"]) else ""
    X = f[f"{root}/X"] if root else f["X"]
    genes = gene_names(f, root)
    n_genes = int(X.attrs["shape"][1])
    assert len(genes) == n_genes, (len(genes), n_genes)
    sel = {g: int(np.where(genes == g)[0][0]) for g in PANEL if (genes == g).any()}
    missing = [g for g in PANEL if g not in sel]
    lut = np.full(n_genes, -1, np.int32)
    lut[list(sel.values())] = np.arange(len(sel))
    indptr = X["indptr"][:]
    n = len(indptr) - 1
    counts = np.zeros((n, len(sel)), np.float32)
    total = np.zeros(n)
    ngene = np.zeros(n, np.int32)
    nonint = 0.0
    s = 0
    while s < n:
        # largest e with nnz(s:e) <= TARGET_NNZ
        e = int(np.searchsorted(indptr, indptr[s] + TARGET_NNZ, side="right")) - 1
        e = min(max(e, s + 1), n)
        a, b = int(indptr[s]), int(indptr[e])
        dat = X["data"][a:b].astype(np.float64)
        ind = X["indices"][a:b]
        if s == 0:
            nonint = float(np.mean(dat[:200_000] != np.round(dat[:200_000]))) if len(dat) else 0.0
        rows = np.repeat(np.arange(e - s), np.diff(indptr[s:e + 1]))
        total[s:e] = np.bincount(rows, weights=dat, minlength=e - s)
        ngene[s:e] = np.bincount(rows, minlength=e - s)
        m = lut[ind] >= 0
        if m.any():
            np.add.at(counts, (rows[m] + s, lut[ind[m]]), dat[m])
        s = e
    if nonint > 0.01:
        raise ValueError(f"{name}: {root or ''}/X is not integer counts ({nonint:.1%} non-integer)")
    obs = {}
    for k in f["obs"].keys():
        if k.startswith("_"):
            continue
        v = read_col(f["obs"], k)
        if v is not None and len(v) == n:
            obs[k] = v
    f.close()
    df = pd.DataFrame(obs)
    df["total"] = total.astype(np.float32)
    df["n_genes"] = ngene
    for g, j in zip(sel, range(len(sel))):
        df[g] = counts[:, j]
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT / f"{name}.parquet", index=False)
    meta = dict(atlas=name, path=str(path), dataset_version_id=vid, matrix=f"{root}/X" if root else "X",
                n_cells=n, genes=list(sel), missing=missing, nonint_frac=nonint,
                seconds=round(time.time() - t0))
    (OUT / f"{name}.meta.json").write_text(json.dumps(meta, indent=1))
    print(f"{name}: {n} cells, {len(sel)} genes ({missing} missing), {meta['seconds']}s", flush=True)


if __name__ == "__main__":
    todo = sys.argv[1:] or [k for k, (p, _) in ATLASES.items() if p.exists()]
    for a in todo:
        extract(a)
