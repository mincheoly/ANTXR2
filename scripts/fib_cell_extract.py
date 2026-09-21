"""Stream fibroblast cells x ALL genes out of the fig5 atlas h5ads.

One pass per source file. Selected rows are appended to per-population scratch
binaries (so peak RAM is one nnz-block, not the whole subset), then assembled
into a CSR .npz per population:

    data, indices, indptr, shape, genes, donor, cell_type, total_counts

`total_counts` is the FULL-transcriptome per-cell count sum, carried so that
memento's size factors can be checked against the unsubsetted library size.

Population definitions are copied verbatim from the fig5 pseudobulk builders
(lineage_src/*_pseudobulk.npz.py) so the cell sets match the published panel.
"""
import os
import numpy as np
import h5py

DATA = r"C:\Data\ANTXR2_workspace\data"
OUT = r"C:\Data\ANTXR2_workspace\fib_cells"

FIBL2 = {'fibroblast', 'alveolar adventitial fibroblast', 'fibroblast of cardiac tissue',
         'fibroblast of breast', 'thymic fibroblast type 1', 'thymic fibroblast type 2'}
ORAL_SITES = ["gingiva", "buccal mucosa", "hard palate"]


def rc(g, c):
    o = g[c]
    if isinstance(o, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) if v.dtype.kind in "SO" else v


def genenames(f):
    v = f["var/feature_name"]
    if isinstance(v, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in v["categories"][:]])
        return cats[v["codes"][:]]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v[:]])


def extract(path, mask, groups, obs_cols, tag, target=4_000_000, out_dir=None):
    """mask: bool per cell; groups: str per cell; obs_cols: dict name -> array."""
    OUT = out_dir or globals()["OUT"]
    os.makedirs(OUT, exist_ok=True)
    scratch = os.path.join(OUT, "_scratch")
    os.makedirs(scratch, exist_ok=True)
    f = h5py.File(path, "r")
    X = f["raw/X" if ("raw" in f and "X" in f["raw"]) else "X"]
    ng = int(X.attrs["shape"][1])
    ip = X["indptr"][:]
    n = len(ip) - 1
    assert len(mask) == n, (len(mask), n)

    cats = sorted(set(np.asarray(groups)[mask]))
    cmap = {c: i for i, c in enumerate(cats)}
    code = np.full(n, -1, np.int32)
    idx = np.where(mask)[0]
    code[idx] = [cmap[g] for g in np.asarray(groups)[idx]]

    safe = [c.replace(",", "").replace(" ", "_").replace("(", "").replace(")", "") for c in cats]
    fh = [(open(os.path.join(scratch, f"{tag}_{s}.dat"), "wb"),
           open(os.path.join(scratch, f"{tag}_{s}.idx"), "wb")) for s in safe]
    rownnz = [[] for _ in cats]
    total = np.zeros(n, np.float64)

    s = 0
    while s < n:
        base = ip[s]
        lo, hi = s + 1, n
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if ip[mid] - base <= target:
                lo = mid
            else:
                hi = mid - 1
        e = max(lo, s + 1)
        a, b = int(ip[s]), int(ip[e])
        dat = X["data"][a:b]
        ind = X["indices"][a:b]
        lens = np.diff(ip[s:e + 1])
        rr = np.repeat(np.arange(s, e), lens)
        np.add.at(total, rr, dat.astype(np.float64))
        gnnz = code[rr]
        m = gnnz >= 0
        if m.any():
            gk = gnnz[m]
            order = np.argsort(gk, kind="stable")   # stable -> row order preserved
            gs = gk[order]
            dsel = dat[m][order].astype(np.float32)
            isel = ind[m][order].astype(np.int32)
            bounds = np.searchsorted(gs, np.arange(len(cats) + 1))
            for gi in range(len(cats)):
                lo_, hi_ = bounds[gi], bounds[gi + 1]
                if hi_ > lo_:
                    fh[gi][0].write(dsel[lo_:hi_].tobytes())
                    fh[gi][1].write(isel[lo_:hi_].tobytes())
        blkcode = code[s:e]
        for gi in np.unique(blkcode[blkcode >= 0]):
            rows = np.where(blkcode == gi)[0]
            rownnz[gi].extend(lens[rows].tolist())
        s = e

    gn = genenames(f)
    f.close()
    for d, i in fh:
        d.close()
        i.close()

    written = []
    for gi, c in enumerate(cats):
        rows = np.where(code == gi)[0]
        nz = np.array(rownnz[gi], np.int64)
        assert len(nz) == len(rows), (c, len(nz), len(rows))
        indptr = np.concatenate([[0], np.cumsum(nz)])
        d = np.fromfile(os.path.join(scratch, f"{tag}_{safe[gi]}.dat"), np.float32)
        ix = np.fromfile(os.path.join(scratch, f"{tag}_{safe[gi]}.idx"), np.int32)
        assert len(d) == indptr[-1] == len(ix), (c, len(d), indptr[-1], len(ix))
        out = os.path.join(OUT, f"{tag}_{safe[gi]}.npz")
        np.savez_compressed(out, data=d, indices=ix, indptr=indptr,
                            shape=np.array([len(rows), ng]), genes=gn,
                            total_counts=total[rows],
                            **{k: np.asarray(v)[rows].astype(str) for k, v in obs_cols.items()})
        os.remove(os.path.join(scratch, f"{tag}_{safe[gi]}.dat"))
        os.remove(os.path.join(scratch, f"{tag}_{safe[gi]}.idx"))
        written.append((c, len(rows), int(indptr[-1]), out))
    return written
