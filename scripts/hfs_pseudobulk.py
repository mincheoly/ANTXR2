"""Low-RAM streaming readers and pseudobulk aggregation for CELLxGENE h5ads.

Every routine here walks the CSR matrix in blocks sized by nonzero count rather
than by row count, so peak memory is set by ``target_nnz`` and not by the file.
This is what makes a 10 GB atlas processable on a 16 GB host without an
out-of-core framework.

Two passes are provided:

``extract``      per-cell counts for a small gene panel plus per-cell library
                 size — used for QC, detection rates and cell-level filters.
``pseudobulk``   full-gene summed counts per group — the input to every
                 CPM-based tissue metric in this arm.

Both read ``raw/X`` when present (CELLxGENE convention for raw counts) and fall
back to ``X`` otherwise.
"""
import h5py
import numpy as np

__all__ = ["read_obs_column", "var_names", "block_bounds", "extract", "pseudobulk"]


def read_obs_column(group, key):
    """Read an AnnData obs/var column (categorical or plain) as a str/num array."""
    o = group[key]
    if isinstance(o, h5py.Group):
        if not ("categories" in o and "codes" in o):
            return None
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                         for c in o["categories"][:]])
        codes = o["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
    v = o[:]
    if v.dtype.kind == "S":
        return np.array([x.decode() for x in v])
    return v.astype(str) if v.dtype.kind not in "fiu" else v


def var_names(f, root=""):
    """Gene symbols, preferring CELLxGENE's ``feature_name`` over the index."""
    vg = f[f"{root}/var"] if root else f["var"]
    if "feature_name" in vg:
        return read_obs_column(vg, "feature_name")
    idx = vg.attrs.get("_index", "_index")
    idx = idx.decode() if isinstance(idx, bytes) else idx
    return read_obs_column(vg, idx)


def block_bounds(indptr, start, n_rows, target_nnz):
    """Largest ``end`` such that rows [start, end) hold <= target_nnz nonzeros."""
    base = indptr[start]
    lo, hi = start + 1, n_rows
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if indptr[mid] - base <= target_nnz:
            lo = mid
        else:
            hi = mid - 1
    return max(lo, start + 1)


def _matrix_root(f):
    return "raw" if ("raw" in f and "X" in f["raw"]) else ""


def extract(path, out_npz, panel, obs_columns=None, target_nnz=4_000_000):
    """Stream per-cell panel counts + library size out of an h5ad.

    Writes ``counts`` (n_cells x n_panel_found), ``total`` (library size),
    ``ngene`` (genes detected) and ``obs_<col>`` arrays to ``out_npz``.
    """
    obs_columns = obs_columns or ["cell_type", "tissue", "tissue_general", "disease",
                                  "development_stage", "donor_id", "assay",
                                  "suspension_type", "sex", "is_primary_data"]
    f = h5py.File(path, "r")
    root = _matrix_root(f)
    X = f[f"{root}/X"] if root else f["X"]
    enc = X.attrs.get("encoding-type", "")
    enc = enc.decode() if isinstance(enc, bytes) else enc
    assert enc == "csr_matrix", f"expected csr_matrix, got {enc!r}"

    names = var_names(f, root)
    n_genes = int(X.attrs["shape"][1])
    indptr = X["indptr"][:]
    n_cells = len(indptr) - 1

    sel = {}
    for g in panel:
        hit = np.where(names == g)[0]
        if len(hit):
            sel[g] = int(hit[0])
    cols = np.array(list(sel.values()), dtype=np.int64)
    lut = np.full(n_genes, -1, dtype=np.int32)
    lut[cols] = np.arange(len(cols), dtype=np.int32)

    counts = np.zeros((n_cells, len(cols)), dtype=np.float32)
    total = np.zeros(n_cells, dtype=np.float64)
    ngene = np.zeros(n_cells, dtype=np.int32)

    s = 0
    while s < n_cells:
        e = block_bounds(indptr, s, n_cells, target_nnz)
        a, b = int(indptr[s]), int(indptr[e])
        dat = X["data"][a:b].astype(np.float64)
        ind = X["indices"][a:b]
        rows = np.repeat(np.arange(e - s), np.diff(indptr[s:e + 1]))
        np.add.at(total, np.arange(s, e), np.bincount(rows, weights=dat, minlength=e - s))
        ngene[s:e] = np.bincount(rows, minlength=e - s)
        m = lut[ind] >= 0
        if m.any():
            np.add.at(counts, (rows[m] + s, lut[ind[m]]), dat[m])
        s = e

    obs = f[f"{root}/obs"] if (root and "obs" in f[root]) else f["obs"]
    meta = {}
    for c in obs_columns:
        if c in obs:
            v = read_obs_column(obs, c)
            if v is not None:
                meta[c] = v
    f.close()

    np.savez_compressed(out_npz, counts=counts, total=total, ngene=ngene,
                        genes=np.array(list(sel)),
                        **{f"obs_{k}": v.astype(str) for k, v in meta.items()})
    return dict(n_cells=int(n_cells), genes=list(sel), root=root or "X")


def pseudobulk(path, mask, group_keys, out_npz, target_nnz=4_000_000, count_nnz=False):
    """Sum raw counts over all genes within each group of cells.

    ``mask``        boolean, which cells to include at all
    ``group_keys``  str array, one key per cell; groups are its sorted unique
                    values over the masked cells

    Writes ``sum`` (n_groups x n_genes, float32), ``ncell``, ``groups``,
    ``genes`` and optionally ``nnz`` (per-group detected-cell counts).
    """
    f = h5py.File(path, "r")
    root = _matrix_root(f)
    X = f[f"{root}/X"] if root else f["X"]
    n_genes = int(X.attrs["shape"][1])
    indptr = X["indptr"][:]
    n = len(indptr) - 1

    cats = np.array(sorted(set(np.asarray(group_keys)[mask])))
    cmap = {c: i for i, c in enumerate(cats)}
    code = np.full(n, -1, np.int32)
    idx = np.where(mask)[0]
    code[idx] = [cmap[x] for x in np.asarray(group_keys)[idx]]

    G = len(cats)
    SUM = np.zeros((G, n_genes))
    NNZ = np.zeros((G, n_genes), np.int32) if count_nnz else None
    NC = np.bincount(code[code >= 0], minlength=G)

    s = 0
    while s < n:
        e = block_bounds(indptr, s, n, target_nnz)
        a, b = int(indptr[s]), int(indptr[e])
        dat = X["data"][a:b].astype(np.float64)
        ind = X["indices"][a:b]
        rows = np.repeat(code[s:e], np.diff(indptr[s:e + 1]))
        m = rows >= 0
        if m.any():
            flat = rows[m].astype(np.int64) * n_genes + ind[m]
            SUM += np.bincount(flat, weights=dat[m], minlength=G * n_genes).reshape(G, n_genes)
            if count_nnz:
                NNZ += np.bincount(flat, minlength=G * n_genes).reshape(G, n_genes).astype(np.int32)
        s = e

    genes = var_names(f, root)
    f.close()

    payload = dict(sum=SUM.astype(np.float32), ncell=NC, groups=cats, genes=genes,
                   total=SUM.sum(1))
    if count_nnz:
        payload["nnz"] = NNZ
    np.savez_compressed(out_npz, **payload)
    return {str(c): int(NC[i]) for i, c in enumerate(cats)}
