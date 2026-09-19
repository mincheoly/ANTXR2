"""Stream panel-gene counts + per-cell library size out of a CELLxGENE h5ad
without ever loading the full matrix (low-RAM host)."""
import h5py, numpy as np, json, sys

PANEL = ["ANTXR2", "ANTXR1",
         "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
         "COL4A1", "COL4A2", "COL1A1", "COL3A1",
         "LAMA4", "LAMB1", "LRP6", "MMP2", "MMP14", "HSPG2", "PDGFRB"]


def read_cat(g, key):
    """Read an obs column (categorical or plain) -> (values_as_str_array)."""
    o = g[key]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = o["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        codes = o["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
        return out
    v = o[:]
    if v.dtype.kind == "S":
        v = np.array([x.decode() for x in v])
    return v.astype(str) if v.dtype.kind not in "fiu" else v


def var_names(f, root):
    vg = f[f"{root}/var"]
    if "feature_name" in vg:
        return read_cat(vg, "feature_name")
    idx = vg.attrs.get("_index", "_index")
    idx = idx.decode() if isinstance(idx, bytes) else idx
    return read_cat(vg, idx)


def extract(path, out_npz, target_nnz=4_000_000):
    f = h5py.File(path, "r")
    root = "raw" if "raw" in f and "X" in f["raw"] else ""
    X = f[f"{root}/X"] if root else f["X"]
    enc = X.attrs.get("encoding-type", "")
    if isinstance(enc, bytes):
        enc = enc.decode()
    assert enc == "csr_matrix", f"expected csr, got {enc}"
    names = var_names(f, root) if root else var_names(f, "")
    n_genes = X.attrs["shape"][1]
    indptr = X["indptr"][:]
    n_cells = len(indptr) - 1

    sel = {}
    for gname in PANEL:
        hit = np.where(names == gname)[0]
        if len(hit):
            sel[gname] = int(hit[0])
    cols = np.array([sel[g] for g in sel], dtype=np.int64)
    gene_order = list(sel)
    lut = np.full(n_genes, -1, dtype=np.int32)
    lut[cols] = np.arange(len(cols), dtype=np.int32)

    counts = np.zeros((n_cells, len(cols)), dtype=np.float32)
    total = np.zeros(n_cells, dtype=np.float64)
    ngene = np.zeros(n_cells, dtype=np.int32)

    start = 0
    while start < n_cells:
        end = start + 1
        base = indptr[start]
        # grow block until ~target_nnz
        lo, hi = start + 1, n_cells
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if indptr[mid] - base <= target_nnz:
                lo = mid
            else:
                hi = mid - 1
        end = max(lo, start + 1)
        a, b = int(indptr[start]), int(indptr[end])
        dat = X["data"][a:b].astype(np.float64)
        ind = X["indices"][a:b]
        rows = np.repeat(np.arange(end - start), np.diff(indptr[start:end + 1]))
        np.add.at(total, np.arange(start, end), np.bincount(rows, weights=dat, minlength=end - start))
        ngene[start:end] = np.bincount(rows, minlength=end - start)
        m = lut[ind] >= 0
        if m.any():
            np.add.at(counts, (rows[m] + start, lut[ind[m]]), dat[m])
        start = end

    obs = f[f"{root}/obs"] if (root and "obs" in f[root]) else f["obs"]
    meta = {}
    for c in ["cell_type", "tissue", "tissue_general", "disease", "development_stage",
              "donor_id", "assay", "suspension_type", "sex", "is_primary_data"]:
        if c in obs:
            meta[c] = read_cat(obs, c)
    f.close()
    np.savez_compressed(out_npz, counts=counts, total=total, ngene=ngene,
                        genes=np.array(gene_order), **{f"obs_{k}": v.astype(str) for k, v in meta.items()})
    return dict(n_cells=int(n_cells), genes=gene_order, root=root or "X")


if __name__ == "__main__":
    r = extract(sys.argv[1], sys.argv[2])
    print(json.dumps(r))
