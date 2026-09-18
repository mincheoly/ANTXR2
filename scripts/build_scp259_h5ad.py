"""Build per-compartment AnnData objects from the SCP259 (Smillie 2019) release.

SCP259 ships three separate gene-sorted MatrixMarket files (Epi, Fib, Imm) plus
one shared cell-metadata table. This script parses each into a cells x genes
CSR AnnData with raw integer counts and the metadata joined on barcode.

Provenance note
---------------
Files come from the Broad Single Cell Portal's own signed URLs, not from a
mirror. The HCA DCP copy of this study was found to serve TRUNCATED matrices
for the Epi and Fib compartments; this script therefore asserts that each
matrix's declared dimensions match its genes/barcodes files and refuses to
build on a mismatch. A completed download is not evidence of correct content.

Capture rate
------------
`q` is set to config.DEFAULT_CAPTURE_RATE (0.1). SCP259 is 10x 3' v2 and the
release ships no per-library sequencing-saturation metric, so the empirical
route in `capture_rate.py` cannot be applied here. This is an assumption, and
memento's absolute mean estimates depend on it; correlation estimates are far
less sensitive to it.
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from config import DEFAULT_CAPTURE_RATE

COMPARTMENT_DIR = {
    "Epi": "5cdc540d328cee7a2efc2348",
    "Fib": "5cdc540d328cee7a2efc2349",
    "Imm": "5cdc540d328cee7a2efc234a",
}
CHUNK_ROWS = 20_000_000


def read_lines(path):
    with open(path, encoding="utf-8") as fh:
        return [l.rstrip("\n") for l in fh if l.strip()]


def mtx_header(path):
    """Return (n_rows, n_cols, nnz, n_header_lines) for a MatrixMarket file.

    n_header_lines counts the '%' comment lines plus the dimension line, i.e.
    exactly what `skiprows` must drop before the coordinate entries begin.
    """
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            n += 1
            if line.startswith("%"):
                continue
            r, c, nnz = (int(x) for x in line.split())
            return r, c, nnz, n
    raise ValueError(f"no header in {path}")


def load_compartment(base, comp, meta):
    d = Path(base) / "expression" / COMPARTMENT_DIR[comp]
    mtx = d / f"gene_sorted-{comp}.matrix.mtx"
    genes = read_lines(d / f"{comp}.genes.tsv")
    bcs = read_lines(d / f"{comp}.barcodes2.tsv")
    n_gene, n_cell, nnz, n_header = mtx_header(mtx)
    if (n_gene, n_cell) != (len(genes), len(bcs)):
        raise SystemExit(
            f"{comp}: matrix is {n_gene}x{n_cell} but genes.tsv has "
            f"{len(genes)} and barcodes has {len(bcs)} -- truncated or mismatched file")

    # rows are genes (file is gene-sorted); accumulate then transpose to cells x genes
    rows = np.empty(nnz, dtype=np.int32)
    cols = np.empty(nnz, dtype=np.int32)
    vals = np.empty(nnz, dtype=np.float32)
    pos = 0
    reader = pd.read_csv(mtx, sep=r"\s+", skiprows=n_header, header=None,
                         names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32,
                                                       "v": np.float32},
                         chunksize=CHUNK_ROWS)
    for ch in reader:
        k = len(ch)
        rows[pos:pos + k] = ch["g"].to_numpy() - 1
        cols[pos:pos + k] = ch["c"].to_numpy() - 1
        vals[pos:pos + k] = ch["v"].to_numpy()
        pos += k
        print(f"  {comp}: {pos:,}/{nnz:,} entries", flush=True)
    if pos != nnz:
        raise SystemExit(f"{comp}: header declares {nnz} nonzeros, file holds {pos}")

    X = sp.coo_matrix((vals, (cols, rows)), shape=(n_cell, n_gene)).tocsr()
    del rows, cols, vals

    obs = pd.DataFrame(index=pd.Index(bcs, name="NAME"))
    obs = obs.join(meta.set_index("NAME"))
    missing = int(obs["Cluster"].isna().sum())
    if missing:
        print(f"  {comp}: {missing} barcodes absent from metadata -- dropped")
        keep = obs["Cluster"].notna().to_numpy()
        X, obs = X[keep], obs.loc[keep]
    obs["q"] = DEFAULT_CAPTURE_RATE
    obs["compartment"] = comp
    A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    A.uns["source"] = "SCP259 Broad Single Cell Portal (Smillie et al. Cell 2019)"
    A.uns["capture_rate_q"] = DEFAULT_CAPTURE_RATE
    return A


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default="data/scp259/SCP259")
    ap.add_argument("--out-dir", default="data/scp259")
    ap.add_argument("--compartments", nargs="+", default=["Fib", "Epi", "Imm"])
    args = ap.parse_args()

    meta = pd.read_csv(Path(args.base) / "metadata" / "all.meta2.txt",
                       sep="\t", low_memory=False)
    if str(meta.iloc[0, 0]).upper() == "TYPE":
        meta = meta.drop(index=0)
    for col in ["nGene", "nUMI"]:
        meta[col] = pd.to_numeric(meta[col], errors="coerce")

    out = Path(args.out_dir)
    for comp in args.compartments:
        print(f"{comp} ...", flush=True)
        A = load_compartment(args.base, comp, meta)
        p = out / f"scp259_{comp}.h5ad"
        A.write_h5ad(p, compression="gzip")
        print(f"  -> {p}  {A.shape[0]} cells x {A.shape[1]} genes  "
              f"states={sorted(A.obs.Health.unique())}", flush=True)
        del A


if __name__ == "__main__":
    main()
