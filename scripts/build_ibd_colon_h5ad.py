"""Convert the IBD colon atlas's raw per-compartment data into 3 compute-ready
h5ad files (one per compartment: Epi/Fib/Imm), so compute_ibd_colon_means.py
can reuse compute_means.py's existing h5py-direct-read machinery unchanged.

Full 30-donor cohort, all 3 compartments, as of 2026-09-09: raw counts come
from `gene_sorted-{comp}.matrix.mtx` and visualization coordinates from
`{comp}.tsne.txt` (attached as `tsne_1`/`tsne_2` obs columns), both
SCP259-direct downloads via a user-provided authenticated bulk-download link,
each verified complete against its own mtx header / expected row count before
use (declared nnz == actual data line count for all 3 matrices; `Epi.tsne.txt`
row count == Epi's full cell count) -- see ANALYSIS_SUMMARY.md.

Supersedes an earlier version of this pipeline that had to fall back to a
17-of-30-donor discovery-cohort Seurat-object export for Epi/Imm, because the
HCA Data Coordination Platform's mirror of those two files turned out to be
truncated at the source (also in ANALYSIS_SUMMARY.md, kept as a standing
lesson). That fallback's outputs (`export_ibd_colon_rds.R`,
`raw/ibd_colon_atlas/rds_export/`, `raw/ibd_colon_atlas/train.*.seur.rds`) are
left in place for provenance but are no longer read by this script.

Usage:
    python build_ibd_colon_h5ad.py --dry-run   # inspect + verify only, no h5ad written
    python build_ibd_colon_h5ad.py             # full run
"""
import argparse
import gc
import logging
import os
import time

import anndata
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from config import (
    IBD_COLON_CLUSTER_COL, IBD_COLON_COMPARTMENTS, IBD_COLON_DIR,
    IBD_COLON_H5AD_BY_COMPARTMENT, IBD_COLON_META_FILE, IBD_COLON_META_NAME_COL,
    IBD_COLON_SUBSETS_FILE,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("build_ibd_colon_h5ad")

CELL_TYPE_LINEAGE_COL = "cell_type_lineage"
DATA_SOURCE_COL = "data_source"


def load_cell_subsets():
    """Fine (51-label) -> coarse lineage map. Tab-separated, no header."""
    df = pd.read_csv(
        os.path.join(IBD_COLON_DIR, IBD_COLON_SUBSETS_FILE),
        sep="\t", header=None, names=["fine", "lineage"],
    )
    mapping = dict(zip(df["fine"], df["lineage"]))
    log.info(f"cell_subsets.txt: {len(mapping)} fine->lineage mappings, "
              f"{df['lineage'].nunique()} distinct lineages: {sorted(df['lineage'].unique())}")
    return mapping


def load_meta():
    """all.meta2.txt: tab-separated, row 0 = header, row 1 = an SCP/Seurat 'TYPE'
    metadata-type row (group/numeric), not real data -- must be skipped."""
    meta = pd.read_csv(
        os.path.join(IBD_COLON_DIR, IBD_COLON_META_FILE), sep="\t", skiprows=[1],
    )
    meta = meta.set_index(IBD_COLON_META_NAME_COL)
    if meta.index.duplicated().any():
        n_dup = meta.index.duplicated().sum()
        raise RuntimeError(f"all.meta2.txt has {n_dup} duplicate NAME/barcode values -- "
                            f"cannot join 1:1 onto per-compartment barcodes")
    log.info(f"all.meta2.txt: {meta.shape[0]} rows, columns={list(meta.columns)}")
    return meta


def verify_mtx_complete(matrix_path, tag):
    """Required check after the earlier HCA-corruption episode: a matching
    file size / HTTP success is NOT sufficient evidence a text mtx file is
    complete -- compare its own header's declared nnz against the actual
    number of data lines before trusting it."""
    with open(matrix_path) as f:
        f.readline()  # banner
        declared_nnz = int(f.readline().split()[2])
    with open(matrix_path, "rb") as f:
        actual_lines = sum(1 for _ in f) - 2  # minus banner + dims line
    if actual_lines != declared_nnz:
        raise RuntimeError(
            f"[{tag}] {matrix_path} is INCOMPLETE: header declares {declared_nnz} data "
            f"lines, file actually has {actual_lines} ({actual_lines/declared_nnz:.1%}) -- "
            f"do not proceed (see the 2026-09-09 HCA-corruption lesson in ANALYSIS_SUMMARY.md)"
        )
    log.info(f"  [{tag}] verified complete: {actual_lines} data lines == declared nnz")


def read_mtx_triplet(matrix_path, genes_path, barcodes_path, tag):
    genes = pd.read_csv(genes_path, header=None)[0].astype(str)
    barcodes = pd.read_csv(barcodes_path, header=None)[0].astype(str)
    log.info(f"  [{tag}] {len(genes)} genes, {len(barcodes)} barcodes")

    with open(matrix_path) as f:
        f.readline()  # %%MatrixMarket banner
        dims = f.readline().split()
    n_rows, n_cols = int(dims[0]), int(dims[1])
    if n_rows != len(genes) or n_cols != len(barcodes):
        raise RuntimeError(
            f"[{tag}] mtx dims ({n_rows} x {n_cols}) don't match "
            f"genes ({len(genes)}) x barcodes ({len(barcodes)}) -- orientation assumption wrong"
        )
    verify_mtx_complete(matrix_path, tag)

    t0 = time.time()
    mat = scipy.io.mmread(matrix_path)  # genes x cells, per mtx header
    log.info(f"  [{tag}] mmread done in {time.time()-t0:.0f}s, dtype={mat.dtype}, "
              f"row.dtype={mat.row.dtype}, nnz={mat.nnz}")
    # Downcast before the COO->CSR conversion below, which otherwise briefly
    # holds both the old COO arrays and new CSR arrays at full int64 width --
    # for a ~174M-nonzero matrix that transpose+conversion step was enough to
    # get OOM-killed on this machine at int64 width. Values here are small
    # raw UMI counts (max seen ~24,563) and indices are well under 2^31 (max
    # ~210K cells / 20K genes), so int32 has enormous headroom on both.
    mat.data = mat.data.astype(np.int32, copy=False)
    mat.row = mat.row.astype(np.int32, copy=False)
    mat.col = mat.col.astype(np.int32, copy=False)
    mat = sp.csr_matrix(mat.T)  # -> cells x genes, anndata convention
    gc.collect()
    max_val = mat.data.max()
    sample_n = min(5_000_000, len(mat.data))
    frac_integer = np.mean(np.isclose(mat.data[:sample_n], np.round(mat.data[:sample_n])))
    log.info(f"  [{tag}] post-transpose shape={mat.shape}, max={max_val:.1f}, "
              f"frac_integer(sample)={frac_integer:.4f}")
    return mat, genes, barcodes


def load_tsne(comp, expected_barcodes):
    path = os.path.join(IBD_COLON_DIR, f"{comp}.tsne.txt")
    tsne = pd.read_csv(path, sep="\t", header=None, names=["barcode", "tsne_1", "tsne_2"], index_col="barcode")
    n_expected = len(expected_barcodes)
    if len(tsne) != n_expected or set(tsne.index) != set(expected_barcodes):
        log.warning(f"  [{comp}] tsne.txt has {len(tsne)} rows vs {n_expected} cells in the "
                     f"matrix -- barcode sets {'differ' if set(tsne.index) != set(expected_barcodes) else 'match'}")
    else:
        log.info(f"  [{comp}] tsne.txt: {len(tsne)} rows, full barcode match")
    return tsne


def build_compartment_adata(comp, meta, lineage_map, dry_run=False):
    log.info(f"=== {comp} (full cohort, SCP259 direct) ===")
    matrix_path = os.path.join(IBD_COLON_DIR, f"gene_sorted-{comp}.matrix.mtx")
    if dry_run:
        with open(matrix_path) as f:
            f.readline(); log.info(f"  header: {f.readline().strip()}")
        verify_mtx_complete(matrix_path, comp)
        return None

    mat, genes, barcodes = read_mtx_triplet(
        matrix_path,
        os.path.join(IBD_COLON_DIR, f"{comp}.genes.tsv"),
        os.path.join(IBD_COLON_DIR, f"{comp}.barcodes2.tsv"),
        comp,
    )
    matched = barcodes.isin(meta.index)
    n_unmatched = (~matched).sum()
    if n_unmatched:
        log.warning(f"  [{comp}] {n_unmatched}/{len(barcodes)} barcodes have no all.meta2.txt "
                     f"row -- these cells will be DROPPED")
    else:
        log.info(f"  [{comp}] all {len(barcodes)} barcodes matched all.meta2.txt")

    tsne = load_tsne(comp, barcodes)

    obs = pd.DataFrame(index=pd.Index(barcodes.to_numpy())).join(meta, how="left").join(tsne, how="left")
    obs[DATA_SOURCE_COL] = "full_cohort_mtx"
    var = pd.DataFrame(index=pd.Index(genes.to_numpy()))
    adata = anndata.AnnData(X=mat, obs=obs, var=var)
    keep = matched.to_numpy()
    if not keep.all():
        adata = adata[keep].copy()

    adata.obs["compartment"] = comp
    adata.obs[CELL_TYPE_LINEAGE_COL] = adata.obs[IBD_COLON_CLUSTER_COL].map(lineage_map)
    n_unmapped_lineage = adata.obs[CELL_TYPE_LINEAGE_COL].isna().sum()
    if n_unmapped_lineage:
        log.warning(f"  [{comp}] {n_unmapped_lineage} cells have a Cluster value "
                     f"not present in cell_subsets.txt")

    log.info(f"  [{comp}] final: {adata.shape[0]} cells x {adata.shape[1]} genes, "
              f"{adata.obs[IBD_COLON_CLUSTER_COL].nunique()} cell types, "
              f"{adata.obs['Subject'].nunique()} donors (of 30 full cohort)")
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="verify inputs only, write nothing")
    parser.add_argument("--compartments", nargs="+", choices=IBD_COLON_COMPARTMENTS, default=IBD_COLON_COMPARTMENTS,
                        help="restrict to these compartments (default: all) -- run one per invocation "
                             "to keep peak memory/process lifetime down for the largest matrices")
    args = parser.parse_args()

    lineage_map = load_cell_subsets()
    meta = load_meta()

    for comp in args.compartments:
        adata = build_compartment_adata(comp, meta, lineage_map, dry_run=args.dry_run)
        if args.dry_run:
            continue
        out_path = IBD_COLON_H5AD_BY_COMPARTMENT[comp]
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        adata.write_h5ad(out_path)
        log.info(f"  wrote {out_path}")
        del adata
        gc.collect()

    log.info("done" + (" (dry-run, nothing written)" if args.dry_run else ""))


if __name__ == "__main__":
    main()
