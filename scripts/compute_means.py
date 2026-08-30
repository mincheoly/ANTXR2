"""Compute cell-type-specific mean expression with memento's method-of-moments
estimator, for each downloaded CELLxGene dataset.

For each dataset, cells are processed in per-donor chunks (further split if a donor
exceeds MAX_CHUNK_CELLS) to bound peak memory -- the full raw count matrix for these
datasets (up to ~30-40GB as a sparse CSR in memory) does not fit in this machine's
RAM, so memento's setup_memento/create_groups/compute_1d_moments calls run once per
chunk rather than once for the whole file. Within each chunk we group by
(cell_type x assay); donor identity is the chunk boundary itself. This matches the
plan's requested `dataset_id x donor_id x cell_type` granularity, and keeps assay
identity intact for platforms (e.g. Tabula Sapiens) that mix 10x and Smart-seq2/3.

Data access deliberately bypasses `anndata.read_h5ad(..., backed="r")`: for the
Tabula Sapiens file, opening that (even before touching X) pulled tens of GB into
memory and OOM-killed the process -- `obsp` holds a full 1.14M x 1.14M neighbor
graph (connectivities/distances) that anndata does not appear to keep lazily even in
backed mode. We never need obsp/obsm/varm/uns here, so we read only `obs`, `var` (or
`raw/var`), and `X`/`raw/X` directly via h5py + `anndata.io.sparse_dataset` (which
gives a real lazy, row-sliceable proxy for the on-disk CSR matrix) and
`anndata.io.read_elem` (for the obs/var DataFrames). This keeps peak memory
proportional to one chunk's rows, not the whole file.

Usage:
    python compute_means.py --dry-run [dataset_id ...]   # inspect schema only
    python compute_means.py [dataset_id ...]              # full run (all datasets
                                                            # in manifest.csv if none given)
"""
import argparse
import csv
import gc
import logging
import os
import resource
import sys
import time
import warnings

import anndata
import anndata.io as anndata_io
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

import memento

from config import (
    ASSAY_COL, CAPTURE_RATE_BY_ASSAY, CELL_TYPE_COL, CELL_TYPE_ID_COL,
    DEFAULT_CAPTURE_RATE, DISEASE_COL, DONOR_COL, GENE_ID_COL, GENE_NAME_COL,
    MAX_CHUNK_CELLS, MEMENTO_MIN_CELL_COUNT, OUTPUT_DIR, RAW_DIR, TISSUE_COL,
    TISSUE_GENERAL_COL,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("compute_means")

LOW_N_THRESHOLD = 20

OUTPUT_COLS = [
    "source_collection", "collection_name", "dataset_id", "donor_id", "tissue",
    "cell_type", "gene", "gene_ontology_id", "n_cells", "mean_expression",
    "disease", "assay", "low_confidence", "pipeline_run_date",
]


class ParquetWriterRegistry:
    """Lazily-opened, append-only Parquet writers keyed by output path.

    Each dataset's results can be tens of millions of rows once expanded to
    per-gene x per-group; accumulating those across donors (or datasets) in a
    single in-memory DataFrame before writing is what OOM-killed the first attempt
    at the Tabula Sapiens file (peak per-chunk memory was fine at ~21GB; the
    post-loop pd.concat + column construction over the whole dataset is what spiked
    to ~95GB). Writing one donor's table at a time, immediately, keeps peak memory
    proportional to a single donor's rows.
    """

    def __init__(self):
        self._writers = {}

    def write(self, path, df):
        table = pa.Table.from_pandas(df, preserve_index=False)
        writer = self._writers.get(path)
        if writer is None:
            writer = pq.ParquetWriter(path, table.schema)
            self._writers[path] = writer
        writer.write_table(table)

    def close_all(self):
        for writer in self._writers.values():
            writer.close()
        self._writers = {}


def peak_rss_gb():
    """Peak resident set size of this process so far, in GB (Linux: ru_maxrss is KB)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def available_mem_gb():
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable"):
                return int(line.split()[1]) / 1e6
    return float("nan")


def read_manifest():
    path = os.path.join(RAW_DIR, "manifest.csv")
    with open(path) as f:
        return list(csv.DictReader(f))


def resolve_raw_slot(f, raw_data_location, dataset_id):
    """Pick the on-disk group holding raw integer counts, verifying with a value check.

    Returns (group_path, var_df) where group_path is "raw/X" or "X" and var_df is the
    matching gene metadata (raw/var or var). Only ever reads obs/var/X/raw.X groups --
    never obsp/obsm/varm/uns.
    """
    candidates = []
    loc = (raw_data_location or "").strip()
    has_raw = "raw" in f
    if loc == "raw.X" and has_raw:
        candidates.append("raw/X")
    elif loc == "X":
        candidates.append("X")
    if has_raw:
        candidates.append("raw/X")
    candidates.append("X")

    seen = set()
    for group_path in candidates:
        if group_path in seen:
            continue
        seen.add(group_path)
        try:
            ds = anndata_io.sparse_dataset(f[group_path])
        except Exception as e:
            log.warning(f"  [{dataset_id}] could not open slot {group_path}: {e}")
            continue
        n = ds.shape[0]
        sample_idx = np.linspace(0, n - 1, num=min(200, n), dtype=int)
        sample = ds[sample_idx]
        sample = sample.toarray() if sp.issparse(sample) else np.asarray(sample)
        if sample.size == 0:
            continue
        max_val = float(sample.max())
        frac_integer = np.mean(np.isclose(sample, np.round(sample)))
        looks_raw = max_val > 20 and frac_integer > 0.99
        log.info(f"  [{dataset_id}] slot {group_path}: sample max={max_val:.2f}, frac_integer={frac_integer:.3f}, looks_raw={looks_raw}")
        if looks_raw:
            var_path = "raw/var" if group_path == "raw/X" else "var"
            var_df = anndata_io.read_elem(f[var_path])
            return group_path, var_df
    raise RuntimeError(f"[{dataset_id}] could not identify a raw-counts slot (checked {seen})")


def build_donor_chunks(obs, donor_col):
    """Yield (chunk_label, donor_id, row_positions) tuples, splitting oversized donors."""
    if donor_col not in obs.columns:
        yield ("all", None, np.arange(len(obs)))
        return
    donor_values = obs[donor_col].astype(str)
    groups = list(donor_values.groupby(donor_values).groups.items())
    # Largest donor first: surfaces an OOM/slowdown on the riskiest chunk early
    # rather than after many small chunks have already succeeded.
    groups.sort(key=lambda kv: -len(kv[1]))
    for donor_id, idx in groups:
        positions = obs.index.get_indexer(idx)
        positions = positions[positions >= 0]
        if len(positions) == 0:
            continue
        if len(positions) <= MAX_CHUNK_CELLS:
            yield (str(donor_id), str(donor_id), positions)
        else:
            n_parts = int(np.ceil(len(positions) / MAX_CHUNK_CELLS))
            for i, part in enumerate(np.array_split(positions, n_parts)):
                yield (f"{donor_id}__part{i}", str(donor_id), part)


def run_memento_on_chunk(adata_chunk, cell_type_col, assay_col):
    adata_chunk = adata_chunk.copy()
    if not sp.issparse(adata_chunk.X):
        adata_chunk.X = sp.csr_matrix(adata_chunk.X)
    else:
        adata_chunk.X = adata_chunk.X.tocsr()
    adata_chunk.X = adata_chunk.X.astype(np.float64)

    if assay_col in adata_chunk.obs.columns:
        assay_series = adata_chunk.obs[assay_col].astype(str)
    else:
        assay_series = pd.Series("unknown", index=adata_chunk.obs.index)
    adata_chunk.obs["capture_rate"] = assay_series.map(
        lambda a: CAPTURE_RATE_BY_ASSAY.get(a, DEFAULT_CAPTURE_RATE)
    ).astype(float)
    adata_chunk.obs["_assay_for_grouping"] = assay_series.values

    if cell_type_col not in adata_chunk.obs.columns:
        raise KeyError(f"missing expected column {cell_type_col!r}; have {list(adata_chunk.obs.columns)}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        memento.setup_memento(
            adata_chunk, q_column="capture_rate", min_cell_count=MEMENTO_MIN_CELL_COUNT
        )
        memento.create_groups(adata_chunk, label_columns=[cell_type_col, "_assay_for_grouping"])
        memento.compute_1d_moments(adata_chunk, min_perc_group=0.0, filter_genes=False)

    groups_df = memento.get_groups(adata_chunk)
    genes = adata_chunk.var.index.to_numpy()

    rows = []
    for group_key in adata_chunk.uns["memento"]["groups"]:
        mean_arr = np.asarray(adata_chunk.uns["memento"]["1d_moments"][group_key][0])
        n_cells = adata_chunk.uns["memento"]["group_cells"][group_key].shape[0]
        labels = groups_df.loc[group_key]
        rows.append({
            "cell_type": labels[cell_type_col],
            "assay": labels["_assay_for_grouping"],
            "n_cells": n_cells,
            "mean_expression": mean_arr,
            "gene": genes,
        })
    return rows


def combine_weighted(group_rows):
    """Combine multiple sub-chunk results for the same (donor, cell_type, assay, gene)
    group into one row, weighting the mean by n_cells."""
    df = pd.concat(group_rows, ignore_index=True)
    weight = df["n_cells"].astype(float)
    df["_weighted_mean"] = df["mean_expression"] * weight
    agg = df.groupby(["cell_type", "assay", "gene"], as_index=False).agg(
        n_cells=("n_cells", "sum"),
        _weighted_mean=("_weighted_mean", "sum"),
    )
    agg["mean_expression"] = agg["_weighted_mean"] / agg["n_cells"]
    return agg.drop(columns="_weighted_mean")


STRING_COLS = [
    "source_collection", "collection_name", "dataset_id", "donor_id", "tissue",
    "cell_type", "gene", "gene_ontology_id", "disease", "assay", "pipeline_run_date",
]


def finalize_donor_df(df, true_donor_id, dataset_id, gene_name_by_id, manifest_row):
    df = df.copy()
    df["donor_id"] = true_donor_id
    df["gene_ontology_id"] = df["gene"]
    df["gene"] = df["gene_ontology_id"].map(gene_name_by_id).fillna(df["gene_ontology_id"])
    df["dataset_id"] = dataset_id
    df["collection_name"] = manifest_row["collection_name"]
    df["source_collection"] = manifest_row["collection_id"]
    df["tissue"] = manifest_row.get("tissue")
    df["disease"] = manifest_row.get("disease")
    df["pipeline_run_date"] = time.strftime("%Y-%m-%d")
    df["low_confidence"] = df["n_cells"] < LOW_N_THRESHOLD
    df = df[OUTPUT_COLS]
    # Pin dtypes so every batch produces an identical pyarrow schema -- the
    # ParquetWriter is opened once per output path and every subsequent
    # write_table() must match that first batch's schema exactly.
    for c in STRING_COLS:
        df[c] = df[c].astype(str)
    df["n_cells"] = df["n_cells"].astype("int64")
    df["mean_expression"] = df["mean_expression"].astype("float64")
    df["low_confidence"] = df["low_confidence"].astype("bool")
    return df


def process_dataset(manifest_row, writers, dry_run=False):
    dataset_id = manifest_row["dataset_id"]
    collection_name = manifest_row["collection_name"]
    local_path = manifest_row["local_path"]
    log.info(f"=== {collection_name} / {dataset_id} ===")

    f = h5py.File(local_path, "r")
    log.info(f"  top-level groups: {list(f.keys())}")

    slot_path, var_df = resolve_raw_slot(f, manifest_row.get("raw_data_location"), dataset_id)
    log.info(f"  using raw-counts slot: {slot_path}")

    if dry_run:
        f.close()
        return 0

    raw_X = anndata_io.sparse_dataset(f[slot_path])
    obs = anndata_io.read_elem(f["obs"])
    log.info(f"  obs shape={obs.shape}, columns={list(obs.columns)}")

    gene_ids = var_df.index.to_numpy()
    gene_names = (
        var_df[GENE_NAME_COL].to_numpy()
        if GENE_NAME_COL in var_df.columns
        else gene_ids
    )
    gene_name_by_id = dict(zip(gene_ids, gene_names))

    collection_path = os.path.join(OUTPUT_DIR, f"{collection_name}_celltype_means.parquet")
    combined_path = os.path.join(OUTPUT_DIR, "combined_celltype_means.parquet")

    def flush_donor(true_donor_id, dfs):
        if not dfs:
            return 0
        combined = combine_weighted(dfs) if len(dfs) > 1 else dfs[0]
        final_df = finalize_donor_df(combined, true_donor_id, dataset_id, gene_name_by_id, manifest_row)
        writers.write(collection_path, final_df)
        writers.write(combined_path, final_df)
        return len(final_df)

    t0 = time.time()
    n_chunks = 0
    n_rows_written = 0
    pending_donor = None
    pending_dfs = []
    for chunk_label, true_donor_id, positions in build_donor_chunks(obs, DONOR_COL):
        n_chunks += 1
        if true_donor_id != pending_donor:
            n_rows_written += flush_donor(pending_donor, pending_dfs)
            pending_donor, pending_dfs = true_donor_id, []
        try:
            counts = raw_X[positions]
            sub = anndata.AnnData(
                X=counts, obs=obs.iloc[positions].copy(), var=var_df.copy()
            )

            if CELL_TYPE_COL not in sub.obs.columns:
                log.warning(f"  [{dataset_id}] chunk {chunk_label}: missing {CELL_TYPE_COL}, skipping")
                continue

            rows = run_memento_on_chunk(sub, CELL_TYPE_COL, ASSAY_COL)
            for r in rows:
                pending_dfs.append(pd.DataFrame({
                    "cell_type": r["cell_type"],
                    "assay": r["assay"],
                    "n_cells": r["n_cells"],
                    "mean_expression": r["mean_expression"],
                    "gene": r["gene"],
                }))
            del sub
        except Exception as e:
            log.error(f"  [{dataset_id}] chunk {chunk_label} FAILED: {e}")
        if n_chunks % 5 == 0:
            gc.collect()
            log.info(
                f"  ... {n_chunks} chunks done, {time.time()-t0:.0f}s elapsed, "
                f"peak_rss={peak_rss_gb():.1f}GB, mem_available={available_mem_gb():.1f}GB, "
                f"rows_written={n_rows_written}"
            )

    n_rows_written += flush_donor(pending_donor, pending_dfs)
    log.info(f"  processed {n_chunks} chunks in {time.time()-t0:.0f}s, wrote {n_rows_written} rows")
    f.close()
    return n_rows_written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_ids", nargs="*", help="restrict to these dataset_id(s)")
    parser.add_argument("--dry-run", action="store_true", help="only check schema/raw-slot per dataset")
    args = parser.parse_args()

    manifest = read_manifest()
    if args.dataset_ids:
        manifest = [r for r in manifest if r["dataset_id"] in args.dataset_ids]
        if not manifest:
            print("no matching dataset_id in manifest.csv", file=sys.stderr)
            sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    writers = ParquetWriterRegistry()
    total_rows = 0
    try:
        for row in manifest:
            total_rows += process_dataset(row, writers, dry_run=args.dry_run) or 0
    finally:
        writers.close_all()
    log.info(f"done, {total_rows} total rows written across all datasets")


if __name__ == "__main__":
    main()
