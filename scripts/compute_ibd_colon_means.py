"""Compute memento method-of-moments mean expression per (donor x cell type x
health status) for the IBD colon atlas (Smillie et al., Cell 2019, SCP259 via
HCA), genome-wide, across all 3 compartments (Epi/Fib/Imm).

Architecture mirrors compute_means.py (chunked-by-donor processing to bound
peak memory, streaming Parquet writes) but loops over the 3 per-compartment
h5ad files built by build_ibd_colon_h5ad.py the same way compute_means.py
loops over manifest.csv rows -- all 3 compartments write to the SAME output
parquet (IBD_COLON_OUTPUT_PATH), with `compartment` as a row-level column, per
config.py's documented reason the 3 compartments were never merged into one
raw object (divergent per-compartment gene panels). Kept as a separate,
self-contained script rather than modifying compute_means.py, following the
precedent set by compute_skin_means.py, for the same reason: this schema
(Health as a within-donor grouping axis, no assay column, a 51->7 lineage
label) differs enough that threading extra parameters through the CELLxGene-
specific function would only serve this one caller.

Donor identity: `Subject`, used directly -- verified clean at inspection time
(see config.py's long comment above IBD_COLON_DONOR_COL): exactly 30 distinct
values, each cleanly mapping to 1 or 2 Health values with no collision/alias
pattern found, unlike the Gut Cell Atlas's donor_id. No *_UNIFIED_COL
correction was needed, but the check itself (not just the raw column) is
re-logged here at runtime for the record.

`Health` (Healthy / Inflamed / Non-inflamed) varies WITHIN a donor (18 UC
donors contribute both an Inflamed and a Non-inflamed sample) so it is a
memento grouping label alongside cell type, not a chunk boundary -- the same
role `assay` plays in compute_means.py's cell_type x assay grouping.

Usage:
    python compute_ibd_colon_means.py --dry-run   # inspect schema only
    python compute_ibd_colon_means.py             # full run
"""
import argparse
import gc
import logging
import os
import time
import warnings

import anndata
import anndata.io as anndata_io
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import memento

from compute_means import (
    ParquetWriterRegistry, available_mem_gb, build_donor_chunks, peak_rss_gb,
)
from config import (
    DEFAULT_CAPTURE_RATE, IBD_COLON_CLUSTER_COL, IBD_COLON_COMPARTMENTS,
    IBD_COLON_DONOR_COL, IBD_COLON_H5AD_BY_COMPARTMENT, IBD_COLON_HEALTH_COL,
    IBD_COLON_OUTPUT_PATH, IBD_COLON_RAW_SLOT, MEMENTO_MIN_CELL_COUNT,
)
from build_ibd_colon_h5ad import CELL_TYPE_LINEAGE_COL

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("compute_ibd_colon_means")

LOW_N_THRESHOLD = 20
LABEL_COLS = [IBD_COLON_CLUSTER_COL, IBD_COLON_HEALTH_COL]

OUTPUT_COLS = [
    "dataset", "compartment", "donor_id", "health", "cell_type_fine",
    "cell_type_lineage", "gene", "n_cells", "mean_expression", "assay",
    "low_confidence", "pipeline_run_date",
]
STRING_COLS = [
    "dataset", "compartment", "donor_id", "health", "cell_type_fine",
    "cell_type_lineage", "gene", "assay", "pipeline_run_date",
]


def verify_donor_identity(obs):
    """Required check (project convention): confirm Subject doesn't silently
    collide/alias across people before trusting it as the chunk boundary --
    same class of check that caught the Gut Cell Atlas's donor_id problem."""
    n_subjects = obs[IBD_COLON_DONOR_COL].nunique()
    health_per_subject = obs.groupby(IBD_COLON_DONOR_COL, observed=True)[IBD_COLON_HEALTH_COL].nunique()
    bad = health_per_subject[health_per_subject > 2]
    log.info(f"  donor identity check: {n_subjects} distinct {IBD_COLON_DONOR_COL!r} values, "
              f"health-states-per-subject distribution={health_per_subject.value_counts().to_dict()}")
    if len(bad):
        log.warning(f"  {len(bad)} subjects have >2 distinct Health values -- "
                     f"possible id collision, inspect before trusting: {bad.index.tolist()}")
    else:
        log.info(f"  no subject exceeds 2 Health values -- consistent with clean donor identity "
                  f"(12 single-Health healthy donors + 18 paired-biopsy UC donors)")


def run_memento_on_chunk(adata_chunk, label_cols):
    adata_chunk = adata_chunk.copy()
    if not sp.issparse(adata_chunk.X):
        adata_chunk.X = sp.csr_matrix(adata_chunk.X)
    else:
        adata_chunk.X = adata_chunk.X.tocsr()
    adata_chunk.X = adata_chunk.X.astype(np.float64)

    adata_chunk.obs["capture_rate"] = DEFAULT_CAPTURE_RATE

    for col in label_cols:
        if col not in adata_chunk.obs.columns:
            raise KeyError(f"missing expected column {col!r}; have {list(adata_chunk.obs.columns)}")
        adata_chunk.obs[col] = adata_chunk.obs[col].astype(str)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        memento.setup_memento(
            adata_chunk, q_column="capture_rate", min_cell_count=MEMENTO_MIN_CELL_COUNT
        )
        memento.create_groups(adata_chunk, label_columns=label_cols)
        memento.compute_1d_moments(adata_chunk, min_perc_group=0.0, filter_genes=False)

    # NOT memento.get_groups(): its own implementation (`get_groups` in
    # memento's main.py) deliberately re-encodes any label column with exactly
    # 2 distinct values in this chunk into 0/1 float codes ("so that
    # downstream regression ... receives a numeric design matrix", intended
    # for memento's own binary_test_1d/2d treatment-column convention) --
    # confirmed by reading memento's source after finding "Health" silently
    # turned into 0.0/1.0 for every 2-level (Inflamed/Non-inflamed, no
    # Healthy) donor chunk, i.e. all 18 UC donors. The group KEY string itself
    # is always the real label text (memento builds group_key by joining the
    # raw label values with label_delimiter, get_groups() just re-derives
    # columns from it and then "helpfully" re-codes 2-level ones) -- parsed
    # directly here instead, the same way memento's own get_groups() does
    # internally minus the lossy re-coding step.
    label_delimiter = adata_chunk.uns["memento"]["label_delimiter"]
    genes = adata_chunk.var.index.to_numpy()

    rows = []
    for group_key in adata_chunk.uns["memento"]["groups"]:
        mean_arr = np.asarray(adata_chunk.uns["memento"]["1d_moments"][group_key][0])
        n_cells = adata_chunk.uns["memento"]["group_cells"][group_key].shape[0]
        label_values = group_key.split(label_delimiter)[1:]
        row = dict(zip(label_cols, label_values))
        row.update(n_cells=n_cells, mean_expression=mean_arr, gene=genes)
        rows.append(row)
    return rows


def finalize_chunk_df(rows, donor_id, compartment, lineage_by_cluster, pipeline_run_date):
    dfs = []
    for r in rows:
        df = pd.DataFrame({
            "cell_type_fine": r[IBD_COLON_CLUSTER_COL],
            "health": r[IBD_COLON_HEALTH_COL],
            "n_cells": r["n_cells"],
            "mean_expression": r["mean_expression"],
            "gene": r["gene"],
        })
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True) if len(dfs) > 1 else dfs[0]

    df["dataset"] = "ibd_colon_atlas"
    df["compartment"] = compartment
    df["donor_id"] = donor_id
    df["cell_type_lineage"] = df["cell_type_fine"].map(lineage_by_cluster)
    df["assay"] = "unknown"
    df["pipeline_run_date"] = pipeline_run_date
    df["low_confidence"] = df["n_cells"] < LOW_N_THRESHOLD

    df = df[OUTPUT_COLS]
    for c in STRING_COLS:
        df[c] = df[c].astype(str)
    df["n_cells"] = df["n_cells"].astype("int64")
    df["mean_expression"] = df["mean_expression"].astype("float64")
    df["low_confidence"] = df["low_confidence"].astype("bool")
    return df


def process_compartment(comp, writers, dry_run=False):
    h5ad_path = IBD_COLON_H5AD_BY_COMPARTMENT[comp]
    log.info(f"=== {comp} ({h5ad_path}) ===")
    f = h5py.File(h5ad_path, "r")
    log.info(f"  top-level groups: {list(f.keys())}")

    raw_X = anndata_io.sparse_dataset(f[IBD_COLON_RAW_SLOT])
    obs = anndata_io.read_elem(f["obs"])
    var_df = anndata_io.read_elem(f["var"])
    log.info(f"  obs shape={obs.shape}, columns={list(obs.columns)}")
    log.info(f"  raw slot {IBD_COLON_RAW_SLOT} shape={raw_X.shape}")

    verify_donor_identity(obs)
    lineage_by_cluster = (
        obs.groupby(IBD_COLON_CLUSTER_COL, observed=True)[CELL_TYPE_LINEAGE_COL]
        .agg(lambda s: next(iter(s.dropna()), None)).to_dict()
    )

    if dry_run:
        f.close()
        return 0

    pipeline_run_date = time.strftime("%Y-%m-%d")
    t0 = time.time()
    n_chunks = 0
    n_rows_written = 0
    for chunk_label, donor_id, positions in build_donor_chunks(obs, IBD_COLON_DONOR_COL):
        n_chunks += 1
        try:
            counts = raw_X[positions]
            sub = anndata.AnnData(
                X=counts, obs=obs.iloc[positions].copy(), var=var_df.copy()
            )
            rows = run_memento_on_chunk(sub, LABEL_COLS)
            del sub
            final_df = finalize_chunk_df(rows, donor_id, comp, lineage_by_cluster, pipeline_run_date)
            writers.write(IBD_COLON_OUTPUT_PATH, final_df)
            n_rows_written += len(final_df)
        except Exception as e:
            log.error(f"  [{comp}] chunk {chunk_label} FAILED: {e}")
            continue

        if n_chunks % 5 == 0:
            gc.collect()
            log.info(
                f"  ... {n_chunks} chunks done ({chunk_label}, {len(positions)} cells), "
                f"{time.time()-t0:.0f}s elapsed, peak_rss={peak_rss_gb():.1f}GB, "
                f"mem_available={available_mem_gb():.1f}GB, rows_written={n_rows_written}"
            )

    log.info(f"  [{comp}] processed {n_chunks} chunks in {time.time()-t0:.0f}s, "
              f"wrote {n_rows_written} rows")
    f.close()
    return n_rows_written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="only check schema/donor identity")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(IBD_COLON_OUTPUT_PATH), exist_ok=True)
    writers = ParquetWriterRegistry()
    total_rows = 0
    try:
        for comp in IBD_COLON_COMPARTMENTS:
            total_rows += process_compartment(comp, writers, dry_run=args.dry_run) or 0
    finally:
        writers.close_all()
    log.info(f"done, {total_rows} total rows written across all compartments")


if __name__ == "__main__":
    main()
