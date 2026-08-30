"""Compute memento method-of-moments mean expression per fibroblast subtype for the
skin fibroblast atlas (Steele et al., Nat. Immunol. 2025). Implements
prompts/skin_means.md steps 2-5, given the download resolved in download_skin_data.py.

Architecture mirrors compute_means.py (chunked processing to bound peak memory,
streaming Parquet writes so no full-dataset table is ever held in memory) but is
kept as a separate, self-contained script rather than modifying compute_means.py,
since that pipeline is already verified/complete (see ANALYSIS_SUMMARY.md) and this
dataset's schema differs enough (no donor_id column, no assay column, 4 grouping
label columns instead of 2) that sharing the memento-running function would mean
threading extra parameters through code with no other caller. Only the
dataset-agnostic utilities (ParquetWriterRegistry, peak_rss_gb, available_mem_gb,
build_donor_chunks) are imported and reused as-is.

Per-sample identity recovery (see config.py's long comment for the full story): the
downloaded object has no donor_id/sample_id *column*, but per-sample identity is
encoded in the cell barcode itself (e.g. "..._GSM6111852"), appended during
integration to keep barcodes unique across the ~30 concatenated source studies. A
first pass concluded this wasn't reliably parseable, based on a naive global
"last underscore token" split; per-study it's almost always a clean, uniform
`<barcode>_<GSM or SRS accession>` pattern. `extract_sample_id` below recovers it,
with two studies (SKIN_MESSY_STUDIES) left at coarse GSE-level per user decision
rather than chasing their messier formats.

See config.py (SKIN_* constants) for the remaining schema decisions this script
assumes -- raw-counts slot and grouping label columns -- and why each was made.

Usage:
    python compute_skin_means.py --dry-run   # inspect schema only
    python compute_skin_means.py             # full run
"""
import argparse
import gc
import logging
import os
import re
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
    DEFAULT_CAPTURE_RATE, MEMENTO_MIN_CELL_COUNT, SKIN_DERIVED_COL,
    SKIN_FIBROBLAST_H5AD, SKIN_FIBROBLAST_OUTPUT_PATH, SKIN_LABEL_COLS,
    SKIN_MESSY_STUDIES, SKIN_RAW_SLOT, SKIN_STUDY_COL,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("compute_skin_means")

LOW_N_THRESHOLD = 20
SAMPLE_ID_COL = "_sample_id"  # derived column added to obs before chunking

OUTPUT_COLS = [
    "dataset", "sample_id", "source_study", "patient_status", "lesional_status",
    "disease_category", "fibroblast_subtype", "fibroblast_subtype_skin_nomenclature",
    "gene", "n_cells", "mean_expression", "assay", "raw_counts_slot",
    "low_confidence", "pipeline_run_date",
]
STRING_COLS = [
    "dataset", "sample_id", "source_study", "patient_status", "lesional_status",
    "disease_category", "fibroblast_subtype", "fibroblast_subtype_skin_nomenclature",
    "gene", "assay", "raw_counts_slot", "pipeline_run_date",
]

_PATIENT_DAY_RE = re.compile(r"(patient\d+_Day\d+)")


def extract_sample_id(barcode, source_study):
    """Recover a per-sample id from the cell barcode, per source study.

    Verified against the full obs table before this function was written (see
    ANALYSIS_SUMMARY.md) -- cardinality checked per study to confirm each rule
    produces multiple distinct, plausible samples rather than one repeated
    constant. Studies in SKIN_MESSY_STUDIES are deliberately left unparsed
    (returns source_study itself, i.e. coarse GSE-level) per user decision.
    """
    if source_study in SKIN_MESSY_STUDIES:
        return source_study
    if source_study == SKIN_STUDY_COL:
        # The unlabeled/blank-GSE cohort (this paper's own newly-generated data, not
        # a public GEO deposit) encodes "patientN_DayX" directly in the barcode
        # alongside other tokens (e.g. "..._1_GSM_patient5_Day0_NonLesional") --
        # lesional status is already captured separately in obs, so only the
        # patient+timepoint token is pulled out here.
        m = _PATIENT_DAY_RE.search(barcode)
        return m.group(1) if m else source_study
    if source_study == "Reynolds":
        # Sample codes like "P1"/"E4"/"s3" with inconsistent case, plus occasional
        # "_1"/"_2" replicate suffixes that are kept (likely distinct
        # libraries/technical samples, not necessarily deduplicatable to the same
        # underlying id without checking the source paper).
        suffix = barcode.split("_", 1)[1] if "_" in barcode else source_study
        return suffix.upper()
    # Default: clean single-underscore "<barcode>_<GSM or SRS accession>" pattern,
    # true for the remaining ~27 of 32 source studies.
    parts = barcode.split("_", 1)
    return parts[1] if len(parts) > 1 else source_study


def run_memento_on_chunk(adata_chunk, label_cols):
    """Compute per-group memento means within one chunk, grouping by all of label_cols.

    Unlike compute_means.py's version (fixed cell_type x assay), this takes an
    arbitrary list of obs columns to group by -- this dataset has no assay column
    (flat capture rate for every cell) and needs 4 label columns, not 2.
    """
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
        # memento's create_groups needs plain strings, not pandas categoricals
        # (unused categories from the full-dataset dtype would otherwise leak
        # through as phantom empty groups within a single-GSE chunk).
        adata_chunk.obs[col] = adata_chunk.obs[col].astype(str)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        memento.setup_memento(
            adata_chunk, q_column="capture_rate", min_cell_count=MEMENTO_MIN_CELL_COUNT
        )
        memento.create_groups(adata_chunk, label_columns=label_cols)
        memento.compute_1d_moments(adata_chunk, min_perc_group=0.0, filter_genes=False)

    groups_df = memento.get_groups(adata_chunk)
    genes = adata_chunk.var.index.to_numpy()

    rows = []
    for group_key in adata_chunk.uns["memento"]["groups"]:
        mean_arr = np.asarray(adata_chunk.uns["memento"]["1d_moments"][group_key][0])
        n_cells = adata_chunk.uns["memento"]["group_cells"][group_key].shape[0]
        labels = groups_df.loc[group_key]
        row = {col: labels[col] for col in label_cols}
        row.update(n_cells=n_cells, mean_expression=mean_arr, gene=genes)
        rows.append(row)
    return rows


def finalize_chunk_df(rows, sample_id, source_study, derived_by_subtype, pipeline_run_date):
    dfs = []
    for r in rows:
        df = pd.DataFrame({
            "fibroblast_subtype": r["celltype"],
            "patient_status": r["Patient_status"],
            "lesional_status": r["lesional_vs_nonlesional"],
            "disease_category": r["disease_category_orig"],
            "n_cells": r["n_cells"],
            "mean_expression": r["mean_expression"],
            "gene": r["gene"],
        })
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True) if len(dfs) > 1 else dfs[0]

    df["dataset"] = "skin_fibroblast_atlas"
    df["sample_id"] = sample_id
    df["source_study"] = source_study
    df["fibroblast_subtype_skin_nomenclature"] = df["fibroblast_subtype"].map(derived_by_subtype)
    df["assay"] = "unknown"
    df["raw_counts_slot"] = SKIN_RAW_SLOT
    df["pipeline_run_date"] = pipeline_run_date
    df["low_confidence"] = df["n_cells"] < LOW_N_THRESHOLD

    df = df[OUTPUT_COLS]
    for c in STRING_COLS:
        df[c] = df[c].astype(str)
    df["n_cells"] = df["n_cells"].astype("int64")
    df["mean_expression"] = df["mean_expression"].astype("float64")
    df["low_confidence"] = df["low_confidence"].astype("bool")
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="only check schema/raw-slot")
    args = parser.parse_args()

    f = h5py.File(SKIN_FIBROBLAST_H5AD, "r")
    log.info(f"top-level groups: {list(f.keys())}")

    raw_X = anndata_io.sparse_dataset(f[SKIN_RAW_SLOT])
    obs = anndata_io.read_elem(f["obs"])
    var_df = anndata_io.read_elem(f["var"])
    log.info(f"obs shape={obs.shape}, columns={list(obs.columns)}")
    log.info(f"raw slot {SKIN_RAW_SLOT} shape={raw_X.shape}")

    barcodes = obs.index.astype(str)
    source_studies = obs[SKIN_STUDY_COL].astype(str)
    obs[SAMPLE_ID_COL] = [
        extract_sample_id(bc, gse) for bc, gse in zip(barcodes, source_studies)
    ]
    sample_to_study = dict(zip(obs[SAMPLE_ID_COL], source_studies))
    log.info(
        f"recovered {obs[SAMPLE_ID_COL].nunique()} distinct sample_id chunks "
        f"from {source_studies.nunique()} source studies "
        f"({len(SKIN_MESSY_STUDIES)} left at coarse GSE-level: {sorted(SKIN_MESSY_STUDIES)})"
    )

    derived_by_subtype = (
        obs.groupby("celltype", observed=True)[SKIN_DERIVED_COL]
        .agg(lambda s: next(iter(s.dropna()), None))
        .to_dict()
    )
    log.info(f"fibroblast_subtype -> {SKIN_DERIVED_COL} mapping: {derived_by_subtype}")

    if args.dry_run:
        f.close()
        return

    os.makedirs(os.path.dirname(SKIN_FIBROBLAST_OUTPUT_PATH), exist_ok=True)
    writers = ParquetWriterRegistry()
    pipeline_run_date = time.strftime("%Y-%m-%d")

    t0 = time.time()
    n_chunks = 0
    n_rows_written = 0
    try:
        for chunk_label, sample_id, positions in build_donor_chunks(obs, SAMPLE_ID_COL):
            n_chunks += 1
            try:
                counts = raw_X[positions]
                sub = anndata.AnnData(
                    X=counts, obs=obs.iloc[positions].copy(), var=var_df.copy()
                )
                rows = run_memento_on_chunk(sub, SKIN_LABEL_COLS)
                del sub
                source_study = sample_to_study[sample_id]
                final_df = finalize_chunk_df(
                    rows, sample_id, source_study, derived_by_subtype, pipeline_run_date
                )
                writers.write(SKIN_FIBROBLAST_OUTPUT_PATH, final_df)
                n_rows_written += len(final_df)
            except Exception as e:
                log.error(f"chunk {chunk_label} FAILED: {e}")
                continue

            gc.collect()
            log.info(
                f"... {n_chunks} chunks done ({chunk_label}, {len(positions)} cells), "
                f"{time.time()-t0:.0f}s elapsed, peak_rss={peak_rss_gb():.1f}GB, "
                f"mem_available={available_mem_gb():.1f}GB, rows_written={n_rows_written}"
            )
    finally:
        writers.close_all()
        f.close()

    log.info(f"done, {n_rows_written} total rows written in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
