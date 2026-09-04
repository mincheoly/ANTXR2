"""Multi-gene query/aggregation library, generalizing query_antxr2_table.py's
single-gene pattern (row-group-filtered parquet scan + two-step donor-equal-
weighted aggregation) to an arbitrary gene panel. query_antxr2_table.py itself
is left untouched -- it's a working, already-delivered script; this module
generalizes its pattern for the new heatmap script rather than editing it.
"""
import os
import resource

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def scan_row_groups_for_genes(path, genes, gene_col="gene", cache_path=None):
    """Row-group-filtered scan of a parquet file for a small gene panel.

    Avoids loading the full file into memory -- only row groups that contain at
    least one of `genes` are read in full. `value_set` for pc.is_in must be a
    pyarrow.array, not a pandas.array (the latter raises TypeError).

    If cache_path is given and exists, the cached (unaggregated, gene-filtered)
    frame is loaded instead of rescanning -- lets repeated figure-layout runs
    skip the ~1-2 min full-file scan; the aggregation/curation logic downstream
    of this function is expected to change often during layout iteration, this
    scan itself is not.
    """
    if cache_path and os.path.exists(cache_path):
        return pd.read_parquet(cache_path)

    value_set = pa.array(genes, type=pa.string())
    pf = pq.ParquetFile(path)
    tables = []
    for rg in range(pf.num_row_groups):
        col = pf.read_row_group(rg, columns=[gene_col])[gene_col]
        mask = pc.is_in(col, value_set=value_set)
        if pc.any(mask).as_py():
            t = pf.read_row_group(rg)
            tables.append(t.filter(pc.is_in(t[gene_col], value_set=value_set)))
    df = pa.concat_tables(tables).to_pandas() if tables else pd.DataFrame()
    print(f"  scanned {path}: {len(df)} rows, peak_rss={peak_rss_gb():.2f}GB")

    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        df.to_parquet(cache_path)
    return df


def aggregate_whole_body(df, cell_types=None):
    """Two-step donor-equal-weighted aggregation (binding methodology for this
    project, see query_antxr2_table.py and ANALYSIS_SUMMARY.md -- do not replace
    with a single-step n_cells-weighted pool):
    1) per-donor mean, n_cells-weighted across that donor's own
       collection_name x dataset_id x assay rows, per cell_type x gene.
    2) unweighted average of those donor-level means across donors, per
       cell_type x gene.
    """
    if cell_types is not None:
        df = df[df["cell_type"].isin(cell_types)]

    donor_key = ["cell_type", "gene", "collection_name", "dataset_id", "donor_id"]
    df = df.copy()
    df["_w"] = df["mean_expression"] * df["n_cells"]
    donor_level = df.groupby(donor_key, as_index=False).agg(
        donor_mean=("_w", "sum"), donor_n_cells=("n_cells", "sum")
    )
    donor_level["donor_mean"] = donor_level["donor_mean"] / donor_level["donor_n_cells"]

    agg = donor_level.groupby(["cell_type", "gene"], as_index=False).agg(
        mean_expression=("donor_mean", "mean"),
        n_donors=("donor_id", "nunique"),
        n_cells=("donor_n_cells", "sum"),
    )
    return agg


def aggregate_skin(df, subtypes=None):
    """Healthy/nonlesional filter + simple mean per fibroblast_subtype x gene.
    Each input row is already at sample_id granularity, so this directly
    implements sample-equal-weighting (no further donor-pooling step needed) --
    the same pattern used repeatedly earlier this session.
    """
    healthy = df[(df["patient_status"] == "Healthy") & (df["lesional_status"] == "Nonlesional")]
    if subtypes is not None:
        healthy = healthy[healthy["fibroblast_subtype"].isin(subtypes)]

    agg = healthy.groupby(["fibroblast_subtype", "gene"], as_index=False).agg(
        mean_expression=("mean_expression", "mean"),
        n_samples=("sample_id", "nunique"),
        n_cells=("n_cells", "sum"),
    )
    return agg
