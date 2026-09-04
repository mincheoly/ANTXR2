"""Rebuild combined_celltype_means.parquet by streaming-concatenating the
per-collection parquets.

Needed because compute_means.py writes the per-collection file and the combined
file in the same pass: re-running a subset of datasets with --skip-combined
leaves `combined` stale, and re-running it *without* --skip-combined would
truncate `combined` to only the datasets in that run. This script rebuilds it
from the three per-collection files, which are the authoritative per-dataset
outputs.

Streams row group by row group -- the combined file is ~3 GB / 351M rows and
must never be materialised in memory (a whole-table concat is what OOM-killed
the original pipeline run; see compute_means.py's ParquetWriterRegistry).

Usage: python rebuild_combined.py [--dry-run]
"""
import argparse
import os

import pyarrow.parquet as pq

from config import OUTPUT_DIR

COLLECTIONS = ["tabula_sapiens", "cross_tissue_immune", "gut_cell_atlas"]
COMBINED = os.path.join(OUTPUT_DIR, "combined_celltype_means.parquet")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    parts = [os.path.join(OUTPUT_DIR, f"{c}_celltype_means.parquet") for c in COLLECTIONS]
    missing = [p for p in parts if not os.path.exists(p)]
    if missing:
        raise SystemExit(f"missing per-collection parquet(s): {missing}")

    pfs = [pq.ParquetFile(p) for p in parts]
    ref = pfs[0].schema_arrow
    for p, pf in zip(parts, pfs):
        if not pf.schema_arrow.equals(ref):
            raise SystemExit(
                f"schema mismatch in {p} -- all collections must share one schema "
                "before they can be concatenated into combined"
            )
    total = sum(pf.metadata.num_rows for pf in pfs)
    for p, pf in zip(parts, pfs):
        print(f"  {os.path.basename(p):46s} {pf.metadata.num_rows:>12,} rows")
    print(f"  {'-> combined':46s} {total:>12,} rows")

    if args.dry_run:
        print("[dry run] not writing")
        return

    tmp = COMBINED + ".tmp"
    writer = pq.ParquetWriter(tmp, ref)
    written = 0
    try:
        for p, pf in zip(parts, pfs):
            for rg in range(pf.num_row_groups):
                t = pf.read_row_group(rg)
                writer.write_table(t)
                written += t.num_rows
            print(f"  wrote {os.path.basename(p)} ({written:,} rows so far)")
    finally:
        writer.close()

    if written != total:
        os.remove(tmp)
        raise SystemExit(f"row count mismatch: wrote {written:,}, expected {total:,}")
    os.replace(tmp, COMBINED)
    print(f"\nwrote {COMBINED} ({written:,} rows, "
          f"{os.path.getsize(COMBINED)/1e9:.2f} GB)")


if __name__ == "__main__":
    main()
