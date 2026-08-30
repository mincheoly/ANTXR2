"""Build the per-cell-type ANTXR2 table requested by the user:
mean (memento, donor-averaged), n_cells, n_donors, presence per collection.

Methodology (per user's earlier correction, must persist): two-step aggregation --
1) donor-level mean first (n_cells-weighted across a donor's own
   collection/dataset/assay entries for that cell_type)
2) simple (unweighted) average of those donor-level means, per cell_type

Reads only the ANTXR2 rows via row-group filtered parquet scan (gene column) to
avoid loading the full 351M-row / 3GB file into memory.
"""
import resource

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

PATH = "/data/ANTXR2/celltype_expression/combined_celltype_means.parquet"
GENE = "ANTXR2"

pf = pq.ParquetFile(PATH)
gene_col_idx = pf.schema_arrow.get_field_index("gene")

tables = []
for rg in range(pf.num_row_groups):
    col = pf.read_row_group(rg, columns=["gene"])["gene"]
    mask = pc.equal(col, GENE)
    if pc.any(mask).as_py():
        t = pf.read_row_group(rg)
        t = t.filter(pc.equal(t["gene"], GENE))
        tables.append(t)

table = pa.concat_tables(tables)
df = table.to_pandas()
print(f"ANTXR2 rows: {len(df)}")
print(f"peak RSS so far: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")

# --- step 1: donor-level mean (n_cells-weighted across that donor's own
# collection/dataset/assay entries for a given cell_type) ---
donor_key = ["cell_type", "collection_name", "dataset_id", "donor_id"]
df["_w"] = df["mean_expression"] * df["n_cells"]
donor_level = (
    df.groupby(donor_key, as_index=False)
    .agg(
        donor_mean=("_w", "sum"),
        donor_n_cells=("n_cells", "sum"),
    )
)
donor_level["donor_mean"] = donor_level["donor_mean"] / donor_level["donor_n_cells"]

# --- step 2: simple (unweighted) average of donor-level means, per cell_type ---
agg = (
    donor_level.groupby("cell_type", as_index=False)
    .agg(
        mean_expression=("donor_mean", "mean"),
        n_donors=("donor_id", "nunique"),
        n_cells=("donor_n_cells", "sum"),
    )
)

# --- presence per collection (boolean flags, from the raw per-row table) ---
presence = (
    df.groupby("cell_type")["collection_name"]
    .apply(lambda s: set(s.unique()))
    .rename("collections")
    .reset_index()
)
for coll in ["tabula_sapiens", "cross_tissue_immune", "gut_cell_atlas"]:
    presence[f"in_{coll}"] = presence["collections"].apply(lambda s, c=coll: c in s)
presence = presence.drop(columns=["collections"])

out = agg.merge(presence, on="cell_type", how="left")
out = out.sort_values("mean_expression", ascending=False).reset_index(drop=True)
out.insert(0, "rank", out.index + 1)

out_path = "/tmp/claude-1000/-home-ubuntu-Github-ANTXR2/0d096a79-bc2d-474d-8436-d6d77a55e88a/scratchpad/antxr2_celltype_table.csv"
out.to_csv(out_path, index=False)
print(f"wrote {out_path} ({len(out)} rows)")
print(out.head(15).to_string(index=False))
print(f"final peak RSS: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")
