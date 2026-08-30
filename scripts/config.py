"""Shared configuration for the ANTXR2 cell-type expression pipeline (Phase 1).

Dataset selection rationale (see /home/ubuntu/Github/ANTXR2/prompts/celltype_means.md,
"Design decisions" section, and the clarification round that preceded this pipeline):

Each CELLxGene collection's Curation API listing returns many redundant sub-datasets
(organ/lineage/compartment slices) alongside one canonical dataset covering the same
cells (typically the one with is_primary_data=True). We fetch only the canonical
dataset(s) per collection to avoid double-counting the same physical cells under many
different dataset_id groups. For the Gut Cell Atlas specifically, "Healthy reference"
(genome-wide-ish, healthy only) and "Extended+ - 18485 genes" (reduced panel, but the
only one with the UC/Crohn's/celiac disease samples) are both kept as separate
dataset_id rows, per user decision -- donors likely overlap between the two, which is
acceptable since dataset_id stays a distinct column throughout and rows are never
silently collapsed across datasets.
"""

import os

RAW_DIR = "/data/ANTXR2/raw"
OUTPUT_DIR = "/data/ANTXR2/celltype_expression"

CURATION_API_BASE = "https://api.cellxgene.cziscience.com/curation/v1"

# --- Skin fibroblast atlas (Steele et al., Nat. Immunol. 2025) ---
# See prompts/skin_means.md. Not on CELLxGene Discover -- the plan flagged this as
# a download blocker requiring a stable accession rather than a guessed/scraped URL.
# Resolution: the paper's own web portal (https://cellatlas.io/studies/skin-fibroblast)
# is a JS-rendered SPA backed by a public Strapi CMS API
# (https://strapi-api-dot-haniffa-lab.nw.r.appspot.com/api/). Querying that API for the
# "skin-fibroblast" study (GET /api/studies?filters[slug][$eq]=skin-fibroblast) surfaces
# its "Integrated atlas" dataset record (357,276 cells -- matches the paper abstract's
# stated n exactly), which lists a `resources` entry named "AnnData h5ad" pointing at
# this URL: a public (unsigned, CORS-open, no-auth) GCS object, verified via HEAD
# (27,231,223,312 bytes, last-modified 2024-12-04). This is the literal file the
# official portal itself serves for this study, not a guessed or scraped URL -- the
# scrape was of the CMS API, not the rendered page. A Zarr version of the same object
# also exists at .../skin-fibroblast/zarr/adata_webportal.zarr (what the live portal's
# viewer actually streams from) but h5ad matches the rest of this pipeline's h5py-based
# reading.
SKIN_FIBROBLAST_URL = "https://storage.googleapis.com/haniffalab/skin-fibroblast/adata_webportal.h5ad"
SKIN_FIBROBLAST_DIR = os.path.join(RAW_DIR, "skin_fibroblast_atlas")
SKIN_FIBROBLAST_H5AD = os.path.join(SKIN_FIBROBLAST_DIR, "adata_webportal.h5ad")
SKIN_FIBROBLAST_EXPECTED_SIZE = 27231223312
SKIN_FIBROBLAST_EXPECTED_MD5 = "e2cf64a7d04d6afe87f9278d046c3a46"  # from GCS ETag
SKIN_FIBROBLAST_OUTPUT_PATH = os.path.join(OUTPUT_DIR, "skin_fibroblast_celltype_means.parquet")

# obs schema of the downloaded object, verified at inspection time (see
# ANALYSIS_SUMMARY.md for the full inspection notes and why each choice was made):
#
# - Raw integer counts live at raw/X (equivalently layers/counts; X is
#   log-normalized). No feature_id/feature_name split -- var's index is the gene
#   symbol directly.
# - There is no donor_id/patient_id *column* in obs, but per-sample identity is not
#   actually absent -- it's encoded in the cell barcode (obs index) itself, appended
#   during integration to keep barcodes unique across the ~30 concatenated source
#   studies, and never split back out into its own column. A first pass (naive
#   global "last underscore token" split) wrongly concluded this wasn't reliably
#   parseable; per source study the format is almost always a clean, uniform
#   `<barcode>_<GSM or SRS accession>` (one real GEO/SRA sample per suffix). See
#   `extract_sample_id` in compute_skin_means.py for the per-study parsing rules and
#   SKIN_MESSY_STUDIES below for the two studies where recovery wasn't attempted.
#   The recovered id is stored as `sample_id` (the real per-sample/per-donor grouping
#   axis and chunk boundary) alongside `source_study` (=GSE/author label, coarser,
#   kept for provenance).
# - `celltype` (F1-F7, e.g. "F2: Universal", "F7: Fascia-like myofibroblast") is the
#   complete (no-NaN) fine subtype label and is the primary grouping column.
#   `celltype_skinspecific_nomenclature` is a parallel descriptive name that is 1:1
#   with `celltype` (verified) except the two disease-specific F6/F7 subtypes, which
#   have no skin-specific name (NaN) -- safe to carry as a derived, non-grouping
#   column.
# - There is no single "disease_state" column: `Patient_status` (fine diagnosis,
#   e.g. "Healthy", "Psoriasis", "Dupuytren's") and `disease_category_orig` (coarser
#   bucket) are NOT 1:1 (e.g. Psoriasis spans both "Inflammatory without scarring"
#   and "Nonlesional"), and `lesional_vs_nonlesional` (Lesional/Nonlesional/PostRx)
#   is a third, separate axis. Per the plan's "do not pool healthy and lesional"
#   instruction, all three are kept as independent grouping columns rather than
#   collapsed into one.
# - No assay column exists (technology varies across the ~30 integrated source
#   studies but isn't preserved in this export) -- every group uses the flat
#   DEFAULT_CAPTURE_RATE placeholder; the output `assay` column is always "unknown".
SKIN_RAW_SLOT = "raw/X"
SKIN_STUDY_COL = "GSE"
SKIN_LABEL_COLS = [
    "celltype", "Patient_status", "lesional_vs_nonlesional", "disease_category_orig",
]
SKIN_DERIVED_COL = "celltype_skinspecific_nomenclature"  # 1:1 with celltype, attached post-hoc

# Studies where per-sample barcode-suffix recovery was not attempted, per user
# decision -- coarse GSE-level grouping used instead for these two only:
# - "Ganier": ~90% of cells have a clean `WS_SKN_KCLxxxxxxx` sample code, but ~10%
#   have a different (barcode-order-reversed) format that would need bespoke
#   handling to salvage.
# - "Sole-Boldo": the barcode suffix is a body-site label ("body_solebordo"), not a
#   per-donor id -- not recoverable from the barcode at all.
SKIN_MESSY_STUDIES = {"Ganier", "Sole-Boldo"}

# collection_id -> (subdirectory name, [dataset titles to download])
COLLECTIONS = {
    "e5f58829-1a66-40b5-a624-9046778e74f5": {
        "name": "tabula_sapiens",
        "dataset_titles": ["Tabula Sapiens - All Cells"],
    },
    "62ef75e4-cbea-454e-a0ce-998ec40223d3": {
        "name": "cross_tissue_immune",
        "dataset_titles": ["Global"],
    },
    "f11cb29c-b546-4738-9bd8-66ea621a7bd5": {
        "name": "gut_cell_atlas",
        "dataset_titles": ["Healthy reference", "Extended+ - 18485 genes"],
    },
}

# --- memento parameters ---

# RNA capture efficiency (memento's "q") per assay. memento's own example wrappers
# (run_eqtl, binary_test_1d in the upstream package) use a flat capture_rate=0.1 for
# 10x droplet data without per-chemistry-version tuning, and we couldn't recover a
# more authoritative per-platform table from the paper/docs at pipeline-build time
# (paywalled / not indexed). We use that same 0.1 default for every assay here as a
# literature-typical placeholder, NOT an independently calibrated value -- flagged in
# the output README. Override per-assay below if better estimates become available;
# assay is kept as its own column throughout so re-deriving means with revised q
# values later doesn't require re-grouping.
DEFAULT_CAPTURE_RATE = 0.1
CAPTURE_RATE_BY_ASSAY = {}

# Below this many cells, a donor-chunk read is skipped and its cells fall back to
# the dataset-level (no-donor) grouping bucket instead of being dropped.
MIN_DONOR_CHUNK_CELLS = 1

# Split any single donor's cells into sub-chunks if it exceeds this, to bound peak
# memory during memento's per-chunk setup (sub-chunk results for the same
# donor/cell_type/assay group are recombined via a cell-count-weighted average).
MAX_CHUNK_CELLS = 150_000

# memento setup_memento(min_cell_count=...): memento's own compute_1d_moments fits a
# per-group mean-variance regressor using genes with max count >= 2 *within that
# group*; for groups under ~10 cells this selection is frequently empty (no gene
# reaches count 2), which crashes np.polyfit. memento's own default of 10 exists to
# guard against exactly this, so we keep it rather than pushing lower to "flag, don't
# hide" every group down to n=1 -- groups this small aren't a well-defined
# method-of-moments estimate in the first place. Groups with 10-19 cells still get
# through and are marked low_confidence (see LOW_N_THRESHOLD in compute_means.py).
MEMENTO_MIN_CELL_COUNT = 10

# Standard CZI schema obs column names (verified per-dataset at runtime; see
# compute_means.py's schema check).
CELL_TYPE_COL = "cell_type"
CELL_TYPE_ID_COL = "cell_type_ontology_term_id"
DONOR_COL = "donor_id"
TISSUE_COL = "tissue"
TISSUE_GENERAL_COL = "tissue_general"
ASSAY_COL = "assay"
DISEASE_COL = "disease"

GENE_ID_COL = "feature_id"
GENE_NAME_COL = "feature_name"

