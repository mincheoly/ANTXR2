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

# ANTXR1/ANTXR2 + the ECM-clearance-pathway gene panel used together throughout
# this project's analysis (endocytosis/lysosomal-degradation partners for the
# collagen VI clearance ANTXR1/ANTXR2 are implicated in). Fixed scope -- do not
# fold in the separate activation-axis panel (ACTA2/POSTN/TGFB1), which answers a
# different question and was never part of this list.
GENE_PANEL = ["ANTXR1", "ANTXR2", "MRC2", "CTSB", "CTSK", "MMP14", "TIMP2", "LAMP1"]

# --- Two-arm partner-availability panel (prompts/partner_availability.md, Task 3) ---
#
# ANTXR2/CMG2 has two separable functions with different partner requirements:
# collagen-VI clearance (MRC2 + lysosomal machinery + substrate) and Wnt signal
# transduction (LRP6 + Frizzled). GENE_PANEL above tests only the clearance
# arm's machinery -- not its substrate, and not the Wnt arm at all. Kept as a
# SEPARATE dict of named sub-panels (not folded into GENE_PANEL) so the
# original 8-gene panel and its heatmap stay untouched; EXTENDED_GENE_PANEL
# below is the flat query list derived from it.
GENE_PANEL_ARMS = {
    "RECEPTORS": ["ANTXR1", "ANTXR2"],
    "CLEARANCE_ARM": ["MRC2", "CTSB", "CTSK", "MMP14", "TIMP2", "LAMP1"],
    # the clearance arm's substrate -- the panel tested everything except
    # whether collagen VI itself is present.
    "CLEARANCE_SUBSTRATE": ["COL6A1", "COL6A2", "COL6A3"],
    # LRP6 + Frizzled assembly (Wnt receptor complex), downstream beta-catenin
    # transduction/destruction-complex components, and injury-induced stem-cell
    # renewal markers (Bracq et al. 2025).
    "WNT_ARM": (
        ["LRP5", "LRP6"] + [f"FZD{i}" for i in range(1, 11)]
        + ["CTNNB1", "TCF7L2", "LGR5", "RNF43", "ZNRF3", "AXIN2"]
    ),
}
EXTENDED_GENE_PANEL = [g for arm in GENE_PANEL_ARMS.values() for g in arm]

# --- IBD colon atlas (Smillie et al., Cell 2019, "Intra- and Inter-cellular
# Rewiring of the Human Colon during Ulcerative Colitis") ---
#
# SCP259 on the Broad Single Cell Portal is login-gated: the study page states
# "Please sign in to download data" and its study_files API returns HTTP 401
# (verified at pipeline-build time). Every other dataset in this project had an
# unauthenticated public download URL; this one does not, on its primary host.
#
# Resolution: the identical file set is also hosted on the Human Cell Atlas Data
# Coordination Platform (HCA DCP) as project cd61771b-661a-4e19-b269-6e5d95350de6
# ("HumanColonRewiringUlcerativeColitis"), dataUseRestriction="NRES" (no
# restriction) and no duosId set -- i.e. genuinely open access, not merely listed.
# Verified end-to-end with plain unauthenticated requests: HCA's Azul REST API
# (https://service.azul.data.humancellatlas.org) resolves each file to a signed,
# no-auth-required S3 URL via a 302 redirect from
# GET /repository/files/{file_uuid}?catalog=dcp60&version={version}, and both a
# tiny file (Imm.genes.tsv) and a >1GB file (gene_sorted-Epi.matrix.mtx, with an
# HTTP Range request) round-tripped successfully this way -- no Broad/Google
# login involved anywhere. Catalog is "dcp60" (the current default catalog at
# resolution time; Azul catalogs are versioned/rotated over time by HCA, so
# download_ibd_colon_data.py re-fetches the live default catalog at runtme rather
# than hardcoding "dcp60", the same "IDs stable, resolve URLs at runtime" pattern
# already used by download_data.py for CELLxGene collection IDs).
#
# One file (cell_subsets.txt, the fine->coarse 51-subset lineage map) is not part
# of the HCA file set -- it comes from the paper authors' own analysis repo
# (github.com/cssmillie/ulcerative_colitis), which is also the source of truth
# for all raw file names below (from that repo's run.r).
#
# Not this dataset: GEO GSE114374 ("Composition of the Colonic Mesenchyme...") is
# a related but DISTINCT companion study (Kinchen et al. 2018, 10 samples,
# human+mouse) -- do not substitute it for SCP259/this HCA project.
AZUL_BASE = "https://service.azul.data.humancellatlas.org"
IBD_COLON_HCA_PROJECT_ID = "cd61771b-661a-4e19-b269-6e5d95350de6"
IBD_COLON_DIR = os.path.join(RAW_DIR, "ibd_colon_atlas")
IBD_COLON_COMPARTMENTS = ["Epi", "Fib", "Imm"]

# Raw per-compartment files (10x-like triplet, Matrix Market format -- NOT h5ad,
# unlike every other dataset in this project). Exact names per compartment, from
# cssmillie/ulcerative_colitis/run.r:
IBD_COLON_MATRIX_FILES = {comp: f"gene_sorted-{comp}.matrix.mtx" for comp in IBD_COLON_COMPARTMENTS}
IBD_COLON_GENES_FILES = {comp: f"{comp}.genes.tsv" for comp in IBD_COLON_COMPARTMENTS}
IBD_COLON_BARCODES_FILES = {comp: f"{comp}.barcodes2.tsv" for comp in IBD_COLON_COMPARTMENTS}
# Pre-computed Seurat objects (discovery-cohort subset, NOT the full 30-donor
# cohort -- exact subset size TBD once metadata is inspected), with tSNE already
# computed. Kept as the "processed with visualization coordinates" artifact.
IBD_COLON_SEURAT_RDS_FILES = {comp: f"train.{comp}.seur.rds" for comp in IBD_COLON_COMPARTMENTS}
IBD_COLON_META_FILE = "all.meta2.txt"          # shared per-cell metadata, all compartments
IBD_COLON_SUBSETS_FILE = "cell_subsets.txt"    # fine (51) -> coarse lineage map; from GitHub, not HCA
IBD_COLON_SUBSETS_URL = "https://raw.githubusercontent.com/cssmillie/ulcerative_colitis/master/cell_subsets.txt"

# Kept as 3 SEPARATE per-compartment h5ad files rather than one merged object:
# verified post-download that the 3 compartments' gene panels are NOT identical
# (Epi 20,028 genes, Fib 19,076, Imm 20,529; ~90-95% pairwise overlap, each
# compartment independently gene-filtered by the original authors) and their
# barcodes are disjoint (a cell belongs to exactly one compartment's matrix) --
# so there is no natural shared `var` to merge on without inventing fake zeros
# for genes a compartment's own panel never included. compute_ibd_colon_means.py
# instead loops over these 3 files the same way compute_means.py loops over
# manifest.csv rows, writing all compartments to the ONE shared
# IBD_COLON_OUTPUT_PATH parquet with `compartment` as a row-level column --
# "one dataset = one output parquet" is preserved at the output layer even
# though raw loading stays per-compartment.
IBD_COLON_H5AD_BY_COMPARTMENT = {
    comp: os.path.join(IBD_COLON_DIR, f"ibd_colon_atlas_{comp}.h5ad")
    for comp in IBD_COLON_COMPARTMENTS
}
IBD_COLON_RAW_SLOT = "X"  # these mtx files are already raw counts, single slot, no raw/X duality
IBD_COLON_OUTPUT_PATH = os.path.join(OUTPUT_DIR, "ibd_colon_atlas_celltype_means.parquet")

# obs schema of all.meta2.txt (verified at inspection time): NAME (barcode,
# joins to each compartment's barcodes2.tsv), Cluster (the 51-subset fine cell
# type -- joins 1:1 into cell_subsets.txt's fine label), nGene, nUMI (per-cell
# QC, not used), Subject (donor id -- verified clean: exactly 30 distinct
# values, each mapping to either 1 Health value (12 "Healthy" donors) or 2
# (18 UC donors with paired Inflamed/Non-inflamed biopsies) -- no collision/
# alias pattern found, so no *_UNIFIED_COL correction was needed here, unlike
# the Gut Cell Atlas), Health (Healthy/Inflamed/Non-inflamed), Location
# (Epi/LP -- an anatomical/tissue-layer label, NOT the same axis as
# `compartment`: LP cells are split across the Fib and Imm compartment
# matrices by cell type, so Location and compartment are kept as two separate
# columns, never collapsed). All 365,492 metadata rows matched cleanly against
# the union of all 3 compartments' barcodes (123,006 + 31,872 + 210,614 =
# 365,492 exactly) -- verified, not assumed.
IBD_COLON_META_NAME_COL = "NAME"
IBD_COLON_CLUSTER_COL = "Cluster"
IBD_COLON_DONOR_COL = "Subject"
IBD_COLON_HEALTH_COL = "Health"
IBD_COLON_LOCATION_COL = "Location"
IBD_COLON_SAMPLE_COL = "Sample"

# --- HCA source-data corruption (2026-09-09) -- RESOLVED, kept for provenance ---
#
# Two of the three HCA-hosted raw matrices were found truncated AT THE SOURCE
# (independently re-probed via HTTP Range requests directly against the
# resolved S3 URL, bypassing our own download code entirely -- the S3 object
# itself reported a total size matching the truncated download exactly, not a
# network/resume artifact): gene_sorted-Epi.matrix.mtx had 52.6% of its
# declared data lines, gene_sorted-Imm.matrix.mtx had 4.6%; gene_sorted-Fib
# was independently verified complete. Worked around at the time via
# export_ibd_colon_rds.R, which pulled counts/meta.data/tsne out of the
# discovery-cohort train.{Epi,Fib,Imm}.seur.rds objects (17 of 30 donors) into
# IBD_COLON_RDS_EXPORT_DIR, readable without installing Seurat (a Seurat
# object's S4 slots are plain attributes; only the already-installed `Matrix`
# package is needed for the dgCMatrix counts slot).
#
# SUPERSEDED same day: the user created an SCP259 account and provided a
# bulk-download auth code, giving genuinely complete, full-30-donor-cohort
# copies of all 3 matrices plus real full-cohort tSNE coordinates
# ({Epi,Fib,Imm}.tsne.txt) -- each verified complete against its own mtx
# header / expected row count before use (see
# build_ibd_colon_h5ad.verify_mtx_complete). These now live at the same raw
# paths (gene_sorted-*.matrix.mtx), and build_ibd_colon_h5ad.py reads them
# directly for all 3 compartments -- the RDS-export fallback is no longer
# consumed by that script. The corrupted HCA copies are kept at
# raw/ibd_colon_atlas/corrupted_hca_backup/ (not deleted) and the RDS export
# machinery is left in place (export_ibd_colon_rds.R,
# IBD_COLON_RDS_EXPORT_DIR, the train.*.seur.rds files) as a documented
# fallback if a future dataset needs the same recovery pattern -- see
# ANALYSIS_SUMMARY.md for the full narrative of both the corruption finding
# and the full-cohort resolution.
#
# The RDS metadata's Health column used "Uninflamed" where all.meta2.txt uses
# "Non-inflamed" for the identical concept -- IBD_COLON_HEALTH_MAP was the
# harmonization for that RDS-only path; unused now that all 3 compartments
# read Health from all.meta2.txt directly, kept only for the (still-valid,
# just no longer exercised) RDS export.
IBD_COLON_RDS_EXPORT_DIR = os.path.join(IBD_COLON_DIR, "rds_export")
IBD_COLON_HEALTH_MAP = {"Uninflamed": "Non-inflamed"}  # RDS-export path only, not used by the current full-cohort build

FIGURES_DIR = "/data/ANTXR2/figures"

# --- Phase 2 prep: single-cell co-expression working set -----------------------
#
# Source is the Gut Cell Atlas "Extended+" release (1,596,200 x 18,370). We carve
# out one study (Elmentaite2021) and three cell types into a small standalone
# h5ad -- see /data/ANTXR2/coexpression/README_<variant>.md for the full rationale.
#
# Why Elmentaite2021 specifically: 398,460 cells, 41 donors, 100%
# Non_pathological samples, and gut epithelium + fibroblasts + true myeloid cells
# all drawn from the SAME 36 donors, one lab, one protocol, with 5'/3' chemistry
# roughly balanced. That removes donor, cohort, and protocol confounds from the
# three-way comparison in one move -- which the alternative (gut epithelium here,
# fibroblasts from the skin atlas, myeloid from a third cohort) could not.
GUT_ATLAS_H5AD = os.path.join(
    RAW_DIR, "gut_cell_atlas", "19053a82-9c89-4fb8-bd19-d7b1800b0b7b.h5ad"
)
COEXPR_DIR = "/data/ANTXR2/coexpression"
COEXPR_H5AD = os.path.join(COEXPR_DIR, "elmentaite2021_trio.h5ad")
COEXPR_STUDY = "Elmentaite2021"

# Two working sets are built from the same study, differing only in which
# annotation column defines the groups. Select with
# `subset_coexpression_data.py --variant {celltype,level3}`.
#
# "celltype" uses the CELLxGene-harmonised `cell_type`, which pools author
# subtypes -- e.g. its `macrophage` merges Macrophage + Macrophage_LYVE1 +
# _TREM2 + _MMP9 + _CD5L, and its `fibroblast` merges every fibroblast subtype.
# "level3" uses the authors' own `level_3_annot` for cleaner, narrower
# populations at the cost of cell count. The level3 macrophage is visibly purer
# (COL1A1 0.016 vs 0.030 in the pooled version) and the crypt fibroblast is a
# single defined subtype rather than a mixture.
COEXPR_VARIANTS = {
    "celltype": {
        "label_col": "cell_type",
        "labels": ["fibroblast", "enterocyte", "macrophage"],
        "h5ad": os.path.join(COEXPR_DIR, "elmentaite2021_trio.h5ad"),
    },
    "level3": {
        "label_col": "level_3_annot",
        # NB: the fibroblast label is PI16 (a fibroblast marker gene), not PI1.
        # Plain `Macrophage` deliberately excludes the LYVE1/TREM2/MMP9/CD5L
        # subtypes, which are separate level_3 labels.
        "labels": ["Crypt_fibroblast_PI16", "Enterocyte", "Macrophage"],
        "h5ad": os.path.join(COEXPR_DIR, "elmentaite2021_trio_level3.h5ad"),
    },
}

# The three cell types, each playing a distinct role against the Phase 1 finding
# that gut epithelium has ANTXR2 but neither ANTXR1 nor MRC2:
#   fibroblast  -- canonical case: ANTXR1 + ANTXR2 + MRC2 all present
#   enterocyte  -- the unexplained case
#   macrophage  -- myeloid calibration control (clearance genes canonically high)
#
# NOTE: `macrophage` here is the plain label, NOT "gastrointestinal tract (lamina
# propria) macrophage". The latter is a CELLxGene curation error: all 42,302 cells
# carrying it across the whole Gut Cell Atlas are authored as
# level_3_annot=Lamina_propria_fibroblast_ADAMDEC1 (Mesenchymal/Fibroblast), and
# marker checks agree -- PTPRC 0.4% positive, CD68 2.5%, C1QA 0.5%, but COL1A1
# 89% positive. They are fibroblasts. Do not substitute that label back in.
COEXPR_CELL_TYPES = ["fibroblast", "enterocyte", "macrophage"]

# --- capture rate (memento's `q`) ----------------------------------------------
#
# Two independently-derived estimates are attached as per-cell obs columns, so
# downstream work can test sensitivity to the choice rather than inheriting one
# unexamined number. memento accepts per-cell q via `q_column` and asserts
# max(q) < 1 (memento/main.py:88).
#
# q_chem: memento's own broad droplet rate, from their CELLxGene pipeline
# (publication/cellxgene/make_cube.py: `Q = 0.07  # RNA capture efficiency
# depending on technology`) -- the same value used in their PBMC *co-expression*
# analysis (publication/other/ifn_pbmc/interferon_2d.py). Their wrappers.py
# default is 0.1; 0.07 is the more specific choice and the one used for 2D work.
# Kept as a per-assay dict so the two chemistries can diverge if a better source
# appears -- but both are 0.07 today because NO published per-chemistry q exists.
# This is deliberately not invented: the 5'-vs-3' difference is instead measured
# empirically by the PBMC-scaled column.
Q_CHEM_DEFAULT = 0.07
Q_CHEM_BY_ASSAY = {
    "10x 5' v2": 0.07,
    "10x 3' v2": 0.07,
}

# ASSUMPTION, not a measurement. The paper never states sequencing saturation;
# ENA has read_count=0 for all 89 runs of E-MTAB-9543 (the adult GEX accession),
# so it cannot be recovered for these samples. Only E-MTAB-8901 (developing gut,
# HiSeq 4000, median 352M read pairs/run at ~8,000 cells/reaction -> ~44K
# reads/cell) supports an estimate, and it lands at ~80-85%. Change this one
# constant and re-run to test sensitivity.
ASSUMED_SATURATION = 0.85

# 10x public PBMC references, chemistry-matched to our two assays, used to derive
# the second capture-rate column by comparing UMI depth in matched immune gates.
# pbmc8k is the canonical 3' v2 healthy-donor PBMC run (8,381 cells, 93,552 mean
# reads/cell); sc5p_v2_hs_PBMC_10k is its 5' v2 counterpart.
PBMC_REFERENCE_DIR = os.path.join(RAW_DIR, "pbmc_reference")
PBMC_REFERENCES = {
    "10x 3' v2": {
        "name": "pbmc8k",
        "url": "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz",
        "kind": "mtx_tar",
        "filename": "pbmc8k_filtered_gene_bc_matrices.tar.gz",
    },
    "10x 5' v2": {
        "name": "sc5p_v2_hs_PBMC_10k",
        "url": "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5",
        "kind": "h5",
        "filename": "sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5",
    },
}

# Marker gates used identically on both our data and the PBMC references, so the
# comparison can't drift between them. Deliberately simple positive gating rather
# than clustering -- we only need a population whose UMI depth is comparable, not
# a definitive annotation.
#
# Gate markers must work in BOTH blood and gut tissue, which rules out some
# obvious choices. Notably CD68 was tried and rejected: it is 52% positive with
# mean 3.10 in gut *enterocytes* (while LYZ/C1QA/AIF1/PTPRC are all <1% there),
# because CD68/macrosialin is a lysosomal glycoprotein rather than a
# macrophage-specific marker, and enterocytes have very active endolysosomal
# compartments. LYZ + AIF1 is pan-myeloid in both compartments; EPCAM- guards
# against epithelial carry-in on the tissue side (a no-op for PBMC).
PBMC_GATES = {
    "T": {"pos": ["CD3E"], "neg": []},
    "B": {"pos": ["MS4A1"], "neg": []},
    "NK": {"pos": ["NKG7"], "neg": ["CD3E"]},
    "Mono/Mac": {"pos": ["LYZ", "AIF1"], "neg": ["EPCAM"]},
}

# The gates above must be evaluated on this study's IMMUNE compartment, not on
# the three-cell-type working set -- that set is fibroblast/enterocyte/macrophage
# and contains essentially no lymphocytes, so a T or B gate applied to it would
# match only doublets and ambient RNA. Capture rate is a property of the
# assay/experiment, so it is estimated from immune cells in the same study and
# then applied to the working set's cells.
COEXPR_IMMUNE_LEVEL1 = ["T and NK cells", "B and B plasma", "Myeloid"]

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

# Preferred donor identity, used instead of DONOR_COL wherever the source file
# provides it (currently only the two Gut Cell Atlas files; Tabula Sapiens and
# the Cross-tissue Immune Atlas carry `donor_id` only).
#
# This is not cosmetic. In the gut atlas `donor_id` is wrong in BOTH directions:
#   - collision: `A25` and `A34 (417C)` each map to TWO different people
#     (e.g. A34 (417C) -> D11, 55-74y, 31,370 cells AND D12, 18-34y, 5,814
#     cells -- different age brackets, unambiguously two donors);
#   - aliases: 34 people appear under multiple `donor_id` strings in Extended+
#     (9 in Healthy reference), e.g. D2 -> T036 / T036NEG / T036POS, which are
#     sorted fractions of one donor.
# Net effect in Elmentaite2021: 41 `donor_id` values for 38 actual people.
#
# Grouping by `donor_id` therefore both merges two donors into one group and
# splits single donors across several -- the latter is pseudo-replication, and
# it would propagate straight into memento's donor-level bootstrapping.
DONOR_UNIFIED_COL = "donorID_unified"
TISSUE_COL = "tissue"
TISSUE_GENERAL_COL = "tissue_general"
ASSAY_COL = "assay"
DISEASE_COL = "disease"

GENE_ID_COL = "feature_id"
GENE_NAME_COL = "feature_name"

# --- Phase 2: co-expression analysis (memento point-estimate correlations) -----
#
# Test run on the "level3" trio working set (see COEXPR_VARIANTS above), ahead of
# scaling to the full dataset. Point estimates only -- no bootstrap/hypothesis
# testing (memento's ht_*/binary_test_* functions are never called here).
#
# Groups are donor x cell_type, using the trio h5ad's `donor` obs column, which
# already equals the corrected DONOR_UNIFIED_COL identity (not `donor_id`, which
# has the collision/alias problems described above).
COEXPR_DONOR_COL = "donor"
COEXPR_TARGET_GENE = "ANTXR2"

# Minimum cells required per donor x cell_type group. Below this, the group is
# dropped and reported rather than silently run -- macrophage is small in this
# variant (2,953 cells / 38 donors), so expect a meaningful fraction of
# donor x macrophage groups to fall below this floor.
COEXPR_MIN_GROUP_CELLS = 100

# memento's default gene filter (memento/main.py: setup_memento's
# filter_mean_thresh, compute_1d_moments's min_perc_group). Recorded explicitly
# here rather than left implicit, per the analysis spec.
COEXPR_FILTER_MEAN_THRESH = 0.07
COEXPR_MIN_PERC_GROUP = 0.7

COEXPR_TOP_N_GENES = 50

COEXPR_OUTPUT_H5AD = os.path.join(COEXPR_DIR, "elmentaite2021_trio_level3_coexpr.h5ad")
COEXPR_FIGURES_DIR = os.path.join(FIGURES_DIR, "coexpression")

