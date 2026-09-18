"""Configuration for the HFS / ANTXR2 cross-tissue atlas arm.

This arm asks a different question from the co-expression pipeline in
``config.py``: instead of *which partners are available in a cell type*, it asks
*which tissue's fibroblasts carry the highest collagen VI burden per unit of
ANTXR2 receptor*, i.e. the substrate:receptor ratio that the clearance function
would have to service.

Kept separate from ``config.py`` on purpose — different datasets (CELLxGENE
whole-atlas h5ads rather than the curated co-expression object), different
aggregation (library-size-normalised pseudobulk CPM, not memento moments), and
different output directory.
"""
from pathlib import Path

# ---------------------------------------------------------------- paths
HFS_DIR = Path("/data/ANTXR2/hfs_atlas")
HFS_RAW_DIR = HFS_DIR / "raw"           # downloaded .h5ad files
HFS_PB_DIR = HFS_DIR / "pseudobulk"     # *_pseudobulk.npz
HFS_OUTPUT_DIR = HFS_DIR / "output"     # csv tables
HFS_FIGURES_DIR = HFS_DIR / "figures"

for _d in (HFS_RAW_DIR, HFS_PB_DIR, HFS_OUTPUT_DIR, HFS_FIGURES_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- datasets
# CELLxGENE Discover dataset-version UUIDs, read back out of each downloaded
# file's own uns["citation"] field (not from the portal UI), so the registry and
# the files on disk cannot silently disagree.
#   download URL = https://datasets.cellxgene.cziscience.com/<version_id>.h5ad
CXG_BASE = "https://datasets.cellxgene.cziscience.com/{version_id}.h5ad"

DATASETS = {
    "ts_stromal": dict(
        version_id="2776ea76-6919-458a-b82a-bba8bb6ce5a9",
        collection_id="e5f58829-1a66-40b5-a624-9046778e74f5",
        title="Tabula Sapiens - Stromal", gb=10.06),
    "bmarrow": dict(
        version_id="4aedb645-7487-4ab0-ad95-ee4db1b8eccd",
        collection_id="e5f58829-1a66-40b5-a624-9046778e74f5",
        title="Tabula Sapiens - Bone_Marrow", gb=1.03),
    "gut": dict(
        version_id="c54fbb37-2a1f-451a-815b-3b1fc034b62e",
        collection_id="e33ffcd3-7cbf-4b8c-b0f4-85587ad5019a",
        title="Total - Cells of the human intestinal tract mapped across space and time", gb=5.73),
    "skin": dict(
        version_id="f8b5c4f4-7ece-4314-9f49-fc254bbd5c71",
        collection_id="73f82ac8-15cc-4fcd-87f8-5683723fce7f",
        title="Human healthy adult skin scRNA-seq data", gb=1.03),
    "oral": dict(
        version_id="072148e6-b563-44cf-a1b2-8162fa4b107c",
        collection_id="065ad318-59fd-4f8c-b4b1-66caa7665409",
        title="Mucosal Atlas (Human Oral & Craniofacial Cell Atlas)", gb=0.54),
    "synovium": dict(
        version_id="f758894c-14bc-4bfe-94dd-16dd9945f7d3",
        collection_id="10eb236d-d42d-45b8-8363-c2dcf865f388",
        title="Integrated global cells (JIA synovium)", gb=2.33),
    "tendon_ach": dict(
        version_id="accd7d2e-d67d-4b91-9ea1-f2d1eee78ebf",
        collection_id="a7e81820-297d-4086-8123-cc7dde64e495",
        title="snRNA-seq of the human Achilles tendon", gb=1.12),
    "tendon_quad": dict(
        version_id="2df978a2-4a88-4880-8a7c-cf495bf250d9",
        collection_id="579203e2-182f-47bc-8230-7aa47247e2a4",
        title="snRNA-seq of the human quadriceps tendon", gb=0.32),
    # suspension-type control datasets: both cell and nucleus in one study
    "kidney": dict(
        version_id="4cd166f1-ef51-4137-869d-0a3688bc2bc8",
        collection_id="bcb61471-2a44-4d00-a0af-ff085512674c",
        title="Integrated Single-nucleus and Single-cell RNA-seq of the Adult Human Kidney", gb=2.96),
    "heart": dict(
        version_id="1e1e07c3-bfcb-4a0d-91de-c2614a891409",
        collection_id="3116d060-0a8e-4767-99bb-e866badea1ed",
        title="Combined single cell and single nuclei RNA-Seq data - Heart Global", gb=4.97),
    "muscle": dict(
        version_id="1cb67ec8-c9e3-4d4a-ba23-a1a8ef3f8450",
        collection_id="854c0855-23ad-4362-8b77-6b1639e7a9fc",
        title="Skeletal_muscle", gb=1.44),
}

# ---------------------------------------------------------------- genes
# Panel streamed per-cell by hfs_pseudobulk.extract(); the full-gene pseudobulk
# pass keeps every gene, so collagen groups below are resolved at analysis time.
PANEL = ["ANTXR2", "ANTXR1",
         "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
         "COL4A1", "COL4A2", "COL1A1", "COL3A1",
         "LAMA4", "LAMB1", "LRP6", "MMP2", "MMP14", "HSPG2", "PDGFRB"]

COLLAGEN_GROUPS = {"COL6": ["COL6A1", "COL6A2", "COL6A3"],
                   "COL1": ["COL1A1", "COL1A2"],
                   "COL3": ["COL3A1"],
                   "COL4": ["COL4A1", "COL4A2"],
                   "COL5": ["COL5A1", "COL5A2"]}

# "all collagen" denominator = any gene matching this pattern
COLLAGEN_RE = r"COL\d+A\d+"

# Tabula Sapiens fibroblast labels treated as one population
FIBROBLAST_LABELS = {"fibroblast", "alveolar adventitial fibroblast",
                     "fibroblast of cardiac tissue", "fibroblast of breast",
                     "thymic fibroblast type 1", "thymic fibroblast type 2"}

ORAL_SITES = ["gingiva", "buccal mucosa", "hard palate"]

# ---------------------------------------------------------------- thresholds
MIN_CELLS_PER_DONOR_GROUP = 50   # Tabula Sapiens / oral donor-level groups
MIN_CELLS_SMALL_ATLAS = 25       # synovium / tendon, where donors are few
MIN_CELLS_POOLED = 150           # pooled per-tissue fibroblast floor
MIN_CPM_FOR_RATIO = 1.0          # both arms must clear this in bias log-ratios
PSEUDOBULK_TARGET_NNZ = 4_000_000  # streaming block size (low-RAM hosts)

# ---------------------------------------------------------------- phenotype
# HFS organ-involvement classes, from classification of 45 clinical reports.
# Drives the "affected (core) / occasional / rare / not reported" encoding.
# See HFS_ATLAS_NOTES.md for how the literature set was assembled.
TISSUE_INVOLVEMENT = {
    "Skin": "affected (core)",
    "Large_Intestine": "affected (core)",
    "Small_Intestine": "affected (core)",
    "Tongue": "occasional", "Eye": "occasional", "Lung": "occasional",
    "Liver": "occasional", "Muscle": "occasional", "Spleen": "occasional",
    "Heart": "rare",
    "Fat": "not reported", "Bladder": "not reported", "Vasculature": "not reported",
    "Stomach": "not reported", "Uterus": "not reported", "Salivary_Gland": "not reported",
    "Trachea": "not reported", "Prostate": "not reported", "Pancreas": "not reported",
    "Mammary": "not reported", "Thymus": "not reported", "Ovary": "not reported",
    "Bone_Marrow": "not reported",
}

ORGAN_REPORT_FREQ = {"skin": 41, "gingiva": 29, "joint": 28, "bone": 18, "intestine": 15,
                     "oral mucosa": 4, "ear": 3, "lung": 3, "muscle": 3, "spleen": 3,
                     "eye": 3, "liver": 2, "nose": 2, "thyroid": 2, "heart": 1,
                     "kidney": 1, "brain": 1, "breast": 1, "lymph node": 1, "adrenal gland": 1}
N_CLINICAL_REPORTS = 45

# Organs frequently involved in HFS with no fibroblast single-cell data in any
# atlas surveyed here. This is the standing data gap for this arm.
ABSENT_FROM_ATLASES = {"gingiva": 29, "joint": 28, "bone": 18}
