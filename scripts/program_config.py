"""Configuration for the ANTXR2 gene-program co-expression pilot.

Asks how ANTXR2 co-varies with predefined transcriptional programs in gut
fibroblasts vs enterocytes, using ONE fixed weight vector per program shared
by both cell types. Design decisions are recorded in
planning_summaries/program_coexpression.md; the constants below implement them
and are not to be changed without amending that note.

Kept separate from ``config.py`` (the /data/ANTXR2 Phase 1/2 pipeline, which no
longer exists on this host) and ``hfs_config.py`` (cross-tissue pseudobulk).
"""
import os
from pathlib import Path

# ---------------------------------------------------------------- paths
# Data root; set ANTXR2_ROOT on a host where it lives elsewhere.
ROOT = Path(os.environ.get("ANTXR2_ROOT", "/data/antxr2"))
RAW_DIR = ROOT / "raw" / "gut_cell_atlas"
GENESET_DIR = Path(__file__).resolve().parents[1] / "genesets"   # in the repo: the library is an input
WORK_DIR = ROOT / "program_coexpr"
CELLS_DIR = WORK_DIR / "cells"          # extracted fibroblast + enterocyte counts
OUT_DIR = WORK_DIR / "output"           # csv tables
FIG_DIR = WORK_DIR / "figures"

for _d in (GENESET_DIR, CELLS_DIR, OUT_DIR, FIG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- datasets
# Gut Cell Atlas (Oliver et al. 2024) on CELLxGENE Discover. Identity is
# verified against each file's own uns["citation"] and its matrix nnz, not
# transfer success.
#   healthy   "Healthy reference", 35,461 genes. Holds Elmentaite2021's fetal
#             and most pediatric/adult samples.
#   extended  "Extended+ - 18485 genes", 18,370 genes. The only release with
#             Elmentaite2021's adolescent (13-17 y) and two further 9-12 y
#             donors, which is what lifts the postnatal slice to 12 donors.
SOURCES = {
    "healthy": dict(version_id="b801706d-6345-4af2-b3fe-8e30e0e6150b",
                    dataset_id="bb0818d3-474e-44f4-a82a-169f3c62c94c",
                    h5ad=RAW_DIR / "healthy_reference.h5ad", size=9_399_907_175),
    "extended": dict(version_id="757945c8-a916-431d-aceb-1afbc80a7c55",
                     dataset_id="19053a82-9c89-4fb8-bd19-d7b1800b0b7b",
                     h5ad=RAW_DIR / "extended_plus.h5ad", size=11_576_820_874),
    # Replication cohort (Sanger IBDverse multitissue atlasing cohort): adult
    # terminal-ileum and rectal biopsies + blood, 10x 3' v3, 422 donors,
    # Crohn's and normal. Independent of Elmentaite2021 in donors, lab and
    # chemistry. Standard CELLxGENE obs schema (program_extract adapter).
    "ibdverse": dict(version_id="e3460364-c300-4616-9923-3a0224894998",
                     dataset_id="0e5c75da-7bc9-4a2c-83fc-ef1e378bd5f9",
                     h5ad=ROOT / "raw" / "ibdverse" / "ibdverse.h5ad", size=32_260_433_453),
}
# Same file, crypt (stem/progenitor) epithelium only, for the state contrast;
# a separate npz so the replication superset stays as it was.
SOURCES["ibdverse_crypt"] = SOURCES["ibdverse"]

STUDY_COL = "study"
STUDY = "Elmentaite2021"
# donor_id is wrong in both directions in this atlas (collisions and aliases;
# see config.DONOR_UNIFIED_COL); donorID_unified is the real person.
DONOR_COL = "donorID_unified"
FINE_COL = "level_3_annot"

# Decision 2 (as revised 2026-10-30): the agreed adult-SI slice holds 2 usable
# donors, so the pilot runs two separate donor sets instead -- postnatal SI
# (primary) and fetal SI (developmental contrast). Both cell types always come
# from the same slice. One superset per source (Elmentaite2021 SI, all ages) is
# extracted to CELLS_DIR/fib_ent_<source>.npz; the age slice is applied at
# scoring time.
ADULT = ["18-34", "35-54", "55-74"]
ADOLESCENT = ["13-17"]
PEDIATRIC = ["4-7", "9-12"]
FETAL = ["6-13_first_trimester", "14-20_second_trimester"]
REGION = "Small_intestine"
# A slice filters the source superset on any of: ages (age column), disease,
# regions (region column); a missing key means no filter.
SLICES = {
    "postnatal_si": dict(source="extended", ages=PEDIATRIC + ADOLESCENT + ADULT),
    "fetal_si": dict(source="healthy", ages=FETAL),
    "adult_si": dict(source="extended", ages=ADULT),       # 2 donors; kept for the record
    # IBDverse disease values are "<disease>|<author_inflammation_status>".
    # Replication slice (amendment 2): healthy + Crohn's uninflamed ileum.
    "ibdverse_ileum_uninflamed": dict(source="ibdverse", regions=["ileum"],
                                      disease=["normal|uninflamed", "Crohn disease|uninflamed"]),
    "ibdverse_ileum_inflamed": dict(source="ibdverse", regions=["ileum"],
                                    disease=["Crohn disease|inflamed"]),
    # state contrasts: every ileum disease state; populations pick theirs
    "ibdverse_ileum": dict(source="ibdverse", regions=["ileum"]),
    "ibdverse_ileum_crypt": dict(source="ibdverse_crypt", regions=["ileum"]),
}


def cells_npz(source):
    return CELLS_DIR / f"cells_{source}.npz"


# Each contrast is ANTXR2's program association in population A minus that in
# population B, within one slice. A population is a cell_type from
# program_extract's adapters, optionally narrowed to some subtypes (`fine`);
# `maturation` turns on the enterocyte maturation bins below,
# `min_group_cells` overrides MIN_GROUP_CELLS for that population and
# `disease` narrows it to some disease states (state contrasts: same cell
# type, two states).
#   fib_vs_ent        the original question (postnatal; fetal enterocytes are
#                     below the ANTXR2 floor)
#   fib_vs_smc_fetal  ECM-expressing contrast (user, 2026-10-30): ECM programs
#                     cannot enter fib_vs_ent because enterocytes barely
#                     express them. SMC_PPLP2 is the fetal SI smooth-muscle
#                     subtype with the most donors (15 with >= 100 cells).
CONTRASTS = {
    "fib_vs_ent": dict(slice="postnatal_si", pops={
        "fibroblast": dict(cell_type="fibroblast"),
        "enterocyte": dict(cell_type="enterocyte", maturation=True)}),
    "fib_vs_smc_fetal": dict(slice="fetal_si", pops={
        "fibroblast": dict(cell_type="fibroblast"),
        "smc": dict(cell_type="smooth_muscle", fine=["SMC_PPLP2"])}),
    # replication of fib_vs_ent (pre-registered in the planning note; biopsy
    # stroma is sparse, so fibroblast groups need only 20 cells -- amendment 2)
    "fib_vs_ent_ibdverse": dict(slice="ibdverse_ileum_uninflamed", pops={
        "fibroblast": dict(cell_type="fibroblast", min_group_cells=20),
        "enterocyte": dict(cell_type="enterocyte", maturation=True)}),
    "fib_vs_ent_ibdverse_inflamed": dict(slice="ibdverse_ileum_inflamed", pops={
        "fibroblast": dict(cell_type="fibroblast", min_group_cells=20),
        "enterocyte": dict(cell_type="enterocyte", maturation=True)}),
    # state contrasts (planning note, "State contrast plan"): inflamed vs
    # uninflamed Crohn's (primary) or vs normal (secondary), no donor in both
    **{f"{ct}_infl_vs_{ref}": dict(slice=sl, pops={
        "inflamed": dict(cell_type=ct, disease=["Crohn disease|inflamed"], **extra),
        ref: dict(cell_type=ct, disease=[state], **extra)})
       for ct, sl, extra in [("crypt", "ibdverse_ileum_crypt", {}),
                             ("enterocyte", "ibdverse_ileum", dict(maturation=True))]
       for ref, state in [("uninflamed", "Crohn disease|uninflamed"),
                          ("normal", "normal|uninflamed")]},
}

# Programs derived from the discovery run, appended to the library at scoring
# time (collection "discovery", equal positive weights). IEG_core = genes in
# the top-15 PC1 loadings of >= 3 of the 9 enterocyte-leaning discovery
# programs (fib_vs_ent_fine_mat3_cov25). Fixed before the replication data was
# read; they need MIN_GENES_DISCOVERY genes above the floor in both
# populations instead of MIN_GENES.
DISCOVERY_PROGRAMS = {
    "IEG_core": ["CDKN1A", "JUN", "JUNB", "ZFP36", "MAFF", "CCND1", "KLF6", "DUSP1",
                 "PPP1R15A", "NR4A1"],
}
DISCOVERY_SET = ["cytosig:EGF", "cytosig:TRAIL", "cytosig:FGF2", "cytosig:CD40L",
                 "collectri:JUND", "collectri:TCF4", "collectri:SIRT1", "collectri:PML",
                 "hallmark:TNFA_SIGNALING_VIA_NFKB"]
MIN_GENES_DISCOVERY = 6

# Fixed a-priori modules for the state contrast, from Bracq et al. 2025 (Cmg2
# KO: Ly6a+ fetal-like cells form after DSS but fail to become Lgr5+ ISCs;
# lower nuclear beta-catenin and Lgr5/Axin2/Ascl2/Ccnd1; Cyr61/Ctgf up). Same
# rules as DISCOVERY_PROGRAMS, collection "bracq".
HYPOTHESIS_PROGRAMS = {
    "WNT_output": ["LGR5", "AXIN2", "ASCL2", "CCND1", "RNF43", "ZNRF3", "SMOC2", "EPHB2",
                   "EPHB3", "SOX9"],
    "YAP_fetal": ["CCN1", "CCN2", "AMOTL2", "ANKRD1", "AXL", "ANXA1", "CLU", "TACSTD2",
                  "LAMC2"],
}


# ---------------------------------------------------------------- memento
# Decision 8: the 2026-09-04 corrected size-factor settings (memento's own
# defaults shifted fibroblast/enterocyte correlations positive), q from the
# 2026-09-19 revision, the house raw-mean floor, and the q=0.15 group floor.
Q = 0.15
SHRINKAGE = 0.0
TRIM_PERCENT = 0.5
FILTER_MEAN_THRESH = 0.07     # unnormalised raw count mean, per cell type
MIN_GROUP_CELLS = 100         # donor x fine-subtype groups below this are dropped
TARGET = "ANTXR2"

# ---------------------------------------------------------------- programs
# Decision 6: a gene enters a program only if it clears the floor in BOTH cell
# types; the program needs >= MIN_GENES such genes covering >= MIN_COVERAGE of
# its members measured in the dataset. Programs failing only because of one
# cell type are kept in the tables, flagged not_comparable.
MIN_GENES = 10
MIN_COVERAGE = 0.5
CYTOSIG_TOP = 100             # top |weight| genes per cytokine signature
PROGENY_TOP = 500             # decoupler's default footprint size

# ---------------------------------------------------------------- maturation
# ANTXR2 rises with enterocyte differentiation, and in postnatal enterocytes its
# excess coupling over matched genes sits on the maturation axis (planning note,
# 2026-10-30 follow-up). Enterocytes are therefore binned on a FIXED marker
# score -- canonical mature-absorptive genes minus the immature/crypt marker
# OLFM4, each z-scored within donor -- and grouped donor x bin, so program
# associations are measured within a maturation state. Fixed markers rather
# than the data-derived PC2 so the conditioning is not learned from the genes
# being tested.
MATURATION_MARKERS = {
    "mature": ["APOB", "APOA4", "ALPI", "SI", "LCT", "MTTP", "SLC5A1"],
    "immature": ["OLFM4"],
}
N_MATURATION_BINS = 3         # within-donor quantile bins; 0 disables

# ---------------------------------------------------------------- inference
# The question is the fibroblast - enterocyte DIFFERENCE (user, 2026-10-30), so
# one-sample machinery (within-cell-type permutation, per-cell-type anchor
# percentiles) is gone. The null for the difference is 50 genes matched to
# ANTXR2 on raw mean + detection in BOTH cell types: a program can differ
# between cell types for any gene, and only ANTXR2's excess over that is about
# ANTXR2.
N_ANCHORS = 50
N_BOOT = 1000                 # donor bootstrap, paired across cell types
SEED = 5

# ---------------------------------------------------------------- landscape
# ANTXR2 co-expression landscape (planning note, "Reframing ... landscape"):
# IBDverse ileum, every author cell type with >= 10 donors holding a
# >= 100-cell group in some state (rule on cell counts only, fixed
# 2026-11-03), plus pooled fibroblasts at the replication's floor of 20.
# One npz per population in LANDSCAPE_CELLS (program_extract --source
# ibdverse_landscape); cells of the matrix are population x state.
SOURCES["ibdverse_landscape"] = SOURCES["ibdverse"]
LANDSCAPE_CELLS = CELLS_DIR / "landscape"
LANDSCAPE_OUT = OUT_DIR / "landscape"
LANDSCAPE_POPS = [
    "B germinal centre/plasmablast", "Naive B", "Memory B PTPRJ", "Memory B PDE4D TEX9",
    "Enterocyte progenitor/precursor crypt OLFM4lo", "Enterocyte top villus",
    "Monocyte S100A8/9", "Macrophage CD163lo", "Macrophage CD163hi",
    "Plasma B CD38lo", "Plasma B CD38int", "Plasma B CD38hi IgG IgA",
    "Ileum goblet MUC2lo BCAS1int", "Ileum goblet MUC2hi SPINK4hi",
    "Ileum goblet MUC2hi BCAS1hi", "Ileum goblet MUC2int",
    "Enterocyte stem OLFM4lo LGR5hi", "Enterocyte stem MKI67",
    "Enterocyte progenitor crypt OLFM4hi",
    "Tissue CD4 Treg", "CD8 TRM TRGC2", "CD4 memory CXCR6", "CD4 memory", "CD4 PASK",
    "CD4 CM/naive", "Tissue resident GD T",
]
LANDSCAPE_FIB = "fibroblast (pooled)"          # ontology classes below, floor 20
LANDSCAPE_FIB_TYPES = ["fibroblast", "subepithelial intestinal fibroblast"]
LANDSCAPE_FIB_FLOOR = 20
LANDSCAPE_STATES = {"normal": "normal|uninflamed",
                    "cd_uninflamed": "Crohn disease|uninflamed",
                    "cd_inflamed": "Crohn disease|inflamed"}
LANDSCAPE_MIN_DONORS = 10     # donor groups for a population x state cell
CORE_PASS_FRAC = 0.5          # core gene: clears the floor in >= half the filled cells
CELL_CORE_COVERAGE = 0.8      # a program is scored where >= 80% of its core passes


def pop_slug(name):
    import re
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")

# Anchors are also matched on genomic span (planning note, "Correction: a
# gene-length axis"): candidates within ANCHOR_SPAN_WINDOW (log10 bp) of
# ANTXR2's span, then the nearest on raw mean + detection.
GENE_SPAN = ROOT / "annotation" / "gene_span.csv"
ANCHOR_SPAN_WINDOW = 0.15


def gene_spans():
    """Series gene symbol / unversioned Ensembl ID -> genomic span (bp)."""
    import pandas as pd
    a = pd.read_csv(GENE_SPAN)
    s = pd.concat([a.drop_duplicates("gene_name").set_index("gene_name").span,
                   a.set_index("gene_id").span])
    return s[~s.index.duplicated()]

# Lineages for the lineage-wide consistency null: a second set of N_ANCHORS
# genes matched to ANTXR2 (span window + raw mean + detection) jointly across
# every filled cell of the lineage and scored in all of them, so "is a
# lineage-wide ANTXR2 profile unusual for one gene" has same-gene nulls. The
# per-cell excess still uses each population's own anchors.
LANDSCAPE_LINEAGES = {
    "epithelial": [p for p in LANDSCAPE_POPS if "Enterocyte" in p or "goblet" in p],
    "lymphocyte": [p for p in LANDSCAPE_POPS if p.endswith(" B") or " B " in p or p.startswith("B ")
                   or "CD4" in p or "CD8" in p or "Treg" in p or p.startswith("Tissue resident")],
    "myeloid": ["Monocyte S100A8/9", "Macrophage CD163lo", "Macrophage CD163hi"],
    "fibroblast": [LANDSCAPE_FIB],
}
IBDVERSE_LANDSCAPE_CELLS = LANDSCAPE_CELLS

# Replication atlas (planning note, "Step 4 plan"): GCA Extended+ postnatal
# small intestine, normal donors, every level_3_annot type with >= 10 donors
# holding a >= 100-cell group. ANTXR2_LANDSCAPE=gca switches every landscape
# script to it; it is scored with the IBDverse landscape's weights
# (LANDSCAPE_WEIGHTS_FROM), never re-fitted.
SOURCES["gca_landscape"] = SOURCES["extended"]
GCA_LANDSCAPE_CELLS = CELLS_DIR / "landscape_gca"
GCA_LINEAGES = {
    "epithelial": ["Enterocyte", "TA", "Goblet", "Goblet_progenitor", "Goblet_cycling",
                   "Paneth", "Epithelial_stem"],
    "lymphocyte": ["Trm_CD8", "Trm_CD4", "Trm_Th17", "B_plasma_IgA1", "B_memory", "B_naive"],
}
GCA_LANDSCAPE_POPS = [p for ps in GCA_LINEAGES.values() for p in ps]
GCA_LINEAGES["fibroblast"] = [LANDSCAPE_FIB]     # level_2_annot == GCA_FIB_LEVEL2, pooled
GCA_FIB_LEVEL2 = "Fibroblast"
# Fetal SI fibroblasts (user, 2026-11-10: boost the fibroblast arm; Extended+
# holds 59,192 from 26 donors). Scored in the ref run as two extra
# populations: fetal alone, and postnatal + fetal pooled (within donor x
# subtype, so the age difference in means does not enter the correlations).
GCA_FIB_FETAL = "fibroblast (fetal)"
GCA_FIB_ALL = "fibroblast (postnatal + fetal)"
SOURCES["gca_fetal_fib"] = SOURCES["extended"]
# Skin fibroblasts (user, 2026-11-15: replicate the fibroblast coupling where
# fibroblasts are plentiful). Human dermal fibroblast atlas (CELLxGENE
# collection 3c4f0970, 25 studies integrated, fibroblasts only): healthy skin,
# 10x 3' v2 / v3. Two populations over the same cells: pooled within donor,
# and within donor x author cluster (cluster_robust).
SOURCES["skin_fib"] = dict(version_id="8ae4471a-2e70-47b5-bc42-4f0df26a9996",
                           dataset_id="a19d1667-a7b5-4556-9e5f-f9bfa690c0f1",
                           h5ad=ROOT / "raw" / "skin" / "dermal_fibroblast_atlas.h5ad",
                           size=1_831_628_432)
SKIN_LANDSCAPE_CELLS = CELLS_DIR / "landscape_skin"
SKIN_ASSAYS = ["10x 3' v2", "10x 3' v3"]
SKIN_FIB_CLUSTER = "fibroblast (by cluster)"
# Single most reliable subtype per tissue (user, 2026-11-18): donor groups
# (one sample per donor in both atlases, so donor x batch = donor).
IBD_FIB_MAIN = "fibroblast ADAMDEC1"     # IBDverse author_cell_type LP fibroblast ADAMDEC1
# Same cells, normal + CD uninflamed pooled (each alone has 9-10 donors with
# >= 20 cells). Its disease labels are rewritten to "normal|uninflamed", so the
# cell id reads "@normal" but holds both uninflamed states.
IBD_FIB_MAIN_UNINF = "fibroblast ADAMDEC1 uninflamed"
SKIN_FIB_MAIN = "fibroblast A3"          # dermal atlas cluster_robust A3 (B1 = stress state, excluded)
LANDSCAPE_EXTRA_POPS = []
LANDSCAPE_ATLAS = os.environ.get("ANTXR2_LANDSCAPE", "ibdverse")
LANDSCAPE_WEIGHTS_FROM = None
if LANDSCAPE_ATLAS == "gca":
    LANDSCAPE_WEIGHTS_FROM = LANDSCAPE_OUT
    LANDSCAPE_CELLS = GCA_LANDSCAPE_CELLS
    # ANTXR2_CELL_COVERAGE lowers the share of a program's IBDverse core that
    # must clear the floor in a GCA cell (exon-only, shallower); its own output
    CELL_CORE_COVERAGE = float(os.environ.get("ANTXR2_CELL_COVERAGE", CELL_CORE_COVERAGE))
    LANDSCAPE_OUT = OUT_DIR / ("landscape_gca" if CELL_CORE_COVERAGE == 0.8
                               else f"landscape_gca_cov{round(CELL_CORE_COVERAGE * 100)}")
    LANDSCAPE_POPS = GCA_LANDSCAPE_POPS
    LANDSCAPE_STATES = {"normal": "normal"}
    LANDSCAPE_LINEAGES = GCA_LINEAGES
elif LANDSCAPE_ATLAS == "skin":
    LANDSCAPE_WEIGHTS_FROM = LANDSCAPE_OUT
    LANDSCAPE_CELLS = SKIN_LANDSCAPE_CELLS
    LANDSCAPE_OUT = OUT_DIR / "landscape_skin"
    LANDSCAPE_POPS = []
    LANDSCAPE_STATES = {"normal": "normal"}
    LANDSCAPE_LINEAGES = {"fibroblast": [LANDSCAPE_FIB]}
else:
    assert LANDSCAPE_ATLAS == "ibdverse", LANDSCAPE_ATLAS

# Reference-scaled run (planning note, "Reference-scaled coupling";
# program_landscape_ref.py), ANTXR2_LANDSCAPE_RUN=ref: weights from every
# gene that clears the floor in >= REF_MIN_CELLS filled cells, each gene on a
# fixed reference SD, 200 anchors, its own output folder. In GCA, ANTXR2's
# own floor drops to 0.05 (exploratory; the gene floor stays 0.07).
LANDSCAPE_RUN = os.environ.get("ANTXR2_LANDSCAPE_RUN", "")
LANDSCAPE_TARGET_FLOOR = FILTER_MEAN_THRESH
REF_MIN_CELLS = 3
if LANDSCAPE_RUN == "ref":
    N_ANCHORS = 200
    if LANDSCAPE_ATLAS == "gca":
        LANDSCAPE_WEIGHTS_FROM = OUT_DIR / "landscape_ref"
        LANDSCAPE_OUT = OUT_DIR / "landscape_gca_ref"
        LANDSCAPE_TARGET_FLOOR = 0.05
        LANDSCAPE_EXTRA_POPS = [GCA_FIB_FETAL, GCA_FIB_ALL]
    elif LANDSCAPE_ATLAS == "skin":
        LANDSCAPE_WEIGHTS_FROM = OUT_DIR / "landscape_ref"
        LANDSCAPE_OUT = OUT_DIR / "landscape_skin_ref"
        LANDSCAPE_EXTRA_POPS = [SKIN_FIB_CLUSTER, SKIN_FIB_MAIN]
    else:
        LANDSCAPE_OUT = OUT_DIR / "landscape_ref"
        LANDSCAPE_EXTRA_POPS = [IBD_FIB_MAIN, IBD_FIB_MAIN_UNINF]
else:
    assert LANDSCAPE_RUN == "", LANDSCAPE_RUN
