# Phase 1 (skin): ANTXR2 mean expression in the skin fibroblast atlas

## Purpose

Compute memento method-of-moments **mean expression** per fibroblast subtype in
human skin, for ANTXR2 and a small companion panel.

Framing note (this drives the gene panel, so don't drop it): HFS is caused by
coding loss-of-function in ANTXR2, so this is **not** a study of how ANTXR2
transcription is regulated. The motivating question is the tissue/cell-type
restriction of the phenotype — the mutation is germline and present in every
cell, but only some tissues manifest. Means are step one toward two questions:
(1) does ANTXR1/ANTXR2 paralog co-presence vary by cell type, and (2) is the
ECM-clearance machinery assembled where ANTXR2 is expressed. Correlation
(Phase 2) is what actually resolves those; means only establish what's present.

## Dataset

**Steele et al., *Nature Immunology* (2025) — single-cell and spatial atlas of
human skin fibroblasts.** Healthy skin plus 23 skin diseases; defines 6 healthy
fibroblast subtypes (papillary, reticular, FRC-like/perivascular, Schwann-like,
etc.) plus 3 disease-specific subtypes including inflammatory myofibroblasts.

- Portal: https://collections.cellatlas.io/skin-fibroblast
- Code: https://github.com/haniffalab/skin_fibroblast_atlas

Chosen over Tabula Sapiens' skin compartment because TS annotates a single
coarse "fibroblast of skin" bucket, which cannot address subtype restriction.

**Download caveat — resolve before scripting.** This atlas does **not** appear
to be on CELLxGene Discover (no collection ID found), so the CELLxGene Curation
API path used for the other Phase 1 datasets does not apply. The portal page is
JS-rendered. Check the GitHub repo and the paper's data-availability statement
for a direct `.h5ad` URL or an ArrayExpress / HCA DCP accession, and prefer a
stable accession over scraping the portal. Do not guess a URL.

## Steps

1. **Download** to `/data/raw/skin_fibroblast_atlas/`, recording source URL /
   accession, file size, checksum, and date in a `manifest.csv`.
2. **Verify raw counts.** memento models the count-sampling process and needs
   raw integer counts, not normalized/log values. Check `adata.X` vs
   `adata.raw.X` (raw should be large and integer-like) and log which slot was
   used. Do not assume.
3. **Inspect annotations before grouping.** Identify the column holding the
   fine-grained fibroblast subtype labels (F1–F6 + disease-specific), and note
   any healthy/lesional and disease-label columns. Report the subtype label
   vocabulary before computing anything.
4. **Compute means with memento** (`pip install memento-de`; standalone local
   pipeline, not the CELLxGene precomputed backend) grouped by
   `donor_id x disease_state x fibroblast_subtype`, falling back to
   `disease_state x fibroblast_subtype` if donor counts are too low.
   Report `n_cells` per group alongside every estimate.
5. **Save** to `/data/celltype_expression/skin_fibroblast_celltype_means.parquet`
   (long format) with columns: `dataset`, `donor_id`, `disease_state`,
   `fibroblast_subtype`, `gene`, `n_cells`, `mean_expression`, `assay`,
   `raw_counts_slot`, `pipeline_run_date`. Add a short README noting provenance.

## Gene panel

Compute genome-wide if tractable (Phase 2 needs many genes from the same file;
avoid a second pass). Otherwise at minimum:

- **Receptors/paralog:** `ANTXR2`, `ANTXR1`
- **Collagen VI + ECM:** `COL6A1`, `COL6A2`, `COL6A3`, `COL1A1`, `COL4A1`
- **Clearance machinery:** `MRC2`, `CTSB`, `CTSK`, `MMP14`, `TIMP2`, `LAMP1`
- **Synthesis/activation axis (confound proxy):** `ACTA2`, `POSTN`, `TGFB1`
- **Senescence:** `CDKN1A`, `CDKN2A`
- **Subtype markers for QC:** whichever the paper uses for F1–F6
- **Cell-identity QC:** `PDGFRA`, `PTPRC`, `PECAM1`, `KRT14`

`ANTXR1` is not optional — the paralog-partition question is the main reason
this dataset is being pulled.

## Notes / expected pitfalls

- Means cannot distinguish "both paralogs in every cell" from "two disjoint
  subpopulations." That is a Phase 2 correlation question; do not over-read
  co-presence of nonzero means here.
- Flag any subtype with fewer than ~20 cells in a group as low-confidence;
  these will likely be unusable for Phase 2 correlations even if reportable now.
- Keep `disease_state` as a column throughout; do not pool healthy and lesional.
- Do **not** report fraction-of-cells-expressing as a primary metric — it is
  confounded by per-cell RNA content across fibroblast subtypes of differing
  size and activation.

## Done when

- [ ] `.h5ad` in `/data/raw/skin_fibroblast_atlas/` + `manifest.csv`
- [ ] Raw-counts slot verified and logged
- [ ] Subtype label vocabulary reported
- [ ] Parquet + README in `/data/celltype_expression/`
- [ ] `n_cells` present for every row
