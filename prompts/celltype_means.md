# Phase 1: ANTXR2 Cell-Type Expression Atlas

## Objective

Produce rigorous, per-cell-type **mean expression** estimates for ANTXR2 (and a
validation panel of related genes) across curated single-cell reference atlases,
using memento's method-of-moments estimator run locally on real count matrices.

This replaces an earlier informal pass that used memento's CELLxGene-precomputed
backend (`memento-cxg`) and a naive "fraction of cells expressing" metric. Two
decisions carry over from that discussion and should be treated as fixed
constraints on this phase, not just background:

1. **Don't use the CELLxGene-precomputed memento backend.** It's still immature.
   Instead, download real datasets to `/data/ANTXR2/raw/` and run the standalone
   `memento` package locally.
2. **Don't use "fraction of cells expressing" as the primary metric.** It's
   confounded by per-cell RNA content / sequencing depth — large mesenchymal
   cells (smooth muscle, stromal, fibroblasts) will look inflated relative to
   small lymphocytes independent of true relative expression. memento's
   method-of-moments mean estimator corrects for this by explicitly modeling
   the sampling process, so use that instead.

Scope note: this phase computes **means only** (as requested). Variance and
gene-gene correlation are Phase 2 work, but see the note in "Design decisions"
below about why it's worth computing means genome-wide (or close to it) now
rather than just for ANTXR2.

---

## 1. Datasets

Three CELLxGene Discover collections, chosen to cover the three biological
threads from the lit review (skin/ECM, gut/Wnt regeneration, immune/AS) with
consistently-processed, well-annotated cell type labels (Cell Ontology terms,
harmonized by CELLxGene's schema — no manual ontology reconciliation needed
across collections).

### 1a. Tabula Sapiens — broad multi-tissue baseline

- **Collection ID:** `e5f58829-1a66-40b5-a624-9046778e74f5`
- **Link:** https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
- **Size:** ~1.1M cells (v1 + v2 combined), 28 organs, 24 donors, ~475+ annotated
  cell types. Includes skin, blood, small intestine, large intestine, uterus,
  vasculature, and more, all processed with a single consistent pipeline.
- **Why:** direct successor to the informal blood-vs-intestine plot, but across
  many more tissues at once and without cross-study batch confounds (same
  pipeline, same consortium). This is the broad-survey / sanity-check dataset.
- **Caveats:** mixes Smart-seq2 and 10x droplet data — check `assay` metadata
  and consider computing separately or at least flagging platform per cell
  type, since the two have different capture efficiency (relevant to the
  memento model). Per-organ, per-cell-type sample sizes can be small for rare
  populations.

### 1b. Cross-tissue Immune Cell Atlas (Domínguez Conde et al., 2022) — deep immune reference

- **Collection ID:** `62ef75e4-cbea-454e-a0ce-998ec40223d3`
- **Link:** https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3
- **Size:** ~330K immune cells, 12 donors, tissues spanning blood, thymus, bone
  marrow, spleen, lymph node, gut, lung, liver, muscle. ~45 fine-grained immune
  cell types/states, with paired BCR/TCR data (not needed here, but present).
- **Why:** needed for the immune/AS thread. Resolves immune lineages at much
  finer granularity than Tabula Sapiens' immune compartment, and — importantly
  — profiles the **same donors' immune cells across multiple tissues**, so a
  macrophage-in-gut vs. monocyte-in-blood comparison (like the one in the
  informal plot) is a within-donor comparison here, not cross-study.
- **Caveats:** healthy donors only — no AS or IBD samples, so this is a mapping
  resource for Phase 1, not a disease-comparison resource. No skin.

### 1c. Gut Cell Atlas — integrated pan-GI resource (Oliver et al., 2024/2025)

- **Collection ID:** `f11cb29c-b546-4738-9bd8-66ea621a7bd5`
- **Link:** https://cellxgene.cziscience.com/collections/f11cb29c-b546-4738-9bd8-66ea621a7bd5
- **Paper:** "Single-cell integration reveals metaplasia in inflammatory gut
  diseases," *Nature* (2024/2025)
- **Size:** ~1.6M cells total — a healthy reference of ~1.1M cells / 385 samples
  / 189 donors, plus ~500K cells from 12 disease datasets (celiac disease,
  ulcerative colitis, Crohn's disease, GI cancers). 136 fine-grained cell
  states, spanning development through adulthood.
- **Why:** needed for the gut/Wnt thread specifically. This is the only
  collection with fine-grained enough epithelial annotation (stem/TA/crypt vs.
  differentiated enterocyte, rather than coarse "epithelial cell" bins) to
  directly test whether ANTXR2 skews toward immature/stem epithelial states —
  the hypothesis raised by the "epithelial cell" vs. "enterocyte of colon"
  gradient in the informal plot. It's also the only one with matched UC/Crohn's
  samples, which sets up a real (human, chronic-disease) differential
  comparison for a later phase — a more defensible one than trying to map onto
  the mouse DSS biology (see prior discussion on why that's discussion-only).
- **Caveats:** still human chronic IBD, not acute injury/recovery — keep the
  distinction from the mouse DSS mechanism explicit in any writeup that cites
  both.

### Not included in Phase 1 (revisit if these threads become priorities)

- A dedicated skin-specific atlas (Tabula Sapiens' skin compartment should be
  enough for an initial pass on the collagen VI / ECM story; go deeper only if
  skin becomes a focus).
- The liver scRNA-seq datasets used in the Huang et al. ANTXR2-liver-fibrosis
  paper (GSE136103, GSE181483) — relevant if the endothelial/MMP2 thread gets
  picked up later.
- A joint/synovium atlas — relevant if the AS GWAS thread gets pursued directly
  rather than through the immune cross-tissue atlas.

---

## 2. Downloading to `/data/ANTXR2/raw/`

- Create `/data/ANTXR2/raw/` if it doesn't already exist.
- **Resolve download URLs at runtime via the CELLxGene Discover Curation API**
  rather than hardcoding dataset-level URLs (those rotate; collection IDs are
  stable):
  `GET https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}`
  returns the current datasets in a collection and their asset download URLs.
  Verify the exact current endpoint/response shape when implementing — CZI has
  changed this API before.
- For each of the 3 collection IDs above: list datasets, download each
  dataset's `.h5ad` asset.
- **Folder structure:**
  ```
  /data/ANTXR2/raw/
    tabula_sapiens/{dataset_id}.h5ad
    cross_tissue_immune/{dataset_id}.h5ad
    gut_cell_atlas/{dataset_id}.h5ad
    manifest.csv   # collection_id, dataset_id, title, tissue, download_url,
                    # file size, checksum, download timestamp
  ```
- **Raw counts check (important, easy to get wrong):** memento's estimator
  models the count-sampling process and requires actual raw integer counts,
  not normalized/log values. CELLxGene schema conventions put raw counts in
  `adata.raw.X` when `adata.X` holds normalized values — but verify this per
  dataset (`adata.raw.X.max()` should be a large integer-ish value; normalized
  `adata.X` will look like small floats). Don't assume; check and log which
  slot was used for each dataset in the manifest.
- Some collections may include non-RNA modalities (ATAC, spatial) or QC-failed
  cells — filter to standard scRNA-seq assays and cells passing the
  collection's own QC flags if present.

---

## 3. Computing cell-type-specific means with memento

- Install the standalone package: `pip install memento-de`
  (repo: https://github.com/yelabucsf/scrna-parameter-estimation) — this is
  the local pipeline, not the CELLxGene-precomputed backend.
- For each downloaded dataset:
  1. Load raw counts (verified slot from step 2).
  2. Group cells by the harmonized `cell_type` field (Cell Ontology label +
     ID — already consistent across all three collections thanks to
     CELLxGene's schema).
  3. Within each `dataset_id x cell_type` group, run memento's method-of-moments
     mean estimator.
  4. Record `n_cells` per group alongside the mean — needed to flag unstable
     estimates from rare cell types (this came up directly in review of the
     informal plot: ICC/neuroblast-type rare populations need their sample
     size visible, not hidden).
- **Grouping granularity:** group by `dataset_id x donor_id x cell_type` if
  donor metadata is available and cell counts support it, so donor/batch
  effects can be assessed later; otherwise fall back to `dataset_id x
  cell_type`. Don't collapse straight to a single cross-collection cell-type
  mean — keep collection/dataset as a column throughout, since technical
  characteristics differ across the three collections and cross-collection
  comparability shouldn't be assumed silently.

---

## 4. Saving to `/data/ANTXR2/celltype_expression/`

- Create `/data/ANTXR2/celltype_expression/` if it doesn't already exist.
- One output table per source collection (long format, parquet):
  ```
  /data/ANTXR2/celltype_expression/
    tabula_sapiens_celltype_means.parquet
    cross_tissue_immune_celltype_means.parquet
    gut_cell_atlas_celltype_means.parquet
    combined_celltype_means.parquet   # all three, tagged by source_collection
    README.md                          # column definitions, provenance, date run
  ```
- **Columns:** `source_collection`, `dataset_id`, `donor_id` (if available),
  `tissue`, `tissue_general`, `cell_type`, `cell_type_ontology_id`, `gene`,
  `n_cells`, `mean_expression`, `disease` (healthy/celiac/UC/Crohn's etc. —
  relevant for the Gut Cell Atlas collection), `assay`, `pipeline_run_date`.

---

## Design decisions to confirm before running

1. **Gene scope: genome-wide (or a large panel) vs. ANTXR2-only.** Recommend
   computing means for **all genes, or at least a large panel**, not just
   ANTXR2 — Phase 2 (coexpression) will need many genes' worth of data from
   these same raw files, and memento's mean estimator is cheap. Doing it now
   avoids a second full pass over the raw `.h5ad` files later. If genome-wide
   is too slow/large to store, at minimum include: `ANTXR2`, ECM genes
   (`COL6A1`, `COL6A2`, `COL6A3`, `COL4A1`, `LAMA1`), Wnt pathway genes
   (`LRP6`, `AXIN2`, `ASCL2`, `LGR5`, `WNT2B`, `WNT5A`, `WNT5B`), proliferation
   (`MKI67`), and canonical markers for cell types flagged in the prior review
   (`ACTA2`/`MYH11` for smooth muscle, `KIT`/`ANO1` for ICC, `EPCAM`, `PTPRC`,
   `PECAM1`) for QC sanity-checking.
2. **Smart-seq2 vs. 10x in Tabula Sapiens** — compute separately, or pool with
   an assay flag? Recommend keeping `assay` as a column and not pooling
   silently.
3. Confirm the CELLxGene Curation API response shape before writing the
   download script — it's the one piece of this spec most likely to have
   drifted since this doc was written.

## Definition of done

- [ ] `/data/ANTXR2/raw/` populated with `.h5ad` files for all datasets in the 3
      collections, plus `manifest.csv`
- [ ] Raw-counts slot verified and logged per dataset
- [ ] `/data/ANTXR2/celltype_expression/` populated with the 4 parquet files + README
- [ ] Spot check: ANTXR2 mean in Tabula Sapiens skin fibroblasts / collagen
      genes correlate directionally as expected (positive control from the
      Bürgi et al. collagen VI literature)
- [ ] Spot check: `n_cells` column present and used to flag any cell type
      with fewer than ~20 cells in a given group as low-confidence
