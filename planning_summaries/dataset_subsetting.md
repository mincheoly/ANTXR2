# Plan: prepare a subset h5ad + capture-rate columns for single-cell co-expression

## Context

Phase 1 (cell-type mean expression) is complete for both atlases. The next
phase is genuinely different: **single-cell co-expression** — whether the *same
individual cells* express ANTXR2 together with its paralog and clearance
partners, as opposed to the cell-type-mean co-presence we've measured so far.
Means cannot distinguish "both genes in every cell" from "two disjoint
subpopulations averaging out", and that distinction is the whole question.

This plan covers **preparation only** — carving out a small, well-characterised
working dataset and attaching the capture-rate metadata memento needs. No
correlation is computed here.

Why this specific slice: investigating the enterocyte label showed the
whole-body "218,750 cells" figure was inflated by double-counting between two
overlapping Gut Cell Atlas releases (all 75 "Healthy reference" donors are also
in "Extended+"), and that the label is really only two deeply-sampled studies.
**Elmentaite2021** turned out to be the ideal first-pass slice: 398,460 cells,
41 donors, 100% `Non_pathological`, and — critically — it contains gut
epithelium, fibroblasts, and true myeloid cells from **the same 36 donors**,
one lab, one protocol, balanced 5′/3′ chemistry. That removes essentially every
technical confound from the three-way comparison at once.

## Correction that changes the cell-type trio — read first

**`gastrointestinal tract (lamina propria) macrophage` is not a macrophage.**
All 42,302 cells carrying that `cell_type` across the entire Gut Cell Atlas are
authored as `level_1_annot=Mesenchymal`, `level_2_annot=Fibroblast`,
`level_3_annot=Lamina_propria_fibroblast_ADAMDEC1`. Markers confirm the authors:
PTPRC 0.4% positive, CD68 2.5%, LYZ 0.5%, C1QA 0.5%, but COL1A1 **89%**
positive. The CELLxGene ontology mapping (CL:0000865) is wrong, systematically,
across all 10+ constituent studies.

Consequences to handle:
1. The trio uses plain **`macrophage`** (4,694 cells, 41 donors; C1QA 96%,
   CD68 82%, AIF1 90% — correct) as the myeloid member.
2. This label is currently in the published heatmap's HFS-diffuse/systemic
   category *described as a macrophage*, and appears in the ANTXR1 top-20 table
   in `ANALYSIS_SUMMARY.md` also described as a macrophage (its ANTXR1-high
   rank is explained by it being a fibroblast). **Fixing those is out of scope
   for this plan but must be tracked** — add a note to `ANALYSIS_SUMMARY.md`'s
   open-items list as part of this work, since leaving it uncorrected would
   propagate a wrong claim.

## Target cells

Source: `/data/ANTXR2/raw/gut_cell_atlas/19053a82-9c89-4fb8-bd19-d7b1800b0b7b.h5ad`
(11.6 GB, 1,596,200 × 18,370). Filter `study == "Elmentaite2021"` AND
`cell_type` in:

| cell type | cells | donors | role | ANTXR2 raw mean |
|---|---:|---:|---|---:|
| fibroblast | 42,278 | 41 | canonical: ANTXR1 + ANTXR2 + MRC2 all present | 0.224 |
| enterocyte | 35,062 | 36 | the unexplained case: ANTXR2 present, ANTXR1/MRC2 absent | 0.078 |
| macrophage | 4,694 | 41 | myeloid calibration control | 0.116 |

**~82,034 cells total.** All three clear memento's 0.07 `filter_mean_thresh`;
enterocyte only just, so treat a null result there as underpowered rather than
as evidence of absence.

## Output

`/data/ANTXR2/coexpression/elmentaite2021_trio.h5ad`, containing:

- **`X` = raw integer counts** copied from the source's `raw/X` (verified
  integer; the source's `X` is log-normalized and must not be used — memento
  models the count-sampling process). Keep float32 as-is.
- **`var`** = the source `var`, unmodified (18,370 genes, Ensembl index with
  `feature_name` carrying symbols — note this differs from the skin atlas,
  where the index *is* the symbol, so `GENE_PANEL` must be mapped through
  `feature_name`). Do not subset genes: Phase 2 needs an open gene scope, and
  the file is small enough regardless.
- **`obs`** = all 39 source columns for the selected cells, **with unused
  categoricals dropped** (`.cat.remove_unused_categories()` on every categorical
  column). This is not cosmetic: `compute_skin_means.py:133` documents that
  leftover unused categories leak into memento's `create_groups` as phantom
  empty groups.
- **Plus the two new capture-rate columns** (below).
- **Deliberately excluded**: `obsp` (the 1.6M × 1.6M neighbour graph — this is
  the exact structure that OOM-killed the Tabula Sapiens run, see
  `compute_means.py`'s docstring), `uns` colour arrays (their lengths are tied
  to full-category counts and would be wrong after subsetting), `varm`/`varp`
  (empty anyway). **Keep `obsm`** (`X_umap`, `X_scANVI`) — trivially
  row-sliceable, 22 floats/cell, useful for sanity plots.
- Written with `write_h5ad(..., compression="gzip", compression_opts=4)` to
  match the source's compression profile.

Expected size: ~82k cells × ~1,453 nnz/cell ≈ 1.4 GB in memory, well under the
29 GB available; on disk roughly 0.5–1 GB gzipped. No chunking needed — this
fits comfortably in one pass, unlike the Phase 1 pipelines.

## The two capture-rate columns

Both are written as per-cell `obs` columns (memento's `q_column` takes per-cell
values and asserts `max() < 1` — see `memento/main.py:88`). Values are constant
within an `assay`; a per-cell refinement scaled by each cell's own `n_counts`
is noted as an option but not the default.

### Column 1 — `capture_rate_chem` (chemistry × saturation)

```
capture_rate_chem = q_chem[assay] × detected_fraction(saturation)
```

- **`q_chem`**: memento's own broad rate, `0.07`, from
  `publication/cellxgene/make_cube.py:84` (`Q = 0.07  # RNA capture efficiency
  depending on technology`) — the same value their CELLxGene-wide pipeline and
  their PBMC *co-expression* analysis (`publication/other/ifn_pbmc/interferon_2d.py`)
  use. Stored as a per-assay dict so 3′ v2 and 5′ v2 can diverge later; both
  default to 0.07, since **no published per-chemistry q exists** — documented
  explicitly rather than invented. (memento's `wrappers.py` uses 0.1; 0.07 is
  the more specific choice and the one used for their 2D work.)
- **`saturation`**: a single named constant `ASSUMED_SATURATION = 0.85` in
  `config.py`, clearly commented as an **assumption, not a measurement**.
  Justification and its limits: the paper never states saturation; ENA has
  `read_count = 0` for all 89 runs of E-MTAB-9543 (the adult GEX accession);
  only E-MTAB-8901 (developing gut, HiSeq 4000, median 352M read pairs/run,
  ~8,000 cells/reaction → ~44K reads/cell) supports a ~80–85% estimate.
- **`detected_fraction`**: the Poisson correction, per your call — *not* the
  literal product. 10x's "sequencing saturation" is the read-duplicate
  fraction, so it must be inverted through the read→molecule relation:
  ```
  saturation s = 1 - (1 - e^-λ)/λ      # solve numerically for λ
  detected_fraction = 1 - e^-λ
  ```
  At s = 0.85 this gives λ ≈ 6.65 and detected ≈ 0.9987, so
  `capture_rate_chem ≈ 0.0699`. **Expect the saturation term to be nearly
  inert** — that is the honest result at high saturation, and the code should
  print λ and the detected fraction so it's visible rather than buried. (The
  literal product would have given 0.0595, ~15% lower.)
- Implement the solve with `scipy.optimize.brentq` over λ ∈ (1e-6, 500);
  store λ and `detected_fraction` in `uns` for auditability.

### Column 2 — `capture_rate_pbmc` (public-PBMC comparison)

Logic: for a *given cell type*, UMIs detected scale with capture efficiency, so
the ratio of UMI depth between our data and a reference of known q transfers
that q:

```
capture_rate_pbmc[assay] = q_ref × median( medUMI_ours[gate, assay] / medUMI_ref[gate] )
```
taken over matched immune gates.

- **Reference**: 10x's own public PBMC datasets, chemistry-matched. `pbmc8k`
  (Chromium 3′ v2, 8,381 cells, 93,552 mean reads/cell) is the confirmed 3′ v2
  match — filtered matrix and `web_summary.html` both public at
  `cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/`. The 5′ v2 counterpart
  must be selected at implementation time from 10x's dataset catalogue (a 5′ v2
  healthy-donor PBMC GEX run); if no clean 5′ v2 PBMC dataset exists, fall back
  to applying the 3′-derived scaling to both and say so in the column's
  documentation.
- **`q_ref` = 0.07**, memento's broad rate for droplet data (same source as
  above), i.e. we assume the 10x reference sits at memento's nominal droplet
  capture and scale our data relative to it.
- **Matched gates**, chosen to avoid needing full PBMC annotation — simple
  marker-positive gating applied identically to both datasets: T (`CD3E`+),
  B (`MS4A1`+), NK (`NKG7`+ `CD3E`−), monocyte/macrophage (`LYZ`+ `CD68`+).
  Our side already has curated labels for these, so the gate is a
  cross-check against the label rather than the only evidence.
- **Document the main limitation prominently**: this assumes a given immune
  cell type has the same absolute mRNA content in gut tissue as in blood.
  Tissue-resident lymphocytes are plausibly more activated and RNA-richer than
  circulating ones, which would inflate our apparent q. This is why it's a
  second, independent column rather than a replacement for column 1 — the two
  disagreeing is informative.
- Reassuring prior: our immune populations already sit in the normal 10x v2
  range (naive B 3,917 / 1,894 medUMI for 5′ / 3′ v2; CD8 memory T 3,361 /
  2,101), so a scaling factor far from 1 would itself be a red flag.

**Expect the two chemistries to differ.** 5′ v2 yields ~1.7–2× the UMIs of 3′ v2
in the *same* cell types from the *same* study — enterocyte splits 21,371 (5′) /
13,691 (3′), so this is not a corner case. Column 2 will capture that
difference; column 1 will not (both get 0.07). That divergence is expected and
worth surfacing in the summary output, not smoothing over.

## Implementation

New files in `scripts/`:

- **`capture_rate.py`** — pure functions, no I/O, so the math is testable in
  isolation:
  - `saturation_to_detected_fraction(saturation) -> (lambda, detected_fraction)`
    (brentq solve, described above)
  - `chem_capture_rate(assay, saturation) -> float`
  - `pbmc_scaled_capture_rate(our_umi_by_gate, ref_umi_by_gate, q_ref) -> float`
  - `gate_cells(adata_or_counts, var_names) -> dict[gate, bool mask]` — the
    marker gating, shared by both datasets so the definition can't drift.
- **`download_pbmc_reference.py`** — fetches the 10x public PBMC filtered
  matrices + `web_summary.html` into `/data/ANTXR2/raw/pbmc_reference/`,
  following the resumable-download + manifest pattern already in
  `download_skin_data.py` (md5, size, URL, timestamp recorded).
- **`subset_coexpression_data.py`** — the main entry point: opens the source
  with the established `h5py` + `anndata.io.sparse_dataset()` /
  `read_elem()` idiom (**never** `read_h5ad(backed="r")` — see
  `compute_means.py`'s docstring for why), selects rows, builds the AnnData,
  attaches both capture-rate columns, writes the h5ad and a README.

Additions to `scripts/config.py` (no rewrites):
`GUT_ATLAS_H5AD`, `COEXPR_DIR = "/data/ANTXR2/coexpression"`,
`COEXPR_STUDY = "Elmentaite2021"`, `COEXPR_CELL_TYPES` (the three),
`Q_CHEM_BY_ASSAY` (defaulting both to 0.07), `ASSUMED_SATURATION = 0.85`,
`PBMC_REFERENCE_URLS`. Leave `DEFAULT_CAPTURE_RATE = 0.1` alone — it belongs to
the Phase 1 pipeline and changing it would silently alter that pipeline's
meaning.

Reuse rather than reimplement: `peak_rss_gb()` / `available_mem_gb()` from
`compute_means.py` for the progress logging convention. `build_donor_chunks` is
**not** needed here (82k cells fit in one pass) — and note it yields
largest-group-first, which would scramble row order.

A short `README.md` alongside the output h5ad, documenting: source file and
filter, the mislabelled-macrophage correction, both capture-rate derivations
with their assumptions stated as assumptions, and the tissue/chemistry
composition — following the precedent of
`/data/ANTXR2/celltype_expression/README.md`.

## Verification

Run `python subset_coexpression_data.py` and confirm from its printed report:

- **Shape**: 82,034 ± small cells × 18,370 genes; per-cell-type counts match
  the table above (42,278 / 35,062 / 4,694).
- **Counts are raw**: `X` values integer-valued, max > 20, min non-zero = 1;
  row sums match the source `obs['n_counts']` (verified upstream: median 7,031
  vs 7,033 — should match exactly per cell now).
- **Marker sanity, which is what caught the curation error**: fibroblast
  COL1A1 ~90% positive / PTPRC ~0%; macrophage C1QA ~96% / CD68 ~82% /
  COL1A1 ~2.5%; enterocyte EPCAM ~86%. If the macrophage rows come back
  COL1A1-high, the wrong label was selected.
- **ANTXR2 raw means** reproduce 0.224 / 0.078 / 0.116 (fibroblast /
  enterocyte / macrophage) — direct check against this session's numbers.
- **Capture-rate columns**: both present, `0 < q < 1` everywhere (memento
  asserts `max() < 1`); `capture_rate_chem` ≈ 0.0699 for both assays;
  `capture_rate_pbmc` differs between 5′ v2 and 3′ v2 by roughly the observed
  ~1.7–2× UMI ratio. Print both, plus λ and detected-fraction.
- **Categoricals**: no categorical column retains unused categories
  (`all(len(col.cat.categories) == col.nunique())`).
- **Round-trip**: re-open the written file with `anndata.read_h5ad` (safe now —
  no `obsp`) and confirm shape, `obs` columns, and that `obsp` is absent.
- **Donor overlap preserved**: the 36 donors shared across all three cell types
  are still shared in the output.
- **memento smoke test** (cheap, catches integration problems early): run
  `setup_memento(adata, q_column="capture_rate_chem")` +
  `create_groups(adata, label_columns=["cell_type"])` +
  `compute_1d_moments(...)` on the subset and confirm it completes and returns
  three groups. Do **not** run `compute_2d_moments` yet — gene-pair selection
  is the next planning step, not this one.

## Explicitly not in scope

- No correlation / 2D moments / `binary_test_2d` — this is preparation only.
  (For reference, the API is `compute_2d_moments(adata, gene_pairs)`, which
  requires explicit pairs — there is no all-pairs default — plus
  `ht_2d_moments` / `get_corr_matrix`; gene-pair scoping deserves its own plan.)
- No Kong2023, and no disease contrast. Elmentaite2021 has 29,594
  Crohn-disease cells but *every* sample is `Non_pathological` — Crohn donors'
  non-inflamed tissue, thin per cell type (enterocyte 3,020 / macrophage 1,269 /
  fibroblast 598). Kong2023 is the place for inflamed-vs-healthy, within-study,
  once the method is calibrated on clean tissue.
- No colonocyte: 11 donors, only 9 shared with enterocyte, and 92% 5′ v2 vs
  enterocyte's 61% — a small-vs-large-intestine contrast there would be partly
  a chemistry contrast.
- No re-alignment of FASTQs.
- No change to `DEFAULT_CAPTURE_RATE` or any Phase 1 output.
