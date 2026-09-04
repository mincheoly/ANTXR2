# Plan: ECM-clearance-pathway heatmap with GAPO/HFS tissue-involvement annotation

## Context

Over this session we built a chain of analysis on ANTXR2 (CMG2, causes hyaline
fibromatosis syndrome / HFS) and its paralog ANTXR1 (TEM8, causes GAPO syndrome):
a co-expression scatter of the two paralogs across cell types, a check of the
ECM-clearance-machinery gene MRC2 (which turned out to be essentially absent in
gut epithelium despite robust ANTXR2 there), and a clinically-grounded correction
of the initial "ANTXR1 backs up ANTXR2" framing — redundancy requires *shared
function*, which we only have weak evidence for, and skin fibroblasts (the
tissue with HFS's most visible pathology) actually have *robust* ANTXR1
co-expression, which is evidence against naive redundancy conferring protection,
not for it.

The user now wants this tied together into one durable figure: cell types /
tissues, annotated with which disease (GAPO, HFS) clinically affects them,
against expression of ANTXR1, ANTXR2, and the ECM-clearance panel
(MRC2/CTSB/CTSK/MMP14/TIMP2/LAMP1). The data (both mean-expression atlases) is
already downloaded and computed — this is a pure query + visualization task, no
new pipeline run needed. Given the corrected understanding from this session,
**the figure must not smuggle the redundancy/vulnerability framing back in** as
a computed score — it shows raw ingredients side by side and lets the reader
compare, consistent with what we agreed after the correction.

## Data sources (already computed, verified — see `ANALYSIS_SUMMARY.md` and
`/data/ANTXR2/celltype_expression/README.md` for full provenance/caveats)

- **Whole-body atlas**: `/data/ANTXR2/celltype_expression/combined_celltype_means.parquet`
  (351M rows; columns `collection_name, dataset_id, donor_id, cell_type, gene,
  n_cells, mean_expression, ...`). 238 `cell_type` values. Binding aggregation
  method (already implemented in `scripts/query_antxr2_table.py`, reuse this
  exact pattern): **two-step donor-equal-weighted** — per-donor
  n_cells-weighted mean across that donor's own rows for a `cell_type`, then
  unweighted average across donors. Never replace with a single-step
  n_cells-weighted pool.
- **Skin fibroblast atlas**: `/data/ANTXR2/celltype_expression/skin_fibroblast_celltype_means.parquet`
  (64M rows; columns `sample_id, patient_status, lesional_status,
  fibroblast_subtype, gene, n_cells, mean_expression, ...`). 13
  `fibroblast_subtype` values. Aggregation (pattern already used repeatedly
  this session): filter `patient_status=="Healthy" & lesional_status=="Nonlesional"`,
  then simple mean of `mean_expression` grouped by `fibroblast_subtype x gene`
  (rows are already at `sample_id` granularity, so this is directly
  sample-equal-weighted).
- Query pattern for both (reuse, don't reinvent): `pyarrow.parquet.ParquetFile`,
  scan row groups filtering the `gene` column via `pyarrow.compute.is_in(col,
  value_set=pa.array(genes))` — **must be `pa.array`, not `pd.array`**, the
  latter raised `TypeError` earlier this session — concatenate only matching
  row groups, aggregate in pandas. ~1-2 min per full-file scan for an 8-gene
  panel.
- Environment: conda env `antxr2` has matplotlib 3.11.1 and seaborn 0.13.2
  installed; **plotly is not installed**. No plotting code exists anywhere in
  this repo yet — this is the first. Output goes under
  `/data/ANTXR2/figures/ecm_clearance_heatmap/`, matching the existing
  convention of generated data artifacts living under `/data/ANTXR2/...`
  rather than in the git-tracked repo.

## Gene panel (8 genes, fixed scope — do not add ACTA2/POSTN/TGFB1)

`ANTXR1, ANTXR2, MRC2, CTSB, CTSK, MMP14, TIMP2, LAMP1` — exactly the
"ECM clearance pathway" panel already used together this session. Add as
`GENE_PANEL` in `scripts/config.py`.

## Row curation

### Whole-body panel (curated subset of the 238 cell types, not all 238)

Per your answer, thin site-specific GI splits (duodenum/ileum/jejunum
enterocyte, each 1-4 donors) are **excluded** — only well-powered generic
labels shown. Categories (row order = category block order below, then
descending ANTXR2 within block):

- **HFS — nodular/fibrotic (skin/gingiva/perianal)**: `fibroblast`,
  `fibroblast of gingiva`, `myofibroblast cell`, `adventitial cell`,
  `keratinocyte`, `melanocyte`, `sebocyte`
  *(myofibroblast and adventitial cell moved here from a generic "other
  organs" bucket — myofibroblasts are the classic fibrosis-effector cell,
  directly relevant to nodule formation, not muscle)*
- **HFS — GI tract (diffuse hyaline deposition + functional phenotype, NOT
  nodular — see caption note below)**: `enterocyte`, `colonocyte`,
  `intestine goblet cell`, `intestinal crypt stem cell`, `paneth cell`,
  `intestinal tuft cell`, `M cell of gut`,
  `gastrointestinal tract (lamina propria) macrophage`,
  `interstitial cell of Cajal`, `enteric neuron`, `enteroglial cell`
- **HFS — other organs (severe/ISH-form hyaline deposition: muscle, lymphoid
  stroma)**: `smooth muscle cell`, `skeletal muscle satellite stem cell`,
  `fibroblastic reticular cell`, `follicular dendritic cell`
- **GAPO-relevant**: `corneal epithelial cell`, `retinal pigment epithelial cell`,
  `retinal blood vessel endothelial cell`, `endothelial cell`
  *(eye + vasculature proxies; `fibroblast` above is also the GAPO mechanistic
  proxy — actin cytoskeleton/ECM-turnover defect reported in GAPO fibroblasts
  — noted in caption rather than duplicated as a second row)*
- **Comparison / neither disease (baseline)**: `tissue-resident macrophage`
  *(CTSB/CTSK/LAMP1/MRC2 are canonically myeloid-high — this row shows that
  baseline explicitly, so high clearance-gene signal there isn't misread as
  disease-specific)*, `hepatocyte`, `pancreatic acinar cell`, `naive B cell`,
  `neuron`

~30 rows total. List lives in a new `scripts/cell_type_curation.py`, kept
separate from query/plot code so it's easy to review/diff on its own.

Category assignment above still governs the **row color-strip annotation**,
but per your feedback, no longer dictates **row order** — see "Row ordering:
hierarchical clustering" below, which replaces the earlier
category-block-then-magnitude-sort plan for both panels.

### Skin panel (all 13 subtypes)

Disease-emergent annotation — **resolved this session** via a targeted check
against the published paper (not the original naming-guess): the paper states
three fibroblast populations have **no healthy-skin counterpart** —
inflammatory myofibroblasts, (plain) myofibroblasts, and fascia-like
myofibroblasts. Mapped onto our object's actual labels (which use a different
F-number scheme than the published F1-F8):

- **Disease-emergent**: `F6: Myofibroblast`, `F6: Inflammatory myofibroblast`,
  `F7: Fascia-like myofibroblast`
- **Healthy (paper's F1-F5, our F1-F5 + F2/3 + three F4 variants)**:
  everything else — `F1: Superficial`, `F2: Universal`,
  `F2/3: Stroma_PPARG+`, `F3: FRC-like`, `F4: DP_HHIP+`, `F4: DS_DPEP1+`,
  `F4: TNN+COCH+`, `F5: NGFR+`, `F5: RAMP1+`
- **`F_Fascia` — genuinely unresolved, flag distinctly (not coerced to either
  category)**: the paper states healthy fascial fibroblasts were merged into
  `F2: Universal`, not kept as their own group — yet our downloaded object has
  a separate `F_Fascia` label with its own cells (2,978 cells, 3 donors). This
  is a real discrepancy between the object's labeling and the secondary source
  used to check it (a PMC full-text extraction, not the primary supplementary
  table, which wasn't accessible). Render `F_Fascia`'s annotation as
  "unresolved" (e.g. a third, visually distinct hatch/gray marker), not silently
  as healthy or disease.

## Row ordering: hierarchical clustering (per your feedback)

Replaces the earlier plan of category-block order + within-category magnitude
sort (whole-body) / fixed F1→F7 paper order (skin). Applied independently to
each panel, on **rows only** — column order stays fixed (`ANTXR1, ANTXR2`
first, bolded, then the clearance panel in a set order) since that pairing is
itself part of what the figure is meant to make legible, not something to let
clustering scramble.

- Distance/linkage: `scipy.cluster.hierarchy.linkage` on each row's
  log10-transformed 8-gene vector (the same log matrix used for the heatmap
  cells, not raw values — puts all genes on a comparable scale before
  clustering), **average linkage, correlation distance** (clusters by
  *shape* of the expression profile across genes, not raw magnitude — two
  cell types that are both "clearance-gene-high, paralog-low" cluster
  together even if their absolute expression levels differ) — call out this
  choice explicitly in a caption/code comment since it's a real methodological
  decision, not the only reasonable one (euclidean distance would instead
  group primarily by overall magnitude).
- Leaf order from `scipy.cluster.hierarchy.dendrogram(linkage_matrix,
  no_plot=True)["leaves"]` reorders the matrix rows — used only to determine
  row order, per your feedback no dendrogram is drawn on the figure itself
  (no extra `GridSpec` column for it).
- The category color strip and the disease-emergent skin annotation are
  **still computed and shown** (as the adjacent color-strip column) — they
  just no longer determine position. Where clustering pulls a "HFS-GI" row
  next to a "comparison/baseline" row, that's real signal (their expression
  *profiles* are similar) worth seeing, not something to suppress by forcing
  category blocks.
- New function in `plot_ecm_clearance_heatmap.py`:
  `cluster_row_order(log_matrix) -> list[str]` (row labels in leaf order),
  called by both `build_whole_body_matrix()` and `build_skin_matrix()` after
  the log transform, before drawing.
- This also resolves a minor tension in the original plan: the skin panel's
  "keep it in the paper's own F1→F7 order" reasoning was about aiding
  interpretation with only 13 rows — clustering is a strictly more informative
  version of that same goal (data-driven similarity instead of a label-number
  ordering) and is now applied the same way in both panels, for consistency
  of method.

## Coverage-gap handling

**Footnote/caption text, not placeholder rows** (avoids the figure implying
"measured near-zero" where the truth is "not measurable at all"). Exact text,
rendered on the figure via `fig.text` and echoed in both output CSVs' header
comments:

> Not represented in this atlas (no matching Cell Ontology label): bone/
> cartilage/joint/synovium (GAPO metaphyseal dysplasia; HFS joint contractures
> and osteopenia/osteolysis), teeth/odontogenic tissue (GAPO pseudoanodontia),
> thyroid and adrenal gland (HFS-ISH severe-form hyaline deposition sites),
> gonadal tissue (GAPO hypogonadism). Eye/vascular rows shown are proxies, not
> optic-nerve- or scalp-vein-specific.

## Figure design

- One `matplotlib` figure, **two panels stacked vertically** (whole-body on
  top, ~30 rows; skin below, 13 rows) — not side-by-side, since the two
  panels' very different row counts don't share a natural aspect ratio the way
  the earlier interactive scatter's two roughly-square panels did. Each
  panel's height sized proportionally to its own row count so labels stay
  legible in both.
- **Two independent color scales**, one per panel, each with its own colorbar
  and its own log10 pseudocount floor (`min(nonzero in that panel)/2`, same
  technique validated for the earlier scatter plot — note: that plot's exact
  file no longer exists on disk from this session, so the floor is freshly
  derived here, not reused verbatim). Rationale: the README's own caveat
  (flat, uncalibrated capture rate across all assays) means cross-panel
  absolute-scale comparison isn't well supported by the data — a shared scale
  would visually imply a comparability that doesn't hold. Caption states this
  explicitly: "color scales are independent per panel; compare patterns within
  a panel, not color intensity across panels."
- Each panel's row axis, left to right: categorical color strip (whole-body:
  5 categories; skin: healthy / disease-emergent / unresolved) → main heatmap
  cells — no dendrogram drawn. Legend for the color strip; row order in both
  is the clustering leaf order (see "Row ordering" above, computed but not
  visualized), not category or a fixed reading order.
- Column order: `ANTXR1, ANTXR2` first (bolded tick labels, visually anchoring
  the pair being compared), then `MRC2, CTSB, CTSK, MMP14, TIMP2, LAMP1`.
- Row tick labels include sample size, e.g. `"fibroblast of gingiva (n=62)"` —
  this is a static image, so sample size has to be visible directly rather
  than in a hover tooltip.
- Skin panel only (13 rows, room for it): numeric value annotated directly in
  each cell. Whole-body panel (~30 rows): no in-cell numbers, color + colorbar
  carries the signal.
- **Hard constraint, carried through every function in the implementation, not
  just stated in prose**: no function computes or stores an ANTXR2:ANTXR1
  ratio, "redundancy score," "vulnerability score," or any other composite
  cross-gene metric for display or in the output CSVs. The figure shows raw
  per-gene expression and independently-sourced disease annotations only. A
  caption line states this directly: *"Descriptive juxtaposition of raw
  expression values and independently-sourced disease-tissue annotations; not
  a composite risk, redundancy, or predictive score."* (The one narrow
  exception: the verification step below computes the ratio transiently,
  printed to stdout only, purely to cross-check against this session's earlier
  scatter-plot finding — never written to a file or plotted.)
- Save `heatmap.png` + `heatmap.svg`, plus `whole_body_curated.csv` and
  `skin_curated.csv` (long format: cell type/subtype, category/annotation,
  gene, mean_expression, log10_mean_expression, n_donors or n_samples,
  n_cells) — kept as two separate CSVs, matching the source parquets' own
  "kept in its own file" convention.

## Implementation

New files in `scripts/`:
- **`cell_type_curation.py`** — pure data: `WHOLE_BODY_CATEGORIES` dict,
  `SKIN_SUBTYPE_ORDER` list, `SKIN_DISEASE_EMERGENT` dict (with the
  `F_Fascia: None` unresolved marker), `COVERAGE_GAP_NOTE` string. No I/O.
- **`gene_panel_query.py`** — library module:
  `scan_row_groups_for_genes(path, genes, gene_col="gene") -> pd.DataFrame`
  (generalizes `query_antxr2_table.py`'s scan loop to an arbitrary gene list;
  optional `cache_path` to skip re-scanning the 3GB/400MB files on repeat runs
  during figure-layout iteration);
  `aggregate_whole_body(df, cell_types=None) -> pd.DataFrame` (the two-step
  donor-equal-weighted method, generalized to multiple genes, pre-filterable
  to the curated list);
  `aggregate_skin(df, subtypes=None) -> pd.DataFrame` (healthy/nonlesional
  filter + simple mean, generalized to multiple genes).
- **`plot_ecm_clearance_heatmap.py`** — `__main__` entry point:
  `build_whole_body_matrix()`, `build_skin_matrix()` (query + aggregate +
  pivot to wide `cell_type x gene`, return `(matrix, row_meta)`);
  `log_transform(matrix) -> (log_matrix, floor)`;
  `cluster_row_order(log_matrix) -> list[str]` (average-linkage, correlation
  distance, see "Row ordering" above — returns leaf order only, no dendrogram
  object retained for drawing); `draw_heatmap_panel(...)` (shared low-level
  drawing routine for both panels — color strip, main cells, colorbar, tick
  labels; no dendrogram axis); `main()` wires it together, writes the two
  CSVs, builds the figure, prints the verification checks below.
- `config.py` gets one addition: `GENE_PANEL = [...]`. `query_antxr2_table.py`
  is left untouched (it's a working, already-delivered single-gene script;
  the new module generalizes its pattern rather than editing it).

**Before finalizing colors** (implementation-time step, not solved in this
plan): reload the `dataviz` skill, follow `choosing-a-form.md` (heatmap =
sequential single-hue for the magnitude cells — this is a magnitude-by-category
matrix, the textbook case), pull the sequential ramp and the categorical
palette from `color-formula.md`/`palette.md` for the row color strips, and run
`scripts/validate_palette.js` on the final categorical hex values before
locking them into `cell_type_curation.py`.

## Verification (printed at the end of `main()`, not silently trusted)

- Exactly 8 distinct `gene` values returned from each parquet scan (catches a
  silent symbol mismatch).
- MRC2 near-zero and comparable in magnitude across all GI-category
  whole-body rows and across all 13 skin subtypes — should match this
  session's finding; print MRC2's min/max within each panel.
- ANTXR2:ANTXR1 ratio per skin subtype, printed only, spot-checked against
  "< 1 in 9 of 13 subtypes" from the earlier scatter (never written to a file
  — see the hard constraint above).
- Every row has non-zero, non-NaN `n_donors`/`n_samples` — print any row that
  doesn't (would indicate a cell-type-name typo against the atlas's actual
  238-value vocabulary).
- Pre-/post-log-transform min/max per panel, to confirm the floor isn't
  compressing real dynamic range to invisibility.
- Programmatic re-check of the coverage-gap claim: substring-match the full
  238-value vocabulary against `chondro|osteo|thyroid|adrenal|dental|tooth|
  synov|gonad` to confirm nothing relevant was missed (cheap insurance, since
  this was already checked by hand this session).
- `n_donors` should be constant across all 8 genes within a given whole-body
  row (same donor cohort contributes every gene for a fixed cell_type) — print
  a warning if not, since it would mean a gene is silently missing from some
  donor's data.

## Not in scope for this figure (explicitly, so it isn't added by accident)

- No composite redundancy/vulnerability score (see hard constraint above).
- No ACTA2/POSTN/TGFB1 (activation-axis genes) — different panel, different
  question, not requested here.
- No placeholder rows for tissues the atlas doesn't cover — footnote only.
- No Phase 2 (correlation/significance testing) — this is a means-level
  descriptive figure, consistent with everything computed so far.
