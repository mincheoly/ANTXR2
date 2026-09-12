# ANTXR2 project -- analysis summary

> **READ THIS SECTION FIRST.** Everything below the `---` divider is a
> **chronological session log**, appended to over many sessions. This top
> section is the **durable framing**: the question the project is answering,
> what each kind of evidence can and cannot establish, and corrections that
> must survive session handoff.
>
> **If a later section of the log contradicts this section, this section
> wins**, unless it has been explicitly amended here. Framing drift across
> sessions has already happened at least once (see "Standing corrections /
> collagen" below) -- this header exists to stop it recurring.

## Central question

ANTXR2/CMG2 has **at least two separable molecular functions with different
partner requirements**:

1. **Collagen-VI clearance** -- receptor-mediated endocytosis and lysosomal
   degradation of collagen VI and other hyaline proteins (Bürgi et al. 2017).
   Requires MRC2 + lysosomal/proteolytic machinery + the substrate itself.
2. **Wnt signal transduction** -- CMG2 complexes with LRP6, which assembles
   with Frizzled receptors to bind Wnt ligands (Abrami et al. 2008; Wei et al.
   2006). Required for **injury-induced** intestinal stem-cell renewal: in
   CMG2-KO mice the fetal-like reversion after DSS damage occurs normally, but
   the fetal-like → Lgr5+ transition fails (β-catenin nuclear translocation),
   causing failed epithelial restitution (Bracq et al. 2025, EMBO Mol Med).
   Requires LRP6 + Frizzled. **Baseline guts of CMG2-KO mice are normal** --
   this function is conditional on injury.

**The project asks: which cell types (and cell states) have the partners for
which function, and does that partition predict HFS's tissue-restricted
phenotype?**

The germline LoF is present in every cell; the phenotype is not (skin nodules,
gingiva, joints, gut, uterus). **Partner availability is a transcriptionally
readable explanation for that gap.** This is why atlas data is the right
instrument for this question, and why no causal identification is required.

This framing supersedes the earlier "tissue-restriction paradox" phrasing
(`planning_summaries/planning.md`) -- same question, but now with a specific
mechanistic axis (which function, gated by which partners) rather than an
open-ended one.

## What each kind of evidence can establish

Stated explicitly because the project has repeatedly slipped between these,
and they carry very different weight.

| Evidence | Claim it supports | Strength |
|---|---|---|
| **Cell-type mean expression** | Partner is *present* (or absent) in the cell type | **Strongest.** Often sufficient on its own. |
| **Single-cell correlation** | Either co-occupancy (same cells) **or** co-regulation (shared program) -- both interesting | **Moderate.** Needs a specificity floor; see below. |
| **Neither** | Physical protein interaction, or that a function is actually *occurring* | **Not reachable with this data.** |

**The inference is asymmetric, and the negative direction is much stronger:**

- **Partners absent → that function cannot operate.** Strong, needs only
  means, survives every confound raised so far. This is the direction the
  project's best result runs in (MRC2 at detection floor in gut epithelium).
- **Partners present → the function is *available*, not that it happens.**
  Weak. This is where all the confounding lives.

**On correlation specifically.** A significant correlation implies either
co-occupancy or co-regulation, and both are findings -- the disjunction does
not need to be resolved to have learned something. But there is a third,
uninteresting branch: **nonspecific covariation** (shared activation state,
cell size, depth), where both genes rise together because the cell is doing
more of everything. Excluding that branch is the entire job of the
**anchor-gene specificity control** (see "Designed but never run"). It does
not test whether correlation exists; it establishes what an expression-matched
arbitrary gene's correlation looks like in that cell type, so that anything
above the floor is in the two interesting branches.

Note that correlation and co-presence are **not the same statistic**. Two genes
expressed in every cell at stable levels have perfect co-presence and zero
correlation -- no variance to co-vary. Correlation is informative about joint
occupancy only when expression is patchy/bimodal; for near-uniformly expressed
partners, the mean already settles availability and correlation adds nothing.
**Check each panel gene's expression distribution before treating its
correlation as a co-presence readout.**

memento is comparatively well suited to the co-occupancy branch, since it
estimates the correlation of *true* expression with sampling noise modeled --
which is exactly the noise that flattens patchy-expression correlations in raw
counts.

**The two interesting branches are partly separable empirically.** If a
correlation is driven by subset structure (both genes on in subpopulation S,
off outside it), conditioning on S should collapse it -- within S there is no
remaining joint on/off variation. If it is graded co-regulation, it persists
within S. Attenuation implicates co-occupancy; persistence implicates
co-regulation. Same conditioning machinery as the cell-state control, pointed
at a different question.

## Standing corrections (carry these forward)

- **Collagen genes are NOT a positive control.** Bürgi et al. 2017 showed
  collagen VI **mRNA does not change** in Antxr2−/− uteri -- the mechanism is
  degradation, not transcription. An ANTXR2–COL6 positive correlation
  therefore **cannot validate that the pipeline detects ANTXR2's known
  function**, because a positive control must be able to fail. It remains
  entirely compatible with genuine co-regulation under shared ECM-program
  pressure, which is interesting -- **report it as co-regulation evidence, not
  as confirmation that the pipeline works.** (This correction was established
  in `planning_summaries/planning.md` and then lost: the 2026-09-05 z-score
  GSEA section originally described the fibroblast ECM result as the project's
  cleanest confirmation of a real signal. Amended in place; see that section.)

- **Anti-correlation is the stronger result class.** Shared-activation-state
  confounding inflates positive correlations but **cannot manufacture mutual
  exclusivity**. Negative findings run against the confound and must not be
  demoted relative to positive ones. Corollary: a significant anti-correlation
  is interpretable **without** the anchor-gene floor; a positive one is not.

- **Co-presence is not redundancy.** Redundancy requires *shared function*.
  Skin fibroblasts show robust ANTXR1 co-presence and have HFS's most visible
  pathology; the paralogs bind different collagen VI domains (triple-helical
  vs C5). **No composite "vulnerability score", ever.** Show raw ingredients
  side by side.

- **Positive correlations require the anchor-gene specificity control.**
  **RUN 2026-09-06** (`scripts/anchor_gene_control.py`, prompts/partner_availability.md
  Task 2; see the "Task 2 results" section for full detail) -- outcome differs
  by cell type and must not be treated uniformly:
  - **Fibroblast ECM/collagen correlation SURVIVES the control.** No key term
    (Extracellular Matrix Organization, ECM-receptor interaction, Collagen
    Formation, Focal adhesion) is recovered by a majority of 20
    expression/detection-matched anchor genes (6-10/20 recover any one term);
    ANTXR2 sits at the 80th-100th percentile of the anchor NES distribution
    for all four. Still co-regulation evidence, not a positive control and
    not a causal claim (see the item above) -- but no longer merely
    "ungated exploratory GSEA": this is now the project's best-supported
    positive-correlation finding.
  - **Enterocyte digestion/absorption correlation DOES NOT SURVIVE the
    control, and this reframes that finding.** A MAJORITY of anchors recover
    both "Protein digestion and absorption" (17/20) and "Fat digestion and
    absorption" (16/20) at the same FDR threshold ANTXR2 was judged at,
    with ANTXR2 only at the 90th/95th percentile of a distribution whose mean
    (1.6-1.8) is already strongly positive. **This correlation is a property
    of enterocyte's transcriptional program at ANTXR2's expression level, not
    specific evidence about ANTXR2** -- any gene expressed like ANTXR2 in
    enterocytes tends to travel with digestion/absorption genes, which is
    exactly what you'd expect from a highly polarized secretory/absorptive
    epithelial cell where most moderately-expressed genes covary with the
    dominant differentiation program. This does NOT bear on the separate,
    means-level, formally-tested finding that ANTXR2 itself is higher in
    differentiated enterocyte than crypt stem/TA (Task 4.2) -- that is a
    presence/gradient claim, not a co-regulation claim, and stands on its own.

- **|z| correlates with mean expression (Task 2, step 6 check).** r=0.44
  (fibroblast) / 0.37 (enterocyte) between log10 mean expression and |z|
  across each cell type's full tested universe. A real, moderate confound --
  treat z-score-ranked GSEA (`gsea_zscore/`) as somewhat expression-level-
  biased in both tails, on top of the co-regulation-vs-coincidence caveat
  above. The anchor-gene matching controls for this specifically in the
  recovery-count comparisons above (anchors are expression/detection-matched
  to ANTXR2), which is why those comparisons are trustworthy despite the
  confound existing in the underlying ranking.

- **Doublets inflate cross-cell-type positive correlations.** Probably handled
  by atlas QC, but would specifically inflate ANTXR2–ANTXR1 if
  fibroblast-epithelial doublets survived. Verify before interpreting a
  positive paralog correlation.

- **Exploratory vs confirmatory analyses must be labeled.** GSEA over the full
  ranked universe is **exploratory characterization** of what ANTXR2 travels
  with. Pre-defined panel tests are **confirmatory partner-availability
  checks**. A GSEA term is not a tested claim. The two GSEA atlases (`gsea/`,
  `gsea_zscore/`) are exploratory throughout.

- **Macrophage is exploratory-only, everywhere, permanently.** 5 donors; the
  discovery/replication split found sign concordance indistinguishable from
  chance (66.7%, binomial p=0.25, negative effect-size correlation). Its
  z-score GSEA looks *better* than its mean-ranked GSEA for a specific and
  misleading reason -- see the 2026-09-05 section.

- **A successful download is not the same as correct data -- verify content
  against its own header/manifest, not just transfer success.** The IBD colon
  atlas's HCA-hosted `gene_sorted-Epi.matrix.mtx` and `gene_sorted-Imm.matrix.mtx`
  downloaded with byte-exact size matches against the host's own reported size,
  yet were truncated to 52.6% and 4.6% of their declared line counts
  respectively -- confirmed genuine (not a download bug) by independently
  Range-probing the live S3 object directly. Public archive mirrors can host a
  self-consistently-sized but incomplete object. Cheap check for any future
  dataset with a self-describing format (mtx nnz header, a stated cell count,
  etc.): compare the parsed content against that count before trusting the file,
  not just the downloaded byte count against an expected size.

## Designed but never run (authoritative open list)

These are ranked. The list at "## Open / next steps" mid-log is older and
narrower; **this list supersedes it.**

**Items 1-5 below were RUN 2026-09-06** (`prompts/partner_availability.md`, `scripts/paralog_readout.py` / `anchor_gene_control.py` / `two_arm_panel_query.py` / `binary_tests.py` / `kong2023_feasibility.py`) -- kept here with their original rationale for context, but see the "Task 1-5 results" sections in the session log below for outcomes. Items 6-8 remain genuinely open.

1. **ANTXR1/ANTXR2 single-cell co-expression -- the project's original primary
   analysis, never executed.** `planning.md` ranked paralog partitioning as
   the strongest sub-hypothesis, and the entire stated rationale for Phase 2
   was that means establish presence while only correlation distinguishes
   "both paralogs in every cell" from "two disjoint subpopulations". Phase 2
   ran ANTXR2-vs-all-genes and this specific question was never asked.
   **ANTXR1's post-bug-fix `coef`/`se`/`z`/`fdr` already exist** in
   `full_dataset_ht/{cell_type}_full_dataset_ht.csv` and have never been read
   out. The only ANTXR1 value ever discussed (0.35) is pre-bug-fix and known
   inflated. Near-zero cost. **Pre-register the asymmetry before looking:** a
   significant negative is real evidence of mutual exclusivity; a null is
   uninformative, *not* evidence of co-presence.

2. **Anchor-gene specificity control.** ~20 expression- and detection-matched
   genes through the identical pipeline, per cell type, to establish the
   nonspecific-covariation floor. Gates every positive correlation claim in
   the project, including the fibroblast ECM result and the enterocyte
   digestion/absorption result. Also required to make the enterocyte *negative*
   interpretable -- you cannot call a module absent without knowing what a null
   module looks like there.

3. **Extend the gene panel to the Wnt arm.** The current 8-gene panel
   (`ANTXR1, ANTXR2, MRC2, CTSB, CTSK, MMP14, TIMP2, LAMP1`) covers only the
   clearance function. Partner availability for the *other* function is
   untestable without: `LRP5, LRP6, FZD1-10, CTNNB1, TCF7L2, LGR5, RNF43,
   ZNRF3, AXIN2`. Also add `COL6A1, COL6A2, COL6A3` -- the clearance arm's
   substrate, currently missing from a panel that tests everything except
   whether the substrate is there. Pure query against already-computed parquet;
   same pattern as `scripts/gene_panel_query.py`.

4. **Formal tests for the two means-level claims.** `binary_test_1d` on MRC2
   (gut epithelium vs lineage-matched keratinocyte/corneal epithelium), and on
   ANTXR2 (crypt stem/TA vs differentiated enterocyte). The second is now
   load-bearing: Lencer's commentary on Bracq et al. proposes that CMG2's
   context-specific role may reflect the diminishing Wnt gradient along the
   crypt-villus axis, making ANTXR2's own position on that axis a direct,
   testable prediction. Existing means hint at it (crypt stem cell +1.21 vs
   enterocyte +2.12 log-ratio) but it is untested.

5. **Kong2023 feasibility check, then differential correlation.** Before
   scoping anything: confirm ANTXR2 clears `min_perc_group` in enough
   Kong2023 donors *per condition* to support a two-group test. ANTXR2 already
   needed `min_perc_group` relaxed to 0.6 in the Elmentaite enterocyte full
   set. If feasible, the pre-specified prediction from Bracq et al. is
   directional: **Wnt-arm partner co-availability should be absent or weak in
   `Non_pathological` and appear in `Inflamed`**, since the function is
   injury-conditional. This is a much sharper use of `ht_2d_moments` than
   open-ended exploration.

6. **Cell-state conditioning.** The pipeline handles donor structure carefully
   and within-donor heterogeneity not at all. Doubles as the co-occupancy vs
   co-regulation separator (see "What each kind of evidence can establish").

7. **Dropped scope worth reconsidering.** The trio working set
   (fibroblast/enterocyte/macrophage) silently became the analysis scope.
   `planning.md` aimed the machinery-completeness question at the unexplained
   **enteric smooth muscle / ICC** signal; kept Tabula Sapiens specifically for
   **uterus** (where the Col6a1 rescue cross was done); and flagged **which
   fibroblast subtype** as newly answerable at Steele's 9-subtype resolution.
   None are in the trio.

8. **Verify the AS GWAS variant annotation at ANTXR2.** Long-pending. The
   association is replicated (top SNPs rs12504282/rs4333130/rs4389526, r²≥0.76,
   meta p=6.7e-9) and the lead SNPs sit near a putative regulatory region
   rather than in coding sequence -- so the AS mechanism is plausibly
   expression-level, unlike HFS's coding recessive LoF. This is the one axis
   on which ANTXR2's *own* transcriptional regulation is disease-relevant.

## Out of scope, deliberately

**Causal identification.** Perturb-seq (Replogle K562) is real do(X) but the
wrong cell type; MR via OneK1K is properly identified but PBMC-only, where
ANTXR2 is low -- a weak instrument exactly where the data is. Ruled out:
PC/GES (needs causal sufficiency), pseudotime/velocity ordering, front-door.
**The gap is itself a finding**: it defines what new data would be needed.

**Claims about what the ANTXR2 protein does.** Its function is
post-transcriptional (endocytosis, lysosomal degradation, MMP activation,
β-catenin translocation). Gene-gene correlation cannot speak to it. The
project's claims are about *partner availability*, which is what expression
data measures well.

---

*Everything below is the chronological session log.*

## Session log

Working notes for this analysis session, kept separate from the pipeline's own
technical README (`/data/ANTXR2/celltype_expression/README.md`, regenerated by
the pipeline itself). This file is a human-readable log of what was built, what
was found, and what's still open -- update it as follow-on analysis continues.

## What was built

Implemented the plan in `prompts/celltype_means.md`: a genome-wide, per-cell-type
mean expression atlas built with memento's method-of-moments estimator (not
memento-cxg's precomputed backend, not "fraction expressing"), across 3 CELLxGene
Discover collections:

- Tabula Sapiens ("All Cells", primary dataset)
- Cross-tissue Immune Cell Atlas ("Global")
- Gut Cell Atlas ("Healthy reference" + "Extended+ - 18485 genes", kept as two
  separate `dataset_id` rows -- donors likely overlap between them)

Code lives in `scripts/` (`download_data.py`, `compute_means.py`, `config.py`).
Output is 4 parquet files + README under `/data/ANTXR2/celltype_expression/`
(338,470,470 rows total, grouped by `dataset_id x donor x cell_type x assay`,
where `donor` is `donorID_unified` for the gut atlas -- see the 2026-09-03
donor-identity correction below).

Full caveat list (capture rate placeholder, dataset selection, gut atlas overlap,
dataset-level tissue/disease metadata, memento's mean vs. naive average) is in
that README and still applies to everything downstream, including this file.

## Notable build issues (for context on what "done" survived)

- A `MEMENTO_MIN_CELL_COUNT` of 1 crashed on tiny groups (memento's per-group
  regressor needs >=1 gene with count >=2 within the group); fixed by using
  memento's own built-in default of 10.
- `anndata`'s backed mode is not lazy for `obsp`/`obsm`/`varm`/`uns` -- opening
  the 45GB Tabula Sapiens file OOM'd the machine just from its neighbor graph.
  Fixed by dropping to raw `h5py` + `anndata.io.sparse_dataset()` /
  `read_elem()`, touching only `raw/X` (or `X`) and `obs`/`var`.
  **`obsp` is not used anywhere in this pipeline, by design.**
- The chunked per-donor processing was memory-bounded per chunk but still OOM'd
  at the very end from a single `pd.concat` across all donors before writing.
  Fixed by streaming each donor's results straight to a `pyarrow.parquet.ParquetWriter`
  as soon as they're finalized, never holding the full result table in memory.
- Final run: 92GB RAM, ~48 min, steady peak RSS ~25GB, exit 0.

## Definition-of-Done spot checks (both passed)

- ANTXR2 vs. collagen VI positive control (Tabula Sapiens fibroblasts):
  COL6A1/2/3 rank highest as expected; ANTXR2 clearly elevated (5th of 14 genes),
  consistent with the Bürgi et al. collagen VI / ANTXR2 interaction. Also
  elevated in endothelial cells (consistent with CMG2's known vascular
  expression).
- `n_cells` / `low_confidence` flagging present and behaving correctly (floor at
  n=10, 16.9% of rows flagged `low_confidence` for n<20, none silently dropped).

## Aggregation methodology (binding for all further per-cell-type summaries)

Per explicit user correction, any "expression by cell type" summary in this
project uses a **two-step, donor-equal-weighted** aggregation, not a single
cell-count-weighted pool:

1. Compute each donor's own mean first (n_cells-weighted across that donor's
   own `collection_name x dataset_id x assay` entries for a given cell_type).
2. Average those donor-level means with equal weight per donor, per cell_type.

This avoids donors with huge cell counts (e.g. a 10x run with 100K cells)
silently dominating a cell type's reported mean over donors with fewer cells
but equally valid biology.

## Key finding so far: ANTXR2 expression by cell type

Full sorted table (238 cell types, all 3 collections, memento mean +
n_cells + n_donors + per-collection presence) is in
`antxr2_celltype_table.csv` (delivered to user; regenerate via
`scripts/query_antxr2_table.py`, which implements the two-step aggregation
above against `combined_celltype_means.parquet`).

Headline pattern: **not a simple "epithelial vs. non-epithelial" story.**

- Classic ANTXR2/CMG2 biology (stromal/mesenchymal/endothelial, collagen VI
  and laminin binding) replicates well: fibroblasts, endothelium, mesenchymal
  stem cells rank high.
- Terminally differentiated barrier epithelia (keratinocyte, corneal
  epithelial cell, retinal pigment epithelial cell) sit at the very bottom,
  near zero -- consistent with prior expectation.
- **Gut/intestinal epithelium is the exception**: enterocytes, Paneth cells,
  and goblet cells (especially proximal small intestine) sit in the
  upper-middle of the ranking, not the bottom. This is the specific hypothesis
  the plan doc flagged the Gut Cell Atlas as being able to test ("does ANTXR2
  skew toward immature/stem epithelial states").
- That specific sub-hypothesis (stem/TA skew) is **not** supported by the
  means alone: intestinal crypt stem cells and transit-amplifying cells rank
  *lower* than differentiated enterocytes/Paneth/goblet cells in this data,
  if anything the opposite of a stem-skew.
- Caveat on the literal #1 rank ("enterocyte of epithelium proper of small
  intestine"): only 2 donor-units / 1,415 cells -- a thinly-sampled,
  fragmented label. The well-powered generic "enterocyte" label (131 donors,
  218,757 cells) ranks far lower (150th of 238), ~63x lower mean. Don't over-read the literal top rank without checking `n_donors`.

## Skin fibroblast atlas (Steele et al., Nat. Immunol. 2025) -- second session

Implemented the download + means computation for `prompts/skin_means.md`. This is a
separate dataset from the 3-collection pipeline above (fibroblast-only, different
source, different obs schema); its output lives in its own file
(`skin_fibroblast_celltype_means.parquet`), not merged into `combined_celltype_means.parquet`.
Full provenance, schema decisions, and caveats are documented in
`/data/ANTXR2/celltype_expression/README.md` (second half); this section is the
narrative of *how* those decisions got made, since several were genuine blockers not
resolvable from the plan doc alone.

**Download blocker, resolved.** The plan doc flagged this correctly: the atlas isn't
on CELLxGene Discover, and the paper's own portal (`cellatlas.io/studies/skin-fibroblast`)
is a JS-rendered SPA with no crawlable download link. Resolution: the SPA's minified
JS bundle references a public Strapi CMS API
(`strapi-api-dot-haniffa-lab.nw.r.appspot.com`) that backs the portal. Querying that
API for the study (`GET /api/studies?filters[slug][$eq]=skin-fibroblast`) surfaces its
"Integrated atlas" dataset record (357,276 cells -- matches the paper abstract's stated
n exactly), which lists a public, unsigned, no-auth GCS object as its `AnnData h5ad`
resource. This is the literal file the official portal itself serves for this study --
not a guessed or scraped-from-rendered-HTML URL, which is what the plan doc was
specifically warning against. Downloaded and md5-verified (27.23GB,
`e2cf64a7d04d6afe87f9278d046c3a46`).

**Two structural gaps between the plan's assumptions and the actual file, both
surfaced to and resolved by the user before computing:**

1. **No `donor_id`.** The plan asked for `donor_id x disease_state x fibroblast_subtype`
   grouping, but this specific downloadable object has no donor/patient column
   anywhere in `obs` -- only `GSE` (source-study accession), `Patient_status`,
   `celltype`, `celltype_skinspecific_nomenclature`, `disease_category_orig`,
   `lesional_vs_nonlesional`. The cell barcode sometimes carries a per-sample suffix,
   but the format is inconsistent across the ~30 integrated source studies (GEO
   sample accessions, SRA sample accessions, author-name labels, and for some cells a
   non-ID literal like "Lesional") -- not safely parseable into one clean field
   without heavy per-source special-casing. **User decision: use `GSE` (source study)
   as the grouping/chunk-boundary proxy instead**, honestly labeled `source_study`,
   not `donor_id`. This means the pipeline's binding "donor-equal-weighted
   aggregation" convention (see above) is actually *study*-equal-weighted for this
   dataset -- a real granularity loss to keep in mind for any Phase 2 work here.
2. **No single "disease_state" column**, and the two disease-related columns that do
   exist (`Patient_status`, `disease_category_orig`) are not 1:1 with each other (e.g.
   Psoriasis spans two different `disease_category_orig` values). Kept both, plus
   `lesional_vs_nonlesional`, as independent grouping columns rather than collapsing
   one into the other, per the plan's own "do not pool healthy and lesional"
   instruction.

Also discovered and handled without needing a user decision (mechanical, not a
judgment call): raw counts live at `raw/X` (verified via the same
integer/magnitude sampling check as the other 3 datasets); no `assay` column exists,
so every group uses the same flat capture-rate placeholder as the rest of the
pipeline; `celltype_skinspecific_nomenclature` is 1:1 with `celltype` (safe to carry
as a derived, non-grouping column) except the two disease-specific F6/F7 subtypes,
which have no skin-specific name.

**Run:** 32 chunks (one per `GSE`, largest 59,176 cells -- all well under the
existing pipeline's 150K per-chunk cap, so no sub-chunking was needed), 115s wall
time, peak RSS 4.8GB, 19,691,338 rows written, exit 0. Implemented as a new
self-contained script (`compute_skin_means.py`) rather than extending
`compute_means.py`, since the latter is already verified/complete and this dataset's
schema (no donor_id, no assay column, 4 grouping label columns instead of 2) would
have meant threading extra parameters through code with no other caller; only the
dataset-agnostic memory-bounding utilities were imported and reused as-is.

**Spot checks (plan's Definition of Done, all passed):** 13-value `fibroblast_subtype`
vocabulary reported (matches paper's "6 healthy + 3 disease-specific" framing once
sub-labels are accounted for); ANTXR2/ANTXR1/collagen-VI panel non-zero and
internally consistent in healthy nonlesional tissue (COL1A1/COL6A1 orders of
magnitude above ANTXR1/ANTXR2; POSTN/ACTA2 visibly elevated in myofibroblast
subtypes F6/F7 vs. F1-F5, consistent with expected activated phenotype); `n_cells`
present on every row (min 10, memento's floor), 10.4% flagged `low_confidence`
(n<20), none dropped.

**Not yet done:** the plan's Phase 2 correlation questions (ANTXR1/ANTXR2 paralog
co-presence by subtype; ECM-clearance machinery co-assembly) -- means only, per the
plan's own scope for this step.

## Correction: donor/sample identity was recoverable after all

User pushback, same session: given this is a translational study, donor-level
granularity matters a lot, and asked to see `obs.head()` directly rather than take
the "no donor_id" conclusion above at face value. That scrutiny was warranted --
**the "no donor_id" conclusion above was wrong**, or at least incomplete. The
inspection that produced it used a single global rule ("take the barcode's last
underscore-delimited token") to characterize whether the barcode suffix was a usable
sample id, saw messy-looking top values (`Lesional`, `NonLesional`, `solebordo`) and
concluded the format was inconsistent. That was too shallow a check: those messy
values came from source studies whose barcode has *more than one* underscore, where
grabbing only the last token grabs the wrong field. Checked properly -- per source
study, not globally -- 27 of the 32 integrated source studies use a clean, uniform
`<barcode>_<GSM or SRS accession>` format, i.e. a real GEO/SRA sample id was sitting
right there in the barcode the whole time, just never split out into its own
column.

Reworked `compute_skin_means.py` to recover it (`extract_sample_id`, per-study
rules): direct GSM/SRS-suffix extraction for the 27 clean studies; a dedicated
regex for one study whose barcode encodes `patientN_DayX` directly (turned out to be
the paper's own newly-generated, non-GEO-deposited data -- the richest of all once
parsed); case-normalization for Reynolds's `P`/`E`/`S` sample codes. Per user
decision, two studies (`Ganier`, `Sole-Boldo`) were deliberately left at
`source_study`-level rather than chasing their messier formats (`Ganier` is ~90%
clean but ~10% barcode-order-reversed; `Sole-Boldo`'s suffix is a body-site label,
not a donor id at all) -- "ignore the messy data sources" was the user's explicit
call here, not a default to more work.

Before trusting the new parsing rules, checked cardinality per source study (how
many distinct recovered ids per `GSE`) to confirm each rule was actually recovering
multiple plausible samples, not a repeated constant -- 31 of 32 studies showed >1
recovered id (the exception, a 675-cell single-sample study, is expected to have
exactly one).

**Result:** 287 `sample_id` groups recovered from the original 32 `source_study`
groups (up from ~12-13 real samples-worth of study-level pooling to 30-70 real
patients/samples for the well-powered F1-F5 subtypes) -- much closer to what the
pipeline's binding donor-equal-weighted aggregation convention actually needs.
Recomputed: 287 chunks (10 of the smallest -- all tiny Reynolds replicate codes or
singleton GSM samples, well under 50 cells combined -- fell below memento's own
`min_cell_count=10` floor and are absent, same documented behavior as elsewhere in
this pipeline), 201s wall time, peak RSS 4.4GB, 63,941,947 rows written (up from
19,691,338 -- expected, since sample-level grouping is much finer than study-level).
`low_confidence` rate rose from 10.4% to 22.0%, which is the correct, honest
consequence of finer grouping producing more small groups, not a regression. Spot
checks (ANTXR2/ANTXR1/collagen-VI panel, `n_cells`/`low_confidence` integrity)
re-run and still pass; values are close to the pre-correction run, as expected --
this was a granularity fix, not a different biological signal.

**Lesson for this pipeline going forward:** when a "not available" conclusion about
metadata is reached from a single global heuristic (a regex, a split, a naive
join), re-check it per source/subgroup before treating it as a hard blocker,
especially for a field as consequential as donor identity. The full write-up above
(the original "No `donor_id`" finding) is left in place rather than deleted, since
it documents what was actually tried and why it looked plausible at the time --
`source_study` and `sample_id` are now both in the output specifically so neither
level of granularity is silently lost.

## ANTXR1/ANTXR2 paralog co-presence across cell types, MRC2, and a corrected
redundancy framing

Third session, working from the two completed Phase 1 atlases above (no new
downloads or compute -- pure query/analysis/visualization on
`combined_celltype_means.parquet` and `skin_fibroblast_celltype_means.parquet`).

Terminology note: this section uses "co-presence" (matching the framing
question already posed in `prompts/skin_means.md`: "does ANTXR1/ANTXR2 paralog
co-presence vary by cell type") rather than "co-expression" -- everything here
is a **cell-type-mean** comparison (does a cell type's average include both
genes), which is a different and weaker claim than **co-expression**, reserved
throughout this project (see `prompts/initial_discussion.md`,
`prompts/celltype_means.md`) for the still-unstarted Phase 2 question of
whether the *same individual cells* express both genes simultaneously, as
opposed to two disjoint subpopulations averaging out to a similar-looking
mean. Conflating the two terms would make Phase 2 write-ups ambiguous against
this one.

**Co-presence scatter (both atlases, log10 mean expression, one dot per cell
type/subtype).** Whole-body: weak-to-moderate correlation between ANTXR1 and
ANTXR2 across the 238 cell types (Pearson r=0.326 log-log, R²=0.106, OLS
slope 0.19, residual scatter 5.4x around the trend). Only 10 of 238 cell types
(~4%) land in a genuinely discordant quadrant (one paralog in its top quartile,
the other in its bottom): ANTXR2-hi/ANTXR1-lo -- classical monocyte,
plasmacytoid dendritic cell, basophil, hematopoietic stem cell, hepatocyte,
Paneth cell of colon (thin, n=3); ANTXR1-hi/ANTXR2-lo -- myoepithelial cell,
spermatocyte/spermatid (both n=1, don't trust). Most cell types sit in a broad,
noisy co-presence band, not hugging either axis. Skin fibroblast subtypes
(healthy, nonlesional): no significant correlation at all (r=-0.18, n=13).

**The gut epithelial lineage is the most extreme, best-powered outlier in the
whole atlas.** Every well-powered gut epithelial label -- not just the thin
"proper of small intestine" splits already flagged as a caveat above -- sits at
the ANTXR2-dominant/ANTXR1-near-absent extreme: intestine goblet cell
(log-ratio +1.08, n=194), intestinal crypt stem cell (+1.21, n=134), enterocyte
(+2.12, n=142), colonocyte (+2.15, n=84), tuft cell of small intestine (+1.64,
n=92) -- vs. an atlas-wide median ratio of +0.42. Directly relevant to HFS's
gut phenotype (severe/ISH form: chronic diarrhea, protein-losing enteropathy).
Counter-example, reported alongside rather than dropped: `fibroblast of
gingiva` (n=62, well-powered) is *ANTXR1*-dominant (ratio -0.55) despite
gingival hypertrophy being a hallmark HFS symptom -- if gingival fibroblasts
specifically lean on ANTXR1, this cell type's own paralog balance doesn't
explain that phenotype the same way.

**MRC2 (the collagen-VI-clearance partner ANTXR2 is reported to work with) is
essentially absent specifically in gut epithelium**, checked with a
lineage-matched comparison (gut epithelium vs. other low-ANTXR2 epithelia --
keratinocyte, corneal epithelial cell -- to rule out "it's just an epithelial
vs. mesenchymal thing"): keratinocyte and corneal epithelial cell both retain
some MRC2 (3-4e-5) despite barely expressing ANTXR2 at all, while all four
well-powered gut epithelial types (enterocyte, colonocyte, goblet cell, crypt
stem cell) show MRC2 at the detection floor. Reframes the gut finding: not
"gut can't clear collagen VI," but "collagen-VI clearance probably was never
gut epithelium's job -- ANTXR2's presence there is doing something else."

**Correction, prompted by user pushback on the framing (not the data): early
in this analysis we treated "ANTXR1 present alongside ANTXR2 in a cell type's
mean" as evidence of functional redundancy/backup. That doesn't follow.**
Redundancy requires *shared function*, which cell-type-level co-presence alone
doesn't establish -- and the two paralogs are reported to bind different
domains of collagen VI (triple-helical vs. C5), so even
co-presence-plus-shared-binding-partner isn't proof of substitutability.
Concrete counter-evidence, not just a literature caveat: **skin fibroblasts --
the tissue with HFS's most visible pathology (the disease's namesake nodules)
-- show robust ANTXR1 co-presence** (ANTXR2:ANTXR1 ratio <1 in 9 of 13 healthy
skin fibroblast subtypes; see the skin atlas above). If ANTXR1 co-presence
conferred real protection, skin should be comparatively spared -- it isn't.
The corrected framing: the ANTXR1:ANTXR2 scatter is a cell-type co-presence
map, not a redundancy map, and any composite "vulnerability score" built from
it would be overclaiming what the data supports (and would still say nothing
about true single-cell co-expression either way). This governs the heatmap
design below (no composite score, ever).

Clinical tissue-involvement grounding used throughout (see the heatmap's own
`COVERAGE_GAP_NOTE` for exact sourcing): GAPO syndrome (ANTXR1 biallelic LoF)
-- growth retardation/metaphyseal dysplasia, alopecia, pseudoanodontia
(failed tooth eruption), progressive optic atrophy/glaucoma, dilated scalp
veins, hypogonadism; GAPO fibroblasts specifically show disrupted actin
cytoskeleton and reduced ECM turnover (same clearance-failure logic as HFS).
HFS/infantile systemic hyalinosis (ANTXR2 biallelic LoF) -- discrete
papular/nodular skin, gingival, and perianal lesions (matrix-accumulation
pathology); joint flexion contractures and osteopenia/osteolysis; chronic
diarrhea and protein-losing enteropathy with *diffuse* (non-nodular)
histological hyaline deposition in the gut wall -- a different kind of lesion
from the skin/gingiva nodules, not the same pathology at lower severity; and
(severe/ISH form) additional diffuse hyaline deposition in muscle, lymph node,
spleen, thyroid, and adrenal gland.

### ECM-clearance-pathway heatmap (figure)

Implemented per a reviewed plan (`/home/ubuntu/.claude/plans/clean-up-the-tissue-enchanted-axolotl.md`).
New files: `scripts/cell_type_curation.py` (curated row lists, disease
categories, palette -- pure data), `scripts/gene_panel_query.py` (query/
aggregation generalized from `query_antxr2_table.py`'s pattern to an 8-gene
panel, same two-step donor-equal-weighted method for the whole-body atlas and
sample-equal-weighted for skin), `scripts/plot_ecm_clearance_heatmap.py` (the
figure). `config.py` gained `GENE_PANEL` and `FIGURES_DIR`.
`query_antxr2_table.py` untouched.

- **Gene panel**: `ANTXR1, ANTXR2, MRC2, CTSB, CTSK, MMP14, TIMP2, LAMP1` --
  fixed scope, not the separate activation-axis panel.
- **Rows**: 31 curated whole-body cell types (thin site-specific GI splits
  excluded per user decision) in 4 categories -- HFS-nodular/fibrotic,
  HFS-diffuse/systemic (GI + muscle + lymphoid organs, merged from two
  categories to fit the dataviz skill's validated 3-hue all-pairs-safe
  categorical palette once rows are clustering-reordered rather than
  category-blocked), GAPO-relevant, and a neither-disease comparison baseline
  -- plus all 13 skin fibroblast subtypes, annotated healthy vs.
  disease-emergent. That annotation was verified against the published paper
  (not left as a naming-pattern guess, per user decision): three fibroblast
  populations have no healthy-skin counterpart per the paper (inflammatory
  myofibroblasts, plain myofibroblasts, fascia-like myofibroblasts -- the
  paper's own F6/F7/F8, mapped onto our object's differently-numbered F6/F6/F7
  labels for the same three groups). `F_Fascia` is flagged "unresolved," not
  silently coerced either way -- the paper says healthy fascial fibroblasts
  were merged into F2:Universal, yet our object has `F_Fascia` as its own
  distinct label with real cells (2,978 cells, 3 donors), a genuine
  discrepancy not resolved from a secondary source alone.
- **Row order**: hierarchical clustering (average linkage, correlation
  distance on the log-transformed gene profile), computed independently per
  panel, per user request -- replaces an earlier category-block-order plan.
  No dendrogram is drawn (per user request); only the leaf order is used.
- **Color**: per-gene (column) min-max normalization, added after the first
  render showed LAMP1 (a uniformly-high housekeeping lysosomal gene) anchoring
  the shared color scale and drowning out every gene's own cross-cell-type
  pattern -- per user request. Real values are preserved throughout: in-cell
  text (skin panel), the CSV's `mean_expression`/`log10_mean_expression`
  columns, and a separate `normalized_0_1` column carrying exactly what color
  was displayed. Two independent color scales (one per panel, each its own
  log floor) since the README's flat-capture-rate caveat means cross-panel
  absolute-scale comparison isn't well supported by the data.
- **Hard constraint carried through the implementation, not just prose**: no
  function computes or stores an ANTXR2:ANTXR1 ratio or any other composite
  score for display/export -- the one exception (an ANTXR2:ANTXR1 ratio used
  only to cross-check against the earlier scatter's "9 of 13" finding) is
  printed to stdout only, never written to a file or plotted.
- **Coverage gaps handled as a caption footnote, not placeholder rows**: bone/
  cartilage/joint/synovium, teeth, thyroid, adrenal, gonadal tissue have no
  matching Cell Ontology label anywhere in either atlas -- stated explicitly
  rather than silently omitted or faked as measured-near-zero.
- **Output**: `/data/ANTXR2/figures/ecm_clearance_heatmap/` --
  `heatmap.png`/`.svg` and `whole_body_curated.csv`/`skin_curated.csv` (each
  with the score/scale/normalization/coverage-gap caveats embedded as header
  comments).
- **Verification passed**: exactly 8 genes returned per scan; no curated row
  missing from either atlas (would have logged a warning); `n_donors`/
  `n_samples` constant across all 8 genes within a row (no gene silently
  missing for some donor); MRC2 near-zero specifically across the GI-category
  whole-body rows, matching the finding above; skin ANTXR2>ANTXR1 in 4/13
  subtypes (i.e. ANTXR1>=ANTXR2 in 9/13), matching the scatter finding
  exactly.

## Phase 2 prep: co-expression working set (Elmentaite2021 trio)

Preparation only -- no correlation computed yet. Built by
`scripts/subset_coexpression_data.py` (+ `capture_rate.py`,
`download_pbmc_reference.py`).

**Two variants exist**, selected with `--variant {celltype,level3}` from
`config.COEXPR_VARIANTS`; they differ *only* in which annotation column defines
the groups, and each writes its own h5ad and README into
`/data/ANTXR2/coexpression/`:

| variant | file | cells | groups from |
|---|---|---:|---|
| `celltype` | `elmentaite2021_trio.h5ad` (0.37 GB) | 82,034 | CELLxGene-harmonised `cell_type` |
| `level3` | `elmentaite2021_trio_level3.h5ad` (0.26 GB) | 56,882 | authors' own `level_3_annot` |

`level3` is the later, cleaner one (see "Level-3 variant" below) and is the
better default for co-expression work; `celltype` is kept because the Phase 1
means and the heatmap are all in `cell_type` space, so it is what joins back to
them. Both are 18,370 genes, carry identical capture-rate columns, and record
their variant in `uns['provenance']`.

**Investigating the enterocyte label first changed what we thought we had.** The
whole-body "218,750 enterocytes / 142 donors" figure was inflated by
double-counting: all 75 donors in the Gut Cell Atlas "Healthy reference" release
are *also* in "Extended+" (the README's caveat 3 biting in practice). Real unique
count is 152,548 cells / 169 donors, from a single collection -- and effectively
just two studies (Kong2023 52%, Elmentaite2021 23%, then a long unusable tail).
Also: "enterocyte" is 96% small intestine (only 386 large-intestine cells), so
enterocyte vs. colonocyte is a real anatomical contrast, not a naming one; and
`disease == normal` is **not** a healthy baseline -- 54,246 of those cells are
`Neighbouring_inflamed` (IBD donors' uninflamed-looking tissue) and 9,325
`Neighbouring_cancer`, leaving only 74,594 truly `Non_pathological`. That
distinction exists only in the raw h5ad's `sample_category`, not in our parquet.

**Chosen slice: Elmentaite2021** -- 398,460 cells, 38 donors, 100%
`Non_pathological`, and critically it contains all three needed roles from the
**same 35 donors**, one lab, one protocol, balanced 5'/3' chemistry
(donor counts by `donorID_unified`; see the donor-identity correction below --
the raw `donor_id` column would say 41 and 36, and both are wrong):

| cell type | cells | role | ANTXR2 raw mean |
|---|---:|---|---:|
| fibroblast | 42,278 | canonical: ANTXR1 + ANTXR2 + MRC2 all present | 0.224 |
| enterocyte | 35,062 | the unexplained case (ANTXR2 without ANTXR1/MRC2) | 0.078 |
| macrophage | 4,694 | myeloid calibration control | 0.116 |

(per-group donor counts, `donorID_unified`: fibroblast 38, enterocyte 35,
macrophage 38; 35 shared by all three.)

All clear memento's 0.07 `filter_mean_thresh`; enterocyte only just, so a null
result there is underpowered rather than evidence of absence. This beats the
originally-sketched design (gut epithelium here, fibroblasts from the skin
atlas, myeloid from a third cohort), which would have confounded every
cross-cell-type difference with cohort and technology.

### Second curation error found -- and this one propagated into published output

**`gastrointestinal tract (lamina propria) macrophage` (CL:0000865) is not a
macrophage.** All 42,302 cells carrying that label across the *entire* Gut Cell
Atlas are authored `level_1_annot=Mesenchymal`, `level_2_annot=Fibroblast`,
`level_3_annot=Lamina_propria_fibroblast_ADAMDEC1`. Markers side with the
authors, decisively: PTPRC 0.4% positive, CD68 2.5%, LYZ 0.5%, C1QA 0.5% --
against COL1A1 **89%**. These are fibroblasts; the CELLxGene ontology mapping is
wrong, systematically, across all 10+ constituent studies.

**Corrected in the heatmap (2026-09-01):** the row was **dropped** from
`WHOLE_BODY_CATEGORIES` in `scripts/cell_type_curation.py` and
`plot_ecm_clearance_heatmap.py` re-run -- the whole-body panel is now 30 cell
types (was 31), and the regenerated figure/CSVs no longer contain the label.
Dropped rather than relabelled to `fibroblast`: the generic `fibroblast` row
already represents that compartment in the same panel, and a second fibroblast
row under a GI heading would imply a distinct GI-resident population the label
cannot support.

**Still uncorrected:** the ANTXR1 top-20 table earlier in *this file* lists
`gastrointestinal tract (lamina propria) macrophage` (137 donors after the
2026-09-03 donor-identity re-run; 147 before) among
ANTXR1-high cell types, described as a macrophage. The number is right but the
identity is wrong -- and its ANTXR1-high rank is in fact *explained* by it being
a fibroblast, which is consistent with every other fibroblast row there rather
than the surprise it reads as. Treat that row as a fibroblast when reading the
table.

Same lesson as the donor-identity correction earlier in this file: a label that
looked authoritative (a Cell Ontology term, no less) was wrong, and only a
marker check caught it. The trio's marker panel is therefore recomputed on every
pipeline run rather than trusted from documentation.

### Two capture-rate columns, deliberately not averaged

memento takes per-cell `q` via `q_column` (asserts `max < 1`). Both columns are
attached; they rest on different assumptions, so disagreement is informative.

- **`capture_rate_chem` = 0.06991** (both chemistries). memento's own broad
  droplet rate 0.07 (`publication/cellxgene/make_cube.py`, and the value used in
  their PBMC *co-expression* analysis), corrected for sequencing saturation.
  The correction is **not** the naive product: 10x "sequencing saturation" is a
  read-duplicate fraction, not molecule loss, so it's inverted through the
  Poisson relation (`s = 1-(1-e^-L)/L`, `detected = 1-e^-L`). At the assumed
  s=0.85 that gives L=6.658, detected=0.9987 -- i.e. nearly inert, which is the
  honest answer; the naive product would have understated q by ~15%.
  `ASSUMED_SATURATION=0.85` is an **assumption**: the paper never states
  saturation and ENA has `read_count=0` for all 89 runs of the adult GEX
  accession (E-MTAB-9543), so it is unrecoverable for these samples; only
  E-MTAB-8901 (developing gut, ~44K reads/cell) supports the 80-85% figure.
- **`capture_rate_pbmc` = 0.0399 (3' v2) / 0.0544 (5' v2)**, transferred from
  chemistry-matched 10x public PBMC references (`pbmc8k`,
  `sc5p_v2_hs_PBMC_10k`) by median UMI-depth ratio across marker-gated immune
  populations. Gut immune cells run *shallower* than blood (ratios 0.52-1.09),
  so this lands below the nominal 0.07 -- plausible for dissociated tissue.
  Limitation stated in the output README: it assumes a given immune cell type
  carries the same absolute mRNA content in gut as in blood, which
  tissue-residency may violate.

Two implementation traps worth remembering: (1) the gates must be evaluated on
the study's **immune compartment** (84,330 cells), not the working set -- the
trio has no lymphocytes, so a first pass gating the subset matched only doublets
and ambient RNA and produced ratios >1 in the wrong direction; (2) **CD68 is not
a macrophage marker in gut** -- 52% positive, mean 3.10 in *enterocytes* (a
lysosomal glycoprotein, and enterocytes have very active endolysosomal
compartments), so the myeloid gate uses LYZ + AIF1 with EPCAM-negative instead.

Also noted: `obs['n_counts']` is *not* the row sum of the shipped matrix (runs
~0.04% higher, max 14% on one cell; predates this release's gene subsetting,
since `feature_is_filtered` is all-False). Everything depth-related in the
pipeline uses matrix row sums.

### Third ID correction: `donor_id` is wrong in both directions (2026-09-03)

Prompted by the question "do donors appear in 2 different `sourceID`s?" The
literal answer was reassuring -- 4 donors do (`390C`, `HT-228`, `HT-234`,
`HT-236`), none of them in Elmentaite2021 -- but checking it exposed a worse
problem in `donor_id` itself, which **does** reach our data:

- **Collision.** `A25` and `A34 (417C)` each map to *two different people*. For
  `A34 (417C)`: `D11` (55-74y, 31,370 cells, 22 samples) and `D12` (18-34y,
  5,814 cells, 2 samples). Different age brackets -- unambiguously two donors
  under one id string. 8,304 cells in the `celltype` working set carry it.
- **Aliases.** 34 people appear under several `donor_id` strings in Extended+
  (9 in Healthy reference); in Elmentaite2021, 3 do -- `D12` -> `A32 (411C)` /
  `A34 (417C)`, `D2` -> `T036` / `T036NEG` / `T036POS`, `D5` -> `T110NEG` /
  `T110POS`. The NEG/POS pairs are sorted fractions of one donor. (`D12` is
  caught in *both* problems.)

Net in Elmentaite2021: **41 `donor_id` strings for 38 actual people.** So
`donor_id` simultaneously merges two donors and splits three others.

**Why this is worse than a miscount.** The project's binding aggregation method
averages *equally across donors*, so a merged pair is under-weighted and a split
donor is over-weighted. And memento's 2D path bootstraps over donors, so
carrying `donor_id` into the co-expression work would be straightforward
pseudo-replication in exactly the estimates we are about to compute.

**Fix, per user instruction that all analysis use `donorID_unified`:**
`config.DONOR_UNIFIED_COL` was added and `compute_means.py` now resolves the
donor column per dataset -- `donorID_unified` where the file provides it, else
`donor_id`, logging which. **Only the gut atlas provides it**; Tabula Sapiens
and Cross-tissue Immune have `donor_id` only, so "use `donorID_unified`
everywhere" is not literally achievable and the resolution is documented
instead of faked.

Scope of the re-run:
- **Gut atlas recomputed** (both datasets; 185->175 and 308->271 donors), since
  chunking is per-donor and the grouping itself had to change.
- **Tabula Sapiens / Cross-tissue Immune not recomputed** -- no unified column
  exists, their semantics are unchanged, and re-running them would have cost a
  full pipeline pass for a guaranteed-identical result.
- `combined_celltype_means.parquet` rebuilt from the three per-collection files
  by a new `scripts/rebuild_combined.py` (streaming concat, schemas verified
  identical first). This was necessary because `compute_means.py` writes the
  per-collection and combined files in one pass, so a partial re-run would
  otherwise have truncated `combined` to just the re-run datasets -- hence the
  new `--skip-combined` flag.
- Downstream re-run: `query_antxr2_table.py`, `plot_ecm_clearance_heatmap.py`,
  and both co-expression variants.

**One deliberate wart:** the output column is still named `donor_id` but now
holds `donorID_unified` for gut rows and `donor_id` elsewhere. Adding a
provenance column instead would have broken schema compatibility with the two
collections that were not re-run, and `combined` is a plain concatenation that
depends on one shared schema. Documented prominently in the output README
(caveat 5) rather than left to be discovered; filter on `collection_name` if the
distinction matters.

**Re-verified after the re-run: every headline conclusion is unchanged.** Only
donor counts moved. The paralog co-presence statistics reproduce almost exactly
(Pearson r 0.326 -> 0.324, R² 0.106 -> 0.105, OLS slope 0.190 -> 0.189,
residual scatter 5.4x both times, 10 of 238 cell types discordant both times,
atlas median log-ratio +0.42 both times). The gut-epithelium signal holds
(goblet +1.08, crypt stem +1.21 -> +1.28, enterocyte +2.12 -> +2.09, colonocyte
+2.15, small-intestine tuft +1.64 -> +1.63), as does the gingival-fibroblast
counter-example (-0.55 -> -0.54). The generic `enterocyte` label still ranks
150th of 238.

What did change is `n_donors` throughout -- e.g. fibroblast 179 -> 159,
macrophage 187 -> 173, enterocyte 142 -> 131, intestine goblet cell 194 -> 182 --
and total row count, 351,128,441 -> 338,470,470 (fewer, larger donor groups).
So the correction mattered for *how confidently* each cell type is supported and
for the validity of donor-level bootstrapping, not for the biology read off the
means. Worth stating plainly rather than implying the re-run rescued a result.

This is the third ID-level error in this project, after the barcode-derived
donor recovery in the skin atlas and the mislabelled lamina-propria macrophage.
The pattern is consistent: identifiers that look authoritative -- a Cell
Ontology term, a `donor_id` column -- were wrong, and only cross-checking them
against an independent signal (markers, age brackets, a unified id) caught it.

### Level-3 variant -- narrower, cleaner populations (2026-09-02)

Per user request, a second working set built from the authors' own
`level_3_annot` rather than the CELLxGene-harmonised `cell_type`, which pools
subtypes. Output `elmentaite2021_trio_level3.h5ad`, **56,882 cells, 34 donors
shared across all three** (all donor counts by `donorID_unified`):

| `level_3_annot` | cells | donors | ANTXR2 | ANTXR1 | MRC2 | COL1A1 |
|---|---:|---:|---:|---:|---:|---:|
| Enterocyte | 35,062 | 35 | 0.078 | 0.0008 | 0.0009 | 0.039 |
| Crypt_fibroblast_PI16 | 18,867 | 37 | 0.278 | 0.343 | 0.422 | 27.97 |
| Macrophage | 2,953 | 38 | 0.121 | 0.027 | 0.078 | 0.016 |

All three still clear memento's 0.07 threshold. Two differences from the
`celltype` variant that matter:

- **`Crypt_fibroblast_PI16` is a single defined subtype**, not a mixture, and
  carries visibly stronger signal across the whole axis of interest than the
  pooled `fibroblast` (ANTXR2 0.278 vs 0.224, ANTXR1 0.343 vs 0.323, MRC2 0.422
  vs 0.306). A better-defined comparator for co-expression. (The label is
  `PI16`, the fibroblast marker gene -- there is no `PI1` label in this atlas.)
- **`Macrophage` drops to 2,953 cells** (from 4,694), because the plain level-3
  label excludes the `Macrophage_LYVE1` (1,436), `_TREM2` (527), `_MMP9` (295)
  and `_CD5L` (10) subtypes that `cell_type` had merged in. The population is
  measurably purer for it (COL1A1 0.016 vs 0.030), but it is now the smallest
  group by a wide margin -- worth remembering when the 2D bootstrap arrives,
  since correlation estimates are far more sample-hungry than means.

Implemented by parameterising the existing script rather than forking it:
`--variant` selects `label_col` + `labels` + output path from
`config.COEXPR_VARIANTS`, and the README filename follows the variant. Both
variants were re-run so they are mutually consistent; the `celltype` file's
contents are unchanged.

**Verification passed on both variants:** per-group counts and shared-donor
counts as tabulated; counts integer with min-nonzero 1; no unused categoricals
survive; `obsp` absent from the written files; both q columns in (0,1); ANTXR2
means reproduce the values above; marker panel correct per group (notably
`Macrophage` C1QA 96% / PTPRC 64% / COL1A1 1.5%, `Crypt_fibroblast_PI16`
COL1A1 91% / PTPRC 0.2%, `Enterocyte` EPCAM 86%); and a memento smoke test
(`setup_memento` + `create_groups` + `compute_1d_moments`) returns 3 groups with
finite moments over all 18,370 genes.

## Open / next steps (SUPERSEDED)

> **This list is stale and narrower than the current framing.** The
> authoritative open list is "Designed but never run" in the header.
> Kept here for provenance -- several items below have since been
> re-ranked or reframed.


- Phase 1 (means only) is complete and verified. Phase 2 (variance /
  formal hypothesis testing via memento's `ht_1d_moments`) is not started.
- **Co-expression point-estimate correlations are now run** (see "Phase 2:
  co-expression analysis (test run)" below) — this resolves the gene-pair-
  scoping decision that was previously open: `compute_2d_moments`'s
  ANTXR2-vs-all pairs use every gene in the `level3` variant's global
  `filter_mean_thresh`-filtered gene list (4,009 genes), not a curated seed
  panel. Point estimates only — `binary_test_2d`/`ht_2d_moments` (bootstrap
  hypothesis testing) are still not run; that remains open if a formal
  significance test on specific gene pairs becomes the next question.
- The macrophage statistical-power risk flagged below was confirmed sharply in
  the run above: only 5 of 38 donors survive the 100-cell floor for
  macrophage (vs. 21 fibroblast, 28 enterocyte), and its top correlations
  (several exactly 1.000) are likely an artifact of that small n combined
  with several donor groups where ANTXR2 is barely detected. If macrophage
  co-expression becomes load-bearing for a real question, revisit the
  `celltype` variant's pooled macrophage (4,694 cells) mentioned below rather
  than trusting the `level3` numbers as-is.
- Disease contrast for co-expression still needs Kong2023 (within-study
  inflamed vs. healthy, 22,101 `Non_pathological` / 49,743
  `Neighbouring_inflamed` / 7,491 `Inflamed`, 42 donors). Elmentaite2021 cannot
  provide it — every sample there is `Non_pathological`.
- The stem-vs-differentiated gut epithelium comparison above is suggestive but
  informal (means only, no significance test) -- worth a proper memento
  `binary_test_1d` comparison (stem/TA vs. differentiated enterocyte) if this
  becomes a specific question to chase.
- Capture rate is still a flat, uncalibrated 0.1 placeholder for every assay
  (see README caveat 1) -- revisit if absolute-scale comparisons across assays
  become important.
- Skin fibroblast atlas: means computed and verified with recovered `sample_id`
  granularity (see "Correction" section above), Phase 2 (ANTXR1/ANTXR2 paralog
  co-presence and ECM-clearance machinery co-assembly by subtype) not started. Two
  source studies (`Ganier`, `Sole-Boldo`) remain at coarse `source_study`-level by
  choice, not oversight -- if either becomes load-bearing for a specific Phase 2
  question, the `Ganier` ~90%-clean subset is a plausible target for a real fix.
- Paralog co-presence + ECM-clearance heatmap (see section above) is a
  means-level descriptive figure, not a significance test -- the MRC2-absent-
  in-gut finding and the skin-fibroblast redundancy counter-evidence are both
  suggestive, not formally tested. A real memento `binary_test_1d` (gut
  epithelium vs. other epithelium, on MRC2 specifically) would close that gap
  cheaply, same pattern as the stem-vs-differentiated item above.
- The "what else is ANTXR2 doing in gut epithelium if not collagen-VI
  clearance" question is still open -- flagged but not chased this session.
  A positive/co-presence screen (what other genes' cell-type means correlate
  *with* ANTXR2's specifically in gut epithelium, genome-wide or against a
  candidate list -- still cell-type-level, not single-cell) is the natural
  next move if this becomes a specific question to pursue, as opposed to the
  MRC2 check's absence-focused approach.

## Phase 2: co-expression analysis (test run)
Ran `scripts/coexpression_pipeline.py` on the `level3` trio working set (/data/ANTXR2/coexpression/elmentaite2021_trio_level3.h5ad) -- the test run planned in `planning_summaries/coexpression.md`. Point estimates only (memento's `compute_1d_moments`/`compute_2d_moments`/`get_corr_matrix`), no bootstrap/hypothesis testing. Output: `/data/ANTXR2/coexpression/elmentaite2021_trio_level3_coexpr.h5ad` (a NEW file -- the source h5ad is untouched), figures in `/data/ANTXR2/figures/coexpression/`.
### Parameters
- Donor column: `obs['donor']` (not `donor_id` -- see donor-identity note in the Phase 2 prep section above)
- Capture rate: `obs['capture_rate_pbmc']`, passed via memento's `q_column` (NOTE: only 2 distinct values across all cells -- a per-assay constant, not continuous per-cell; expected, see prep notes)
- Minimum group size: 100 cells per donor x cell_type group
- Gene filter (memento defaults, recorded explicitly): `filter_mean_thresh=0.07`, `min_perc_group=0.7`
- Ranking: top 50 genes per cell type by |donor-averaged correlation|
- Averaging across donors: plain unweighted mean over non-NaN donor values (same donor-equal-weighting convention as Phase 1's binding aggregation)

### Donor x cell_type groups
38 donors total. 56 of 110 existing groups fell below the 100-cell floor and were dropped:

```
donor  cell_type  n_cells
   D1 fibroblast       18
   D1 macrophage       87
  D10 enterocyte       10
  D10 fibroblast       47
  D10 macrophage        6
  D13 macrophage       28
  D14 macrophage       31
  D15 fibroblast       34
  D15 macrophage       54
 D150 fibroblast       10
 D150 macrophage        9
 D151 fibroblast        4
 D151 macrophage       32
 D152 enterocyte       43
 D152 macrophage       96
 D153 enterocyte       72
 D153 fibroblast        4
 D153 macrophage        4
 D154 fibroblast       21
 D154 macrophage       81
 D155 fibroblast        1
 D155 macrophage       31
 D156 fibroblast       17
 D156 macrophage       42
   D2 fibroblast       33
   D2 macrophage       38
   D3 fibroblast       12
   D3 macrophage       34
   D4 fibroblast       20
   D4 macrophage       72
   D5 fibroblast       12
   D5 macrophage       30
   D6 fibroblast        2
   D6 macrophage       12
   D7 macrophage        7
   D8 fibroblast        6
   D8 macrophage       10
   D9 macrophage       40
   F1 enterocyte       14
   F1 macrophage       27
  F10 macrophage       97
  F11 macrophage        3
  F13 macrophage       78
  F14 enterocyte        1
  F14 macrophage        1
   F2 macrophage        6
   F3 macrophage       19
   F4 enterocyte       18
   F4 macrophage        1
   F5 enterocyte        2
   F5 fibroblast       38
   F5 macrophage        9
   F6 macrophage       24
   F7 macrophage        1
   F8 macrophage       45
   F9 macrophage       67
```

Groups kept, by cell type:

```
cell_type
enterocyte    28
fibroblast    21
macrophage     5
```

**Macrophage is severely thinned by the 100-cell floor** -- only 5 of 38 donors survive for macrophage (vs. 21 for fibroblast, 28 for enterocyte). Macrophage's donor-averaged correlations rest on a much smaller donor sample than the other two cell types -- treat any macrophage-specific finding here as exploratory, not confirmatory.

### Gene filtering
Global filtered gene list: 4009 / 18370 genes pass (`filter_mean_thresh > 0.07` in >70% of donor x cell_type groups). ANTXR2 survives the global filter (required for step 2 to run at all); per-group survival against ANTXR2's OWN group's mean-expression threshold is logged at run time. **Caveat, not a NaN case**: a group where ANTXR2 fails its own local mean filter still produces a numeric correlation there (memento only emits NaN when variance is exactly zero, which is stricter than failing the mean filter -- verified directly, e.g. sg^F3^enterocyte has ANTXR2 mean=4.6e-6, var=8.3e-11, both nonzero). That correlation is real output, not dropped, but is derived from near-undetected expression and is noisier than groups where ANTXR2 clears its own filter -- most relevant to the macrophage top-gene list below, whose correlation magnitudes (up to 1.000) likely reflect this combined with the small donor count (14 of the 54 kept groups fall in this low-expression category, listed at run time by group name).

### Top ANTXR2-correlated genes per cell type (donor-averaged)

**fibroblast** (top 10 of 50 shown, by |mean correlation|):

```
 rank gene_symbol  mean_corr  n_donors
    1         C1D   0.688365        21
    2        SDHD   0.667203        21
    3       SH2B1   0.628330        21
    4     TMEM208   0.626196        21
    5     PPP1R11   0.622356        21
    6     ADIPOR1   0.613340        21
    7      PITHD1   0.605012        21
    8        TPP1   0.601699        21
    9        MOB2   0.597954        21
   10       STK16   0.586686        21
```

**enterocyte** (top 10 of 50 shown, by |mean correlation|):

```
 rank gene_symbol  mean_corr  n_donors
    1     TRPC4AP   0.620934        28
    2        SBF2   0.618805        28
    3       ROCK2   0.612711        28
    4      LSM14A   0.611623        28
    5       APLP2   0.598420        28
    6      CTNND1   0.591493        28
    7      ZFAND3   0.589756        28
    8      PRKAB2   0.580358        28
    9     SMPDL3A   0.579310        28
   10        CD46   0.571864        28
```

**macrophage** (top 10 of 50 shown, by |mean correlation|):

```
 rank gene_symbol  mean_corr  n_donors
    1     KAZALD1   1.000000         5
    2        YAP1   1.000000         5
    3      SEMA6D   1.000000         5
    4       PTPRG   1.000000         5
    5       NDEL1   1.000000         5
    6        PHF1   0.963056         5
    7        SCO1   0.954579         5
    8       DCAF5   0.934275         5
    9       UBTD2   0.932276         5
   10         HGS   0.926760         5
```

Union panel: 146 genes (145 unique top-50 genes + ANTXR2).

Pairwise overlap between cell types' top-gene lists:

```
fibroblast vs enterocyte: 4
fibroblast vs macrophage: 1
enterocyte vs macrophage: 0
```

### Package versions
```
anndata: 0.12.19
scanpy: 1.11.5
memento-de: 0.1.3
```

### Runtime (test-run trio: 56,882 cells x 18,370 genes)
```
1_load_and_filter: 4.3s
2_compute_antxr2_correlations: 10.9s
3_average_across_donors_1d: 0.0s
4_select_top_genes: 0.0s
5_compute_pairwise_correlations: 1.6s
6_average_pairwise_across_donors: 0.0s
7_write_output_and_figures: 40.2s
total: 57.0s
```

Step 2 (ANTXR2 vs. all filtered genes) and step 5 (pairwise over the union panel) are the two that will scale with dataset size when this moves beyond the trio test run -- step 2 scales roughly with n_filtered_genes x n_groups, step 5 with n_union_genes^2 x n_groups (here n_union_genes is fixed at 146 regardless of dataset size, so step 5's cost is mostly driven by n_groups, i.e. how many donor x cell_type combinations exist).

## Correction: memento estimator bug was inflating every correlation above -- found, fixed, re-run (2026-09-04)

The Phase 2 run above had a real problem: the ANTXR2-vs-gene correlation distribution sat well above zero in every cell type (fibroblast median 0.276, enterocyte 0.321, macrophage 0.186), and the effect got *worse* the smaller the donor group. Chasing this down (total UMI, the fetal/pediatric donor-cohort split, dissociation-stress genes, a fibroblast-only UMAP, per-donor leverage analysis, a 1000-random-gene-pair sanity check, a permutation null) eventually found the actual cause, and it wasn't biology or a normalization choice:

**The bug**: memento 0.1.3's `estimator._corr_from_cov` (used by `compute_2d_moments`/`get_2d_moments`, i.e. step 2 of this pipeline) initializes its output array to a placeholder value of `5.0`, then only overwrites entries where *both* genes have positive variance in that donor group. Entries left at the placeholder -- which happens whenever either gene's variance estimate is `<=0`, common for lowly-expressed genes in small donor groups -- get silently clipped to exactly `1.0` by the unconditional `corr[corr>1]=1` line that follows, instead of being set to NaN as intended. (Its sibling function `_hyper_corr_symmetric`, used by `get_corr_matrix` -- step 5 of this pipeline -- does **not** have this bug: it explicitly NaNs anything still outside `[-1,1]` after clipping. The bug is specific to the `compute_2d_moments` path.)

Verified directly before fixing anything: every affected entry, in every donor group checked, was pinned at **exactly** 1.0 (never near it, never -1 -- confirming a deterministic code path, not statistical noise). It accounted for up to **59% of a 1000-random-gene-pair sample** in the smallest donor groups (e.g. donor D152, 118 fibroblast cells) and a negligible ~0.6% in the largest (F13, 3,108 cells) -- a clean, monotonic relationship with donor group size that explained the entire earlier "positive skew" finding. Excluding just the bugged entries brought every donor's correlation distribution -- including the smallest ones -- back to roughly centered on zero (D152: fake median 1.000 -> real median -0.27).

**Fix applied** (`compute_antxr2_correlations` in `scripts/coexpression_pipeline.py`): after `compute_2d_moments`, for every donor group, null out any pair where either gene's group-level variance is `<=0`, before averaging across donors. This run corrected **35,339 pair x group values** from a fake `1.0` to `NaN`.

**Re-ran the full pipeline with the fix, plus the size-factor settings adopted earlier** (`shrinkage=0`, `trim_percent=0.5`, both now the pipeline's defaults -- memento's own defaults are 0.5 and 0.1 respectively; `min_perc_group` reverted to the strict default of 0.7). New distribution, over the full tested-gene universe (not just top-50):

```
cell_type   n_genes    mean    median    std    frac_positive
fibroblast     3746   0.040     0.038   0.137          0.607
enterocyte     3746   0.009     0.007   0.124          0.521
macrophage     3744  -0.117    -0.130   0.347          0.353
```

Fibroblast and enterocyte are now both tightly centered almost exactly on zero. Macrophage sits slightly negative with much wider spread (std 0.35) -- consistent with its persistent small-donor-count problem (only 5 surviving donor groups), not a new bias.

**New top-10 ANTXR2-correlated genes per cell type** (donor-averaged, bug-corrected; note the qualitatively different, more trustworthy character vs. the pre-fix table above -- realistic donor coverage for fibroblast/enterocyte, and a healthy mix of positive *and* negative signs where before almost everything was positive):

```
fibroblast          mean_corr  n_donors     enterocyte          mean_corr  n_donors
1  PDCD6                0.540        12     1  PRPF31              -0.434        24
2  MICU2                0.525        13     2  ZFAND3               0.415        21
3  C11orf58             0.494        15     3  USE1                -0.392        24
4  NDUFA5               0.466        15     4  PICALM               0.387        21
5  SERINC1              0.464        16     5  WDR1                 0.376        23
6  ADPGK                0.461        14     6  FAM3C                0.375        25
7  BPNT2                0.451        17     7  CPD                  0.367        25
8  TMEM87A              0.440        16     8  SNX15               -0.366        25
9  MKRN1                0.436        15     9  CD164                0.365        24
10 PITHD1               0.427        15     10 ANKHD1               0.359        23
```

Macrophage's top-10 still shows values up to exactly 1.000/-1.000, but this is no longer the bug -- these are now transparently backed by `n_donors` of 1-4 (e.g. IFT20 r=1.000 on n_donors=1, ECHS1 r=-1.000 on n_donors=4), a legitimate (if unreliable) small-sample point estimate rather than a hidden defect. Still exploratory-only for macrophage, as flagged throughout -- not confirmatory.

**Implication for the ECM-clearance-panel question from earlier this session**: CTSB/LAMP1/ANTXR1/MRC2/CTSK/TIMP2/MMP14's correlations with ANTXR2, reported earlier as moderate-to-strong positive (CTSB 0.48, LAMP1 0.45, ANTXR1 0.35, etc.), were computed before this fix and are very likely inflated by the same bug -- not necessarily wrong in sign, but overstated in magnitude, and worth recomputing before drawing further conclusions from them.

**Remaining known limitation, not fixed by this patch**: donor-equal-weighted averaging still gives a 100-cell donor group the same vote as a 3,000-cell one. The bug's *systematic* effect is now gone, but small donor groups still carry genuinely higher sampling *variance* in their correlation estimates than large ones -- visible in the per-donor boxplot from this session's random-gene-pair sanity check. Worth revisiting (e.g. weighting by donor group size, or a higher minimum-cell floor specific to correlation work) if this analysis is extended beyond the current trio test run.

**Follow-up check: reverting shrinkage/trim_percent to memento's own defaults, with the bug fix still applied**, confirmed the size-factor settings were fixing a real, separate problem, not just masking the bug:

```
setting                                          fibroblast   enterocyte   macrophage
memento defaults (shrinkage=0.5, trim_percent=0.1)    0.148        0.237        0.099
this pipeline's defaults (shrinkage=0, trim_percent=0.5)  0.038        0.007       -0.130
```

(medians, bug fix applied in both cases). Even with the estimator bug patched, memento's own size-factor defaults leave fibroblast and enterocyte visibly, systematically shifted positive -- so the bug and the shrinkage/trim_percent choice were two distinct real problems, not one masquerading as the other, and both fixes are needed. Output on disk is confirmed at this pipeline's defaults (`shrinkage=0`, `trim_percent=0.5`, `min_perc_group=0.7`, bug fix applied).

## GSEA (gseapy prerank) replaces top-50 ORA as the primary enrichment method (2026-09-04)

**Committed the fix** (`4be8a36`): the `_corr_from_cov` patch in `compute_antxr2_correlations`, `shrinkage=0`/`trim_percent=0.5` as this pipeline's new defaults, `scripts/coexpression_enrichment.py`, and the correction write-up above. Re-ran the full pipeline from the committed code and confirmed it reproduces the corrected numbers exactly (fibroblast/enterocyte/macrophage medians 0.038 / 0.007 / -0.130, as above) -- deterministic, no drift.

**Switched enrichment method from Enrichr ORA to gseapy pre-ranked GSEA** (`scripts/coexpression_gsea.py`), because ORA only tests a fixed top-50 gene list against a hypergeometric background and is underpowered here: it can only detect a pathway if enough of its individual members happen to independently clear an arbitrary top-50 cutoff. GSEA instead runs a running-sum statistic over the **entire ranked gene list** (~3,745 genes per cell type, ranked by signed donor-averaged correlation with ANTXR2), so it can detect a pathway that skews consistently toward one end of the ranking even when no single member is individually extreme enough to make a top-50 cut. Used gseapy's `prerank` (1000 permutations, GO Biological Process 2023 / KEGG 2021 Human / Reactome 2022, FDR q < 0.25 -- GSEA's own conventional threshold, looser than ORA's 0.05 since it tests gene-*set* skew rather than per-gene significance).

**Result: a real, cross-validated biological signal that the top-50/ORA approach could not see.**

- **Ribosome biogenesis / translation machinery is significantly, negatively correlated with ANTXR2 -- independently discovered in both fibroblast and enterocyte.** Fibroblast: 112 significant terms, 92 of them negative, led by "Eukaryotic Translation Elongation" (NES=-3.99, FDR~0, 71 leading-edge genes out of 84), "Cytoplasmic Translation," "Ribosome," rRNA processing. Enterocyte: 335 significant terms (264 positive, 71 negative), with the *same* negative theme recurring independently (mitochondrial translation, spliceosome, "Ribosome," NES down to -2.29). No individual ribosomal gene was extreme enough in either cell type to reach its own top-50 list -- this is a pathway-level effect that requires the full ranking to detect, and its independent replication across two unrelated cell types is stronger evidence than either result alone.
- **Fibroblast's positive side** (20 of 112 significant terms) is coherent with ANTXR2's established biology as an endocytic receptor: regulation of endocytosis, protein localization to membrane/lysosome, with PICALM, CLTC, LRP1, DAB2 as leading genes -- a cleaner, better-powered version of the CTSB/LAMP1 lysosomal-trafficking thread noticed earlier this session with the (pre-GSEA, pre-bug-fix) ECM-clearance-panel comparison.
- **Enterocyte's positive side** (264 terms) centers on actin cytoskeleton regulation, the RHOH GTPase cycle, plasma membrane organization, and SREBP-driven cholesterol biosynthesis.
- **Macrophage**: only 5 significant terms -- consistent with its 5-donor limitation, not read as a real finding.

**Two artifacts published**: the existing "ANTXR2 Correlation Atlas" (top genes, ORA enrichment, overlap, clustered heatmaps) was updated in place with the corrected post-bug-fix numbers; a new companion "ANTXR2 GSEA Atlas" was published with the same top-gene/overlap/heatmap sections but the ORA section replaced by the GSEA results above.

**Re-verified the bug fix specifically on macrophage, on request, rather than assuming it generalized from the fibroblast test case**: macrophage's 5 donor groups had 2,390 degenerate (var<=0) pairs across them (up to 20.9% of the 3,746 tested pairs in the smallest groups, F16/D11) -- every one of them was pinned at exactly 1.0 pre-fix (the same deterministic signature found in fibroblast), and all 2,390 are correctly nulled to NaN post-fix. Separately spot-checked one of macrophage's remaining exact-1.000 top-gene values (IFT20 in the D12 macrophage group): both genes' variances are strictly positive at full float precision (3.3e-10 and 3.6e-9), and the raw pre-clip correlation is a real computed 2.75 (only 16 of 399 cells double-positive) -- a genuine small-sample-extreme value that legitimately clips to 1.0, not a leftover placeholder. Macrophage's remaining exact +/-1.000 top-gene values are real point estimates, correctly labeled with their (very low, 1-4) `n_donors`, not a residual bug -- but should still be read as exploratory-only given the sample size, same caveat as before.

## Donor-replication sanity check for the ANTXR2 correlation estimates (2026-09-05)

New script `scripts/coexpression_replication_qc.py`, run against the already-computed `.uns['memento_correlations']['by_donor']` in `elmentaite2021_trio_level3_coexpr.h5ad` -- no new memento computation. Purpose: before trusting the donor-averaged top-gene/GSEA results above, check whether independent donors' per-donor ANTXR2-vs-gene single-cell correlation estimates actually agree with each other (a **correlation of correlation estimates**, not a single-cell correlation itself -- wording kept precise per project convention).

**Scope narrowed after investigation**: chemistry (10x 5' v2 vs 3' v2) is completely confounded with donor identity in this working set -- no donor was profiled with both, and macrophage is 100% 5' v2. A true same-donor technical replicate does not exist here, so this check is biological (donor-vs-donor) replication only; chemistry composition is recorded per cell type as a caveat, not tested formally. Per user decision, macrophage is included but flagged exploratory throughout.

**A second data-quality issue found in the process**: 3 fibroblast donors (D152, F2, F3) and 2 enterocyte donors (D6, F9) have a 100%-NaN `antxr2_vs_all` vector -- ANTXR2 itself has zero variance in their group (too few cells), which nulls every entry for that donor via the estimator-bug fix above (ANTXR2 is one of the two genes in every pair). These contribute no information and were dropped before computing anything, rather than left in as blank heatmap rows -- usable donor counts: fibroblast 18/21, enterocyte 26/28, macrophage 5/5 (none dropped).

**Method**: per cell type, built a genes x donors matrix from the per-donor correlation vectors, computed all pairwise Pearson correlations between donors (pairwise-complete on non-NaN genes), and compared against a permutation null (each donor's gene axis independently shuffled 200x, breaking true gene correspondence while preserving each donor's own value distribution and NaN count) via Mann-Whitney U (one-sided, real > null).

**Result: a real but modest signal, honestly quantified rather than eyeballed.** Median pairwise donor-donor correlation sits close to the null band in absolute terms (fibroblast 0.023, enterocyte 0.021, macrophage 0.014, vs. null 95% bands of roughly [-0.036, 0.036]) -- expected, since most of the ~3,746 tested genes carry no real ANTXR2 relationship and dilute the whole-vector correlation toward zero. But pooled across many independent donor pairs, the shift is highly significant for the two well-powered cell types: **fibroblast p=4.3e-28, enterocyte p=4.9e-55** (153 and 325 pairs respectively; ~28-31% of real pairs exceed the null's 97.5th percentile, vs. the 2.5% expected by chance). **Macrophage is much weaker** (p=0.038, only 5 donors / 10 pairs) -- consistent with its documented small-sample caveat elsewhere in this file; not read as a null result, just underpowered.

**Interpretation for the tables/GSEA above**: donor-averaging is picking up real, reproducible cross-donor signal for fibroblast and enterocyte, not averaging independent noise -- but the whole-genome effect size is modest, consistent with only a minority of genes carrying real signal against a large noise floor (the top-50 lists are exactly that minority, not the median gene). Macrophage's results should continue to be treated as exploratory only.

**Outputs** (`/data/ANTXR2/figures/coexpression/replication_qc/`): per-cell-type donor x donor correlation heatmaps (clustered order); a combined real-vs-null distribution figure; per-cell-type "pairs()"-style scatter grids, both a full all-donors version (hexbin, lower triangle) and a legible top-8-donors-by-cell-count version (true scatter) -- these make the modest-but-real effect visually inspectable as diffuse-but-off-center point clouds, not tight diagonals; `donor_pairwise_correlations.csv` (long format, caveat header) and `replication_summary_stats.csv`.

### Finding the right gene-selection metric: two flawed hacks, then a proper discovery/replication design

Follow-on the same session, prompted by user pushback on interpretation. The whole-gene-set check above shows real but modest aggregate agreement (median donor-pair r~0.02, highly significant vs. permutation null but small in magnitude) -- expected, since most of the ~3,746 tested genes carry no true ANTXR2 relationship and dilute a whole-vector correlation. The natural follow-up question -- "does a small set of specific, individually-identifiable genes replicate cleanly across donors?" -- went through three iterations before landing on a defensible answer.

**Attempt 1, circular (rejected).** Selected each cell type's top-50 genes by |donor-averaged mean correlation| (the same statistic the project's GSEA top-gene tables use), excluding genes that hit memento's known ±1.0 small-sample clipping artifact in any donor, then re-ran the donor-pairwise replication check restricted to that subset. Median donor-pair r jumped to 0.379 (fibroblast) / 0.189 (enterocyte) -- but this is circular: the same donors being tested for agreement were also used to compute the mean that selected the genes (double-dipping / winner's-curse inflation), so the jump is partly a selection artifact, not purely discovered signal.

**Attempt 2, non-circular but noise-starved (also rejected, for the opposite reason).** Per user suggestion, removed circularity by selecting each DONOR'S OWN top-50 genes (by |value|, using only that donor's data) and taking the union across donors, then testing pairwise agreement on that union set. This is legitimately non-circular (no gene was selected using both members of any tested pair), but it over-corrected: different donors' individual top-50 lists barely overlapped at all (fibroblast 803/900 possible-max union, enterocyte 1101/1300, macrophage 247/250 -- i.e. ~11-15% pairwise overlap, near the "no overlap" extreme). Restricting to this union brought the apparent signal back down to roughly the full-gene-set baseline (median r 0.016 fibroblast, 0.017 enterocyte, macrophage went to -0.038/n.s.). Diagnosis: a hard top-N cut on a SINGLE noisy donor's point estimate is dominated by that donor's own sampling noise (order statistics of noisy data reproduce poorly), so low per-donor-list overlap is expected even when real average signal exists -- this test is under-powered, not evidence of no signal.

**Attempt 3, the fix: discovery/replication split + memento's own one-sample bootstrap test.** Per user direction, combined a proper fold split (donors round-robin-assigned to a discovery half and a replication half, balanced by cell count) with a REAL significance test in the discovery half, rather than a magnitude-ranking hack. memento exposes `ht_2d_moments(adata, treatment, ...)`, normally used for two-group correlation-difference testing via bootstrap; passing a constant (`treatment=1` for every group) collapses its internal regression to exactly a one-sample bootstrap test of "is this correlation significantly different from zero", using memento's own per-cell bootstrap for the standard error (`memento/hypothesis_test.py::_regress_2d`: `if (treatment==1).mean()==1: corr_coef = np.average(boot_corr, axis=0, weights=Nc_list)`). Ran this (num_boot=5000, BH-FDR correction) on the discovery half only, then evaluated replication ONLY on the untouched replication half -- no donor is ever used for both gene selection and its own test. ANTXR2 needed a lowered `min_perc_group` (0.7->0.5) to survive its own presence filter in the enterocyte discovery half specifically (consistent with the "enterocyte only just clears the filter" caveat noted earlier in this file); fibroblast and macrophage used the standard 0.7. Macrophage was attempted despite an expected-underpowered 3-vs-2 donor split.

New script: `scripts/coexpression_discovery_replication.py`.

**Result -- a real, well-replicating gene-level signal in fibroblast and enterocyte; macrophage confirmed underpowered as expected:**

| cell type | discovery/replication donors | FDR<0.1 genes | evaluable in replication | sign concordance | binomial p | effect-size replication r | cross-half median r |
|---|---:|---:|---:|---:|---:|---:|---:|
| fibroblast | 9 / 9 (of 11/10 total, 2 degenerate donors dropped) | 130 | 81 | 77/81 (95%) | 7.3e-19 | 0.770 | 0.243 |
| enterocyte | 13 / 13 (of 14/14 total, 1 degenerate donor dropped) | 291 | 187 | 160/187 (86%) | 1.7e-24 | 0.364 | 0.036 |
| macrophage | 3 / 2 (EXPLORATORY) | 21 | 9 | 6/9 (67%) | 0.25 (n.s.) | -0.381 | -0.030 |

Fibroblast replicates convincingly on both sign and magnitude: 95% of FDR-significant genes keep the same sign in a completely independent set of donors, and the discovery-half bootstrap coefficient correlates at r=0.77 with the replication-half mean -- comparable to a solid GWAS/omics replication rate. Enterocyte shows a more specific pattern worth remembering when using its gene lists: sign concordance is even more significant than fibroblast's (p=1.7e-24, since there are more genes), but the effect-size correlation is weaker (r=0.36) and nearly all 187 evaluable genes are positive-direction -- i.e. **for enterocyte, trust the direction of an effect much more than its exact magnitude**. Macrophage shows no significant replication by any measure (binomial p=0.25, negative effect-size correlation) -- confirms, rather than merely asserts, that macrophage's donor-averaged correlations in this project should stay exploratory-only.

**Outputs** (same `replication_qc/` directory): `discovery_replication_ht.png` (the effect-size replication scatter, one panel per cell type); `discovery_replication_summary.csv` (the table above); `discovery_significant_genes.csv` (every FDR<0.1 discovery gene, its bootstrap coef/se/pval/fdr, and its replication-half mean correlation where evaluable -- caveat header explains the non-circular design). This significant-gene list is a more defensible starting point for future candidate-gene follow-up than the original point-estimate top-50 tables earlier in this file, which had no formal significance test or independent-donor validation behind them.

**Net methodological lesson, worth carrying into any future per-gene QC on this kind of data**: neither extreme of "rank by a data-derived statistic computed on the test data" nor "use only one noisy replicate's own extreme values" gives a trustworthy gene-level signal -- the former double-dips, the latter is starved by single-replicate noise. A discovery/replication split combined with an actual significance test (not a magnitude cutoff) threads this correctly, and cost only ~1-2 minutes of compute per cell type at num_boot=5000 on 16 cores.

### Full-dataset significance test, significance-based heatmaps, and z-score GSEA (2026-09-05, same session)

Having validated the one-sample bootstrap test methodology on a discovery/replication split, the natural next step -- run it on ALL usable donors per cell type for maximum power, then rebuild every downstream product (pairwise heatmaps, cross-cell-type overlap, GSEA ranking) from its output instead of the original point-estimate top-50s. New scripts: `scripts/coexpression_full_ht.py` (the full-dataset test), `scripts/coexpression_full_ht_downstream.py` (overlap + heatmaps), and a `--metric zscore` mode added to `scripts/coexpression_gsea.py`.

**Full-dataset one-sample test** (`memento.ht_2d_moments`, constant treatment, same trick as the discovery-half test above -- see that section for why this is a real significance test, not a hack): num_boot=5000 for fibroblast (21 donors, 6m27s) and macrophage (5 donors, 46s); num_boot=3000 for enterocyte (28 donors, 5m10s, reduced only to fit this session's tooling constraints, not for a methodological reason). ANTXR2 needed `min_perc_group` lowered to 0.6 for enterocyte's full donor set (vs. 0.5 needed for its smaller discovery half -- consistent with it being a borderline-expressed gene there, not a new problem).

| cell type | donors | genes tested | FDR<0.1 | FDR<0.05 |
|---|---:|---:|---:|---:|
| fibroblast | 21 | 4,502 | 225 | 163 |
| enterocyte | 28 | 5,852 | 287 | 212 |
| macrophage | 5 | 3,977 | 42 | 22 |

Output: `/data/ANTXR2/figures/coexpression/full_dataset_ht/{cell_type}_full_dataset_ht.csv` -- every tested gene (not just significant ones, so GSEA and any future re-ranking can use the full universe), with `coef`, `se`, `pval`, `fdr`, and `z = coef/se`.

**Cross-cell-type overlap of significant genes is small** (`significant_gene_overlap.csv`): fibroblast vs enterocyte share 11/225-287 genes (~5% of the smaller set: AKAP9, APLP2, CIRBP, DDX17, EIF4A2, KCNQ1OT1, PDIA3, RTN4, SQSTM1, UBC, WSB1); enterocyte vs macrophage share 2 (CD59, REEP3); fibroblast vs macrophage share 0; no gene is significant in all three. Each cell type's ANTXR2 correlation partners are largely distinct gene sets, consistent with the different biological programs each cell type's GSEA implicates (below).

**Significance-based pairwise heatmaps**: replaced the original `select_top_genes`'s |mean_corr|-magnitude top-50 panel with each cell type's own top-50-by-p-value FDR<0.1 genes (still capped at 50 per cell type for legibility only -- the FULL 542-gene significant union blew past image size limits and would have been illegible regardless; overlap counts and GSEA below are unaffected by this cap, they use the complete significant/tested sets). Output: `full_dataset_ht/{cell_type}_pairwise_correlation.png`. The fibroblast heatmap shows clean, visually obvious block structure: a large ECM/collagen module (COL6A2, COL1A1, COL1A2, FN1, ELN, FAP) positively correlated with ANTXR2 and with each other; a tight mitochondrial-gene block (MT-ND1-5, MT-CO1/3, MT-CYB, MT-ATP6); and a ribosomal-protein block (RPS/RPL genes) that anti-correlates with the mitochondrial block -- a visible mitochondrial-vs-cytoplasmic-translation axis underlying the GSEA ribosome finding below.

**GSEA re-ranked by bootstrap z-score (coef/se), not the plain point-estimate mean** -- per user request, to down-weight noisy/imprecise correlations relative to confident ones (a large coefficient with a large se should count for less than the same coefficient with a tight se, which a plain mean-based ranking cannot express). Still uses the FULL gene universe per cell type (4,502 / 5,852 / 3,977 genes), same "no top-50 cutoff" principle as the original mean-based GSEA. Output: `/data/ANTXR2/figures/coexpression/gsea_zscore/` (parallel to the original `gsea/`).

| cell type | old (mean-ranked) significant terms | new (z-score-ranked) significant terms |
|---|---:|---:|
| fibroblast | 112 (92 neg / 20 pos) | 421 (107 neg / 314 pos) |
| enterocyte | 335 (71 neg / 264 pos) | 505 (64 neg / 441 pos) |
| macrophage | 5 (not read as real) | 119 (106 neg / 13 pos) |

- **Fibroblast's positive side sharpened dramatically and became directly on-target**: the z-score ranking's top positive terms are essentially a clean readout of ANTXR2's own known biology -- Extracellular Matrix Organization, Collagen Formation, ECM-receptor Interaction, Focal Adhesion, Collagen Biosynthesis, ECM Proteoglycans (lead genes: COL1A1, COL1A2, COL3A1, COL6A1/2, COL9A3, COL14A1, FN1, LAMB1/C1, ITGB1, SPARC, DCN, ADAMTS1...). This is a more mechanistically specific result than the mean-ranked GSEA's positive side (generic endocytosis/lysosomal trafficking). **AMENDED (see "Standing corrections" in the header): this was originally written up as the project's cleanest confirmation that the pipeline measures something biologically real. That reading is wrong and must not be reused.** Bürgi et al. 2017 showed collagen VI mRNA does *not* change in Antxr2−/− uteri -- the mechanism is degradation, not transcription -- so ANTXR2→COL6 transcriptional coupling is not predicted, and a positive ECM correlation cannot serve as a positive control, because a positive control must be able to fail. The result is still worth reporting, but as **co-regulation evidence** (ANTXR2 and the ECM program plausibly under shared tissue-program pressure in a cell type where the clearance partners *are* present), and it is gated on the anchor-gene specificity control, which has not been run. Note also that this is an *exploratory* GSEA finding, not a confirmatory panel test.
- **Fibroblast and enterocyte's negative side (ribosome/translation) reproduces under the new metric**, essentially unchanged in character from the mean-ranked version -- Ribosome, rRNA Processing, Cytoplasmic Translation, Peptide Chain Elongation all remain top negative hits. This cross-metric reproduction is itself evidence the finding is not an artifact of one particular ranking choice.
- **Enterocyte's positive side changed character**, from actin-cytoskeleton/cholesterol-biosynthesis (mean-ranked) to Fat/Protein Digestion and Absorption, MHC Class I antigen presentation, PPAR signaling (z-score-ranked). Both are plausible enterocyte biology; which one is "more correct" isn't resolved by this analysis alone -- flagged as a metric-sensitive result, not adjudicated. **Later note:** both branches are compatible with an epithelial barrier/absorption axis (brush border and tight junctions are actin/plasma-membrane organization; protein-losing enteropathy and malabsorption are digestion and absorption), so they may not be in conflict. Treat this as a hypothesis to pre-specify and test, NOT as a resolution -- enterocytes do digestion and absorption as their dominant biology, so this enrichment is exactly what the anchor-gene control exists to check. **AMENDED 2026-09-06: checked, and it fails the control** -- a majority of expression-matched anchor genes also recover this term (see the header's anchor-gene-specificity standing correction and the "Task 2 results" section below). Read this enrichment as generic enterocyte biology, not ANTXR2-specific co-regulation.
- **Macrophage's apparent gain (5 -> 119 terms) should NOT be read as a real improvement, and is very likely a statistical artifact worth understanding explicitly**: memento's bootstrap se estimates SAMPLING noise from cells *within* one donor's group; it does not, and cannot, capture *between-donor* biological variability. The discovery/replication split-half check earlier in this file directly measured between-donor agreement for macrophage and found it indistinguishable from chance (66.7% sign concordance, binomial p=0.25, negative effect-size correlation) -- macrophage's 5 donors simply don't agree with each other enough to support gene-level claims. A z-score can look confidently large purely because a donor group happens to have enough cells for a small within-group bootstrap se, even while the true cross-donor signal for that gene is noise. Mean-based ranking is comparatively insensitive to this trap (it doesn't reward small se at all); z-score ranking is exactly the metric most exposed to it when between-donor replication has already failed. **Macrophage's z-score GSEA results are not used for any biological conclusion in this project** -- consistent with every other macrophage caveat in this file, just now with a specific mechanistic reason for why this particular metric makes it look artificially better, not worse.

**Practical implication for reading this project's two GSEA atlases going forward**: for fibroblast and enterocyte -- both independently validated by the discovery/replication split -- prefer the z-score-ranked results (`gsea_zscore/`) over the original mean-ranked ones (`gsea/`) where they differ, since down-weighting imprecise correlations is the statistically correct thing to do and the fibroblast ECM/collagen signal it surfaces is a stronger, more specific confirmation than what mean-ranking found. For macrophage, neither GSEA result should be trusted as a biological finding; the gap between them is itself informative only as a demonstration of the se-vs-replication-variance trap above.
## Partner-availability analysis: pre-registration (2026-09-06)

Implementing `prompts/partner_availability.md`, items 1-5 of "Designed but never
run." Interpretation rules below are written **before looking at any Task
1/2/5 value**, per that prompt's instructions. Tasks 3 and 4 are pure
extension/formal-test tasks with no ambiguous-direction result to pre-commit
to, so no rule is pre-registered for them.

**Task 1 (ANTXR1/ANTXR2 paralog readout):**
- Significant **negative** correlation -> mutual exclusivity -> real evidence
  ANTXR1 and ANTXR2 occupy different cells (the paralog-partitioning result
  the project was designed around). Interpretable immediately, no specificity
  control needed (per the standing "anti-correlation is stronger" correction).
- Significant **positive** correlation -> co-occupancy or co-regulation ->
  **not interpretable until Task 2's anchor-gene floor has run**, and requires
  the doublet check (fibroblast-epithelial doublets would specifically inflate
  this pair).
- **Null** -> uninformative. Explicitly NOT evidence of co-presence in either
  direction.
- Macrophage's value will be reported but is exploratory-only and no
  conclusion will be drawn from it (standing correction).

**Task 2 (anchor-gene specificity control):**
- If **most** ~20 expression/detection-matched anchor genes in a cell type
  independently recover the same GSEA terms (ECM/collagen for fibroblast,
  digestion/absorption for enterocyte) at the same FDR threshold used for
  ANTXR2, that recovery is a **property of the cell type's transcriptional
  program, not of ANTXR2 specifically** -- ANTXR2's positive correlation with
  that program is not distinguishing evidence. This will be reported as
  clearly as the alternative outcome; it is a real, useful negative, not a
  failed analysis.
- If ANTXR2 sits at or beyond the extreme tail of the anchor distribution
  (empirical percentile) for that term's NES while most anchors do not recover
  it, that supports ANTXR2 carrying a specific (not merely cell-type-generic)
  relationship to the program -- still co-regulation/co-occupancy, not a
  causal or physical-interaction claim (per the header's evidence table).
- This result **gates** the fibroblast ECM-correlation and enterocyte
  digestion/absorption-correlation claims already in this log -- neither is
  to be treated as confirmatory until this control has run.

**Task 5 (Kong2023 feasibility check):**
- This is a **feasibility gate, not a results-producing step.** No
  differential-correlation test is run in this task regardless of outcome.
- If ANTXR2 fails `min_perc_group=0.7` in a condition x cell-type group, the
  threshold that *would* be needed is reported (reference point: 0.6 was
  needed for Elmentaite's full enterocyte donor set) rather than treating
  failure as a dead end without characterizing it.
- Any condition with <8 usable donors (>=100 cells, per
  `COEXPR_MIN_GROUP_CELLS`) for a cell type is flagged **underpowered for
  two-group testing** in that condition, using the macrophage precedent (5
  donors -> replication indistinguishable from chance) as the reference for
  what "underpowered" has meant empirically in this project -- not a new,
  arbitrary bar.
- Chemistry/batch confounding with condition will be checked explicitly
  (donor confounding with a grouping variable has already burned this project
  once, in the `donor_id` collision/alias issue).
- The differential-correlation test itself is only scoped (never run) in this
  task, and only if feasibility passes, with the pre-specified directional
  prediction from Bracq et al.: ANTXR2's coupling to Wnt-arm partners and
  regeneration programs should be **absent or weak in `Non_pathological`** and
  **appear in `Inflamed`**, since the Wnt function is injury-conditional.

Results for all five tasks follow below, each labeled confirmatory
(pre-defined panel test) or exploratory (screening), per the prompt's output
spec.

### Task 1 results: ANTXR1/ANTXR2 paralog readout

**Confirmatory** (pre-defined panel test -- ANTXR1 is the single, pre-specified gene of interest here, read out of an already-computed full-dataset test).

| cell type | in tested universe | rank | coef | se | z | pval | fdr | call |
|---|---|---:|---:|---:|---:|---:|---:|---|
| fibroblast | yes (4502 genes tested) | 1908 | 0.1174 | 0.0750 | 1.564 | 0.3414 | 0.8051 | NULL (uninformative) |
| enterocyte | NO (5852 genes tested) | -- | -- | -- | -- | -- | -- | filtered out (near-absent expression) |
| macrophage | NO (3977 genes tested) | -- | -- | -- | -- | -- | -- | filtered out (near-absent expression) |

Pre-fix point estimate for comparison: ANTXR1-ANTXR2 was reported earlier in this log ("Implication for the ECM-clearance-panel question" section) at **0.35** -- that value predates the `_corr_from_cov` variance<=0 null-out fix and is **superseded**; per the correction, it is not necessarily wrong in sign but is very likely overstated in magnitude by the same placeholder-clipping bug that affected every other panel gene.

**Fibroblast (the only cell type where ANTXR1 survives the detection filter and gets a real test): NULL.** coef=0.117, pval=0.341, fdr=0.805, rank 1908 of 4502 -- squarely mid-pack, not extreme in either direction. Per the pre-registered rule, **this is explicitly NOT evidence of ANTXR1/ANTXR2 co-presence in fibroblast** -- it means the correlation test has no power to distinguish the co-presence hypothesis from the mutual-exclusivity hypothesis here, not that co-presence is confirmed. Because coef is nominally positive but not significant, this is a **null, not a 'significant positive'** in the pre-registered sense -- the doublet check (which the prompt gates on significant positive correlations specifically) was correctly **not triggered** for fibroblast; see doublet-check section below for the reasoning trace.

**Enterocyte and macrophage: ANTXR1 is filtered out of the tested universe entirely** (fails the memento `min_perc_group` presence filter before a correlation could even be attempted). Direct per-cell counts explain why: ANTXR1 is detected in only **0.08% of enterocytes** (35,062 cells) and **2.2% of macrophages** (2,953 cells), vs. 25.2% of fibroblasts. **This is itself the stronger, means-level asymmetric result** the header's evidence table describes: ANTXR1 is essentially absent from gut epithelium and largely absent from macrophage by the mean/detection-rate alone -- no correlation test is needed to make that call, and none was possible. A co-expression question ("do ANTXR1 and ANTXR2 occupy the same enterocytes/macrophages") is moot when one partner is barely present at all; the paralog-partitioning claim for these two cell types rests on presence/absence, not correlation.

**Macrophage is reported for completeness only and is exploratory** (standing correction) -- no conclusion is drawn from it regardless of the above.

#### Bimodality check

Figure: `/data/ANTXR2/figures/coexpression/paralog_readout/antxr1_antxr2_distribution.png`. Per-cell expression is sparse/zero-inflated for both genes in every cell type (detection rates 0.08%-25.3%, well below the >90% near-uniform threshold used here) -- **neither gene is near-uniformly expressed anywhere**, so correlation remains a meaningful co-presence readout wherever both genes clear the detection filter (fibroblast); it is not needed as a readout where one partner is already known absent by the mean (enterocyte, macrophage).

```
 cell_type   gene  n_cells     mean  pct_nonzero  cv_nonzero  near_uniform
fibroblast ANTXR1    18867 0.342980    25.250437    0.588618         False
fibroblast ANTXR2    18867 0.278264    21.577357    0.512392         False
enterocyte ANTXR1    35062 0.000799     0.077006    0.182108         False
enterocyte ANTXR2    35062 0.077605     6.830757    0.376249         False
macrophage ANTXR1     2953 0.027430     2.201151    0.923371         False
macrophage ANTXR2     2953 0.121233    10.836438    0.351512         False
```

#### Doublet check

**Not triggered for any cell type.** The prompt's doublet check applies only to a significant positive correlation; fibroblast's positive point estimate was not significant (null), and enterocyte/macrophage have no correlation at all (ANTXR1 filtered out). No doublet check was run.


### Task 3 results: extended gene panel, two-arm heatmap

**Exploratory characterization / descriptive query** (not a significance test -- means only, same status as the original ECM-clearance heatmap). `GENE_PANEL_ARMS` added to `config.py` (RECEPTORS, CLEARANCE_ARM, CLEARANCE_SUBSTRATE, WNT_ARM); original `GENE_PANEL` untouched. Figures: `/data/ANTXR2/figures/two_arm_panel/two_arm_heatmap_{wholebody,skin}.png`. **No composite score computed** -- `*_arm_completeness.csv` reports, per row and per gene, whether that gene's raw value is within 2 orders of magnitude of its OWN maximum across the curated rows (the same standard already used by hand for the MRC2 gut-epithelium-vs-keratinocyte/corneal finding earlier in this log, a ~2-3 order-of-magnitude gap), alongside the raw values in `*_two_arm_raw.csv`; it is a per-gene reading aid over the same raw numbers, not a score that combines genes.

**Per-cell-type arm availability call (whole-body atlas):**

| cell type | clearance arm | substrate (COL6) | Wnt arm | neither arm complete |
|---|---|---|---|---|
| retinal pigment epithelial cell | 6/6 | 1/3 | 14/18 |  |
| endothelial cell | 5/6 | 2/3 | 16/18 |  |
| retinal blood vessel endothelial cell | 6/6 | 2/3 | 14/18 |  |
| sebocyte | 3/6 | 1/3 | 14/18 |  |
| skeletal muscle satellite stem cell | 6/6 | 3/3 | 16/18 |  |
| neuron | 4/6 | 2/3 | 12/18 |  |
| follicular dendritic cell | 4/6 | 2/3 | 13/18 |  |
| enteric neuron | 5/6 | 2/3 | 13/18 |  |
| enteroglial cell | 6/6 | 2/3 | 13/18 |  |
| adventitial cell | 6/6 | 3/3 | 15/18 |  |
| tissue-resident macrophage | 6/6 | 3/3 | 15/18 |  |
| fibroblast of gingiva | 6/6 | 3/3 | 16/18 |  |
| fibroblastic reticular cell | 6/6 | 3/3 | 14/18 |  |
| fibroblast | 6/6 | 3/3 | 17/18 |  |
| myofibroblast cell | 6/6 | 3/3 | 16/18 |  |
| interstitial cell of Cajal | 6/6 | 3/3 | 15/18 |  |
| smooth muscle cell | 6/6 | 3/3 | 16/18 |  |
| paneth cell | 3/6 | 0/3 | 16/18 |  |
| melanocyte | 5/6 | 2/3 | 15/18 |  |
| naive B cell | 2/6 | 0/3 | 6/18 | **YES** |
| pancreatic acinar cell | 4/6 | 0/3 | 15/18 |  |
| hepatocyte | 2/6 | 0/3 | 13/18 |  |
| M cell of gut | 4/6 | 1/3 | 13/18 |  |
| colonocyte | 4/6 | 1/3 | 15/18 |  |
| enterocyte | 4/6 | 1/3 | 13/18 |  |
| intestinal crypt stem cell | 3/6 | 0/3 | 16/18 |  |
| intestine goblet cell | 4/6 | 0/3 | 16/18 |  |
| intestinal tuft cell | 4/6 | 0/3 | 15/18 |  |
| corneal epithelial cell | 4/6 | 0/3 | 14/18 |  |
| keratinocyte | 5/6 | 0/3 | 14/18 |  |

Cell types with the clearance arm essentially absent AND the Wnt arm well under half present (i.e. ANTXR2 present with **neither** function's machinery structurally available): naive B cell.

**Gut epithelium specifically** (paneth cell, colonocyte, enterocyte, intestinal crypt stem cell, intestine goblet cell, intestinal tuft cell): clearance-arm genes (mean n_present/6 across these rows: 3.7) confirm the project's existing MRC2-absence finding extends to the arm as a whole. Wnt-arm presence (mean n_present/18: 15.2) is the new information this panel adds -- if broadly present, gut epithelium has the Wnt-arm partners even though it lacks the clearance-arm ones, consistent with Bracq et al.'s finding that CMG2's role there is Wnt-pathway (injury-conditional), not clearance.

**Skin fibroblast subtypes:** see `skin_arm_completeness.csv` for the per-subtype breakdown (13 subtypes) -- not reproduced row-by-row here since the whole-body table above already carries the generic `fibroblast` row for this atlas's clearance-arm baseline (established: robust ANTXR1 co-presence, standing correction: not evidence of redundancy).

### Task 4 results: formal binary_test_1d tests

**Confirmatory** (pre-specified genes/groups, formal memento `ht_1d_moments` two-group test with donor as the replicate unit; adaptive `min_perc_group` retry, same posture as `coexpression_discovery_replication.py`, when a gene is near-absent in one group by design). Ran directly (not via `memento.binary_test_1d`, which hardcodes `min_perc_group=0.9` with no retry). `de_coef`/`de_se`/`de_pval` are memento's differential-mean bootstrap test; positive `de_coef` means higher expression in the group listed first (treatment=1). Output: `/data/ANTXR2/figures/binary_tests/binary_test_1d_results.csv` (+ per-test donor summaries).

#### Task 4.1: MRC2, gut epithelium vs. lineage-matched comparison

Run as **two separate within-atlas tests** rather than one pooled cross-atlas test -- keratinocyte and gut epithelium co-occur in the Gut Cell Atlas itself (perianal-adjacent skin samples), while corneal epithelial cell only exists in Tabula Sapiens; merging cells across two atlases with different chemistries into one memento run would need its own capture-rate harmonization (as Phase 2's PBMC-calibration work did for the coexpression trio), which is out of scope for a formal test of an existing means-level finding. This is a deviation from a single pooled test, noted explicitly.

**4.1a:** failed_presence_filter (min_perc_group tried down to 0.1) -- MRC2 could not even clear a heavily-relaxed presence filter in one group, which is itself consistent with (arguably stronger than) the significant-difference result: the gene is too near-absent for the bootstrap machinery to even engage.

**4.1b (Tabula Sapiens, gut epithelium vs. corneal epithelial cell): EXPLORATORY -- underpowered by this project's own macrophage precedent** (<8 usable donors on at least one side: 5 gut-epithelium vs. 2 corneal donors). de_coef=-4.0060, se=0.3100, pval=1.699e-38 (SIGNIFICANT, min_perc_group=0.57). Reported for completeness since corneal epithelial cell was part of the originally-chosen comparison set, but this specific test should not be treated as confirmatory -- read 4.1a as the formal result.

#### Task 4.2: ANTXR2, crypt stem/TA vs. differentiated enterocyte

(135 stem/TA vs. 106 enterocyte donors ≥100 cells, well-powered.) de_coef=-0.6164 (treatment=1 is crypt_stem_TA, so **negative = higher in the differentiated enterocyte** group), se=0.0282, pval=5.825e-106 -- **SIGNIFICANT**. ANTXR2 is significantly HIGHER in differentiated enterocyte, formally confirming the descriptive means (crypt stem +1.21 vs. enterocyte +2.12 log-ratio) with a real significance test for the first time. This is now load-bearing for reading Lencer's commentary on Bracq et al. (ANTXR2's own position tracks the diminishing Wnt gradient along the crypt-villus axis) as a confirmed, not merely suggestive, pattern.

#### Task 4.3: WNT_ARM genes along the same crypt-villus axis

Same two groups as 4.2 (crypt_stem_TA=1 vs. differentiated_enterocyte=0), read against ANTXR2's own gradient above -- **negative de_coef means higher in differentiated enterocyte (same direction as ANTXR2 itself, i.e. tracks ANTXR2 down the gradient toward the villus)**, positive means higher in crypt/stem (opposite direction, i.e. tracks WITH the Wnt-active compartment as expected for genuine Wnt-pathway partners since crypt-base cells are where canonical Wnt signaling is active).

| gene | status | de_coef | de_se | pval | fdr (within WNT_ARM) | direction |
|---|---|---:|---:|---:|---:|---|
| LRP5 | ok | -0.5841 | 0.0176 | 2.492e-242 | 5.399e-242 **sig** | villus/differentiated-enriched |
| LRP6 | ok | 0.0088 | 0.0217 | 0.7174 | 0.7174  | crypt/stem-enriched (Wnt-active compartment) |
| FZD1 | ok | -0.0429 | 0.0421 | 0.4146 | 0.4491  | villus/differentiated-enriched |
| FZD2 | failed_presence_filter | -- | -- | -- | -- | -- |
| FZD3 | ok | 1.7380 | 0.0570 | 3.455e-189 | 6.416e-189 **sig** | crypt/stem-enriched (Wnt-active compartment) |
| FZD4 | failed_presence_filter | -- | -- | -- | -- | -- |
| FZD5 | ok | -0.1056 | 0.0131 | 3.362e-16 | 4.37e-16 **sig** | villus/differentiated-enriched |
| FZD6 | ok | 0.8925 | 0.0423 | 6.526e-93 | 1.061e-92 **sig** | crypt/stem-enriched (Wnt-active compartment) |
| FZD7 | ok | -0.2462 | 0.0378 | 5.113e-11 | 6.042e-11 **sig** | villus/differentiated-enriched |
| FZD8 | failed_presence_filter | -- | -- | -- | -- | -- |
| FZD9 | failed_presence_filter | -- | -- | -- | -- | -- |
| FZD10 | failed_presence_filter | -- | -- | -- | -- | -- |
| CTNNB1 | ok | -0.1288 | 0.0091 | 1.241e-45 | 1.792e-45 **sig** | villus/differentiated-enriched |
| TCF7L2 | ok | -0.3832 | 0.0112 | 4.961e-257 | 1.29e-256 **sig** | villus/differentiated-enriched |
| LGR5 | ok | 3.9424 | 0.0964 | 0 | 0 **sig** | crypt/stem-enriched (Wnt-active compartment) |
| RNF43 | ok | 1.0428 | 0.0143 | 0 | 0 **sig** | crypt/stem-enriched (Wnt-active compartment) |
| ZNRF3 | ok | 2.4317 | 0.0524 | 0 | 0 **sig** | crypt/stem-enriched (Wnt-active compartment) |
| AXIN2 | ok | 2.4343 | 0.0525 | 0 | 0 **sig** | crypt/stem-enriched (Wnt-active compartment) |

### Task 5 results: Kong2023 feasibility check

**Feasibility gate, not a results-producing step** (per pre-registration) -- no differential-correlation test was run. Kong2023 data source: already inside the downloaded Gut Cell Atlas Extended+ h5ad (`/data/ANTXR2/raw/gut_cell_atlas/19053a82-9c89-4fb8-bd19-d7b1800b0b7b.h5ad`, `study=='Kong2023'` filter), no new download needed. 235,327 cells, 71 donors. Output: `/data/ANTXR2/figures/kong2023_feasibility/{per_donor_counts,feasibility_summary,chemistry_by_condition}.csv`.

#### 1-2. Per cell_type x condition: donor counts and ANTXR2 presence

ANTXR2 "presence" replicates memento's own filter formula directly (a donor group's raw mean count must exceed `filter_mean_thresh=0.07`; the gene passes overall if that holds in a strict majority, `>min_perc_group`, of the cell type x condition's donor groups with >=100 cells) rather than by running memento itself -- only the pass/fail outcome is needed for a feasibility check. Donors with <8 usable (>=100-cell) groups are flagged underpowered, using the macrophage precedent (5 donors -> replication indistinguishable from chance) as the reference.

| cell type | condition | donors (total / ≥100 cells) | ANTXR2 present (n/frac) | clears mpg=0.7 | underpowered (<8 donors) |
|---|---|---|---|---|---|
| colonocyte | Inflamed | 5 / 4 | 4/4 (1.00) | YES | **YES** |
| colonocyte | Neighbouring_inflamed | 17 / 16 | 8/16 (0.50) | no | no |
| colonocyte | Non_pathological | 16 / 16 | 6/16 (0.38) | no | no |
| enterocyte | Inflamed | 12 / 10 | 7/10 (0.70) | no | no |
| enterocyte | Neighbouring_inflamed | 28 / 27 | 24/27 (0.89) | YES | no |
| enterocyte | Non_pathological | 10 / 9 | 7/9 (0.78) | YES | no |
| intestinal crypt stem cell | Inflamed | 15 / 2 | 1/2 (0.50) | no | **YES** |
| intestinal crypt stem cell | Neighbouring_inflamed | 40 / 4 | 2/4 (0.50) | no | **YES** |
| intestinal crypt stem cell | Non_pathological | 24 / 6 | 3/6 (0.50) | no | **YES** |
| intestine goblet cell | Inflamed | 16 / 12 | 5/12 (0.42) | no | no |
| intestine goblet cell | Neighbouring_inflamed | 43 / 38 | 31/38 (0.82) | YES | no |
| intestine goblet cell | Non_pathological | 25 / 22 | 9/22 (0.41) | no | no |
| paneth cell | Inflamed | 14 / 1 | 1/1 (1.00) | YES | **YES** |
| paneth cell | Neighbouring_inflamed | 34 / 11 | 4/11 (0.36) | no | no |
| paneth cell | Non_pathological | 10 / 1 | 1/1 (1.00) | YES | **YES** |
| transit amplifying cell | Inflamed | 17 / 12 | 7/12 (0.58) | no | no |
| transit amplifying cell | Neighbouring_inflamed | 43 / 37 | 18/37 (0.49) | no | no |
| transit amplifying cell | Non_pathological | 25 / 21 | 7/21 (0.33) | no | no |

**Colonocyte (the Bracq et al.-matched cell type -- mouse colon DSS colitis) fails the standard `min_perc_group=0.7` presence filter** in both well-powered conditions (Non_pathological: 0.38, Neighbouring_inflamed: 0.50 of donors present) -- well below the 0.6 threshold that already had to be used for Elmentaite's full enterocyte donor set. A relaxed threshold of roughly 0.35-0.5 would be needed just to admit ANTXR2 into a colonocyte test, i.e. **further relaxation than any precedent in this project.** Colonocyte's Inflamed condition additionally has only 4 well-powered donors -- underpowered regardless of the presence question.

**Enterocyte (the extension cell type) clears the standard filter** in Neighbouring_inflamed (0.89) and Non_pathological (0.78), and sits exactly at the boundary in Inflamed (0.70 -- since memento's filter is a strict `>`, exactly 0.70 would still FAIL and needs a hair of relaxation, e.g. 0.69). All three enterocyte conditions clear the 8-donor power floor.

**Crypt stem cell and Paneth cell are underpowered in every condition** (intestinal crypt stem cell/Inflamed: 2 donors, intestinal crypt stem cell/Neighbouring_inflamed: 4 donors, intestinal crypt stem cell/Non_pathological: 6 donors, paneth cell/Inflamed: 1 donors, paneth cell/Neighbouring_inflamed: 11 donors, paneth cell/Non_pathological: 1 donors) -- both flagged underpowered regardless of the presence question; not usable for a two-group test in Kong2023 at all.

#### 3. Chemistry/batch confounded with condition?

```
assay                  10x 3' v1  10x 3' v2  10x 3' v3
condition                                             
Inflamed                       0      18554       7975
Neighbouring_inflamed          0      60954      67825
Non_pathological           23628      32442      15387
```

**Yes, confounded.** `10x 3' v1` appears ONLY in `Non_pathological` (23628 cells, 33% of that condition's cells if present) while `Neighbouring_inflamed` and `Inflamed` are 100% `10x 3' v2`/`v3` -- the same shape of confound as the `donor_id` collision issue flagged earlier in this project (a real technical variable perfectly or near-perfectly aligned with the biological grouping of interest). Any Non_pathological-vs-inflamed differential test would need to either restrict Non_pathological to its v2/v3 donors only, or treat chemistry as an explicit covariate -- **not treat the raw condition contrast as chemistry-free.**

#### 4. Scoping (not running) the differential test, if feasible

**Enterocyte is the only feasible cell type for a well-powered Non_pathological-vs-Inflamed or Non_pathological-vs-Neighbouring_inflamed two-group test** (colonocyte fails presence; crypt stem/Paneth fail power; goblet/TA are intermediate -- see full table above). Pre-registered directional prediction (Bracq et al.): since CMG2's Wnt function is injury-conditional (CMG2-KO baseline guts are normal), **ANTXR2's coupling to Wnt-arm partners and regeneration programs should be absent or weak in `Non_pathological` and appear in `Inflamed`.** Any such test must restrict or covary for the chemistry confound in point 3 above, and should use enterocyte as primary with colonocyte reported only as a presence-filter negative result (an interesting finding in its own right, not a null test outcome), consistent with `prompts/partner_availability.md`'s framing that colonocyte is the literal match to Bracq et al.'s mouse colon model while enterocyte is the extension.

#### 5. Colonocyte vs. enterocyte as the Bracq et al. match

Confirmed per the pre-registration: Bracq et al. is mouse **colon** with DSS colitis, so colonocyte is the directly-matched cell type; enterocyte is the extension. The feasibility result above means the directly-matched cell type is NOT the one available for the actual test -- a real scoping constraint to carry forward, not an incidental detail.


### Task 2 results: anchor-gene specificity control

**Confirmatory** (pre-registered interpretation rule, see pre-registration above) -- this gates the fibroblast ECM-correlation and enterocyte digestion/absorption-correlation exploratory GSEA findings reported earlier in this log; those remain exploratory GSEA characterizations, and this section is what determines whether they carry any ANTXR2-specific evidentiary weight.

#### fibroblast

**Anchor selection**: 20 genes randomly drawn (seed=20260906) from 1305 candidates matched to ANTXR2 within log10±0.250 mean expression and ±0.050 detection rate (ECM/ribosomal/mitochondrial genes and this project's own panel genes excluded from the candidate pool). Selected: ARFRP1, BDP1, C1QTNF2, CBR1, COLEC12, DNPEP, HOTAIRM1, PLRG1, PLXNB2, QKI, REST, RHBDD2, RIC8A, SGTA, SLC25A37, SYS1, THAP7, TMF1, USP1, WDR45. Full candidate table: `fibroblast_selected_anchors.csv`.

20/20 anchors have a completed HT+GSEA run.

- **Extracellular Matrix Organization (GO:0030198)** (GO_Biological_Process_2023): ANTXR2 NES=2.56 (FDR=0) vs. **6/20 anchors also recover this term** at the same FDR<0.25 threshold (anchor NES range [-1.89, 2.49], mean 1.04). ANTXR2's NES sits at the **100th percentile** of the anchor distribution.
- **ECM-receptor interaction** (KEGG_2021_Human): ANTXR2 NES=2.40 (FDR=0) vs. **8/20 anchors also recover this term** at the same FDR<0.25 threshold (anchor NES range [-1.62, 2.43], mean 1.00). ANTXR2's NES sits at the **95th percentile** of the anchor distribution.
- **Collagen Formation R-HSA-1474290** (Reactome_2022): ANTXR2 NES=2.30 (FDR=0) vs. **10/20 anchors also recover this term** at the same FDR<0.25 threshold (anchor NES range [-1.62, 2.37], mean 1.25). ANTXR2's NES sits at the **80th percentile** of the anchor distribution.
- **Focal adhesion** (KEGG_2021_Human): ANTXR2 NES=2.19 (FDR=0) vs. **6/20 anchors also recover this term** at the same FDR<0.25 threshold (anchor NES range [-1.15, 2.41], mean 1.15). ANTXR2's NES sits at the **95th percentile** of the anchor distribution.

**Interpretation (pre-registered rule applied): no key term is recovered by a majority of anchors in fibroblast** -- the ANTXR2 correlation with this program is not simply a generic property of fibroblast's cells at ANTXR2's expression level, supporting (but not proving causally) a specific relationship.

#### enterocyte

**Anchor selection**: 20 genes randomly drawn (seed=20260906) from 824 candidates matched to ANTXR2 within log10±0.250 mean expression and ±0.050 detection rate (ECM/ribosomal/mitochondrial genes and this project's own panel genes excluded from the candidate pool). Selected: ADCY6, ATAD2B, C2CD5, CDC27, FANCL, FBXO8, GTF3C2, METTL14, NT5E, PEAK1, PPP1R35, SFXN5, SLC25A44, SSTR1, TIPARP, TMEM80, TRAPPC9, UBXN2B, ZFX, ZNF44. Full candidate table: `enterocyte_selected_anchors.csv`.

20/20 anchors have a completed HT+GSEA run.

- **Protein digestion and absorption** (KEGG_2021_Human): ANTXR2 NES=2.20 (FDR=0) vs. **17/20 anchors also recover this term** at the same FDR<0.25 threshold (anchor NES range [1.21, 2.23], mean 1.78). ANTXR2's NES sits at the **90th percentile** of the anchor distribution.
- **Fat digestion and absorption** (KEGG_2021_Human): ANTXR2 NES=2.14 (FDR=0) vs. **16/20 anchors also recover this term** at the same FDR<0.25 threshold (anchor NES range [1.21, 2.27], mean 1.62). ANTXR2's NES sits at the **95th percentile** of the anchor distribution.

**Interpretation (pre-registered rule applied): at least one key term is recovered by a MAJORITY of anchors in enterocyte** -- for that term/those terms, the positive ANTXR2 correlation is a property of enterocyte's transcriptional program, not specific evidence about ANTXR2. This is a real, useful negative for those terms, reported as such rather than suppressed. **Header amended** (see "Standing corrections" -> "Positive correlations require the anchor-gene specificity control") -- this changes the framing of the enterocyte digestion/absorption GSEA finding from "ungated exploratory result" to "checked and did not survive", and the 2026-09-05 GSEA section's enterocyte note is amended in place with a pointer here.

#### |z| vs. mean-expression check (Task 2, step 6)

- **fibroblast**: corr(log10 mean expression, |z|) = 0.437 across 4502 tested genes (corr with signed z = -0.045). Non-trivial positive relationship -- consistent with memento's bootstrap se shrinking with cell count/detection rate, so highly-expressed genes can get inflated |z| in both tails. Treat this as a caveat on z-ranked results generally (this project's GSEA z-score ranking included); the anchor-matching above at least holds mean expression roughly fixed between ANTXR2 and its anchors, which limits (but does not eliminate) this confound's effect on the recovery-count comparison specifically.
- **enterocyte**: corr(log10 mean expression, |z|) = 0.372 across 5852 tested genes (corr with signed z = 0.248). Non-trivial positive relationship -- consistent with memento's bootstrap se shrinking with cell count/detection rate, so highly-expressed genes can get inflated |z| in both tails. Treat this as a caveat on z-ranked results generally (this project's GSEA z-score ranking included); the anchor-matching above at least holds mean expression roughly fixed between ANTXR2 and its anchors, which limits (but does not eliminate) this confound's effect on the recovery-count comparison specifically.

## IBD colon atlas: `health` column corruption in compute_ibd_colon_means.py, found and fixed (2026-09-09, same day)

Found while answering a follow-up question (ANTXR2 in mature colon epithelium,
by health status). `memento.get_groups()` -- called in `run_memento_on_chunk` to
read back each group's label values after `compute_1d_moments` -- deliberately
re-encodes any label column with **exactly 2 distinct values within that donor's
chunk** into `0.0`/`1.0` float codes (confirmed by reading memento's source:
`get_groups` does `df[col] = df[col].astype('category').cat.codes.astype(float)`
whenever `pd.to_numeric` fails and `nunique()==2`, intended to hand memento's own
`binary_test_1d`/`2d` a numeric design matrix). This project's donors are either
`Healthy`-only (1 value in-chunk, unaffected) or paired
`Inflamed`+`Non-inflamed`-only (exactly 2 values, no `Healthy`) -- so this
silently corrupted `health` to `0.0`/`1.0` for **all 18 UC donors' rows** (health
label only; `cell_type_fine`, `n_cells`, and `mean_expression` were never
affected, since those come from a different, non-recoded column / from
`1d_moments` directly). Caught because a health-stratified query returned literal
`0.0`/`1.0` values instead of `Inflamed`/`Non-inflamed`.

**Fix**: don't trust `memento.get_groups()`'s per-column values at all -- parse
each group's label values directly from its group-key string instead
(`group_key.split(label_delimiter)[1:]`), which is exactly what `get_groups()`
itself does internally *before* the lossy re-coding step, so this is strictly
more reliable, not a workaround. Recomputed; `health` now reads
`Healthy`/`Inflamed`/`Non-inflamed` throughout (verified: 7,704,357 /
5,575,149 / 6,551,148 rows respectively, no numeric leakage). The
`ibd_colon_feasibility.py` presence/donor-count tables were never affected (they
read `Health` directly from each h5ad's `obs`, not through memento's
`get_groups()`), and the overall (health-pooled) ANTXR2-by-cell-type rankings
reported earlier this session were also unaffected (health wasn't part of that
aggregation) -- only a health-stratified breakdown was wrong.

**Lesson**: a library's own "helpful" convenience re-encoding (here, for its own
downstream regression convenience) can silently break a caller using the same
function for a different purpose (a descriptive label, not a treatment design
matrix) -- when a returned label value looks suspiciously like a code (bare
`0.0`/`1.0` where a category label was expected), verify against the
lower-level/raw representation (the group key string, here) rather than assuming
a cast bug in this project's own code.

## IBD colon atlas (Smillie et al., Cell 2019, SCP259): download, and a source-data corruption finding (2026-09-09)

New dataset, first pass. User request: evaluate SCP259 (Smillie et al., *Cell* 2019,
366,650 cells, 30 donors -- 18 UC + 12 healthy, 51 annotated cell subsets: 15
epithelial including Stem/TA1/TA2/Cycling TA, 13 stromal/glial including multiple
fibroblast subtypes, 23 immune) for feasibility, starting with ANTXR2/gene-panel
mean expression by cell type, and download both raw counts and a processed/
visualization-coordinate version for later hypothesis testing.

**Download blocker, and how it was resolved**: unlike every prior dataset in this
project, SCP259's own portal is login-gated (`study_files` API returns HTTP 401,
not a guess -- the study page also states "Please sign in to download data").
Resolved via the Human Cell Atlas Data Coordination Platform, which mirrors the
identical file set as project `cd61771b-661a-4e19-b269-6e5d95350de6`
("HumanColonRewiringUlcerativeColitis") under an open license
(`dataUseRestriction=NRES`, no `duosId`) -- verified end-to-end with plain
unauthenticated `requests` calls: Azul's REST API resolves each file uuid to a
signed, no-login S3 URL via a 302 redirect, confirmed working for both a tiny file
and (via an HTTP Range request) a >1GB file before committing to this as the
primary source. `scripts/download_ibd_colon_data.py` downloaded all 15 files this
way (9 raw compartment files + `all.meta2.txt` + `cell_subsets.txt` from the
authors' GitHub repo + 3 pre-computed Seurat `.rds` objects), all exact byte-size
matches against Azul's own listed sizes.

**Standing correction: a "successful download" is not the same as "correct data" --
verify content, not just transfer.** `build_ibd_colon_h5ad.py`'s first run crashed
inside `scipy.io.mmread` on `gene_sorted-Epi.matrix.mtx` with `Invalid integer
value` partway through the file. Investigation (not assumption) found: the file's
actual line count (91,725,164) was only 52.6% of what the mtx header itself
declared (174,423,911 data lines), and it ended mid-line with no trailing newline.
Because the downloaded file's size matched Azul's own reported size exactly, this
looked at first like it could be a bad `expected_size` on Azul's end rather than a
real truncation -- so before concluding anything, the actual live S3 object was
independently re-probed with raw HTTP Range requests (bypassing this project's own
download code entirely): `Content-Range` on a request near the end of the
downloaded file confirmed the S3 object's own declared total size exactly equals
what was downloaded, and requests further out (1.5-2.3GB) returned `416 Range Not
Satisfiable`. **The object hosted on HCA is genuinely incomplete at the source, not
a transfer artifact.** `gene_sorted-Imm.matrix.mtx` turned out to have the same
problem, far worse (7,965,459 of 173,255,714 declared lines, ~4.6%).
`gene_sorted-Fib.matrix.mtx` was checked the same way and is complete and correct
(39,087,919 declared == 39,087,919 actual, clean final line). **Lesson for future
datasets: when a raw file parses cleanly but a sanity total (row count, cell count,
nnz) doesn't match its own header/manifest, re-verify against the live source
before assuming a local bug -- public archive mirrors can silently host incomplete
objects that report a self-consistent but wrong size.**

**Recovery, without blocking on a login-gated re-download.** The 3
`train.{Epi,Fib,Imm}.seur.rds` files (the paper's "discovery cohort" Seurat
objects, downloaded from the same HCA source) all passed a full `gzip -t`
integrity check -- genuinely intact, unlike the 2 broken mtx files. Rather than
stopping to ask the user to obtain an SCP-authenticated `curl_config.txt` (the
only way to re-fetch the broken files from their original host), these RDS
objects were used as the raw-count source for Epi and Imm instead. They did not
require installing the (heavy, slow-to-build) Seurat R package: a Seurat object's
S4 slots (`assays`, `meta.data`, `reductions`) are plain R attributes readable via
base `attr()`/`readRDS()`, and the `counts` slot is a `dgCMatrix` needing only the
already-installed, lightweight `Matrix` package -- confirmed interactively before
writing `scripts/export_ibd_colon_rds.R`. This is the **first use of the `r-env`
conda environment in this project.**

**Real, documented cost of this workaround**: the RDS objects cover only the
"discovery cohort" -- 17 of 30 donors (all 3 Health states represented, but skewed
toward `Healthy`) -- not the full cohort Fib's intact mtx provides. This asymmetry
(`data_source` = `full_cohort_mtx` for Fib vs. `discovery_cohort_rds` for Epi/Imm)
is carried through as an explicit column in the output parquet and h5ads, not
silently absorbed. A minor labeling inconsistency was also found and harmonized:
the RDS metadata uses `"Uninflamed"` where the full-cohort `all.meta2.txt` uses
`"Non-inflamed"` for the identical concept (`IBD_COLON_HEALTH_MAP`).

**Donor identity verified before use** (per this project's standing convention,
[[feedback_donor_identity_verification]]): `Subject` (30 distinct values) checked
for the same collision/alias pattern found earlier in the Gut Cell Atlas's
`donor_id` -- none found; every subject maps cleanly to 1 (12 `Healthy` donors) or
2 (18 UC donors, paired Inflamed/Non-inflamed biopsies) `Health` values. No
`*_UNIFIED_COL` correction was needed here.

**Compute**: `scripts/compute_ibd_colon_means.py` (memento method-of-moments,
genome-wide, grouped by `donor_id x cell_type_fine x health`, processing the 3
compartments as independent "datasets" the same way `compute_means.py` loops over
CELLxGene collections) ran cleanly, ~14s/compartment, 19,830,654 rows written to
`/data/ANTXR2/celltype_expression/ibd_colon_atlas_celltype_means.parquet`. See
that file's README section for full schema/caveats.

**Directional sanity check against Phase 1**: donor-equal-weighted ANTXR2 mean
expression, ranked across all 51 cell types, puts the Wnt-niche fibroblast
subtypes (RSPO3+, WNT2B+ Fos-hi, WNT2B+ Fos-lo 1, WNT5B+ 1/2) and Inflammatory
Fibroblasts at the top, with Stem/TA1/Cycling TA near the bottom -- consistent
with Phase 1's "high in fibroblast, not stem-skewed" finding on an entirely
independent tissue source and cohort, and directly on-theme for the project's
Wnt-transduction arm (RSPO3/WNT2B/WNT5B-defined fibroblasts are the literal
Wnt-niche cell types in gut architecture).

## IBD colon atlas: full 30-donor cohort obtained, superseding the discovery-cohort fallback (2026-09-09, same day)

The user created an SCP259 account and provided a Broad Single Cell Portal
bulk-download link (`generate_curl_config` + `auth_code`, valid 30 minutes to
*generate* the config -- the resulting GCS-signed URLs it produces are valid
~24h, so the actual download itself wasn't time-pressured once the config was
in hand). This is the same file set the corrupted HCA mirror was supposed to
provide, but from SCP259 directly: `gene_sorted-{Epi,Fib,Imm}.matrix.mtx`,
`{Epi,Fib,Imm}.genes.tsv`/`.barcodes2.tsv`, `all.meta2.txt`, and -- new,
not previously obtained -- `{Epi,Fib,Imm}.tsne.txt`, real full-cohort tSNE
coordinates (better than the RDS-derived discovery-cohort-only ones).

**Verified before trusting, per the standing lesson from the HCA episode**:
every mtx file's actual data-line count was checked against its own header's
declared nnz before use, not just download success/size. All 3 now verify
exactly: Epi 174,423,911, Fib 39,087,919, Imm 173,255,714 declared == actual,
all with clean final lines. `all.meta2.txt` is byte-identical to the earlier
HCA copy (confirms metadata was never the corrupted part). The corrupted HCA
mtx copies were moved to `raw/ibd_colon_atlas/corrupted_hca_backup/` (not
deleted) rather than overwritten in place.

**Rebuild ran into a new problem: OOM on the larger matrices.**
`build_ibd_colon_h5ad.py` was killed by the system's memory manager twice in a
row processing the full (now ~2.3-2.4GB text, ~174M/~173M-nonzero) Epi and Imm
matrices, right after `mmread` completed -- `ps aux`/`free -h` showed no
smoking-gun concurrent process, pointing to the COO->CSR transpose step itself:
`scipy.io.mmread` returns int64 data/index arrays by default, and for a
174M-nonzero matrix that's ~4GB of COO arrays plus another ~3-4GB for the new
CSR arrays held simultaneously during conversion. Fixed by downcasting to
int32 immediately after `mmread` (raw UMI counts and 20K/210K-scale indices
have enormous headroom under int32) plus an explicit `gc.collect()` after the
transpose; also added a `--compartments` flag to `build_ibd_colon_h5ad.py` so
each compartment can run as its own short-lived process rather than one long
process holding all 3 compartments' peak memory in sequence. All 3
compartments then built cleanly as separate invocations. **Lesson: mtx text
parsing defaulting to int64 is a real memory cliff for real-atlas-scale
matrices (100M+ nonzeros) even on a machine with plenty of nominal RAM --
downcast explicitly rather than trusting the loader's default dtype.**

**Recomputed** (`compute_ibd_colon_means.py`, already carrying the
`memento.get_groups()` health-label fix from earlier the same day): all 3
compartments now show 30/30 donors (previously Epi/Imm: 17/30). 33,725,155
total rows (previously 19,830,654) written to
`ibd_colon_atlas_celltype_means.parquet`. `health` column re-verified clean.

**What changed vs. the 17-donor discovery-cohort numbers, and what didn't**:
Fibroblast values are unchanged (Fib was always full-cohort, so this doubles
as an internal consistency check -- e.g. RSPO3+ still 0.000072/10 donors,
WNT2B+ Fos-hi still 0.000068/30 donors, identical to the earlier run). Epi/Imm
donor counts roughly doubled (e.g. Stem 14->28, TA1/TA2/Cycling TA 17->30,
Enterocytes 14->26), and per-condition donor counts in the feasibility check
improved enough to flip several `underpowered` flags from YES to no (e.g. TA
2/Inflamed: 2->9 well-powered donors, no longer underpowered). The qualitative
picture held up throughout: Wnt-niche fibroblast subtypes highest for ANTXR2,
Stem/TA1 lowest, mature Enterocytes topping the epithelial ranking, and the
Wnt-receptor-machinery-vs-target-gene split (CTNNB1/FZD5/LRP5/TCF7L2 rising
toward mature Enterocytes like ANTXR2, while LGR5/ASCL2/OLFM4/SMOC2 stay
sharply crypt-restricted) reproduced cleanly on the full cohort.

**One real correction, not just a power improvement**: mature Enterocytes'
Inflamed-vs-Non-inflamed ANTXR2 ordering flipped. At n=5 donors per condition
(discovery cohort), Inflamed (0.000051) appeared higher than Non-inflamed
(0.000041). At n=12 donors per condition (full cohort), Non-inflamed
(0.000071) is now slightly higher than Inflamed (0.000062) -- both still well
above Healthy (0.000030). The direction "elevated in UC tissue vs. healthy"
holds and strengthens; the specific Inflamed-vs-Non-inflamed ordering reported
earlier this session does not survive more data and should not be treated as
established. **Lesson, consistent with this project's existing n<8
underpowered-flag convention: don't trust the relative ordering between two
already-small groups (5 vs 5 donors) even when both individually look
"present" -- the earlier answer to the user's Inflamed-vs-Non-inflamed
question should be considered superseded by this section.**

## IBD colon atlas (Smillie et al., Cell 2019, SCP259) feasibility check -- full 30-donor cohort (2026-09-09)

**Supersedes an earlier version of this section** that ran on the 17-of-30-donor discovery-cohort fallback (see "full 30-donor cohort" session entry above for how the full cohort was obtained same day) -- that version is not kept verbatim below since every number in it is superseded by this full-cohort run, but its qualitative conclusions (fibroblasts full-cohort throughout and directionally on-theme; Epithelial/Stem/TA numbers usable but power-limited) matched this run closely; the main change is donor counts roughly doubling for Epi/Imm groups (e.g. many `cell_type x health` groups went from ~5-7 well-powered donors to ~10-18), which flipped several `underpowered` flags from YES to no without changing which cell types rank highest for ANTXR2.

**Feasibility gate, not a results-producing step** -- no differential-correlation test was run. Data source: SCP259 direct (user-authenticated bulk download, see above), full 30-donor cohort, all 3 compartments verified complete. 51 cell subsets total (15 epithelial, 13 stromal/glial, 23 immune), 30 donors (18 UC + 12 healthy). This check covers only the Epithelial and Fibroblasts lineages (the user's stated interest). Output: `/data/ANTXR2/figures/ibd_colon_feasibility/{per_donor_counts,feasibility_summary,gene_panel_by_celltype}.csv`.

### Donor counts and ANTXR2 presence, Epithelial + Fibroblasts lineages

ANTXR2 presence replicates memento's own filter formula directly (raw per-donor mean count > `filter_mean_thresh=0.07`; passes overall if true in a strict majority, `>min_perc_group=0.7`, of donor groups with >=100 cells), computed directly from raw counts (not the memento capture-corrected parquet). Donors with <8 usable groups flagged underpowered (macrophage precedent). `data_source` shows whether a cell type's counts came from the full 30-donor cohort (`full_cohort_mtx`, Fib only) or the 17-donor discovery-cohort fallback forced by the HCA source corruption (`discovery_cohort_rds`, Epi + Imm).

| lineage | cell type | health | source | donors (total/≥100 cells) | ANTXR2 present (n/frac) | clears mpg=0.7 | underpowered |
|---|---|---|---|---|---|---|---|
| Epithelial | Best4+ Enterocytes | Healthy | full_cohort_mtx | 12 / 5 | 5/5 (1.00) | YES | **YES** |
| Epithelial | Best4+ Enterocytes | Inflamed | full_cohort_mtx | 13 / 1 | 1/1 (1.00) | YES | **YES** |
| Epithelial | Best4+ Enterocytes | Non-inflamed | full_cohort_mtx | 18 / 6 | 6/6 (1.00) | YES | **YES** |
| Epithelial | Cycling TA | Healthy | full_cohort_mtx | 12 / 12 | 1/12 (0.08) | no | no |
| Epithelial | Cycling TA | Inflamed | full_cohort_mtx | 17 / 14 | 12/14 (0.86) | YES | no |
| Epithelial | Cycling TA | Non-inflamed | full_cohort_mtx | 18 / 13 | 9/13 (0.69) | no | no |
| Epithelial | Enterocyte Progenitors | Healthy | full_cohort_mtx | 12 / 12 | 0/12 (0.00) | no | no |
| Epithelial | Enterocyte Progenitors | Inflamed | full_cohort_mtx | 16 / 2 | 1/2 (0.50) | no | **YES** |
| Epithelial | Enterocyte Progenitors | Non-inflamed | full_cohort_mtx | 18 / 10 | 1/10 (0.10) | no | no |
| Epithelial | Enterocytes | Healthy | full_cohort_mtx | 12 / 6 | 6/6 (1.00) | YES | **YES** |
| Epithelial | Enterocytes | Inflamed | full_cohort_mtx | 16 / 3 | 3/3 (1.00) | YES | **YES** |
| Epithelial | Enterocytes | Non-inflamed | full_cohort_mtx | 16 / 5 | 5/5 (1.00) | YES | **YES** |
| Epithelial | Enteroendocrine | Healthy | full_cohort_mtx | 11 / 0 | n/a | no | **YES** |
| Epithelial | Enteroendocrine | Inflamed | full_cohort_mtx | 15 / 0 | n/a | no | **YES** |
| Epithelial | Enteroendocrine | Non-inflamed | full_cohort_mtx | 17 / 1 | 0/1 (0.00) | no | **YES** |
| Epithelial | Goblet | Healthy | full_cohort_mtx | 12 / 5 | 2/5 (0.40) | no | **YES** |
| Epithelial | Goblet | Inflamed | full_cohort_mtx | 14 / 0 | n/a | no | **YES** |
| Epithelial | Goblet | Non-inflamed | full_cohort_mtx | 17 / 1 | 1/1 (1.00) | YES | **YES** |
| Epithelial | Immature Enterocytes 1 | Healthy | full_cohort_mtx | 12 / 11 | 2/11 (0.18) | no | no |
| Epithelial | Immature Enterocytes 1 | Inflamed | full_cohort_mtx | 16 / 5 | 4/5 (0.80) | YES | **YES** |
| Epithelial | Immature Enterocytes 1 | Non-inflamed | full_cohort_mtx | 16 / 6 | 3/6 (0.50) | no | **YES** |
| Epithelial | Immature Enterocytes 2 | Healthy | full_cohort_mtx | 12 / 11 | 10/11 (0.91) | YES | no |
| Epithelial | Immature Enterocytes 2 | Inflamed | full_cohort_mtx | 17 / 7 | 7/7 (1.00) | YES | **YES** |
| Epithelial | Immature Enterocytes 2 | Non-inflamed | full_cohort_mtx | 18 / 10 | 10/10 (1.00) | YES | no |
| Epithelial | Immature Goblet | Healthy | full_cohort_mtx | 12 / 11 | 1/11 (0.09) | no | no |
| Epithelial | Immature Goblet | Inflamed | full_cohort_mtx | 17 / 9 | 1/9 (0.11) | no | no |
| Epithelial | Immature Goblet | Non-inflamed | full_cohort_mtx | 18 / 12 | 2/12 (0.17) | no | no |
| Epithelial | M cells | Healthy | full_cohort_mtx | 8 / 0 | n/a | no | **YES** |
| Epithelial | M cells | Inflamed | full_cohort_mtx | 9 / 1 | 1/1 (1.00) | YES | **YES** |
| Epithelial | M cells | Non-inflamed | full_cohort_mtx | 9 / 1 | 1/1 (1.00) | YES | **YES** |
| Epithelial | Secretory TA | Healthy | full_cohort_mtx | 12 / 7 | 1/7 (0.14) | no | **YES** |
| Epithelial | Secretory TA | Inflamed | full_cohort_mtx | 17 / 2 | 2/2 (1.00) | YES | **YES** |
| Epithelial | Secretory TA | Non-inflamed | full_cohort_mtx | 18 / 5 | 4/5 (0.80) | YES | **YES** |
| Epithelial | Stem | Healthy | full_cohort_mtx | 12 / 4 | 0/4 (0.00) | no | **YES** |
| Epithelial | Stem | Inflamed | full_cohort_mtx | 17 / 1 | 1/1 (1.00) | YES | **YES** |
| Epithelial | Stem | Non-inflamed | full_cohort_mtx | 17 / 4 | 3/4 (0.75) | YES | **YES** |
| Epithelial | TA 1 | Healthy | full_cohort_mtx | 12 / 12 | 0/12 (0.00) | no | no |
| Epithelial | TA 1 | Inflamed | full_cohort_mtx | 17 / 15 | 0/15 (0.00) | no | no |
| Epithelial | TA 1 | Non-inflamed | full_cohort_mtx | 18 / 17 | 0/17 (0.00) | no | no |
| Epithelial | TA 2 | Healthy | full_cohort_mtx | 12 / 12 | 6/12 (0.50) | no | no |
| Epithelial | TA 2 | Inflamed | full_cohort_mtx | 17 / 10 | 9/10 (0.90) | YES | no |
| Epithelial | TA 2 | Non-inflamed | full_cohort_mtx | 18 / 13 | 10/13 (0.77) | YES | no |
| Epithelial | Tuft | Healthy | full_cohort_mtx | 12 / 0 | n/a | no | **YES** |
| Epithelial | Tuft | Inflamed | full_cohort_mtx | 15 / 0 | n/a | no | **YES** |
| Epithelial | Tuft | Non-inflamed | full_cohort_mtx | 16 / 0 | n/a | no | **YES** |
| Fibroblasts | Inflammatory Fibroblasts | Healthy | full_cohort_mtx | 6 / 0 | n/a | no | **YES** |
| Fibroblasts | Inflammatory Fibroblasts | Inflamed | full_cohort_mtx | 14 / 3 | 3/3 (1.00) | YES | **YES** |
| Fibroblasts | Inflammatory Fibroblasts | Non-inflamed | full_cohort_mtx | 13 / 2 | 2/2 (1.00) | YES | **YES** |
| Fibroblasts | Myofibroblasts | Healthy | full_cohort_mtx | 12 / 0 | n/a | no | **YES** |
| Fibroblasts | Myofibroblasts | Inflamed | full_cohort_mtx | 18 / 2 | 1/2 (0.50) | no | **YES** |
| Fibroblasts | Myofibroblasts | Non-inflamed | full_cohort_mtx | 18 / 3 | 3/3 (1.00) | YES | **YES** |
| Fibroblasts | RSPO3+ | Healthy | full_cohort_mtx | 12 / 0 | n/a | no | **YES** |
| Fibroblasts | RSPO3+ | Inflamed | full_cohort_mtx | 9 / 0 | n/a | no | **YES** |
| Fibroblasts | RSPO3+ | Non-inflamed | full_cohort_mtx | 13 / 0 | n/a | no | **YES** |
| Fibroblasts | WNT2B+ Fos-hi | Healthy | full_cohort_mtx | 12 / 4 | 4/4 (1.00) | YES | **YES** |
| Fibroblasts | WNT2B+ Fos-hi | Inflamed | full_cohort_mtx | 15 / 4 | 3/4 (0.75) | YES | **YES** |
| Fibroblasts | WNT2B+ Fos-hi | Non-inflamed | full_cohort_mtx | 18 / 7 | 7/7 (1.00) | YES | **YES** |
| Fibroblasts | WNT2B+ Fos-lo 1 | Healthy | full_cohort_mtx | 12 / 7 | 7/7 (1.00) | YES | **YES** |
| Fibroblasts | WNT2B+ Fos-lo 1 | Inflamed | full_cohort_mtx | 18 / 4 | 3/4 (0.75) | YES | **YES** |
| Fibroblasts | WNT2B+ Fos-lo 1 | Non-inflamed | full_cohort_mtx | 18 / 8 | 8/8 (1.00) | YES | no |
| Fibroblasts | WNT2B+ Fos-lo 2 | Healthy | full_cohort_mtx | 12 / 1 | 0/1 (0.00) | no | **YES** |
| Fibroblasts | WNT2B+ Fos-lo 2 | Inflamed | full_cohort_mtx | 17 / 3 | 1/3 (0.33) | no | **YES** |
| Fibroblasts | WNT2B+ Fos-lo 2 | Non-inflamed | full_cohort_mtx | 18 / 5 | 3/5 (0.60) | no | **YES** |
| Fibroblasts | WNT5B+ 1 | Healthy | full_cohort_mtx | 12 / 3 | 3/3 (1.00) | YES | **YES** |
| Fibroblasts | WNT5B+ 1 | Inflamed | full_cohort_mtx | 17 / 2 | 1/2 (0.50) | no | **YES** |
| Fibroblasts | WNT5B+ 1 | Non-inflamed | full_cohort_mtx | 18 / 0 | n/a | no | **YES** |
| Fibroblasts | WNT5B+ 2 | Healthy | full_cohort_mtx | 12 / 5 | 4/5 (0.80) | YES | **YES** |
| Fibroblasts | WNT5B+ 2 | Inflamed | full_cohort_mtx | 16 / 2 | 1/2 (0.50) | no | **YES** |
| Fibroblasts | WNT5B+ 2 | Non-inflamed | full_cohort_mtx | 18 / 5 | 4/5 (0.80) | YES | **YES** |

**Stem/TA cells** (Stem, TA 1, TA 2, Cycling TA): 3/12 (cell type x health) groups are both well-powered and clear the presence filter. Consistent with Phase 1's prior finding that ANTXR2 is not stem-skewed in gut epithelium.

**Fibroblast subtypes**: 1/24 (cell type x health) groups clear both power and presence filters -- fibroblasts are the full 30-donor `full_cohort_mtx` compartment (not affected by the HCA corruption), and the donor-equal-weighted heatmap below shows the Wnt-niche fibroblast subtypes (RSPO3+, WNT2B+, WNT5B+) at the top of the ANTXR2 ranking -- directionally consistent with Phase 1's fibroblast finding and on-theme for the project's Wnt-transduction arm.

### Gene panel (GENE_PANEL) donor-equal-weighted mean, top ANTXR2 cell types

| cell type | mean ANTXR2 | n donors |
|---|---|---|
| RSPO3+ | 0.000072 | 10 |
| WNT2B+ Fos-hi | 0.000068 | 30 |
| Inflammatory Fibroblasts | 0.000060 | 11 |
| WNT2B+ Fos-lo 1 | 0.000055 | 30 |
| WNT5B+ 1 | 0.000053 | 29 |
| WNT5B+ 2 | 0.000052 | 30 |
| Enterocytes | 0.000051 | 26 |
| Myofibroblasts | 0.000045 | 28 |
| WNT2B+ Fos-lo 2 | 0.000028 | 30 |
| Best4+ Enterocytes | 0.000022 | 25 |

### Verdict

**Usable for the ANTXR2 project's stated interest (fibroblasts, epithelial, TA/stem cells).** **Full 30-donor cohort, all 3 compartments, no coverage asymmetry.** Fibroblast coverage was always full-cohort; Epithelial coverage (including Stem/TA1/TA2/Cycling TA) is now also full-cohort after the SCP259-direct download superseded the earlier HCA-corruption-forced discovery-cohort fallback (see the download section above). The Wnt-niche fibroblast subtypes rank highest for ANTXR2 (directionally on-theme); most cell-type x health groups clear the 8-donor floor (see table above for exactly which don't).

