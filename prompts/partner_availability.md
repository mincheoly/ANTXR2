# Partner-availability analysis: ANTXR2's two functional arms

## Context

Read `ANALYSIS_SUMMARY.md`'s **header section** (everything above the `---`
divider) before starting. It defines the project's framing, what each kind of
evidence can establish, and standing corrections that override anything in the
chronological log below the divider. This prompt implements items 1–5 of that
header's "Designed but never run" list.

**The framing in one line:** ANTXR2/CMG2 has two separable functions with
different partner requirements -- collagen-VI clearance (needs MRC2 +
lysosomal machinery + substrate) and Wnt signal transduction (needs LRP6 +
Frizzled) -- and the project is mapping which cell types have the partners for
which function.

**Two asymmetries govern interpretation throughout. Do not lose them:**

- **Absence is strong, presence is weak.** Partners absent → that function
  cannot operate (needs only means, survives every confound). Partners present
  → the function is *available*, not that it happens.
- **Anti-correlation is strong, positive correlation is weak.** Shared
  activation state inflates positive correlations but cannot manufacture
  mutual exclusivity. A significant negative is interpretable on its own; a
  significant positive is not interpretable until Task 2 has run.

Everything here is **confirmatory panel testing against pre-defined gene
sets**, not exploratory screening. Do not expand a panel mid-analysis because
a result looks interesting.

## Environment

- conda env `antxr2`. Run scripts as
  `conda run -n antxr2 python scripts/<name>.py`.
- Reuse existing patterns rather than reinventing:
  - Parquet queries: `scripts/gene_panel_query.py` (`pyarrow.parquet.ParquetFile`,
    row-group scan filtering `gene` via `pyarrow.compute.is_in` with a
    **`pa.array`, not `pd.array`** value set).
  - Aggregation for the whole-body atlas: **two-step donor-equal-weighted**
    (per-donor `n_cells`-weighted mean within a `cell_type`, then unweighted
    mean across donors). Never a single-step pooled weighted mean.
  - memento correlation pipeline: `scripts/coexpression_pipeline.py`, including
    the `_corr_from_cov` variance<=0 null-out fix and this pipeline's defaults
    (`shrinkage=0`, `trim_percent=0.5`).
- Append findings to `ANALYSIS_SUMMARY.md`'s **session log** (below the
  divider). If a finding changes the framing or a standing correction, amend
  the **header** explicitly and say so in the log entry.

## Tasks, in dependency order

### Task 1 -- ANTXR1/ANTXR2 paralog readout (do this first; near-zero cost)

The project's original primary analysis, never executed. The numbers already
exist -- this is a readout, not a computation.

**Pre-register before looking at any value.** Write the interpretation rule
into the summary *first*:
- Significant **negative** → mutual exclusivity → real evidence that ANTXR1
  and ANTXR2 occupy different cells, which is the paralog-partitioning result
  the project was designed around. Interpretable immediately, no specificity
  control needed.
- Significant **positive** → co-occupancy or co-regulation → **not
  interpretable until Task 2 establishes the nonspecific-covariation floor**,
  and requires the doublet check below.
- **Null → uninformative.** Explicitly NOT evidence of co-presence. Say so.

Steps:
1. Pull ANTXR1's row (`coef`, `se`, `pval`, `fdr`, `z`) from
   `/data/ANTXR2/figures/coexpression/full_dataset_ht/{cell_type}_full_dataset_ht.csv`
   for fibroblast and enterocyte. Report macrophage's value but label it
   exploratory-only and draw no conclusion from it.
2. Report ANTXR1's rank within each cell type's full tested universe, not just
   its raw value.
3. **Bimodality check.** Compute the per-cell expression distribution of
   ANTXR1 and ANTXR2 within each cell type. If either is near-uniformly
   expressed, correlation is a poor co-presence readout there and the mean
   already settles availability -- state this rather than over-reading the
   correlation.
4. **Doublet check** (only if the correlation is positive): confirm the atlas
   QC removed doublets, and check whether ANTXR2-high/ANTXR1-high cells are
   enriched for cross-lineage marker co-expression. Fibroblast-epithelial
   doublets would specifically inflate this pair.
5. Report the ANTXR1 value alongside the pre-fix 0.35 that appears earlier in
   the log, and state plainly that the pre-fix value was inflated by the
   `_corr_from_cov` bug and is superseded.

### Task 2 -- Anchor-gene specificity control

Establishes the nonspecific-covariation floor. **Gates every positive
correlation claim in the project**, including the fibroblast ECM result and
the enterocyte digestion/absorption result.

1. Per cell type, select ~20 anchor genes matched to ANTXR2 on **both** mean
   expression and detection rate (fraction of cells nonzero) within that cell
   type. Match within a tolerance band; record the band and the selected genes.
   Exclude known ECM, ribosomal, and mitochondrial genes from the anchor set
   so the null isn't contaminated by the modules under test.
2. Run each anchor gene through the **identical** pipeline ANTXR2 went through:
   `coexpression_full_ht.py` (one-sample bootstrap, same `num_boot`, same
   filters), producing a full ranked z-vector per anchor.
3. Run `coexpression_gsea.py --metric zscore` on each anchor's ranking, same
   gene-set libraries.
4. Report, per cell type: how many anchor genes recover the ECM/collagen terms
   (fibroblast) or digestion/absorption terms (enterocyte) at the same FDR,
   and where ANTXR2's NES for those terms falls in the anchor distribution
   (empirical percentile).
5. **Interpretation rule, also pre-registered:** if most anchors recover the
   same terms, the ANTXR2 result is a property of the cell type, not of
   ANTXR2. Report that outcome as clearly as the alternative -- it is a real
   and useful negative.
6. Also check whether `|z|` correlates with per-gene mean expression across
   the tested universe. memento's bootstrap `se` shrinks with cell count and
   detection rate, so highly expressed genes may get inflated `|z|` in **both**
   tails. If present, note it as a caveat on all z-ranked results and consider
   expression-decile matching.

### Task 3 -- Extend the gene panel to the Wnt arm

The current 8-gene panel covers only the clearance function, so partner
availability for the Wnt arm is currently untestable. Pure query against
already-computed parquet; no pipeline run.

1. Add to `GENE_PANEL` in `scripts/config.py`, as **named sub-panels** rather
   than one flat list:
   - `CLEARANCE_ARM`: `MRC2, CTSB, CTSK, MMP14, TIMP2, LAMP1` (existing)
   - `CLEARANCE_SUBSTRATE`: `COL6A1, COL6A2, COL6A3` (**new** -- the panel
     currently tests everything except whether the substrate is present)
   - `WNT_ARM`: `LRP5, LRP6, FZD1-10, CTNNB1, TCF7L2, LGR5, RNF43, ZNRF3,
     AXIN2` (**new**)
   - `RECEPTORS`: `ANTXR1, ANTXR2` (existing)
2. Query both atlases (whole-body + skin fibroblast) for the extended panel,
   reusing the existing aggregation exactly.
3. Produce a **two-arm heatmap**: same curated cell-type rows as the existing
   ECM-clearance figure, with the sub-panels as labeled column blocks so the
   two arms can be read against each other per row.
4. **No composite score of any kind** -- standing correction. Raw values side
   by side only.
5. Report explicitly, per cell type: which arm's partners are present, which
   are absent, and which cell types have ANTXR2 with **neither** arm complete.

### Task 4 -- Formal tests for the two means-level claims

Both are currently descriptive means with no significance test.

1. `binary_test_1d` on **MRC2**: gut epithelium (enterocyte, colonocyte,
   goblet, crypt stem cell) vs the lineage-matched comparison already chosen
   (keratinocyte, corneal epithelial cell). This tests the project's strongest
   result -- that the clearance arm is structurally unavailable in gut
   epithelium -- which currently rests on untested means.
2. `binary_test_1d` on **ANTXR2**: crypt stem/TA vs differentiated enterocyte.
   Now load-bearing: Lencer's commentary on Bracq et al. proposes CMG2's
   context-specific role may reflect the diminishing Wnt gradient along the
   crypt-villus axis, making ANTXR2's position on that axis a direct testable
   prediction. Existing means hint at it (crypt stem cell +1.21 vs enterocyte
   +2.12 log-ratio) but it is untested.
3. Run the same test for the `WNT_ARM` genes along the crypt-villus axis, so
   the partner gradient can be read against ANTXR2's own.

### Task 5 -- Kong2023 feasibility check (gate before any Phase 3 scoping)

**Do the feasibility check first and report before proceeding.** The
differential-correlation analysis is only worth scoping if ANTXR2 survives
filtering.

1. Load Kong2023. Report per `cell_type` x condition
   (`Non_pathological` / `Neighbouring_inflamed` / `Inflamed`) x donor: cell
   counts, and whether ANTXR2 clears `min_perc_group` at the standard 0.7.
   Report what threshold *would* be needed if 0.7 fails. (For reference,
   ANTXR2 needed 0.6 for Elmentaite's full enterocyte donor set.)
2. Report how many donor groups per condition would survive the 100-cell
   floor. If any condition has <8 usable donors for a cell type, flag that
   cell type as underpowered for two-group testing in that condition, using
   the macrophage precedent (5 donors → failed replication) as the reference
   for what "underpowered" has meant in this project.
3. Check whether chemistry/batch is confounded with condition, as it was with
   donor in the Elmentaite working set.
4. **Only if feasible**, scope the differential test with this pre-specified
   directional prediction: since CMG2's Wnt function is injury-conditional
   (CMG2-KO baseline guts are normal), **ANTXR2's coupling to Wnt-arm partners
   and regeneration programs should be absent or weak in `Non_pathological`
   and appear in `Inflamed`.** Register this before running.
5. Note for scoping: Bracq et al. is mouse **colon** with DSS colitis, so
   **colonocyte** is the matched cell type for direct comparison; enterocyte
   is the extension.

## Expected outputs

1. **Scripts** in `scripts/`, one per task, each runnable independently and
   logging every parameter/threshold used. Follow existing naming
   (`coexpression_*.py`, `gene_panel_query.py`).
2. **A pre-registration block appended to `ANALYSIS_SUMMARY.md`'s session log
   BEFORE any results**, stating the interpretation rules for Tasks 1, 2 and 5
   as written above. This is the point of the exercise -- do not write it
   after seeing output.
3. **Results appended to the session log**, with each finding labeled
   **confirmatory** (pre-defined panel test) or **exploratory** (screening).
4. **Figures** under `/data/ANTXR2/figures/`, matching the existing convention
   that generated artifacts live under `/data/ANTXR2/` rather than in the
   git-tracked repo.
5. **Header amendments** if and only if a result changes the framing or a
   standing correction -- with the amendment noted in the log entry that
   prompted it.

## Do not

- Expand any panel mid-analysis because something looked interesting. Panels
  are pre-defined; new hypotheses go in a new prompt.
- Build a composite vulnerability/redundancy score. Standing correction.
- Treat a null correlation as evidence of co-presence.
- Interpret any positive correlation before Task 2 has run.
- Draw a biological conclusion from macrophage.
- Report a GSEA term as a tested claim.
