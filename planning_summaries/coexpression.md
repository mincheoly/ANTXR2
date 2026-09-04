## ANTXR2 coexpression analysis (macrophage / enterocyte / fibroblast trio)

## Context
Test run of the ANTXR2 coexpression pipeline on a curated subset containing
only macrophages, enterocytes, and fibroblasts, before scaling to the full
dataset.

Data location: `/data/ANTXR2/coexpression/elmentaite2021_trio_level3.h5ad`

## Environment
- Use the `memento-de` package (PyPI: `pip install memento-de`). Confirm the
  installed version and check its current API before writing code — do not
  assume function names from older memento tutorials, since `memento-de`'s
  interface differs from the earlier `scrna-parameter-estimation` package.
- Record package versions (memento-de, scanpy, anndata) in the analysis
  summary for reproducibility.
- This run uses point estimates only, no bootstrap/hypothesis testing (see
  step 3) — check whether `memento-de`'s estimation functions can be run
  standalone without triggering its full inference/testing pipeline, and use
  that lighter-weight path if available, since it should also make this run
  meaningfully faster.

## Instructions

1. **Read and validate the data**
   - Load the h5ad file.
   - Confirm `cell_type` and `capture_rate_pbmc` columns exist in `obs`, and
     that `ANTXR2` exists in `var_names`.
   - Report cell counts per `cell_type` — flag any cell type below the
     minimum cell count threshold (see step 2) before running memento on it.
   - Confirm `.X` contains raw/unnormalized counts (memento's estimators
     assume raw counts) — if not, locate raw counts (e.g. in `.layers` or
     `.raw`) and note where they were found.
   - Sanity-check `capture_rate_pbmc`: confirm it's populated for all cells
     in the trio, and report its range/distribution per cell type.

2. **Construct groups and filter genes**
   - Construct memento groups by `donor` × `cell_type` (one group per
     donor/cell-type combination). Estimates must exist for every
     donor×cell_type group individually — this donor-level granularity is
     needed before the cross-donor averaging in steps 3b and 5b below.
   - Enforce a minimum group size of **100 cells** per donor×cell_type
     group; flag and exclude any group below this threshold rather than
     silently running on it (report which donor/cell-type combinations were
     dropped and why, since this affects the denominator in the averaging
     steps).
   - Use the PBMC-calibrated per-cell capture rate estimates in
     `obs['capture_rate_pbmc']` as the capture-efficiency input, rather than
     a single global/default value. Check `memento-de`'s API for how it
     expects per-cell capture rate to be supplied (e.g. a column reference
     vs. a precomputed per-group value) and confirm the per-cell values are
     actually being used, not silently collapsed to a group-level constant.
   - Apply the default gene filter, but explicitly record what that default
     threshold is in the summary.
   - Log the number of genes retained per donor×cell_type group after
     filtering, and confirm ANTXR2 survives the filter in every group (if
     it's filtered out of a given group, that group can't be analyzed for
     ANTXR2 — note it and exclude it from the averaging in step 3b rather
     than silently skipping).

3. **Compute ANTXR2 vs. all-genes correlation estimates (per donor×cell_type group)**
   - Compute point-estimate correlations for ANTXR2 against every other
     filtered gene, separately for each donor×cell_type group.
   - No inference needed at this stage — point estimates only, no
     bootstrapping, no confidence intervals, no p-values/FDR correction.
   - Handle and log genes where the estimate fails or returns NaN in a
     given donor group (e.g. due to sparsity) rather than dropping them
     silently.

   **3b. Average across donors, per cell type**
   - For each cell type, average the ANTXR2 correlation estimate for each
     gene across that cell type's donor groups, producing one averaged
     correlation value per gene per cell type.
   - State and use an explicit averaging method — a plain (unweighted) mean
     across donors, unless there's reason to weight by each donor group's
     cell count; pick one and record it.
   - For genes missing/filtered out in some donor groups but not others,
     average only over the donors where the gene passed filtering, and
     record how many donors contributed to each gene's averaged value (a
     gene averaged over 2 donors is less reliable than one averaged over
     all donors for that cell type — surface this rather than hiding it).

4. **Select top 50 genes per cell type**
   - Rank genes by **absolute magnitude** of the donor-averaged correlation
     from step 3b (largest |mean correlation| first), independent of sign.
   - Exclude ANTXR2 itself from its own top-50 list.
   - Preserve the sign of the correlation in the output (positive vs.
     negative coexpression), even though ranking is magnitude-based.
   - Note ties-handling if relevant.

5. **Re-compute pairwise correlations among top genes (per donor×cell_type group)**
   - Take the **union** of the top-50 gene sets across all three cell types
     (≤150 unique genes, likely fewer due to overlap), plus ANTXR2.
   - Compute the pairwise correlation matrix over this unioned gene panel
     **separately for each donor×cell_type group** — do not pool cells
     across donors or across cell types for this step.
   - Point estimates only, consistent with step 3 (no bootstrap/inference).

   **5b. Average across donors, per cell type**
   - For each cell type, average the pairwise correlation matrices from
     step 5 across that cell type's donor groups (element-wise mean over
     the matrix), producing one averaged pairwise matrix per cell type.
   - Use the same averaging method (and same missing-donor handling) as
     step 3b, for consistency.

## Expected outputs

1. **Updated h5ad file** with memento correlation results stored in a
   documented, consistent location (e.g. `.uns['memento_correlations']`) —
   specify the structure, keeping donor-level and donor-averaged results
   distinguishable, e.g.:
   - `.uns['memento_correlations']['by_donor'][cell_type][donor]` → step 3
     ANTXR2-vs-all estimates and step 5 pairwise matrix for that group
   - `.uns['memento_correlations']['donor_averaged'][cell_type]` → step 3b
     averaged ANTXR2-vs-all estimates (with n_donors per gene) and step 5b
     averaged pairwise matrix
2. **Scripts** saved to `scripts/` in the repo, organized by step (e.g.
   `01_load_and_filter.py`, `02_compute_antxr2_correlations.py`,
   `03_average_across_donors.py`, `04_top_genes.py`,
   `05_pairwise_correlations.py`, `06_average_pairwise_across_donors.py`),
   each runnable independently given the prior step's output, with logging
   of parameters/thresholds used at each stage.
3. **Heatmap or clustered correlation plots**, one per cell type, for the
   donor-averaged union-gene pairwise correlation matrix, saved to a
   `figures/` or `plots/` folder — useful for a quick sanity check on this
   test run.
4. **`ANALYSIS_SUMMARY.md`** updated with:
   - Parameters used at each step (gene filter threshold, min cell count of
     100 per donor×cell_type group, capture rate source
     `capture_rate_pbmc`, ranking criterion — top 50 by magnitude of
     donor-averaged correlation, note that this run used point estimates
     only, averaging method used across donors)
   - Donor and cell counts per donor×cell_type group; any groups dropped
     for falling below the minimum cell count
   - Gene counts pre/post filter, per group
   - Top ANTXR2-correlated genes per cell type (donor-averaged, brief
     summary table), including how many donors contributed to each
   - Size of the unioned top-gene set and overlap across cell types
   - Package versions, for reproducibility
   - Runtime, to help estimate cost of scaling to the full dataset