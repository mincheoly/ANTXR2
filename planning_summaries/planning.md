# ANTXR2 project — conversation summary (handoff)

## Corrections to prior assumptions

- **DSS memory was wrong in mechanism.** Bracq et al. 2025 (EMBO Mol Med):
  CMG2-KO mice are normal at baseline, and CMG2 does *not* worsen acute DSS
  injury (identical DAI, weight loss, lipocalin-2). The defect is in
  **recovery** — fetal-like reversion happens normally, but the fetal-like →
  Lgr5+ transition fails (β-catenin nuclear translocation). Not
  "disproportionate injury," not "protective." No colitis-protective effect
  found in the literature.
- **Collagen genes are not a positive control.** Bürgi et al. 2017 showed
  collagen VI mRNA does *not* change in Antxr2−/− uteri — the mechanism is
  degradation, not transcription. ANTXR2→COL6 transcriptional coupling is
  already ruled out in that direction.
- **Pathway scores don't fix causality.** Aggregating targets reduces noise in
  Y; it does not close a backdoor. Conceded.
- **Steele atlas is likely not on CELLxGene Discover** — no collection ID
  found. Hosted at collections.cellatlas.io/skin-fibroblast + GitHub
  haniffalab/skin_fibroblast_atlas. Needs a separate download path.

## Where the framing landed

- ANTXR2's function is **post-transcriptional** (endocytosis, lysosomal
  degradation, MMP activation, β-catenin translocation). Gene-gene correlation
  therefore cannot speak to what the protein *does*.
- HFS is **coding recessive LoF**, so ANTXR2's own transcriptional regulation
  is not the disease axis either. (Possible exception: the AS common-variant
  GWAS signal is probably regulatory — verify the variant annotation.)
- **Central question instead: the tissue-restriction paradox.** Germline
  mutation in every cell, phenotype restricted to skin nodules, gingiva,
  joints, gut, uterus. Why those tissues? This is what atlas data is actually
  suited to, and needs no causal identification.

## Three testable sub-hypotheses (ranked)

1. **Paralog redundancy (strongest).** ANTXR1/TEM8 LoF → GAPO syndrome, also
   ECM accumulation, also fibroblast senescence, but a non-overlapping tissue
   spectrum. Do the paralogs' expression domains partition across cell types in
   a way that predicts which tissues break in which disease? Caveat: the two
   bind different collagen VI domains (triple-helical vs C5), and at least one
   paper disputes anthrax receptors as collagen VI receptors at all.
2. **Machinery completeness.** Where is the full clearance toolkit
   (MRC2/cathepsins/lysosomal) co-assembled vs. absent? Cell types with ANTXR2
   but no clearance machinery → candidate non-ECM roles. Aim this at the
   unexplained enteric smooth muscle / ICC signal.
3. **Substrate burden (weakest).** COL6 production vs clearance capacity ratio;
   mRNA is a poor proxy for protein flux.

Plus, more modestly: **which fibroblast subtype** — never asked; Steele's
9-subtype resolution makes it newly answerable.

## The memento tie-in

Means establish presence; **correlation resolves co-presence**. A nonzero mean
for both paralogs is compatible with "both in every cell" (redundancy holds) or
"two disjoint subpopulations" (no backup at the cell level). Only correlation
distinguishes these.

Key asymmetry: the shared-activation-state confound inflates *positive*
correlation but cannot manufacture mutual exclusivity. So **anti-correlation is
the strong result here** — it runs against the confound. Conversely,
ANTXR2–COL6 positive correlation is exactly what the confound produces, so it
stays weakest and needs a specificity check against a generic activation score.

## Methods decisions

- Use **standalone memento** (`memento-de`), not the CELLxGene precomputed
  backend (immature).
- **Mean expression, not fraction-positive** — the latter is confounded by
  per-cell RNA content (large mesenchymal cells vs small lymphocytes).
- **Don't build genome-wide modules.** With one seed gene: rank by memento
  correlation with FDR control, run GSEA-preranked for interpretation, compare
  gene sets across atlases via preservation/overlap. cNMF only if a broader
  program map is later needed.
- Phase 1 groupings must match Phase 2 strata exactly. `n_cells` on every row —
  rare types (ICC, neuroblast) will have reportable means but unusable
  correlations.

## Datasets

| Dataset | Role |
|---|---|
| Steele skin fibroblast atlas (2025) | Skin; 9 fibroblast subtypes + spatial + 23 diseases. Download path unresolved. |
| Gut Cell Atlas (Oliver 2024/25) | Gut; fine epithelial states (stem/TA vs enterocyte), UC/Crohn's samples |
| Cross-tissue Immune Atlas (Domínguez Conde 2022) | Immune; same donors across tissues → within-donor macrophage vs monocyte |
| Tabula Sapiens | **Repositioned:** not the skin/gut/immune source anymore. Kept for (a) tissues nothing else covers — esp. **uterus**, where the Col6a1 rescue cross was done, (b) unbiased scan capability (source of the ICC/SMC surprise), (c) technically consistent cross-tissue ranking |

## Deferred / discussion-only

- **Mouse DSS regeneration**: no scRNA-seq exists from CMG2-KO; Ayyaz
  GSE123516 is irradiation-based and WT-only. Cite, don't re-derive.
- **Causal identification**: Perturb-seq (Replogle K562 genome-wide, 9,866
  genes — check if ANTXR2 is included) is real do(X) but wrong cell type. MR
  via OneK1K is properly identified but PBMC-only, where ANTXR2 is low → weak
  instrument exactly where the data is. Ruled out: PC/GES (needs causal
  sufficiency), pseudotime/velocity ordering, front-door. **The gap itself is a
  finding** — it defines what new data would be needed.

## Files written

- `docs/phase1_celltype_expression.md` — 3-collection CELLxGene pull (predates
  the reframing; needs revision)
- `docs/phase1_skin_mean_expression.md` — Steele skin, current framing
- `data/raw/`, `data/celltype_expression/` created

## Next decisions pending

- Revise `phase1_celltype_expression.md` for the new framing (tissue-restriction
  question, ANTXR1 in panel, Tabula Sapiens repositioned)?
- Genome-wide means vs. panel-only
- Tabula Sapiens Smart-seq2 vs 10x: separate or flag?
- Verify the AS GWAS variant annotation at ANTXR2
