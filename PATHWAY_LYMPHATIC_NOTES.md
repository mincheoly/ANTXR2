# ANTXR2 -- cross-tissue co-expression pathways, endothelial annotation granularity, and lymphatic endothelium

Companion to `ANALYSIS_SUMMARY.md` (authoritative framing and standing
corrections), `HFS_ATLAS_NOTES.md` (COL6/ANTXR2 load arm) and
`IBD_WNT_NOTES.md` (colitis Wnt arm). Sessions 2026-10-03 to 2026-10-05.

Same convention as the other notes: **the "Corrections" section wins over
anything stated earlier in this file.** Every number below was re-derived
from the saved post-fix outputs listed at the end on 2026-10-05, not restated
from session memory.

Evidence classes used throughout, per `ANALYSIS_SUMMARY.md`: **means-level**
(pseudobulk presence / enrichment; strongest), **co-expression** (memento
correlation; moderate, needs a specificity control), **pathway enrichment**
(exploratory characterization only). No causal claims anywhere.

---

## 0. What the project holds onto (scope set 2026-10-03)

Two results are treated as concrete; everything else in the earlier broad
analyses is considered under-powered or uninterpretable at present:

1. **COL6 / ANTXR2 load across tissues** (`HFS_ATLAS_NOTES.md`): a
   descriptive ranking (gingiva, buccal mucosa, JIA synovium, skin, hard
   palate highest; gut the inverse configuration). The organ-involvement
   enrichment claim is withdrawn there and stays withdrawn.
2. **ANTXR2 opposite the Wnt-output gradient in intestinal epithelium** --
   means-level, crypt -> differentiated, **not** within-cell co-expression
   (that was retracted in `IBD_WNT_NOTES.md`). Holds in the Gut Cell Atlas
   (formal memento test), SCP259 colon, and Tabula Sapiens small / large
   intestine (Section 1). Not a general epithelial rule.

## 1. The intestinal Wnt-axis observation does not generalize across epithelia

`scripts/wnt_axis_cross_tissue.py`, Tabula Sapiens epithelial extracts,
donor-paired pseudobulk, progenitor -> differentiated pool per tissue.
ANTXR2 percentile = rank of its log2FC among genes expressed >= 0.1 CP10K.

| tissue | axis | donors | ANTXR2 log2FC (pct) | donors same sign | Wnt targets |
|---|---|---|---|---|---|
| Small intestine | stem/TA -> enterocytes | 5 | +0.72 (92nd) | 5/5 | crypt-side (OLFM4, ASCL2, LGR5, ZNRF3) |
| Large intestine | stem/TA -> colonocytes | 3 | +0.66 (92nd) | 3/3 | crypt-side (OLFM4, SOX9, MYC) |
| Lung airway | basal -> club/goblet/ciliated | 4 | +0.85 (86th) | 3/4 | basal markers down, see below |
| Tongue | basal -> squamous | 4 | -0.34 (49th) | 4/4 | basal-side |
| Trachea | basal -> ciliated/secretory | 5 | -0.19 (42nd) | 3/5 | basal-side |
| Salivary | basal -> duct | 4 | -0.11 (67th) | 2/4 | weakly basal-side |
| Mammary | basal -> luminal | 2 | +0.01 (51st) | 1/2 | basal-side |
| Prostate | basal -> luminal | 3 | +2.20 (98th) | 3/3 | **luminal-side** |

- Intestine replicates a third time. Elsewhere ANTXR2 is flat while Wnt
  targets are basal-side, and in prostate both rise together.
- Lung airway is a **club-cell** effect (ANTXR2 club > basal in 4/4 donors),
  not an inverse-Wnt relationship: the "Wnt" genes falling basal -> club are
  basal markers (LGR6, LEF1, APCDD1), while the general feedback readouts
  AXIN2 and ZNRF3 are equal or higher in club cells.
- Non-gut Wnt targets (AXIN2, NKD1, NOTUM, SP5) sit mostly below 0.1 CP10K, so
  their non-gut fold changes are near the detection floor; LGR6, LEF1, ZNRF3
  are the better-measured ones. 2-5 donors per tissue.

**Statement:** in intestinal epithelium, across three datasets, ANTXR2 rises
with differentiation (top ~8% of genes) as Wnt output falls; this is not a
general feature of epithelial differentiation.

## 2. Cross-tissue co-expression pathways (exploratory)

Input: the 100 per-population one-sample memento scans from the previous
session (`coexpr_out/`, `coexpr_comp/`; fibroblast 25, endothelial 26,
epithelial 16, immune 33; median 4 donors per population).

### 2.1 Why the raw scans cannot be ranked as-is

- Every population's genome-wide distribution is shifted positive (median r
  +0.04 to +0.19; 70-98% of genes positive) -- consistent with memento's
  default size-factor settings (shrinkage 0.5 / trim 0.1), which
  `antxr2_coexpr.py` uses; the project had earlier found these push
  correlations positive.
- In 98/100 populations a gene's raw mean predicts its correlation with ANTXR2
  (median Spearman +0.25), so unadjusted ranking enriches highly expressed
  gene sets regardless of ANTXR2.
- The one-sample bootstrap SE is within-group cell sampling; with a median of
  4 donors per population, "significant vs zero" is not a usable criterion.

**Ranking statistic used everywhere below:** the correlation minus a running
median (301 genes) of correlations at similar expression
(`coexpr_pathway_gsea.ranking`, "resid"). This removes the shift and the
trend (residual Spearman with expression -0.02).

### 2.2 Recurrence across tissues (GSEA prerank, Hallmark / Reactome / KEGG)

`coexpr_pathway_gsea.py` + `coexpr_pathway_aggregate.py` on the original
scans (`results/pathway_gsea_noGOBP/`). Dominant, recurrent, and **generic**:
ribosome / translation, oxidative phosphorylation, mitochondrial translation
and proteasome negative in every compartment (e.g. KEGG Ribosome negative at
FDR<0.05 in 20/25 fibroblast tissues), and phosphoinositide metabolism /
Rho-GTPase cycles weakly positive. See 2.3 for why these are not ANTXR2-specific.

**Fibroblast ECM is not recurrent across fibroblast tissues.** The project's
earlier best-supported positive correlation (ANTXR2 with ECM/collagen in
Elmentaite gut fibroblasts, 21 donors) does not reproduce across 25
fibroblast tissues: after expression adjustment, Collagen Formation median NES
+0.64 (0 tissues significantly positive, 4 negative) and Hallmark EMT median
NES -1.44 (8 tissues significantly negative, because ECM genes are the most
highly expressed fibroblast genes). Caveat: the Tabula Sapiens populations
have 2-5 donors each.

### 2.3 Specificity against matched anchor genes

`coexpr_anchor_select.py` (5 anchors per population, matched within the
population to ANTXR2 raw mean +/-20% and detection +/-2 pp) ->
`coexpr_point_scan.py` (memento point estimates, ANTXR2 and anchors through an
identical pipeline) -> GSEA -> `coexpr_anchor_compare.py compare --standardize`.
Test per term: ANTXR2 is the max (or min) of the 6 targets with probability
1/6 under exchangeability; 2 x min(one-sided binomial) across populations,
BH within compartment; NES standardized per target to remove ANTXR2's wider
NES spread in immune/epithelial populations.

| compartment | populations | terms FDR<0.05 (+/-) | specific to ANTXR2 | generic (anchors recover it) |
|---|---|---|---|---|
| endothelial | 26 | 57 (41/16) | ECM organization (+, 17/26), scavenger-receptor uptake (+, 16/26), interferon-alpha response (-, 15/26) | ribosome (5/26), OXPHOS (7/26) |
| fibroblast | 25 | 273 (49/224) | unfolded protein response (-, 18/25), Myc targets (-, 16/25), mTORC1 (-, 15/25) | ribosome (9/25, FDR 0.11); collagen formation / EMT not specific |
| immune | 33 | 61 (8/53) | translation (-, 16/33), RNA transport (-, 17/33), long-term potentiation (+, 19/33) | ribosome (11/33, FDR 0.21); 20 of the 61 terms are >30% replication-dependent histone genes -- one cycling-cell signal, not 20 findings |
| epithelial | 16 | 22 (17/5) | weak / mixed | |

The endothelial row is reinterpreted in Section 3: it reflects annotation
granularity (lymphatic cells inside generic "endothelial cell" labels) plus
fibroblast / macrophage contamination, not endothelial co-regulation.

## 3. Endothelial ECM / clearance signal = annotation granularity + contamination

Gene-group check (`endothelial_source_check.py`): mean expression-adjusted
correlation of a gene group with ANTXR2 vs with each of its 5 anchors; cells
are populations where ANTXR2 is the max of the 6 (chance ~ 1/6). Each
contaminant is read on **held-out** genes not used by the filter.

| gene group | original (n=26) | fibroblast cells removed (>=2 DCN/LUM/COL1A1/COL1A2) | blood-vessel EC only |
|---|---|---|---|
| EC scavenger receptors (STAB1/2, SCARF1, MRC1, LYVE1, CD36) | 11 (p=0.002) | 12/25 (p=3e-4) | **4/19 (p=0.39)** |
| basement membrane | 10 (p=0.007) | 8/25 (p=0.045) | 6/19 (p=0.08) |
| fibroblast, held out | 10 (p=0.007) | **6/25 (p=0.23)** | 3/19 (p=0.64) |
| lymphatic, held out (TFF3, RELN, FLT4, MMRN1, PKHD1L1) | 13 (p=1e-4) | 16/25 (p<1e-4) | 7/19 (p=0.03) |
| macrophage, held out (AIF1, CD68, FCER1G, PTPRC, LYZ) | 11/20 (p=1e-4) | 10/19 (p=4e-4) | **8/14 (p=7e-4)** |
| pericyte/SMC | 10/25 (p=0.005) | 5/22 (p=0.30) | 4/18 (p=0.35) |

- **Fibroblast part = doublets / ambient RNA.** Removing fibroblast-marker
  cells abolishes it; in several tissues ANTXR2 itself was concentrated in the
  removed cells (raw mean fat 0.35 -> 0.20, muscle 0.23 -> 0.12, ovary 0.47 ->
  0.34, bladder 0.94 -> 0.74).
- **Scavenger part = lymphatic identity.** It disappears in blood-vessel
  endothelium (lymphatic removal: labelled LEC, or PROX1+CCL21 >= 1 count,
  which recovers 93-100% of labelled LEC in most tissues; spleen and stomach
  excluded because the rule cannot separate them).
- **Macrophage contamination survives a >=2-count filter** in blood-vessel EC
  (ANTXR2 is expressed in myeloid cells).
- **Blood-vessel EC pathway level** (19 tissues, standardized anchor test):
  7 of 1,378 terms at FDR<0.05; ECM organization 3/19 and collagen formation
  at chance. Surviving: interferon-alpha response (-, 11/19), protein
  localization (-, 12/19), xenobiotic / ROS detoxification (-), OXPHOS (-),
  phosphatidylinositol signaling (+, likely residual macrophage).
  Interferon-alpha negative has persisted through every filter; anti-correlation
  is the stronger result class, but ANTXR2's own enrichment is significant in
  only a minority of tissues -- a lead, not a finding.

---

## 4. Lymphatic endothelium (main thread)

Context: HFS, especially the severe form, includes **intestinal
lymphangiectasia** with protein-losing enteropathy, usually attributed to
collagen accumulation compromising lymphatic structure indirectly. Bracq et
al. 2025 (EMBO Mol Med) note this explicitly. The question pursued here:
does ANTXR2 have a lymphatic-intrinsic context?

### 4.1 ANTXR2 is lymphatic-enriched within endothelium (means-level; strong)

`scripts/lymphatic_vs_blood_paired.py` -> `results/lymphatic_vs_blood_paired.csv`.
Donor-paired pseudobulk, labelled lymphatic (LEC) vs other endothelium (BEC)
from the same donor, >= 20 cells each.

| tissue | donors | ANTXR2 CP10K LEC vs BEC | log2FC | percentile |
|---|---|---|---|---|
| Buccal mucosa | 6 | 1.15 vs 0.03 | +3.87 | 99.5 |
| Bladder | 2 | 2.25 vs 0.12 | +3.70 | 99.5 |
| Quadriceps tendon | 2 | 3.63 vs 0.33 | +3.28 | 99.1 |
| Tongue | 3 | 1.30 vs 0.09 | +3.25 | 99.2 |
| Muscle | 2 | 1.56 vs 0.15 | +3.04 | 99.2 |
| Hard palate | 1 | 0.84 vs 0.06 | +3.02 | 98.9 |
| Small intestine | 2 | 1.42 vs 0.14 | +2.98 | 98.6 |
| Synovium (JIA) | 5 | 0.95 vs 0.10 | +2.69 | 98.5 |
| Gingiva | 3 | 0.74 vs 0.08 | +2.42 | 98.2 |
| Vasculature | 2 | 0.93 vs 0.15 | +2.22 | 98.0 |
| Achilles tendon | 4 | 2.10 vs 0.34 | +2.15 | 97.6 |
| Salivary gland | 4 | 0.58 vs 0.18 | +1.49 | 97.1 |
| Thymus | 2 | 0.20 vs 0.07 | +1.06 | 92.5 |
| Lung | 4 | 0.89 vs 0.41 | +0.97 | 89.5 |

**14 tissues; higher in LEC in 41/42 donors; 2-15x; 89.5th-99.5th percentile
of lymphatic-enriched genes.** PROX1, LYVE1, MRC1, CCL21, TFF3 sit at the
98-100th percentile in the same comparison (label check). Several tissues
contribute only 1-2 donors; the consistency across tissues is the evidence,
not any single tissue.

### 4.2 It is not a lymphatic *marker* (cell level, specificity)

`scripts/lymphatic_marker_strength.py` -> `results/lymphatic_marker_{separation,compartments}.csv`.

- **Cell-level separation is the weakest of the genes tested:** AUROC LEC vs
  BEC median 0.65 (range 0.55-0.85, best in bladder), against TFF3 0.94,
  MMRN1 0.93, CCL21 0.92, PROX1 0.92, RELN 0.81, FLT4 0.78, LYVE1 0.76, PDPN
  0.71. ANTXR2 is detected in 39% of LEC vs 9% of BEC; the rest of LEC are
  zero, which caps per-cell discrimination (partly a depth ceiling).
- **No specificity beyond endothelium:** in the same tissue, fibroblasts or
  immune cells usually match or exceed LEC (e.g. Achilles LEC 2.31 /
  fibroblast 7.56 / immune 4.06 CP10K; salivary 0.56 / 2.45 / 1.25). LEC is
  highest only in some oral sites and bladder.

**Statement:** within endothelium, ANTXR2 is concentrated in lymphatic cells;
it is not a lymphatic marker.

### 4.3 Gut lymphatic subtypes: ANTXR2 marks LEC6 (means-level)

Data: `data/gut.h5ad` (Elmentaite 2021 Gut Cell Atlas, CELLxGENE), author
annotation LEC1-LEC6 (~2,650 cells). Panel extracted with
`stream_panel_extract.py` (now takes a custom panel); donor-averaged CP10K in
`results/gut_lymphatic_pseudobulk_by_donor.csv` (donor x subtype groups with
>= 10 cells; fetal = first/second trimester, postnatal = everything else).

| subtype | stage | donors | ANTXR2 | where (cells) | identity |
|---|---|---|---|---|---|
| **LEC6 (ADAMTS4+)** | postnatal | 6 | **1.11** | gut-wall mucosa, mostly colon / rectum / cecum / appendix (373/522) | capillary-like (LYVE1 24.7, CCL21 131) |
| LEC3 (ADGRG3+) | postnatal | 8 | 0.39 | ileal mucosa (358/390) | capillary-like |
| LEC5 (CLDN11+) | postnatal | 3 | 0.29 | largest share mesenteric lymph node (65/136) | valve (FOXC2 1.52, GJA4 5.44, CLDN11 6.41) |
| LEC1 (ACKR4+) | postnatal | 5 | 0.07 | mesenteric lymph node (363/685) | sinus-like |
| LEC4 (STAB2+) / LEC2 (MADCAM1+) | fetal only | 11 / 2 | 0.10 / 0.06 | | |
| blood EC | postnatal | 21 | 0.06 | | |
| fibroblast (Stromal 1-3) | postnatal | 21 | 0.83 | | |

Subtype and stage are partly confounded (LEC2/LEC4 fetal-only; LEC3/5/6
postnatal-only), so within-stage comparisons are the clean ones.

**LEC6 characterization** (`scripts/lec6_characterization.py` ->
`lymphatic/LEC6_vs_otherLEC_paired.csv`, `LEC6_contamination_check.csv`):
LEC6 vs LEC1/LEC3/LEC5 within the same donor and chemistry -- only **3 paired
units (donors A32, A34, A38; all 10x 5' v2; one comparison pool 44 cells)**.

- ANTXR2 log2FC +2.73, up in 3/3, **99.7th percentile of 3,563 genes; 7th
  among genes up in every unit.**
- LEC6 program (top, all up in 3/3): MFAP4, FN1, CSF3, CXCL1, ADAMTS4, IGF1,
  ANTXR2, SLCO2B1, RGS16, RCAN1, MT1E, CCL21, LYVE1, MT2A, COL6A3, HTRA1,
  NFKBIZ, BMP2, PDPN, IL32, TNFAIP3; PIEZO2 also enriched.
  Reading: a gut-wall, capillary-type, matrix-producing (FN1, MFAP4, ADAMTS4,
  some COL6A3), NF-kB-active state.
- **Not fibroblast contamination:** cells with >= 2 fibroblast-marker counts
  0-2.2% in LEC6, same as other LEC; DCN <= 0.11, COL1A1 <= 0.04 CP10K.
  ANTXR2 enrichment holds in cells with zero fibroblast-marker counts
  (detection 28-31% LEC6 vs 2-10% other LEC; CP10K 0.81-1.10 vs 0.05-0.36).
- **Not a stress artefact:** FOS / JUN / HSPA1A are higher in LEC6 in one donor
  (A38) and lower in another (A34).

### 4.4 Partner availability in lymphatic EC (means-level)

From the same gut table (postnatal LEC subtypes vs fibroblasts):

- **Clearance arm absent:** MRC2 0.00-0.11 CP10K in every LEC subtype vs
  0.53-0.76 in fibroblasts; COL6A1/2/3 roughly 0.1-3.8 vs 5-11 in fibroblasts.
  By the project's asymmetric rule (partners absent -> that function cannot
  operate as described), **MRC2-dependent collagen VI clearance by LEC
  themselves is unlikely.** Caveat: MRC2 is one reported partner; ANTXR2 binds
  collagen VI and laminin directly.
- **Wnt arm present:** LRP6, CTNNB1, TCF7L2, FZD4, FZD6 at levels comparable to
  blood EC. Wnt *targets* are too sparse in LEC to test co-expression.

### 4.5 Within-lymphatic co-expression

Population: lymphatic cells pooled from every local source
(`scripts/lymphatic_extract.py`), 6,755 cells x 21,611 shared genes, **23
groups with >= 100 cells from 18 donors**:

- gut atlas: donor x subtype x chemistry (10 groups) -- the only source where
  subtype is fixed within every group;
- Tabula Sapiens endothelium, 10x 3' v3, stomach excluded: donor x tissue
  (5 groups);
- oral / synovium / tendon extracts: donor x site (8 groups).

Groups are never pooled across subtypes *within a group*; Tabula Sapiens and
site groups can still mix capillary / collecting / valve LEC (label is
tissue), which is why the gut-only results are reported separately. One gut
donor-id caveat checked: "A34 (417C)" (a known collision in the Extended+
release) maps to a single donor UUID in this file.

**(a) Memento native one-sample test** (`scripts/memento_onesample_scan.py`,
PR #83 `ht_2d_moments(treatment=None)`, GPU, 5,000 bootstraps).

- **Centering by shrinkage** (`lymphatic/memento_onesample/shrinkage_sweep_summary.csv`):
  median ANTXR2-gene r -0.027 / -0.010 / **+0.004** / +0.030 / +0.131 at
  shrinkage 0 / 0.1 / **0.25** / 0.5 / 1.0 (trim 0.1). 0.25 centers the
  distribution (52% positive) but leaves an expression trend (Spearman +0.31;
  low-expression quintile -0.024, high +0.043). Shrinkage alone could not
  center it until the placeholder bug in Section 5 was fixed.
- Test against zero: 172 genes FDR<0.05, but sorted by expression (56 of the
  negatives in the lowest quintile; 60 of the positives in the highest) --
  the known detection-floor bias and residual trend.
- **Test against the expression-matched expectation** (memento coefficient
  minus running median at that expression, divided by memento's SE):
  175 genes FDR<0.05 (58 +, 117 -); negatives spread across all quintiles.
- memento's SE is within-group cell sampling, not donor-to-donor. Adding
  **donor sign consistency** (binomial p<0.05 across donors): 59 genes;
  **and same sign within the subtype-fixed gut groups: 56 genes (25 +, 31 -)**
  -> `lymphatic/memento_onesample/antxr2_lymphatic_partners_56.csv`.

  - Positive (25): MMRN1, LYVE1, PHKA1, CCL21, HYCC1, SLC38A1, PRMT7, HEXB,
    MGAT5, CAVIN2, TCTN1, EFCC1, ATG16L2, ZNF282, TGFBR2, CD58, PLD3, SYNM,
    SNX2, BBS10, SEMA6A, EPB41L2, EXOC5, NPL, ERMARD.
  - Negative (31): ARHGEF28, DEGS2, OTC, ADM, IL33, STARD10, ZNF687, ARHGAP28,
    TSTD2, NPLOC4, CLPX, ZFYVE19, ZNF608, ECHDC2, AHNAK2, EEF1B2, SMPD1,
    TMEM179B, RAB11B-AS1, EEF1D, RBIS, ATP5MF, SEC22B, PPP1R11, CDKN1C, COG8,
    PIGO, EXOSC6, TTC3, BBLN, C17orf58.
  - Gene-level reading (annotation, not tested programs): capillary-LEC
    identity (MMRN1, LYVE1, CCL21 -- also positive within subtype-fixed gut
    groups); endolysosomal / autophagy / glycan processing (HEXB, PLD3, SNX2,
    ATG16L2, MGAT5, NPL, EXOC5, CAVIN2) -- consistent with an endocytic
    receptor, in cells that lack MRC2; primary cilium (TCTN1, BBS10; thin);
    negative ADM (adrenomedullin), IL33, CDKN1C (IL33 and DEGS2 are strongly
    LEC6-depleted, so partly subtype state); negative translation / ER-Golgi
    housekeeping (EEF1B2, EEF1D, RBIS, SEC22B, COG8, PIGO).
  - Strength: adjusted r roughly 0.1-0.4; **consistent but modest**.

**(b) Donor-level test** (`scripts/lymphatic_gene_test.py`): per group
expression-adjusted r, averaged within donor, t-test across donors (>= 8),
BH. 11,238 genes: **0 at FDR<0.05; 1 at FDR<0.1 (STARD10, -0.23, negative in
15/15 donors, FDR 0.098).** Extreme tail at chance (12 genes p<1e-3 vs 11.2
expected), but a skewed histogram (1,735 genes in the lowest p decile vs
1,124 expected): many weak donor-consistent effects, no strong individual
ones. Power: donor-level SD per gene ~0.2-0.4, so with 18 donors and an
11k-gene burden a gene needs |r| >~ 0.25 in every donor to survive.

**(c) Anchor-z per gene** (`scripts/lymphatic_gene_table.py`): ANTXR2-g
expression-adjusted r vs the same gene's r with 30 matched anchors
(`lymphatic/anchors_lymphatic.csv`), plus same sign in all three sources:
50 genes (68 at |z| >= 3 before the replication filter).

**Convergence (robustness to the control choice, not independent
replication -- same estimates, different null):** the memento-native 56 and
the anchor-z 50 share **11 genes, all same sign** (0.24 expected;
hypergeometric p = 4e-16):

- positive: **HEXB, SNX2** (endolysosomal), **TCTN1** (cilium), HYCC1, PRMT7, ZNF282
- negative: **ADM, IL33**, STARD10, CLPX, ZNF608

These are the most defensible within-lymphatic ANTXR2 partners at present.

**(d) Targeted gene sets vs the 30 anchors, per source**
(`lymphatic/targeted_genesets_vs_anchors.csv`; z of ANTXR2's mean adjusted r
over the set against the anchors; rank among 30 anchors in parentheses):

| set | pooled | gut (subtype fixed) | Tabula Sapiens (mixed) | sites (mixed) |
|---|---|---|---|---|
| capillary (CCL21, LYVE1, PDPN, TFF3) | +2.71 (30) | **+1.63 (29)** | +3.56 (30) | +1.70 (28) |
| valve (FOXC2, GJA4, CLDN11, ITGA9) | -2.22 (1) | **-1.65 (1)** | -0.27 (14) | -1.46 (2) |
| LEC6 program (top 30) | +2.31 (30) | +1.25 (25) | +2.75 (30) | +0.73 (21) |
| Wnt receptors (LRP5/6, FZD4/6, CTNNB1) | -0.23 | +0.34 | +0.12 | -0.90 |
| NF-kB (NFKBIZ, TNFAIP3, CXCL1, CSF3, NFKBIA) | -0.15 | -0.28 | -0.40 | -0.02 |

Capillary-up / valve-down holds **within** the subtype-fixed gut groups as
well as in the mixed sources, which amplify it. So ANTXR2 tracks
capillary-LEC identity both between and, more weakly, within subtypes (gut is
10 groups). No Wnt-receptor or NF-kB co-expression beyond anchors.

### 4.6 Global view (UMAPs, `scripts/lymphatic_umaps.py`)

Figures in `figures/lymphatic_umap/` (A-F). Embeddings use each atlas's own
integrated space (gut: shipped `X_pca`, sub-UMAPs recomputed on the subset;
Tabula Sapiens: shipped scVI compartment UMAP).

- **Why the endothelial signal looked like co-regulation (annotation
  granularity).** In Tabula Sapiens, spleen, fat, ovary, skin, liver, mammary,
  prostate and trachea label 100% of their endothelium generically ("endothelial
  cell") with zero lymphatic-labelled cells, yet each holds tens to hundreds of
  PROX1+ or CCL21+ cells inside that label (fat 592, ovary 533, skin 285;
  `B_ts_granularity_by_tissue.csv`; spleen's 1,085 are likely ambient CCL21).
  On the UMAP the PROX1+CCL21+ lymphatic island contains both labelled-lymphatic
  and generic-labelled cells, and ANTXR2 is enriched across it.
- **Other ANTXR2-high endothelium in Tabula Sapiens:** a PROX1+ but CCL21-
  group (3,998 cells; mostly spleen, fat, heart, muscle; labelled generic)
  with ANTXR2 detection 36%, equal to labelled LEC (32%) -- probably not
  lymphatic (no CCL21); identity unresolved (candidates: splenic sinusoidal,
  venous-valve endothelium). And lung capillary endothelium (8,685 of a
  10,380-cell cluster) at 21%, consistent with lung having the highest BEC
  ANTXR2 and the smallest LEC/BEC difference in 4.1. So "ANTXR2-high
  endothelium" is PROX1+ endothelium broadly, plus lung capillaries.
- **Gut atlas:** LEC form a separate PROX1-high island next to blood-vessel
  endothelium; ANTXR2+ LEC concentrate at one edge of it. Across the whole
  atlas ANTXR2 is highest in mesenchymal and neuronal/glial cells, moderate in
  epithelium and endothelium (figure A, labelled by atlas category).
- **ANTXR2 in gut blood-vessel EC is structured, not uniform** (figure G;
  donor pseudobulk, >= 20 cells per donor x subtype,
  `G_endothelial_subtype_antxr2_summary.csv`). Median CP10K: fetal cycling EC
  0.16 (10 donors), fetal "Mature venous EC" 0.14 (8), fetal venous capillary
  0.11 (12), fetal venous / arterial EC 0.08 (10 / 12); postnatal arterial
  capillary 0.06 (14), mature venous EC 0.03 (11), mature arterial EC 0.008 (8).
  So the higher-ANTXR2 blood-vessel cells are developing (fetal, cycling)
  and venous / capillary endothelium; postnatal large arteries are near zero.
  Stage and segment are partly confounded (venous capillary is fetal-only,
  arterial capillary postnatal-only in this atlas); within postnatal the order
  is capillary > venous > arterial. All of these sit well below postnatal LEC6
  (1.02) and LEC3 (0.38). The same capillary-high pattern appears in Tabula
  Sapiens (lung capillaries above).
- **Gut LEC subtypes form a continuum, not islands.** Developmental stage
  dominates one axis (fetal LEC4 / LEC2 vs postnatal); lymph-node LEC1 and the
  valve markers (ACKR4, CLDN11) sit at one pole, LEC3 (ileal) runs through the
  middle, LEC6 occupies the opposite pole together with ANTXR2, ADAMTS4,
  COL6A3, MFAP4, ADGRG3 and a tight CSF3+ patch. MRC2 is detected in 1% of LEC.
  MKI67+ cycling cells sit in the fetal tip.
- **LEC6 composition:** postnatal only; gut wall (mostly colon / rectum /
  cecum / appendix); 443/522 cells 10x 5' v2 (the paired LEC6 DE in 4.3 is
  within 5' v2, so chemistry is held fixed there).
- **LEC6 quality metrics** (`F_lec_qc_by_donor_subtype.csv`): genes and UMIs
  per cell are mid-range (fetal LEC4 is deeper yet ANTXR2-low, so the LEC6
  ANTXR2 enrichment is not a depth effect); doublet scores are, if anything,
  lower. One LEC6 donor (A26, 11 cells) has median 24% mitochondrial reads; it
  is not one of the three paired DE donors (A32 7.0%, A34 5.6%, A38 3.4%).

### 4.7 Synthesis and hypotheses (none tested)

Supported:
1. Within endothelium ANTXR2 is lymphatic-enriched (14 tissues, 41/42 donors).
2. In gut lymphatics it is a top marker of LEC6, a gut-wall, capillary-type,
   matrix-producing state, and it tracks capillary-LEC identity, also within a
   subtype.
3. LEC lack the clearance co-receptor MRC2 and make little collagen VI; they
   have the Wnt receptors.
4. Within LEC, ANTXR2's cell-to-cell variation has consistent but modest
   partners (11 robust genes; endolysosomal and capillary-identity flavour),
   no strong individual partner, no Wnt-receptor co-expression.

Hypotheses this motivates (all need knockout or functional data):
- **Cell-intrinsic Wnt role in intestinal lymphatics.** Canonical
  Wnt/beta-catenin signalling acts in lymphatic patterning and valve
  formation downstream of flow (Cha et al., Genes Dev 2016 -- *verify*);
  Bracq et al. place CMG2 upstream of beta-catenin nuclear translocation in
  gut epithelium; LEC have the receptors. Note ANTXR2 is valve-*depleted*
  (LEC5 low; valve genes negative), which argues against a valve-specific
  role.
- **Adhesion / anchoring of capillary lymphatics.** ANTXR2's integrin-like vWA
  domain binds laminin and collagen IV; initial lymphatics depend on matrix
  attachment to open under interstitial pressure. Loss of a matrix receptor on
  capillary-type LEC is a route to lymphangiectasia distinct from matrix
  accumulation. Speculative.
- **Endocytic / degradative role without MRC2** (endolysosomal partners above).
- **Collagen VI acts on lymphatics:** Molon et al. 2023 (J Pathol), cited in
  Bracq et al.: collagen VI promotes lymphangiogenesis and drainage in colitis
  recovery -- relevant to how accumulated collagen VI could affect lymphatics
  beyond structure. **ADM** (negative partner) is a lymphatic regulator via
  CALCRL signalling (*verify* primary source).

### 4.8 Data that would test this

- **No Antxr2-knockout transcriptomic or spatial dataset was found.** Bracq et
  al.'s Visium and scRNA-seq panels are re-analyses of **wild-type** data from
  Mayassi et al. 2024 (Nature; GEO GSE245316; DSS injury-recovery course);
  their own knockout data are RNAscope / immunofluorescence / histology
  (data availability: Zenodo, BioStudies; no knockout transcriptomics).
- The Mayassi wild-type data could test whether lymphatic / LEC6-type
  enrichment of Antxr2 holds in mouse (the species of the knockout phenotype).
- Kong2023 LECs (adult gut, inflamed / non-inflamed) for replication: not on
  this host; a CELLxGENE Census query for lymphatic cells only would avoid the
  11.6 GB Extended+ download.
- More donors is the lever for within-LEC co-expression (each doubling lowers
  the detectable effect ~1.4x); a pre-specified panel (~20 genes) is the other.

---

## 5. Corrections (these win over anything above or in session logs)

1. **memento placeholder correlations (library) and my point-estimate
   aggregation (pipeline).** `estimator._corr_from_cov` filled pairs whose
   variance estimate is <= 0 with 5.0 and clipped to exactly +/-1 -- the same
   issue documented in `ANALYSIS_SUMMARY.md` on 2026-09-04 for the older
   pipeline, still present in the 2026-10 release. `ht_2d_moments` skips
   |r| == 1 groups, so memento's own tests were never affected. But
   `coexpr_point_scan.py`, which cell-weight-averages per-group point
   estimates, included them until 2026-10-04. In the pooled lymphatic
   population one 282-cell group (gut A33, LEC1; ANTXR2 raw mean 0) pegged
   every pair at +1 and produced a +0.10 floor that no size-factor setting
   could remove.
   - Fixed in `coexpr_point_scan.py` (drops |r| >= 1 group estimates) and
     upstream in memento (PR #83, second commit: placeholder -> NaN; test
     results unchanged, CPU now matches GPU).
   - All 179 point scans re-run (`scripts/compare_fixed_scans.py`):
     lymphatic rankings changed materially (median Spearman old vs fixed 0.68
     pooled, 0.76 per source); immune 28/33 populations changed (rho
     0.66-0.89); fibroblast 0/25; endothelial and epithelial mostly stable
     (median 0.97).
   - Re-checked: cross-tissue enrichment (31 + 3 populations redone; the rest
     reused because rho >= 0.9 -- `REDONE_populations.csv`), endothelial
     gene-group checks, all lymphatic gene-level analyses, the targeted
     gene-set table. **No substantive conclusion changed**, except that the
     within-gut capillary-up / valve-down association is clearer than
     pre-fix numbers showed (Section 4.5d), and immune histone-heavy terms
     rose (one signal).
2. **First anchor comparison was anti-conservative.** It scored ANTXR2 as
   "top" when most extreme in the direction of its own sign (null probability
   ~1/3, tested against 1/6). Replaced by separate max / min counts, each
   1/6, two-sided (`coexpr_anchor_compare.py`).
3. **Wrong null for "co-expressed with anything".** Counting how many genes
   pass for ANTXR2 vs for each anchor treats anchors as null genes; they have
   real partners of their own, and "is ANTXR2 more connected than a random
   gene" is not a meaningful prior. Withdrawn; replaced by per-gene tests
   (4.5a-c).
4. **Endothelial ECM specificity** (first reported as ANTXR2-specific ECM /
   clearance co-regulation in endothelium) reflects annotation granularity
   (lymphatic cells inside generic endothelial labels; not contamination)
   plus fibroblast / macrophage contamination (Section 3).
5. **Lymphatic "mostly subtype mixing"** (stated 2026-10-04 from pre-fix scans)
   is too strong: the capillary / valve association also holds within
   subtype-fixed gut groups (4.5d).
6. **Donor count**: the lymphatic-vs-blood comparison is 41/42 donors, not
   "39/40" as first reported in session.
7. **memento one-sample test was unavailable in the current release** (#69
   removed the all-ones branch; an all-ones treatment crashed inside
   scikit-learn). Restored explicitly as `treatment=None` in PR #83; the
   original 100 scans (`antxr2_coexpr.py`) were run on a pre-#69 version and
   remain valid.

## 6. Reproducing (outputs under `C:\Data\ANTXR2_workspace`)

Environment (WSL): `~/.venvs/antxr2` (uv; pandas, scipy, gseapy, torch cu124,
memento editable from `~/scrna-parameter-estimation`, branch
`feature/one-sample-correlation-test`). Host has ~13 GB RAM: keep gseapy at
<= 5 workers without GO BP, run memento scans one population at a time;
`coexpr_point_scan.py` chunks pairs (memento's covariance estimator copies the
target column once per pair; unchunked one-vs-all scans peaked at 12 GB).

| step | script | output |
|---|---|---|
| Wnt axis | `wnt_axis_cross_tissue.py comp_cells results/wnt_axis` | `results/wnt_axis/` |
| GSEA on original scans | `coexpr_pathway_gsea.py --libs HALLMARK REACTOME KEGG`; `coexpr_pathway_aggregate.py` | `results/pathway_gsea_noGOBP/` |
| anchors | `coexpr_anchor_select.py` (`--pops` + `--npz-dirs`, or `--npz`, `--n-anchors`) | `results/anchors_*.csv`, `lymphatic/anchors_lymphatic.csv` |
| point scans | `coexpr_point_scan.py <npz> <out> --targets ...` (`--groups-out`, `--shrinkage`) | `results/point_scan/{sf_default,decontam1,decontam2,bloodEC}/`, `lymphatic/scan*/` |
| anchor GSEA | `coexpr_anchor_compare.py split`; `coexpr_pathway_gsea.py`; `coexpr_anchor_compare.py compare --standardize` / `single` | `results/pathway_gsea_{anchors,bloodEC}/` |
| endothelial filters | `endothelial_decontam.py comp_cells <out> [--min-counts 2] [--lymphatic-min 1 --macrophage-min 2 --exclude ts_Spleen ts_Stomach]` | regenerable, not kept |
| endothelial check | `endothelial_source_check.py <out.csv> label=<scan_dir> ...` | `results/endothelial_source_check*.csv` |
| LEC vs BEC | `lymphatic_vs_blood_paired.py comp_cells <out.csv>` | `results/lymphatic_vs_blood_paired.csv` |
| marker strength | `lymphatic_marker_strength.py <ws> <prefix>` | `results/lymphatic_marker_*.csv` |
| gut panel | `stream_panel_extract.extract(..., panel=, obs_cols=)` | regenerable (`results/gut_lymphatic_pseudobulk_by_donor.csv` kept) |
| LEC extracts | `lymphatic_extract.py <ws>` (~2 min) | `lymphatic/*.npz`, regenerable, not kept |
| LEC6 | `lec6_characterization.py lymphatic` | `lymphatic/LEC6_*.csv` |
| memento native | `memento_onesample_scan.py <npz> <out> --shrinkage 0.25 --backend gpu` | `lymphatic/memento_onesample/` |
| gene tables | `lymphatic_gene_table.py`, `lymphatic_gene_test.py` | `lymphatic/antxr2_lymphatic_gene_*.csv` |
| fix audit | `compare_fixed_scans.py` | `results/fixed_vs_old_rankings.csv` |
| UMAPs | `lymphatic_umaps.py <ws> figures/lymphatic_umap` (~3 min; panel cache regenerable, not kept) | `figures/lymphatic_umap/` |
