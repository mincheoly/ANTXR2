# Pilot: cell-type-specific ANTXR2 gene-program co-expression (gut)

Observational co-expression only: which predefined transcriptional programs
ANTXR2 co-varies with in gut fibroblasts vs enterocytes, using ONE fixed weight
vector per program shared by both cell types. Not causal, not regulatory.
Code: `scripts/program_*.py`; data and outputs: `/data/antxr2/` (this host).
The project's standing corrections in `ANALYSIS_SUMMARY.md` apply throughout
(anchor floor for positive correlations; anti-correlation is the stronger
result class; exploratory vs confirmatory labelling).

## Decisions (agreed 2026-10-29/30)

| # | Decision | Choice |
|---|---|---|
| 1 | Dataset | Elmentaite 2021, Gut Cell Atlas on CELLxGENE: "Extended+" (`757945c8-…`, 18,370 genes) for postnatal, "Healthy reference" (`b801706d-…`, 35,461 genes) for fetal |
| 2 | Cells | **Revised 2026-10-30** (adult SI has 2 usable donors, see below): two separate donor sets, small intestine only, both cell types from the same slice. Primary: **postnatal SI** (4–74 y; 12 shared donors). Contrast: **fetal SI** (6–20 wk; 11 shared donors) |
| 3 | Weight scale | Genes standardised within each cell type; v = PC1 of the equal-weight average of the two cell types' memento-corrected correlation matrices; ρ_c = vᵀr_Pg,c / √(vᵀR_PP,c v) |
| 4 | Centering | Pooled within-group moments over donor × `level_3_annot` groups; donor × broad type as sensitivity |
| 5 | Signed sets | Published signs in the equal-weight baseline; published weights (PROGENy, CytoSig) as a third scheme. CollecTRI's published weight is its sign, so no separate third scheme |
| 6 | Membership | Gene must clear raw mean ≥ 0.07 in BOTH cell types; ≥ 10 genes and ≥ 50% of measured members; one-sided failures kept and flagged `not_comparable` |
| 7 | Null | Primary: 50 anchor genes matched to ANTXR2 on raw mean + detection within each cell type (anchor percentile per program × cell type). Secondary: within-group ANTXR2 permutation. Donor bootstrap for stability and Δρ (fib − ent) |
| 8 | memento | q = 0.15, shrinkage 0, trim_percent 0.5, ≥ 100 cells per group; random-pair median-correlation check before the main run |

Smaller defaults: CytoSig = top 100 |weight| genes per cytokine; PC1 sign set
so it correlates positively with the signed equal-weight score; PROGENy top 500.

## Estimator

memento's hypergeometric model corrects only the diagonal of the covariance
(`−(1−q)·mean(X/sf²)`); off-diagonals are plain sample covariances of `X/sf`.
So for a fixed weight vector a (on the `X/sf` scale; a = v / σ_c for the
standardised score) and a target gene g ∉ P, within each group:

    Var_c(s)   = sample var(Σ a_h X_h/sf) − (1−q) Σ a_h² mean(X_h/sf²)
    Cov_c(s,g) = sample cov(Σ a_h X_h/sf, X_g/sf)

Pooling these within-group moments over donor × subtype groups is exactly
"centre on donor × subtype, then correlate". Everything reduces to per-group
sufficient statistics, so a donor bootstrap is a reweighting of groups and no
genome-scale matrix is ever formed. Checked against memento's own
`compute_2d_moments` on single-gene programs.

## Gene-set library (`program_genesets.py`)

Hallmark 50, PROGENy 14 (top 500, signed weights), CollecTRI 1,185 regulons
(signed), CytoSig 43 (top 100), Matrisome 8 Naba categories (umbrella and
tumour-derived NABA sets excluded). ANTXR2 removed before filtering.

## Data finding (2026-10-30): the adult small-intestine slice is too small

Healthy reference, ≥ 100 fibroblasts (`level_2_annot == Fibroblast`) AND
≥ 100 enterocytes (`level_3_annot == Enterocyte`) per donor:

| slice | donors |
|---|---:|
| Elmentaite2021, adult SI | 2 |
| any study, adult SI | 2 |
| Elmentaite2021, pediatric (4–12 y) SI | 6 |
| Elmentaite2021, fetal SI | 11 |
| Elmentaite2021, adult, colon fibroblasts + SI enterocytes | 3 |

The Phase 2 trio's "35 shared donors" were mostly fetal/pediatric and colon.
Adult SI samples in this atlas are largely mucosal/epithelial preparations,
so adult SI fibroblasts are nearly absent. Kong2023 does not help: its SI
cells in Extended+ are epithelium only (no `Fibroblast` label in any
lamina-propria fraction).

Extended+ carries Elmentaite donors the Healthy reference lacks (adolescents
D150–D156; D154/D155 at 9–12 y). Usable donors at the real group rule
(≥ 100 cells per donor × `level_3_annot`, both cell types):

| slice | Healthy reference | Extended+ |
|---|---:|---:|
| adult 18+ | 2 | 2 |
| postnatal 4–74 | 8 | **12** (14 fib / 17 ent donors) |
| fetal 6–20 wk | **11** (15 fib donors) | 11 |

Library coverage: 11,085 of 11,778 library genes in the Healthy reference,
10,481 in Extended+.

**Resolution (user, 2026-10-30):** postnatal SI from Extended+ as the primary
arm, fetal SI from the Healthy reference as a separate, developmentally
distinct donor set. Consequences to carry into every reading of the results:
the postnatal arm is mostly children (8 pediatric, 2 adolescent, 2 adult), and
the fetal arm is a different biological state (immature enterocytes,
developmentally high Wnt), not a replicate. Agreement between arms is
suggestive; disagreement is not a failure to replicate.

## Running it (reproducing every result in this note on a new host)

**Environment.** Python 3.10; numpy 2.2.6, scipy 1.15.3, pandas 2.3.3,
anndata 0.11.4, scanpy 1.11.5, h5py 3.16.0, pyarrow 25.0.1, decoupler 2.1.4,
matplotlib 3.10.9, seaborn 0.13.2, statsmodels 0.15.0, scikit-learn 1.7.2.
memento: `pip install -e` of github.com/yelabucsf/scrna-parameter-estimation
at commit `90566fe` (the original host's checkout also has uncommitted edits
to `ht_2d_moments`' one-sample covariate handling, which this pipeline never
calls; it uses only `setup_memento` for size factors and, under `--validate`,
`compute_1d/2d_moments`). Host used: 8 cores, 61 GB RAM; peak use was well
under that.

**Data root.** Everything lives under `ANTXR2_ROOT` (default `/data/antxr2`;
export it if the new host differs). About 70 GB in all: raw h5ads 51 GB,
extracted cells 9 GB, landscape output 9 GB.

| path under root | what | how to get it |
|---|---|---|
| `raw/gut_cell_atlas/healthy_reference.h5ad` | GCA Healthy reference (9.4 GB) | `https://datasets.cellxgene.cziscience.com/b801706d-6345-4af2-b3fe-8e30e0e6150b.h5ad` |
| `raw/gut_cell_atlas/extended_plus.h5ad` | GCA Extended+ (11.6 GB) | `…/757945c8-a916-431d-aceb-1afbc80a7c55.h5ad` |
| `raw/ibdverse/ibdverse.h5ad` | IBDverse (32.3 GB) | `…/e3460364-c300-4616-9923-3a0224894998.h5ad` |
| `annotation/gencode.v44.basic.annotation.gtf.gz` | gene spans | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz` |
| `genesets/programs_long.csv` | program library | **copy from the old host** (md5 `875a18c2e9107b3e66f947f8456d9832`); regenerating pulls OmniPath's current Hallmark/PROGENy/CollecTRI, which drift |
| `genesets/raw_cytosig/signature.centroid` | CytoSig (needed only to regenerate the library) | `https://raw.githubusercontent.com/data2intelligence/CytoSig/master/CytoSig/signature.centroid` |
| `annotation/gene_span.csv` | spans | `gene_span.py` (md5 `178d0e6fe8c5450e24c0c04104af70cf`) |

`program_extract.py` refuses an h5ad whose byte size, `uns["citation"]`
version or nnz differs from `program_config.SOURCES`.

**Commands** (from `scripts/`; `py` = the environment's python). Run
directories are `output/<contrast>_<grouping>[_mat<N>][_covNN]/`.

```
py program_genesets.py                         # only if not copying programs_long.csv
py gene_span.py
# extraction: one pass per source over its h5ad -> program_coexpr/cells/
py program_extract.py --source extended        # postnatal SI (discovery)
py program_extract.py --source healthy         # fetal SI
py program_extract.py --source ibdverse        # replication, fib + enterocyte
py program_extract.py --source ibdverse_crypt  # state contrast, crypt
py program_extract.py --source ibdverse_landscape   # 27 per-population npz

# iterations 1-3 + discovery ("Iteration 2/3")
py program_scores.py --contrast fib_vs_ent --maturation-bins 3 --min-coverage 0.25 --validate
py program_report.py fib_vs_ent_fine_mat3_cov25
py program_scores.py --contrast fib_vs_smc_fetal --min-coverage 0.25
py program_report.py fib_vs_smc_fetal_fine_cov25
# replication ("Replication result")
py program_scores.py --contrast fib_vs_ent_ibdverse --maturation-bins 3 --min-coverage 0.25 --validate
py program_replicate.py fib_vs_ent_fine_mat3_cov25 fib_vs_ent_ibdverse_fine_mat3_cov25
py program_scores.py --contrast fib_vs_ent_ibdverse_inflamed --maturation-bins 3 --min-coverage 0.25
# state contrasts ("State contrast result")
for r in uninflamed normal; do
  py program_scores.py --contrast crypt_infl_vs_$r --min-coverage 0.25 --validate
  py program_scores.py --contrast enterocyte_infl_vs_$r --maturation-bins 3 --min-coverage 0.25 --validate
done
py program_state.py crypt_infl_vs_uninflamed_fine_cov25 crypt_infl_vs_normal_fine_cov25 \
    enterocyte_infl_vs_uninflamed_fine_mat3_cov25 enterocyte_infl_vs_normal_fine_mat3_cov25
# landscape ("Landscape, first full result")
py program_landscape.py all --validate         # prep -> weights -> score
py program_landscape_report.py                 # -> output/landscape/report/
py program_landscape_fib_subtype.py            # "Fibroblast interferon lean"
```

Each `program_scores.py` run takes 20–40 s; the landscape score stage took
~10 min on the old host. Seeds are fixed (`program_config.SEED`), so bootstrap CIs reproduce
exactly. Earlier `postnatal_si_*` / `fetal_si_*` directories on the old host
come from the pre-`--contrast` interface (Iteration 1) and are superseded by
the `fib_vs_*` runs above.

## Iteration 1 results (2026-10-30). Exploratory, not tested claims

**Estimator check.** On 100+ unclipped group × gene pairs per arm, per-group
single-gene correlations match memento's `compute_2d_moments` to < 1e-14.

**Power and floor.**

| arm | cell type | cells | groups | donors | ANTXR2 raw mean (det) |
|---|---|---:|---:|---:|---|
| postnatal | fibroblast | 4,713 | 21 | 14 | 0.279 (23%) |
| postnatal | enterocyte | 22,100 | 17 | 17 | 0.098 (8.5%) |
| fetal | fibroblast | 39,592 | 48 | 15 | 0.187 (16%) |
| fetal | enterocyte | 12,798 | 11 | 11 | **0.039 (3.7%), below the 0.07 floor** |

Fetal-enterocyte associations, and therefore every fetal Δρ, are below the
house floor and are not interpreted. The fetal arm is fibroblast-only.

**Decision-6 coverage rule is binding.** At ≥ 50% coverage in both cell
types, only 34 (postnatal) / 86 (fetal) programs are comparable; all eight
Matrisome categories, EMT, and most CytoSig sets fail (typically ~60% of
genes pass in each cell type but ~45% pass in both). The `_cov25` runs
(sensitivity only) give 249 / 335.

**Random-pair check (Decision 8).** Median pooled correlation: postnatal
fibroblast +0.007, fetal fibroblast +0.021, fetal enterocyte +0.008, but
**postnatal enterocyte +0.060**: a nonspecific positive covariation floor in
that population.

**Main finding: in postnatal enterocytes ANTXR2's excess coupling is global,
not program-specific.** Across all 110 (391 at cov25) enterocyte units,
ANTXR2's median anchor-z is +1.1 (+1.3), and 36% of units are
anchor-extreme (vs ~5% expected). ANTXR2 sits above its matched anchors on
almost every program (HNF4A, CDX2, heme metabolism, mitotic spindle …),
i.e. it is more tightly tied than matched genes to the dominant shared axis
of enterocyte variation. Per-program anchor percentiles therefore mostly
rediscover that axis. Ranked *relative to ANTXR2's own median z*, the top
enterocyte programs are CytoSig TRAIL, EGF, FGF2, Hallmark TNFα/NF-κB,
IL2/STAT5 (ρ 0.27–0.34). Their PC1 weights are dominated by
immediate-early/stress genes (MCL1, ZFP36, EGR1, DUSP1, JUNB, KLF6,
NFKBIA, CDKN1A), the dissociation-stress signature, so this must be
controlled before it is read as biology. The bottom (most negative relative to
ANTXR2's median) are HDAC1 targets, CytoSig GM-CSF / IL-4, SNAI2/ZEB1 regulons,
PROGENy Estrogen (ρ −0.26 to −0.29). Under the standing correction these are
the stronger result class (shared activation cannot manufacture them). They
sit at the 0th anchor percentile (z −1.9 to −2.8), but their anchor means are
themselves negative (≈ −0.13 to −0.16): most genes expressed like ANTXR2 in
enterocytes anti-correlate with these programs, and ANTXR2 does so more.

**Fibroblasts: no robust program-specific association.**
- Postnatal: 4–6% of units anchor-extreme, i.e. chance level. Tops relative
  to ANTXR2's median: NR4A1, FOXO3, MTF1, STAT2 regulons, LEF1 (ρ 0.16–0.29);
  BACH1 regulon negative (−0.15, 0th pct). PROGENy TGFb 0.19 (96th pct).
- Fetal: ANTXR2 is broadly *above* anchors (median z +0.9) but only 2%
  extreme. **ECM reproduces the earlier trio result in magnitude**:
  Matrisome categories, EMT, PROGENy WNT at ρ 0.29–0.36, 84–96th anchor
  percentile, but z relative to ANTXR2's own median is ≈ 0–0.6, so ECM
  does not stand out from ANTXR2's general coupling in fetal fibroblasts.
  Postnatal fibroblast ECM is weaker (ρ 0.06–0.14, 52–90th pct).
- Fibroblast program profiles across the two independent donor sets
  (434 shared programs): Spearman 0.36 on ρ, **0.14 on anchor-z**. The
  ANTXR2-specific component barely reproduces between postnatal and fetal.

**Fibroblast vs enterocyte (postnatal, the only interpretable Δ).** The
programs with Δρ CI excluding 0 and anchor-extreme Δ are the enterocyte
stress/cytokine set above (fibroblast ρ ≈ 0, enterocyte ρ ≈ 0.25–0.33),
plus GM-CSF (enterocyte −0.29). Given the global enterocyte shift and the
stress-gene loadings, this is not yet evidence of a program ANTXR2 joins in
enterocytes but not fibroblasts.

**Sensitivity.**
- Centering, fine vs broad: postnatal fibroblast ρ r = 0.81 (median
  |Δ| 0.03); enterocyte identical (one subtype).
- Weighting: ρ(PC1) vs ρ(equal) r = 0.59–0.81 by cell type and arm;
  vs published weights r = 0.70–0.86. Program PC1 explains a median 10–12%
  of averaged-correlation variance (IQR ~8–15%), and cos(PC1, equal) ≈ 0.6:
  most programs are weakly coherent in these cells, so PC1 weights partly
  track whichever genes dominate local covariation (the sketch's own caveat).

**Not done in iteration 1 / next steps.**
1. Dissociation-stress control program (van den Brink 2017 list, from the
   source, not from memory) scored alongside, and ANTXR2–program
   associations recomputed conditional on it. This decides whether the
   enterocyte result survives.
2. Characterise the enterocyte global axis directly (top PCs of the pooled
   enterocyte correlation matrix; ANTXR2's loading on them vs anchors'),
   and re-rank programs after removing it.
3. Split-sample weights (learn v on half the donors, score on the other
   half) to remove weight/association reuse.
4. Decide the coverage rule (50% vs 25%) and whether CollecTRI regulons
   under ~15 genes belong in the library.

### Follow-up (2026-10-30): the enterocyte "global axis" is maturation (+ a quality axis)

Enterocytes here are `level_3_annot == "Enterocyte"` only (stem, TA, BEST4
excluded), but that label still spans immature → mature absorptive cells.
Eigen-decomposition of the pooled postnatal enterocyte correlation matrix
(5,954 floor genes); gene–PC correlation = √λ·v:

| PC | share | content | ANTXR2 r | anchors (mean ± sd) | ANTXR2 anchor pct |
|---|---:|---|---:|---|---:|
| 1 | 9.3% | cytoplasmic housekeeping (+) vs MT-* / MALAT1 (−): a cell-quality / mito-fraction axis | 0.19 | 0.155 ± 0.069 | 66 |
| 2 | 3.6% | mature absorptive (APOB, CDHR2, SLC15A1, MTTP, SLC5A1, MME) (+) vs ribosomal (−) | **0.27** | 0.095 ± 0.096 | **98** |
| 3 | 1.6% | OLFM4 / CCND2 / RPL (immature, +) vs APOC3 / GUCA2A / MHC-II | −0.03 | 0.035 ± 0.068 | 16 |

So (i) the positive random-pair floor (+0.06) and the anchors' shared positive
coupling come from PC1, a quality axis ANTXR2 shares no more than matched
genes; (ii) **ANTXR2's excess is on the maturation axis** (98th anchor pct),
consistent with ANTXR2 rising with differentiation (Task 4.2 means result).
ANTXR2 vs markers: APOB +0.25, ADA +0.16, ALPI/LCT/SI +0.12–0.14, OLFM4 −0.01.
The "stress" genes driving the top enterocyte programs load on the mature end
too (CDKN1A +0.43, KLF6 +0.34, ZFP36 +0.33, KLF4 +0.26 on PC2; EGR1/JUNB/FOS
weaker), so those programs are at least partly maturation, not (only)
dissociation stress. Next: condition the enterocyte analysis on maturation
(and on PC1) before reading any enterocyte program association.

## Iteration 2 (2026-10-30): the difference only, with enterocytes conditioned on maturation

**Design changes (user).** (a) Only the fibroblast − enterocyte difference
matters; one-sample machinery (within-cell-type permutation, per-cell-type
anchor percentiles, single-cell-type program units) is removed, and Decision 7
becomes: *the null for the difference is 50 genes matched to ANTXR2 on raw
mean + detection in both cell types*. Headline statistic: **excess = Δρ_ANTXR2
− mean Δρ_anchors** and **excess z = excess / SD(Δρ_anchors)**; one donor
resample per bootstrap draw serves both cell types and carries ANTXR2 and the
anchors together. (b) Enterocytes are grouped donor × maturation bin (3
within-donor quantiles of a fixed score: mean within-donor z of log1p
APOB, APOA4, ALPI, SI, LCT, MTTP, SLC5A1 minus z of OLFM4;
`program_config.MATURATION_MARKERS`). Run directories are now
`<slice>_<grouping>_mat<N>[_covNN]`; `--maturation-bins 0` is the unconditioned
baseline. Fetal is dropped (ANTXR2 below the floor in fetal enterocytes, so no
interpretable difference).

**Maturation bins behave as intended.** ANTXR2 raw mean 0.067 → 0.085 →
0.136 (detection 6.0 → 7.5 → 11.6%) from least to most mature bin. Cost: 13
of 17 enterocyte donors keep all bins at ≥ 100 cells; enterocyte random-pair
median r drops 0.060 → 0.048.

**Results (postnatal, PC1 weights; cov25 shown, cov50 agrees on the 32 shared programs).**
- Conditioning lowers enterocyte ρ broadly (median 0.21 → 0.15) but leaves
  ANTXR2's **excess** nearly unchanged (r = 0.96 between mat0 and mat3, both
  coverages): the jointly matched anchors already absorbed the shared
  maturation/quality axes, so the difference statistic was largely immune.
- **Screen-level calibration: at chance.** 241 programs; |excess z| > 2 in 13
  vs 11.0 expected under N(0,1); max |z| 3.3 vs ~3.0 expected for the largest
  of 241. The bootstrap CI of the excess excludes 0 for 41, but that CI omits
  gene-to-gene variation (the anchor spread), so the anchor z is the right
  yardstick.
- **Enterocyte-leaning (ANTXR2 couples to these in enterocytes beyond matched
  genes):** CytoSig EGF (z −3.2), TRAIL (−2.9), FGF2, CD40L, VEGFA; CollecTRI
  JUND, TCF4, SIRT1, PML, MEF2D; Hallmark TNFα/NF-κB. Fibroblast ρ ≈ 0,
  enterocyte ρ ≈ 0.2–0.3. **One signal, not ten:** their PC1 loadings share an
  immediate-early/stress core (CDKN1A in 8 of 9, JUN, JUNB, ZFP36, MAFF,
  CCND1, KLF6, DUSP1, PPP1R15A, NR4A1). It persists within maturation bins, so
  it is not the maturation axis; whether it is dissociation stress is untested.
- **Fibroblast-leaning:** CollecTRI MTF1 (z +3.3; MT2A, MT1A, SOD1, SLC31A1:
  metal response), JDP2, TP73, MXD1; Hallmark heme metabolism, bile-acid
  metabolism, UV response up; CytoSig IL6, IL27. Individually at |z| 1.8–2.1
  except MTF1, i.e. within the chance band.
- **ECM/collagen programs were NOT tested** (corrected 2026-10-30; an
  earlier version of this line said "no notable excess", which was wrong).
  All 8 Matrisome categories, Hallmark EMT and PROGENy TGFb fail the
  both-cell-types membership rule: e.g. NABA_COLLAGENS 21 genes above floor
  in fibroblasts, 2 in enterocytes, 1 in both; NABA_CORE_MATRISOME 120 / 27 /
  17; EMT 156 / 51 / 44. The TGF-β programs that do pass keep no canonical
  ECM genes. Enterocyte estimation is not the weak side in general (median
  bootstrap CI width of ρ: enterocyte 0.135, fibroblast 0.245). The ECM
  contrast is a presence/absence difference (means), which the fixed
  shared-gene design cannot express as a co-expression difference.
  Fibroblast-only values (iteration 1): postnatal ρ 0.06–0.14 (52–90th anchor
  pct), fetal 0.29–0.36 (84–96th).

**Reading.** Across the library, ANTXR2's fibroblast-vs-enterocyte
program differences are about what a matched gene would show. The one
coherent departure is an immediate-early/stress module that ANTXR2 tracks in
enterocytes (within maturation state) but not fibroblasts. Next: score a
published dissociation-stress signature (van den Brink 2017, from source) as
its own program and as a covariate; if the enterocyte excess survives
conditioning on it, it is a candidate biological difference.

## Iteration 3 (2026-10-30): ECM via an ECM-expressing contrast; estimator fragility

**Why.** ECM programs could not enter fibroblast vs enterocyte (genes not
expressed in enterocytes; see the corrected line above). User chose fetal
fibroblasts vs fetal smooth muscle (SMC_PPLP2; 15 shared donors, Healthy
reference, paired). Code generalised: `program_config.CONTRASTS` names two
populations per contrast (`fib_vs_ent`, `fib_vs_smc_fetal`); extraction now
keeps all mesenchymal lineages (`cells_<source>.npz`); runs are
`OUT_DIR/<contrast>_<grouping>[_mat<N>][_gm][_covNN]`. The generalised code
reproduces the earlier fib-vs-ent numbers exactly.

**ANTXR2 is the weakly measured gene, not the programs.** New diagnostic,
reliability = noise-corrected / raw within-group variance: program scores
0.19–0.75 (median ~0.5–0.6); ANTXR2 only 0.17 (postnatal fib), 0.18
(enterocyte), 0.20 (fetal fib), 0.35 (SMC). With the pooled estimator
(ratio of summed moments) a small corrected ANTXR2 variance in some groups
shrinks the pooled denominator: fetal SMC gave impossible values (SOX2 0.99,
GTF3A 0.99, MECP2 −1.02) while the same programs' per-group SMC ρ had
medians ~0.35 (`group_rho.csv`).

**Second estimator, `--estimator groupmean`:** cell-weighted mean of per-group
ρ clipped to [−1, 1] (memento's own get_2d_moments convention; the house
pipeline's choice). Bounded by construction (max |ρ| 0.55 in SMC).
- fib vs ent (mat3, cov25): stable across estimators (excess r = 0.88);
  still at chance (|z| > 2: 10 vs 11.0 expected). Top: heme metabolism,
  MTF1, FOXO3, TBP lean fibroblast; TCF4, JUND, TRAIL lean enterocyte.
  EGF drops out of the top under groupmean.
- fetal fib vs SMC: estimator-sensitive (excess r = 0.53 between
  estimators), so read only what both agree on.

**ECM in fetal fibroblasts vs SMC.** Comparable at cov25 (collagens 19 genes,
core matrisome 85, glycoproteins 59, basement membranes 17, EMT 118,
PROGENy TGFb 183). Under both estimators ANTXR2–ECM coupling is **weaker in
fibroblasts than in SMC**: groupmean ρ_fib 0.18–0.21 vs ρ_SMC 0.40–0.47,
excess z −2.4 to −3.5; pooled ρ_fib 0.31–0.37 vs ρ_SMC 0.55–0.62, excess z
−1.4 to −1.9. **But this is part of a global shift:** in SMC ANTXR2 couples
more strongly to nearly every program (groupmean median ρ 0.33 vs 0.19,
median excess z −1.96, 205/398 programs |z| > 2; pooled median z −0.52).
Relative to that shift, ECM is ordinary. Top programs by |z| are general
(JUN, CEBPB, STAT3, RELA, hypoxia, HGF), not ECM.

**Reading.** ANTXR2 does not preferentially join ECM programs in fibroblasts
relative to the other ECM-producing mesenchymal lineage; if anything it is
more tightly embedded in SMC transcriptional variation overall (and is
expressed ~3× higher there: raw mean 0.52 vs 0.19). As with postnatal
enterocytes, a population-wide axis needs to be characterised before any
program-specific reading.

**Open decisions.** (1) Primary estimator: groupmean recommended (bounded,
house convention); pooled kept as sensitivity. (2) Characterise the SMC
global axis (eigen-analysis as done for enterocytes). (3) Dissociation-stress
control for the enterocyte IEG signal.

## Replication pre-registration (written 2026-10-30, before any IBDverse data was read)

**Claim under test.** In postnatal SI, ANTXR2 couples to an
immediate-early / inflammatory program more strongly in enterocytes than in
fibroblasts, beyond what expression-matched genes show (discovery:
`fib_vs_ent_fine_mat3_cov25`, pooled estimator, PC1 weights).

**Replication data.** IBDverse multitissue atlasing cohort (CELLxGENE
`e3460364-c300-4616-9923-3a0224894998`; Sanger, 10x 3' v3, adult terminal
ileum biopsies). Independent donors, lab and chemistry. **Primary:
`disease == normal`, ileum.** Fibroblasts = `fibroblast` + `subepithelial
intestinal fibroblast` (myofibroblasts excluded, as in discovery);
enterocytes = `enterocyte of epithelium proper of ileum`, 3 within-donor
maturation bins on the same fixed markers. All other settings identical
(q 0.15, floor 0.07, ≥ 100 cells/group, cov 0.25, 50 jointly matched anchors,
1,000 donor bootstraps). Crohn's ileum is a secondary, exploratory run.

**Pre-specified objects.**
- Discovery set S (9 programs, excess z ≤ −2 in discovery): CytoSig EGF,
  TRAIL, FGF2, CD40L; CollecTRI JUND, TCF4, SIRT1, PML; Hallmark TNFα/NF-κB.
- `discovery:IEG_core`: one new equal-weight program from the genes in the
  top-15 PC1 loadings of ≥ 3 of S: CDKN1A, JUN, JUNB, ZFP36, MAFF, CCND1,
  KLF6, DUSP1, PPP1R15A, NR4A1.

**Pre-specified tests (replication data, pooled estimator primary,
groupmean reported alongside).**
1. **IEG_core:** excess < 0 with donor-bootstrap 95% CI excluding 0 →
   replicated. Excess z reported.
2. **Set S:** mean excess z over S, one-sided against the distribution of
   mean excess z over 10,000 random same-size sets of replication-comparable
   programs (approximate: programs share genes).
3. **Sign concordance:** number of S with excess < 0 (descriptive; overlap
   makes a binomial p optimistic).
4. **Whole profile (descriptive):** Spearman of excess z across programs
   comparable in both datasets.

Replication counts as **supported** if test 1 passes and test 2 has
p < 0.05; **partial** if only one does; **not replicated** otherwise.

*Amendment (same day, still before any IBDverse data was read):* IEG_core
has exactly 10 genes, so the standard MIN_GENES = 10 rule would void test 1 if
any one gene missed the floor in either cell type. Discovery-derived programs
therefore need ≥ 6 of their genes above the floor in both cell types (per-gene
floor unchanged).

*Discovery reference values (circular: IEG_core was defined from these data;
recorded before replication data were read).* IEG_core in
`fib_vs_ent_fine_mat3_cov25`: ρ_fib 0.004, ρ_ent 0.179, excess −0.125
(95% CI −0.238 to +0.022), excess z −2.46; groupmean: excess −0.137
(−0.239 to +0.014), z −2.27. So discovery itself does not meet test 1's CI
criterion with 14 fibroblast / 13 enterocyte donors; the criterion relies on
IBDverse's larger cohort.

*Amendment 2 (2026-10-30, after reading IBDverse cell counts only; no IBDverse
correlation had been computed).* IBDverse ileum biopsies are stroma-poor:
1,847 fibroblasts in 243 normal donors (median 10/donor); 1 donor reaches
100, 30 reach 20 (991 cells). Crohn's ileum: 2,326 fibroblasts, 33 donors
≥ 20. User decision: (a) fibroblast donor × subtype groups need ≥ 20 cells
(enterocyte groups still ≥ 100 per maturation bin); (b) the replication slice
is ileum `normal` + Crohn's `uninflamed` (author_inflammation_status), pooled;
(c) pooled estimator only for the primary verdict (the groupmean estimator
averages per-group ρ and is too noisy at 20-cell groups). Expectation set in
advance: the fibroblast side is underpowered (~1–2k cells vs 4.7k in
discovery), so a failed test is weak evidence against the claim; a pass is
informative. Tests and verdict rule unchanged.

## Replication result (2026-10-30): NOT REPLICATED in healthy/uninflamed ileum

Run: `fib_vs_ent_ibdverse_fine_mat3_cov25` (ileum, normal + Crohn's
uninflamed; 2,006 fibroblasts in 43 donor groups; 124,711 enterocytes, 178
donors, 532 donor × maturation-bin groups; estimator validated against
memento, max error 3.7e-15). `replication.json` / `replication.png` there.

| test | result |
|---|---|
| 1. IEG_core excess | −0.001 (95% CI −0.110 to +0.117), z −0.01; ρ_fib 0.05, ρ_ent 0.07 (discovery: −0.125, ρ_ent 0.18) → **fail** |
| 2. set S mean excess z | +0.54 vs random-set null mean +0.58 (sd 0.31), p = 0.44 → **fail** |
| 3. S with excess < 0 | 1 of 9 |
| 4. whole-profile Spearman | 0.28 (241 shared programs, p = 7e-6) |

**Not a power failure:** the replication CI for IEG_core is narrower than
discovery's and excludes the discovery estimate, and ANTXR2 is measured
better (reliability fib 0.36 / ent 0.47 vs 0.17 / 0.18; enterocyte median UMI
24,158 vs 4,227). The enterocyte-side coupling itself shrinks (ρ_ent IEG_core
0.18 → 0.07). Screen-level calibration in IBDverse is again at chance (|z| > 2:
21 of 399 vs 18.2 expected).

**What does carry over (descriptive, not pre-registered):** the weak
whole-profile agreement comes from the *fibroblast-leaning* side. Same sign,
|z| > 1.5 in both: MTF1 (3.35 → 1.65), MXD1 (1.94 → 2.69), CytoSig IL27,
IL10, Hallmark complement, CollecTRI MYBL2, ELF4, RARA.

**Exploratory, inflamed Crohn's ileum** (`fib_vs_ent_ibdverse_inflamed_*`;
only 532 fibroblasts, 13 donors; fragile): the enterocyte-leaning pattern
reappears. Set S mean excess z −0.92 (random-set p = 0.006), 8 of 9 negative;
IEG_core ρ_ent 0.245, excess −0.11 (CI −0.23 to +0.17). Whole-profile
Spearman vs discovery 0.12.

**Reading.** The enterocyte IEG/inflammatory preference is not a stable
property of healthy ileal enterocytes: it is absent in fresh, deeply
sequenced uninflamed biopsies and present in discovery (Elmentaite surgical
resections / full-thickness tissue, shallower data) and, tentatively, in
inflamed Crohn's tissue. That pattern is what a state- or handling-dependent
stress/inflammation response would look like (dissociation/ischemia time,
or genuine inflammation), not a fixed cell-type difference. Distinguishing
those would need the dissociation-stress control and a better-powered
inflamed cohort.

## State contrast plan (written 2026-10-30, before any state-contrast scoring)

Motivation (user): Bracq et al. 2025 (EMBO Mol Med; PMC12514029) — Cmg2 KO
mice are normal at baseline; after DSS the injury is equally severe but
regeneration fails: Ly6a+ fetal-like cells form, conversion to Lgr5+ ISCs does
not, with less nuclear β-catenin and lower Lgr5/Axin2/Ascl2/Ccnd1. KO is
global; Cmg2 is in enterocytes, ISCs and fibroblasts. Prediction for
wild-type human tissue: ANTXR2's coupling to Wnt output in crypt cells should
be stronger in injured (inflamed) than in uninjured tissue.

**Pairing is impossible.** No IBDverse ileum donor has more than one state
(243 normal|uninflamed, 52 Crohn's uninflamed, 66 Crohn's inflamed). The
contrast is therefore between donors, with the donor bootstrap stratified by
state. In exchange, all three states are biopsies under one protocol, so
tissue handling is matched — an inflamed-vs-uninflamed difference cannot be a
handling artefact, which the discovery-vs-replication difference could be.

**Prior crypt test (IBD_WNT_NOTES.md, SCP259 UC colon) was null with no
power:** ANTXR2 raw mean 0.06–0.14 in stem/TA. In IBDverse ileal stem/crypt
cells (sampled 1,500 cells per type × state) it is 1.1–1.7, detected in
52–63%; median UMI 23–45k. Caveat: 3' v3 data with deep coverage of a long
(~130 kb) gene; anchors are matched on raw mean and detection, so this
affects reliability, not the null.

Populations (ileum):
- `crypt` = author_cell_type ∈ {Enterocyte stem OLFM4lo LGR5hi, Enterocyte
  stem MKI67, Enterocyte progenitor crypt OLFM4hi}; groups donor × subtype,
  floor 100 → 26 inflamed / 35 CD-uninflamed / 99 normal donors.
- `enterocyte` = as in the replication (ontology class), donor × maturation
  bin (3), floor 100.

Contrasts: A = Crohn's inflamed, B = Crohn's uninflamed (primary: isolates
inflammation from disease) or normal (secondary). Excess > 0 = ANTXR2 couples
more to the program in inflamed tissue than matched anchor genes do.

Fixed modules (collection `bracq`, equal positive weights, ≥ 6 genes passing
in both states, fixed now):
- `WNT_output` = LGR5, AXIN2, ASCL2, CCND1, RNF43, ZNRF3, SMOC2, EPHB2, EPHB3,
  SOX9 (Bracq's readouts + canonical intestinal Wnt targets)
- `YAP_fetal` = CCN1, CCN2, AMOTL2, ANKRD1, AXL, ANXA1, CLU, TACSTD2, LAMC2
  (YAP targets incl. Bracq's Cyr61/Ctgf; fetal-like / revival markers)

Hypotheses, primary contrast, PC1 weights:
- **H1 (crypt, Bracq):** `bracq:WNT_output` excess > 0. Supported if the
  donor-bootstrap CI excludes 0.
- **H2 (crypt):** `bracq:YAP_fetal` — no directional prediction (Bracq: CMG2
  is dispensable for YAP reprogramming); reported only.
- **H3 (enterocyte, our discovery pattern):** `discovery:IEG_core` excess > 0,
  CI excluding 0; and set S mean excess z > 0 vs random same-size sets
  (one-sided p < 0.05).
Library-wide screens in both compartments are descriptive and read against
their calibration plot. The secondary (vs normal) contrast is a consistency
check, not a second chance: a hypothesis it alone supports is not supported.
This is an exploratory drill-down to decide where to look next; nothing
here is confirmatory.
Settings as in discovery and replication: fine grouping, 3 maturation bins
(enterocytes), min coverage 0.25, pooled estimator, `--validate`.

## State contrast result (2026-10-30)

Runs (all validated against memento, max error ≤ 1.2e-14; bootstrap
stratified by state since no donor is in both): `crypt_infl_vs_{uninflamed,
normal}_fine_cov25`, `enterocyte_infl_vs_{uninflamed,normal}_fine_mat3_cov25`.
Table: `/data/antxr2/program_coexpr/output/state_contrasts.csv` (`program_state.py`).

| | crypt, vs CD-uninflamed | crypt, vs normal | enterocyte, vs CD-uninflamed | enterocyte, vs normal |
|---|---|---|---|---|
| donors (inflamed / ref) | 26 / 35 | 26 / 99 | 32 / 41 | 32 / 137 |
| ANTXR2 raw mean (inflamed / ref) | 1.33 / 1.16 | 1.33 / 1.06 | 2.87 / 2.17 | 2.87 / 2.03 |
| WNT_output excess [CI] | −0.037 [−0.114, +0.046] | −0.072 [−0.145, +0.005] | −0.009 [−0.060, +0.040] | −0.039 [−0.068, −0.007] |
| IEG_core excess [CI] | −0.019 [−0.147, +0.106] | +0.022 [−0.087, +0.125] | **+0.109 [+0.030, +0.186]** (ρ 0.249 vs 0.125; z 1.54) | **+0.135 [+0.066, +0.195]** (ρ 0.224 vs 0.066; z 2.11) |
| set S mean z (random-set mean; p greater) | +0.39 (+0.35; 0.45) | +0.63 (+0.25; 0.13) | +0.58 (+0.27; 0.20) | +0.90 (+0.44; 0.21) |
| screen \|z\|>2 obs / exp | 23 / 23.4 | 42 / 22.9 | 28 / 22.0 | 156 / 21.9 |

- **H1 (crypt WNT_output up in inflamed): not supported.** The point estimate is
  the opposite sign in both references; CIs include 0. Same for every Wnt
  program in the library (TCF7 is the most negative, z −2.1 / −2.6).
  Baseline crypt coupling to Wnt programs is positive (ρ 0.08–0.39) but only
  0–1.4 anchor SDs above matched genes — not distinctive.
- **H2 (YAP_fetal): untestable.** 3 of 9 genes clear the floor in crypt
  cells; ANXA1, TACSTD2, CCN1, ANKRD1 are essentially absent (raw mean ≤ 0.012)
  even in inflamed Crohn's crypts. The fetal-like regenerative state Bracq
  describes is not represented in these annotated crypt populations, so the
  step their KO blocks (fetal-like → Lgr5+) cannot be examined here.
- **H3 (enterocyte IEG_core up in inflamed): IEG_core passes, set S does not.**
  The IEG_core increase is consistent across all 3 maturation bins (median
  per-group ρ inflamed 0.15–0.21, CD-uninflamed 0.02–0.08, normal ≈ 0) and
  in both references. But z ≈ 1.5–2.1: roughly 1 in 15–30 matched genes
  gains as much coupling, so this is ANTXR2 sitting at the upper edge of a
  general inflamed-state increase rather than a uniquely ANTXR2 effect.
  Library-wide excess z is shifted positive in every state run (median
  +0.2 to +0.5): in inflamed tissue ANTXR2 couples somewhat more to most
  programs than its matched genes do.
- Enterocyte screen, primary contrast at chance; the vs-normal screen is
  strongly overdispersed (156 vs 22) with EGFR/MAPK footprint (PROGENy EGFR,
  z +3.4 / +5.5; ERF, ATF1) gaining and Wnt/estrogen/VEGF footprints losing
  coupling. EGF was in discovery set S, so the MAPK→IEG axis is the common
  thread across discovery, inflamed Crohn's, and this contrast.
- Crypt vs normal is overdispersed too (42 vs 23), led by HSF1/HSF2 (heat
  shock), ECM regulators, coagulation, GATA5, SPDEF — mostly gaining.

**Reading.** In matched-handling biopsies, inflammation raises ANTXR2's
coupling to the immediate-early / MAPK program in enterocytes, which supports
part of the earlier "state-dependent" reading over pure handling artefact.
It does not show the crypt Wnt coupling Bracq's model predicts, and the
regenerative state their model is about is not sampled. Next, if pursued:
locate fetal-like/revival cells in IBDverse (ANXA1/TACSTD2/CLU-high
epithelium may be annotated outside the Stem category) before any further
Wnt test; and ask whether the enterocyte IEG coupling tracks inflammation
severity or the epithelial EGFR ligand (AREG/EREG) environment.

## Reframing (user, 2026-11-03): the ANTXR2 co-expression landscape

The aim is descriptive: what ANTXR2's co-expression landscape looks like
across gut populations and states — not recovering proposed mechanisms (Wnt,
collagen clearance). Those are annotations to read off the map afterwards.
Contrasts (fibroblast − enterocyte, inflamed − uninflamed) become slices of
one matrix rather than separate tests. Agreed scope: IBDverse ileum only,
every population where ANTXR2 is measurable, × 3 states; programs primary,
genome-wide gene-level partners alongside.

### Landscape design (fixed before any landscape expression data was read)

- **Populations** (rule on cell counts only): every ileum `author_cell_type`
  with ≥ 10 donors holding a ≥ 100-cell group in at least one state → 26
  types (B 4, enterocyte 2, myeloid 3, plasma 3, goblet 4, stem/crypt 3,
  T 7). Plus **fibroblast (pooled)** = ontology fibroblast + subepithelial
  intestinal fibroblast, fine = ontology class, floor 20 (as in the
  replication; 14 / 19 / 30 donors), flagged throughout as the one
  exception to the floor.
- **Cells of the matrix** = population × state (normal, CD-uninflamed,
  CD-inflamed). A cell is filled if it has ≥ 10 donor groups and ANTXR2 raw
  mean ≥ 0.07 there; otherwise it is reported as not measurable, never as 0.
- **Groups** = donor × author subtype (one subtype per population, so donor),
  centring as before. Size factors per population (all states jointly).
- **Programs:** the existing library + discovery/bracq modules (signed).
  A program's **core genes** = members clearing the floor in ≥ 50% of filled
  cells; it is scored in a cell if that cell has ≥ 80% of the core and
  ≥ MIN_GENES (6 for fixed modules); the shared weight vector restricted to
  the cell's core genes and renormalised. Coverage is reported per cell.
- **Weights:** one per program for the whole matrix. PC1 of the core-gene
  correlation matrix averaged over filled cells (each entry over the cells
  where both genes pass), oriented to agree with the signed equal-weight
  score; equal weights as sensitivity.
- **Calibration:** per population, 50 anchors matched to ANTXR2 on raw mean
  + detection jointly across its filled states (so state differences within
  a population are calibrated by the same genes). Per cell: excess = ρ_ANTXR2
  − mean ρ_anchor, z = excess / anchor SD.
- **Uncertainty:** one donor bootstrap shared by every population (donors
  resampled within state; IBDverse donors are in one state each), with draws
  saved, so the CI of any difference between two cells comes from the same
  draws.
- **Gene level:** per cell, ANTXR2's corrected correlation with every floor
  gene and the same for the 50 anchors → gene-level z. Used to find partners
  outside the libraries and to annotate programs, not as a separate test.
- Validated against memento's compute_2d_moments on a subset of populations.
- Everything is descriptive. The landscape is read against its own null
  (the anchors' spread, calibration plots per cell), not against hypotheses.

### Correction (2026-11-03, after the first landscape scoring): a gene-length axis

The first landscape pass (anchors matched on raw mean + detection only, as in
every earlier run) showed, in all 59 filled cells, that a gene's correlation
with ANTXR2 rises with its genomic span (Spearman 0.13–0.61), and so does the
anchors' mean correlation (0.24–0.72). Recurrent ANTXR2 "partners" were long
genes (AHNAK, PITPNC1, ATP2B1, CASK, MAP3K7, OSBPL10, TBC1D9). ANTXR2 spans
224 kb, longer than 60–96% of its anchors, so it rides that shared axis harder
than they do: calibrated z still tracked span (0.03–0.39) and the median
program z was positive in almost every cell. IBDverse is whole-cell 10x 3' v3;
the likeliest source is intronic reads, varying per cell with unspliced
content and accumulating in long genes. Not verified at read level.

Fix: anchors are drawn from genes within ±0.15 log10 of ANTXR2's span
(~160–320 kb), then matched on raw mean + detection (`ANCHOR_SPAN_WINDOW`,
`gene_span.py`). Anchors then share the length axis and it cancels in the
excess. Diagnostic kept in the report: span-vs-z per cell should go to ~0.

Also, the per-cell calibration against N(0,1) was the wrong reference (one
gene, correlated programs; observed |z| > 2 counts were mostly *below* the
N(0,1) expectation). It is replaced by ANTXR2's count against the same count
for each anchor (leave-one-out against the other 49).

**Consequence for earlier results:** every run before this (discovery,
replication, state contrasts) matched anchors without span. Their contrasts
difference the length axis between populations only partly; the inflamed
cells show a weaker length axis (span-vs-ρ Spearman drops, e.g. top-villus
enterocytes 0.54 → 0.35), so the enterocyte IEG_core state increase in
particular needs re-checking with span-matched anchors before it is quoted.

### Landscape, first full result (2026-11-03; span-matched anchors)

Output: `/data/antxr2/program_coexpr/output/landscape/` (`cells.csv`,
`landscape.csv`, `genes/`, `boot/`) and `report/` (`program_landscape_report.py`).

**Coverage.** 59 of 81 population × state cells filled. ANTXR2 clears the
floor in all ileal epithelium (stem 1.0–1.6, enterocytes 2.2–3.4 raw mean),
B, plasma and most T populations (0.3–1.7), CD163hi macrophages (1.9–2.0) and
pooled fibroblasts (1.6–2.7). Below floor: B germinal centre, CD163lo
macrophages (very shallow cells), γδ T. Estimator matches memento in every
validated cell (max error 4.5e-14). 437 programs have a matrix-wide core;
lineage-restricted ones do not (collagens 1 core gene, basement membrane,
proteoglycans, the Wnt-output module) — the shared-weights matrix cannot
show ECM programs; that needs a lineage-level panel.

**Calibration.** After span matching, gene-level z is ~unrelated to partner
length (Spearman −0.19 to +0.19, most ~0.05). Per cell, ANTXR2's count of
programs with |z| > 2 sits anywhere from the 0th to the 96th percentile of
its own anchors; the cells where it is most structured are cycling stem cells
in inflamed Crohn's (104 programs vs anchor median 18; 96th), CD4 memory
CXCR6 (92–94th), inflamed enterocyte precursors and top villus (84–88th).

**Lineage-wide profile (same-gene nulls: 50 lineage anchors scored in every
cell).** A first pass that drew a different random anchor per population
gave 125 epithelial programs at p < 0.01 vs 4 expected — an artefact: mixing
genes shrinks the null, while any single real gene carries a consistent
profile. With proper same-gene nulls:
- ANTXR2's lineage-wide profile is **not more structured than a matched
  gene's**: programs with |consensus z| > 2 = 10 (epithelium; 64th
  percentile of 50 matched genes), 6 (lymphocytes; 76th), 10 (fibroblast,
  3 cells; 82nd).
- Its cross-lineage agreement (epithelium vs lymphocyte, Spearman 0.44 over
  124 programs) is also typical: the 8 genes that anchor both lineages give
  0.06–0.51.
- What is specific is *which* programs sit at the extremes, and they are the
  same in both lineages: ANTXR2 couples **less** than every or nearly every
  matched gene to proliferation/growth programs (E2F targets, MYC targets,
  DNA repair, G2M; mTORC1 in lymphocytes; ranks 1–2 of 51) and **more** to
  CytoSig TGF-β3, PROGENy hypoxia and CollecTRI ZNF331 (ranks 46–51), with
  |consensus z| ≈ 2–3 in both. Fibroblasts (3 cells, thin) lean to
  interferon programs. With ~400 programs, a few rank-1/51 values per lineage
  arise by chance; the observation is the coherence of the proliferation set
  across two lineages, not any single program.
- Recurrent gene-level partners (|z| > 3 in ≥ 3 cells, all positive):
  AHNAK (17 cells), KLF6, ANXA1, CD55, MYADM, ITGA2, ITGB6, CEACAM5/1,
  DUOX2/DUOXA2, TXNIP, MCL1, TSC22D3, DUSP6 — mostly epithelial injury /
  inflammation-response genes.

**Reading.** In this map ANTXR2 behaves like an ordinary gene of its length
and expression in how much it co-varies with programs; its identity is a
consistent tilt away from proliferation/biosynthesis and toward TGF-β /
hypoxia / stress-response programs within populations, in both epithelium
and lymphocytes. Natural follow-ups: whether ANTXR2 is lower in cycling cells
within each population (a mean-level check of the same tilt); a lineage-level
panel for ECM/collagen programs in stroma; and re-checking the earlier
contrasts with span-matched anchors.

### Fibroblast interferon lean: what drives it (2026-11-03)

Question (user): why do so many interferon programs top the fibroblast
profile? Checks (`program_landscape_fib_subtype.py` for the subtype rerun and means,
-> `output/landscape/fib_author_subtype/`; the rest from the landscape
tables):
- **One signal, not six.** The top six (IRF1 regulon, CytoSig IFN1/IFNG,
  hallmark IFN-γ/IFN-α, STAT1 regulon) share 20–77% of core genes.
- **Present in every state, normal included.** Per-cell excess z 1.4–3.2;
  bootstrap CIs exclude 0 in normal and CD-uninflamed for IRF1/IFN1/IFNG/
  hallmark IFN-γ, but in inflamed (14 donors) every CI crosses 0.
- **Broad, not driven by a few genes.** Across the ~360 core genes of these
  programs, gene-level z is shifted (median +0.12/+0.15/+0.30 in normal/
  uninflamed/inflamed vs −0.13/−0.12/−0.06 for other genes). Leading genes
  differ by state: inflamed = canonical type-I ISGs (SP110, IFIT2, TRIM22,
  MX1, IRF9, SAMD9L, TNFSF10); uninflamed = IFN-γ antigen-presentation arm
  (GBP4, TAP1, PSMB9, PSME1, SOCS1, ERAP2, CTSS); normal = mixed (EPSTI1,
  PARP9, SAMD9, RNF213 plus non-ISG regulon members AKT1, CDKN1B).
- **Fibroblast-specific.** Same programs: epithelium consensus z 0.5–1.1
  (ranks 39–44/51), lymphocytes −1.4 to +0.1.
- **Partly subtype composition.** Groups are donor × ontology class, and the
  "fibroblast" class pools LP ADAMDEC1, CPM and C3 author subtypes. CPM is
  lower than ADAMDEC1 in both ANTXR2 (0.70 vs 0.90 normal; 0.75 vs 0.98
  uninflamed) and IFN score (−0.14 vs 0.00; −0.10 vs +0.08), in all three
  states — mixing them makes positive coupling. Rerun with donor × author
  subtype groups (floor 20; normal 10 donors/456 cells, uninflamed 10/877,
  inflamed drops out): normal unchanged (IFN1 z 2.7→2.9, IFNG 2.3→2.8,
  IRF1 2.5→2.1, STAT1 2.0→2.3); uninflamed shrinks ~30–45% (IRF1 3.2→2.3,
  hallmark IFN-γ 2.8→1.6, IFN-α 2.4→1.4). Whole-profile Spearman between
  groupings 0.81 (normal) / 0.92 (uninflamed).
- **Gene-level caveat.** Fibroblast gene tables contain ρ > 1 for sparse
  genes (16 in normal, 5 inflamed; e.g. PARP9 0.996, EFHD2 1.32): the
  Poisson-corrected variance is underestimated in small groups. Program
  scores are unaffected in kind but gene-level fibroblast ranks are noisy.

**Reading.** A modest within-fibroblast coupling of ANTXR2 to a broad IFN
response, present in healthy tissue and not produced by subtype mixing
alone; thin (10–25 donors per state, ~20-cell groups). Not yet checked:
immune-doublet/ambient ISG contamination, and whether it holds within
ADAMDEC1 fibroblasts alone.

## Next steps (as of 2026-11-03)

In rough priority order; none started.

1. **Re-check the pre-landscape contrasts with span-matched anchors.** The
   gene-length axis ("Correction") was found after the discovery,
   replication and state-contrast runs, whose anchors are matched on raw
   mean + detection only. Every one of those results, the IEG_core state
   increase (H3) first, is provisional until rerun with anchors inside
   `ANCHOR_SPAN_WINDOW`. Needs `program_scores.py`'s anchor matching to call
   `program_landscape.span_matched`-style candidate filtering.
2. **Mean-level check of the proliferation tilt.** Is ANTXR2 lower in
   cycling (MKI67+/G2M-high) cells within each landscape population? A direct
   expression check of the landscape's strongest lineage-wide signal.
3. **Fibroblast IFN lean.** (a) within LP fibroblast ADAMDEC1 alone, the one
   subtype with enough cells; (b) immune doublet / ambient ISG contamination
   (e.g. coupling to PTPRC/LYZ/CD74 reads in fibroblasts); (c) another
   stromal cohort.
4. **Stromal ECM panel.** Collagen / basement-membrane / proteoglycan
   programs have no core under the landscape's shared weights; score them in
   a stroma-only matrix (IBDverse fibroblast subtypes, myofibroblasts,
   pericytes; Gut Cell Atlas fibroblasts for donors).
5. **Report fixes.** `similarity.png` shows cells with < 100 scored programs
   as blank rows; fibroblast gene tables have rho > 1 for sparse genes
   (flag or drop genes whose corrected variance is below a floor).

