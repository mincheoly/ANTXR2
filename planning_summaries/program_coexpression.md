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
| `annotation/gene_span.csv` | spans | `gene_span.py` (md5 `178d0e6fe8c5450e24c0c04104af70cf`) |

The program library is in the repo, not the data root:
`genesets/programs_long.csv` (md5 `875a18c2e9107b3e66f947f8456d9832`) is the
library every result used, and `genesets/raw_cytosig/signature.centroid` is
the CytoSig input `program_genesets.py` needs. Regenerating the library
fetches Hallmark, PROGENy, CollecTRI and MSigDB matrisome live from OmniPath,
which drift, so rerunning it changes the library; don't, unless that is the
point.

`program_extract.py` refuses an h5ad whose byte size, `uns["citation"]`
version or nnz differs from `program_config.SOURCES`.

**Commands** (from `scripts/`; `py` = the environment's python). Run
directories are `output/<contrast>_<grouping>[_mat<N>][_covNN]/`.

```
# py program_genesets.py                       # rebuilds genesets/programs_long.csv; see above
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


## Proliferation tilt (2026-11-10)

The landscape's clearest lineage-wide signal: ANTXR2 couples less than its
span-matched anchors to E2F / MYC / DNA repair / G2M (epithelium) and
mTORC1 / E2F / MYC (lymphocytes), ranks 1–2 of 51. Taken step by step.

### Step 1: mean level (`program_landscape_cycling.py` -> `output/landscape/cycling/quartiles_cycle.csv`, `stem_pair.csv`)

Within population × state, per donor, Q4 vs Q1 of a Tirosh S + G2M score;
log2 ratio of mean normalised expression for ANTXR2 vs the population's 50
landscape anchors.
- **Epithelium: ANTXR2 is lower in the more-cycling cells in all 26
  cells** (excess < 0 in 26/26; |z| > 1.5 in 11). Not independent tests:
  states share anchors and populations share donors.
- The effect scales with how much Q4 actually enriches cycling cells
  (Spearman of MKI67-detection gain vs excess, all 59 cells: −0.47). Where
  Q4 is ≥ 10 points richer in MKI67+ cells (17 cells), excess < 0 in all 17,
  median z −1.7. Largest: OLFM4hi crypt progenitors (Q4 93% MKI67+ vs 12%):
  ANTXR2 −0.43 log2 while anchors +0.17–0.19; excess −0.6, 91–97% of donors,
  all three states. Goblet SPINK4hi (−0.3 to −0.4, z −1.8 to −2.4) and
  MUC2int (−0.4 to −0.8) the same.
- Within all-cycling "Enterocyte stem MKI67", S/G2M phase barely matters
  (excess −0.1, z −0.4). Between author types, MKI67 vs LGR5hi stem cells
  per donor: ANTXR2 −0.6 log2 vs anchors −0.2 to −0.3; excess −0.27 to
  −0.42 (z −0.8 to −1.4; 78–91% of donors). Same direction, weaker.
- **Lymphocytes: not this.** Q4 holds 0.3–4% MKI67+ cells, i.e. there are
  almost no cycling lymphocytes to compare; excess < 0 in 15/27 cells, with
  CD8 TRM TRGC2 and CD4 memory CXCR6 the exceptions (z −1.1 to −3.1). The
  lymphocyte tilt (mTORC1 / MYC first) must be a growth / activation axis,
  not division.
- Myeloid / fibroblast: no signal (few cycling cells, few donors).

**Reading.** In epithelium the co-expression tilt has a plain mean-level
counterpart: ANTXR2 is ~25% lower in absolute terms (~35% relative to
matched genes) in cycling crypt and goblet-lineage cells than in their
non-cycling neighbours.

### Step 2: does the tilt survive accounting for dividing cells? (`program_landscape_cycle_conditioned.py`)

Landscape prep + score rerun for the 23 epithelial + lymphocyte populations
with the main weights and bootstrap, two ways; `compare` ->
`output/landscape/cycle_compare.csv`, `cycle_compare_cells.csv`.
- **drop** (MKI67-detected cells removed): the tilt is essentially
  unchanged. Over the 44–47 cells scored in both runs, mean excess z
  E2F targets −1.17 → −1.10, MYC V1 −1.02 → −0.95, G2M −0.83 → −0.74;
  epithelial consensus z E2F targets −3.08 → −3.00, MYC V1 −2.24 → −1.99,
  still rank 1–2 of 51. Decisive case: OLFM4hi crypt progenitors (normal),
  where dropping MKI67+ cells removes 70% of cells (16,839 → 5,061; 93 → 30
  donors), E2F targets z −2.03 → −1.99, MYC V1 −2.07 → −2.05.
- **bins** (donor × within-donor tertile of the S + G2M score; tertiles
  separate well, top 86% / bottom 16% MKI67+ in OLFM4hi): too costly. At the
  100-cell group floor only 18 of 53 cells stay filled, and the most
  proliferative populations (goblet SPINK4hi / MUC2int, LGR5hi stem, 2 of 3
  OLFM4hi states) drop out. On the 17 cells left, proliferation z shrinks
  20–50% (G2M most), but the control side shrinks as much (TGF-β3 1.03 →
  0.80, hypoxia 0.90 → 0.50, ZNF331 0.80 → 0.49): precision loss from
  thinner groups, not a cycle-specific explanation. OLFM4hi normal keeps
  E2F targets at −2.1.

**Reading.** Not composition. ANTXR2's negative coupling to E2F / MYC
programs holds among non-dividing cells, i.e. it tracks a graded
proliferative / biosynthetic axis, which also fits the lymphocyte tilt
where almost no cells divide. Caveat: "MKI67 not detected" is not "not
cycling" (dropout; G1/early S); a covariate-adjusted estimator (partial
moments on a cell-level cycle score) would close that gap without losing
donors.

### Step 4 plan: replication in the Gut Cell Atlas (written 2026-11-10, before any GCA landscape expression was read)

- **Data.** GCA Extended+ (`757945c8-…`), postnatal small intestine (all
  non-trimester ages), `disease == normal`, all studies (Kong2023, Elmentaite
  2021, Martin2019, others): donors disjoint from IBDverse; labs other than
  IBDverse's (Kong2023 is Broad). Elmentaite2021 SI enterocytes were read
  earlier in this pilot for the fibroblast/enterocyte question, never for
  the proliferation tilt. Crohn's left out: 19 donors carry both labels,
  which breaks the one-state-per-donor bootstrap.
- **Populations**, by the landscape's rule on cell counts only (≥ 10 donors
  with a ≥ 100-cell `level_3_annot` group): epithelial — Enterocyte (67
  donors), TA (51), Goblet (21), Goblet_progenitor (18), Goblet_cycling (13),
  Paneth (12), Epithelial_stem (12); lymphocyte — Trm_CD8 (22), Trm_CD4 (18),
  Trm_Th17 (17), B_plasma_IgA1 (17), B_memory (13), B_naive (13).
- **Method fixed from IBDverse:** the IBDverse landscape's weights.csv
  (no re-fitting), its core and coverage rules, floors, span-matched
  population and lineage anchors, donor bootstrap — `ANTXR2_LANDSCAPE=gca`
  runs the same scripts on this atlas.
- **Primary tests (epithelium):**
  - P1: lineage consensus z < 0 with ANTXR2 in the bottom 5 of 51 (vs the 50
    lineage anchors) for hallmark E2F_TARGETS and for hallmark
    MYC_TARGETS_V1.
  - P2: mean consensus z over the IBDverse epithelial set (consensus rank
    ≤ 2 there: collectri E2F2, E2F3, E2F4; hallmark E2F_TARGETS,
    MYC_TARGETS_V1, DNA_REPAIR, G2M_CHECKPOINT) below that of 10,000
    random same-size sets of programs scored in GCA, one-sided p < 0.05.
  - P3 (mean level, step 1): in cells whose top cycle-score quartile is
    ≥ 10 points richer in MKI67+ cells, ANTXR2 excess < 0 in most; sign test.
- **Secondary:** lymphocyte P1 (MTORC1_SIGNALING, MYC_TARGETS_V1) and P2
  (IBDverse lymphocyte rank ≤ 3: MTORC1, E2F_TARGETS, MYC_TARGETS_V1/V2);
  the positive side (CytoSig TGFB3, PROGENy Hypoxia, CollecTRI ZNF331)
  consensus z > 0; per-cell tables.
- **Verdict:** replicated if P1 (both programs) and P2 pass; partial if the
  direction holds (both P1 consensus z < 0 and P2 mean z < 0) without the
  thresholds, or only one P1 program passes; not replicated otherwise. P3 is
  reported beside it.
- *Clarification (same day, still before any GCA expression was read):*
  the P2 sets are the named lists above. The rank rule as worded also
  catches CollecTRI SFPQ, PRDM1, ZNF143 (epithelium) and NFYA, CREBBP
  (lymphocytes); they are not in the sets, which are the proliferation /
  growth programs the tilt was described by.

### Step 3: the growth axis (`program_landscape_cycling.py --axis growth [--no-mki67]`)

Same mean-level design, quartiles of a growth score (hallmark MTORC1 + MYC V1
genes minus the cycle genes; 367 genes) -> `cycling/quartiles_growth*.csv`.
- **Lymphocytes: the tilt is a growth axis.** ANTXR2 is lower than its
  anchors in the high-growth quartile in 19/27 lymphocyte cells, 12 with
  z < −1.5, median z −1.17 (cycle axis: −0.07). Strongest in naive B
  (ANTXR2 −0.7 to −1.1 log2 vs anchors ~+0.05; z −2.2 to −3.3; 90–100% of
  donors), memory B PTPRJ (z −2.2 to −3.6), CD4 memory CXCR6 (z −2.1 to
  −3.8), CD8 TRM TRGC2, tissue Treg (z −3.0). Plasma cells are the exception
  (excess +0.1 to +0.3; every gene rises with their growth score).
- Epithelium: 24/26 negative, median z −1.56 (cycle axis −1.29).
- **Among non-dividing cells (MKI67-negative only):** epithelium 20/23
  negative, median z −1.08 (OLFM4hi progenitors −1.2 to −2.3, SPINK4hi
  goblet −1.3 to −2.2); lymphocytes 19/27, median z −0.91, B cells
  unchanged (naive B −2.2 to −3.3). The growth coupling is not division.

**Reading.** Across epithelium and lymphocytes ANTXR2 is lower, at the mean
level, in the cells with the most growth / biosynthesis (mTORC1, MYC)
activity, dividing or not — the mean-level counterpart of the landscape's
negative coupling to MYC / mTORC1 / E2F programs.

### Step 4 result: GCA replication (2026-11-10)

`ANTXR2_LANDSCAPE=gca` landscape (13 populations, 83 normal donors; scorer
matches memento to 6e-15) + `program_landscape_replicate.py` ->
`output/landscape_gca/replication.json`.
- **Pre-registered tests: untestable, not failed.** ANTXR2 clears the raw-mean
  floor in only 5 of 13 cells (enterocyte 0.08, goblet 0.11, goblet
  progenitor 0.07, B memory 0.14, B naive 0.10; IBDverse enterocytes ~2.2),
  and those score 20–46 programs. The floor is the house 0.07, the same one
  the original GCA fib-vs-enterocyte pilot used. Epithelial consensus needs a
  program scored in all 3 filled cells and holds only 5, none named;
  lymphocytes have 2 filled cells (no consensus). P3 has no qualifying cell
  (TA 0.051, stem 0.058, cycling goblet 0.067: just below the floor).
  *Correction (same day):* named programs were scored per cell, just not in
  every cell of a lineage: MYC_TARGETS_V1 in enterocyte (excess z −0.34),
  goblet progenitor (−0.90) and B naive (−1.26); B memory MTORC1 −0.22,
  DNA_REPAIR −0.02, MYC_V1 +0.22, MYC_V2 +0.08. Weak and mostly negative;
  IBDverse top-villus enterocytes are similarly weak (MYC_V1 −0.43), the
  strong cells there being crypt progenitors, whose GCA counterparts sit
  below the floor.
- **Why ANTXR2 is ~25× lower: IBDverse counts intronic reads, this GCA
  release does not.** Gene-by-gene IBDverse / GCA mean ratio vs genomic span,
  ~9,800 genes: Spearman 0.74 (enterocytes) and 0.70 (memory B); genes of
  ANTXR2's length (log10 span 5.3–5.6) are ~4× (log2 +2.0) higher in
  IBDverse relative to the median gene, short genes (< 3 kb) ~2.5× lower.
  ANTXR2 sits on that line (+1.9 / +2.2). This is the same gene-length axis
  found inside IBDverse ("Correction"): most of ANTXR2's IBDverse signal is
  pre-mRNA. The span-matched anchors make the IBDverse comparisons fair, but
  an exon-only atlas measures a much sparser ANTXR2.
- **Exploratory (not pre-registered), mean level with anchors matched per
  cell (`--all-pops`), cycle axis:** TA cells (51 donors; Q4 72% MKI67+ vs
  4%): ANTXR2 −0.21 log2 while anchors +0.75; excess −0.96, z −1.5, 2nd
  percentile of anchors, 71% of donors — the counterpart of IBDverse's OLFM4hi
  progenitors (excess −0.6, 2nd percentile). Not seen in cycling goblet cells
  (13 donors; excess +0.04) or stem cells (12 donors; Q4 only 6% MKI67+).
  Growth axis: goblet cycling −1.12 (z −1.6), B memory −1.44 (z −2.3), but
  B naive +1.38 (z +1.8) and the rest near 0. Anchor SDs are 2–5× IBDverse's
  (ANTXR2 raw mean 0.02–0.14), so most cells cannot resolve the IBDverse
  effect size.

**Reading.** Partial at best: the single best-powered GCA test of the core
finding (TA cells, exon-only counts) agrees in direction and anchor
percentile; everything else is too sparse to say, and the lymphocyte growth
result does not reproduce. A fair replication needs an independent gut
dataset counted with introns (single-nucleus, or Cell Ranger ≥ 7 /
include-introns), where ANTXR2 is measured the way IBDverse measures it.

## Cross-dataset replication (2026-11-10)

Goal (user): co-expression and mean results that hold in both gut datasets,
which differ in platform (IBDverse: ileal biopsies, intron-counted; GCA
Extended+: postnatal small intestine, exon-only). One analysis at a time,
each run in both: (1) the mean result, ANTXR2 rising with differentiation;
(2) co-expression modules in fibroblasts; (3) co-expression modules in
epithelium. Framing: gut only; no "stromal, not epithelial" claim.

### 1. Mean level: ANTXR2 vs epithelial differentiation

`program_landscape_differentiation.py` (both atlases) ->
`output/landscape{,_gca}/differentiation/{ladder,within}.csv`. ANTXR2 against
50 span-matched genes matched on raw mean + detection (candidates >= 0.01).
Ladder: each population vs the stem population (IBDverse LGR5hi stem; GCA
Epithelial_stem), per donor with >= 100 cells of both, anchors matched in the
stem population only. Within: top vs bottom quartile of the fixed maturation
score inside a population.

| ladder vs stem | IBDverse (3 states) | GCA |
|---|---|---|
| dividing (stem MKI67, OLFM4hi progenitor / TA) | excess −0.80 to −0.96 log2, z −2.3 to −3.0, 0–2nd pct, 0% of donors above (9–42 donors) | TA: −0.62, z −1.0, 12th pct, 27% of 11 donors |
| differentiated absorptive (OLFM4lo precursor, top villus / Enterocyte) | ANTXR2 +0.51 to +0.77 while anchors −0.77 to +0.08; excess +0.63 to +1.53, z 0.8–1.5, 82–96th pct, 100% of donors (12–45) | ANTXR2 +0.43, anchors −0.53; excess +0.95, z 0.7, 82nd pct, 91% of 11 donors |
| goblet | ANTXR2 +0.9 to +1.6 but anchors +0.0 to +0.9; excess +0.17 to +1.17, 64–84th pct, 78–100% of donors | 3–6 donors per population; excess −0.4 to +0.2 (no power) |

Within population (maturation quartiles): IBDverse positive in every crypt
and precursor population (excess +0.14 to +0.82, 84–98th pct), flat in top
villus (already mature; 24–70th pct). GCA: TA +0.62 (88th pct, 51 donors),
Enterocyte +0.50 (z 1.4, 92nd pct, 82% of 66 donors); GCA stem −0.97 (8th
pct, 12 donors) disagrees.

**Reading.** Replicates in the absorptive lineage in both datasets, in
direction, donor consistency and anchor percentile, on both tests: ANTXR2
dips in dividing progenitors and rises in differentiated enterocytes more
than ~80–95% of matched genes. It is not an extreme outlier (no z > 2 on the
rising side): about 1 matched gene in 10 rises as much. Goblet: IBDverse
only, and weaker once anchors are subtracted. Note: `stem_pair` in the
cycling script matched anchors in both populations, which biases its excess
toward 0; this script matches in the base population only.

### 2. Fibroblast co-expression: plan (written before any GCA fibroblast correlation was read)

Discovery = IBDverse pooled ileal fibroblasts (landscape, 3 states,
`report/consensus_fibroblast.csv`, 381 programs): the interferon lean.
Replication = GCA Extended+ postnatal small-intestine fibroblasts
(`level_2_annot == Fibroblast`, all studies, `disease == normal`; known so
far only from obs: 13 Elmentaite2021 normal donors + a few others, subtypes
ADAMDEC1 / PI16 / F3). Same pipeline as IBDverse fibroblasts: donor ×
subtype groups with >= 20 cells, >= 10 donors, ANTXR2 raw mean >= 0.07, 50
span-matched anchors, IBDverse PC1 weights (not re-fitted).

One deviation, fixed now: a program is scored in GCA where >= 50% of its
IBDverse core genes clear the floor (>= 10 genes), not 80% — the exon-only
atlas is shallower and at 80% only 20–46 programs scored per cell in step 4.
Output `landscape_gca_cov50/` (`ANTXR2_CELL_COVERAGE=0.5`); the 80% subset is
reported as a sensitivity check.

- **F1 (primary), the IFN set** named on 2026-11-03: collectri:IRF1,
  cytosig:IFN1, cytosig:IFNG, hallmark:INTERFERON_GAMMA_RESPONSE,
  hallmark:INTERFERON_ALPHA_RESPONSE, collectri:STAT1. (a) mean GCA excess z
  vs 10,000 random same-size sets of GCA-scored library programs, one-sided
  greater; (b) same-gene null: ANTXR2's mean ρ over the set ranked among its
  50 anchors' (percentile >= 90 to pass). Untestable if < 3 of 6 are scored.
- **F2, whole profile:** Spearman of excess z, IBDverse fibroblast consensus
  mean z vs GCA fibroblast, programs scored in both (descriptive p; programs
  overlap).
- **F3, other IBDverse extremes** (rank >= 50 of 51, not in the IFN set):
  cytosig:IL27, matrisome:NABA_ECM_REGULATORS, cytosig:FGF2,
  collectri:ZNF91, collectri:FOXO4 — sign and z each. No negative-side set:
  the lowest IBDverse fibroblast program is rank 3.
- Verdict: replicated if F1a p < 0.05 and F1b passes; partial if one does.

### 2. Fibroblast co-expression: result (2026-11-10)

GCA extraction now includes pooled fibroblasts (`gca_landscape` adapter);
`ANTXR2_LANDSCAPE=gca ANTXR2_CELL_COVERAGE=0.5 program_landscape.py all
--pops fibroblast_pooled --validate` -> `output/landscape_gca_cov50/`;
`program_landscape_replicate_fib.py` -> `replication_fibroblast.{json,csv}`.
GCA fibroblast cell: 4,225 cells in 40 donor × subtype groups, 24 normal
donors, ANTXR2 raw mean 0.31 (24% detected), reliability 0.16; 6,244 floor
genes; 380 programs scored (62 at >= 80% coverage); matches memento to 8e-15.

- **F1 (IFN set): partial.** All 6 scored (coverage 0.70–0.80). Excess z:
  IFNG +1.00, hallmark IFN-γ +0.99, IRF1 +0.69, IFN1 +0.33, hallmark IFN-α
  +0.32, STAT1 −0.11 (IBDverse consensus z 2.4–3.0). (a) mean z +0.54 vs
  random 6-program sets: p = 0.022, passes. (b) ANTXR2's set ρ 0.168 vs
  anchors' 0.112: 80th percentile of its 50 anchors, below the 90 required.
- **F2 (whole profile): weak.** Spearman 0.12 over 352 shared programs
  (p = 0.02); 0.14 on the 62 programs at >= 80% coverage.
- **F3:** IL27 +1.05, NABA_ECM_REGULATORS +1.08, FGF2 +0.06, ZNF91 −0.59,
  FOXO4 −1.09: 3 of 5 same sign, none beyond z 1.1.
- GCA's own extremes (fatty-acid metabolism / peroxisome +3.4, SOX9 regulon
  −4.9) are near 0 in IBDverse.

**Reading.** Same direction, a third of the size: the IFN lean is positive in
5 of 6 programs in an independent, exon-only atlas and beats random program
sets, but 1 in 5 matched genes couples to IFN programs as much. By the
pre-set rule this is partial. ANTXR2's fibroblast program profile as a whole
barely carries over (0.12), in line with the earlier postnatal-vs-fetal GCA
comparison (0.14) and with the IBDverse finding that ANTXR2's profile is no
more structured than a matched gene's.

## Reference-scaled coupling (2026-11-10)

Design change (user): shared weights should be learnable for programs some
cell types do not express, and the coupling there should read ~0 instead of
blank. Also: report an activity number beside every coupling, use 200
matched genes, and defer a real testing strategy.

`program_landscape_ref.py`, `ANTXR2_LANDSCAPE_RUN=ref` ->
`output/landscape_ref/` and `output/landscape_gca_ref/` (`coupling.csv`).
- **Weights:** a gene enters a program if it clears the floor in >= 3 filled
  IBDverse cells (was: half of all cells); PC1 of the correlation averaged
  over the cells where each pair is measured. 600 programs (was 437),
  8,760 genes; collagens, basement membrane, proteoglycans now included.
- **Scale:** each gene is divided by a reference SD = 90th percentile of its
  corrected SD (per-10k units) over the cells where it clears the floor,
  i.e. where it varies most. (A median was tried first and rejected:
  lineage-high genes dwarfed their reference, activity up to 69.)
- **Readout:** beta = cov(ANTXR2, score) / (SD ANTXR2 here × program
  reference SD) = rho × activity_sd; excess = beta − mean of 200 span-matched
  anchors; pct = share of anchors below ANTXR2. Reported with
  `activity_mean` (mean per-10k expression of the program's genes),
  `activity_sd` (program SD here / reference) and `frac_pass`.
- **GCA:** IBDverse weights, GCA's own reference SDs, ANTXR2 floor 0.05
  (exploratory): 12 of 14 cells filled (all epithelium incl. TA, stem,
  Paneth; Trm CD4/Th17; B memory/naive; fibroblasts), 600 programs each.
- No bootstrap yet. z and pct are noise where a program is silent (both
  ANTXR2 and anchors ~0); read excess, or filter on activity_sd.
- 200 anchors is a looser match: e.g. IBDverse fibroblasts raw mean 0.9–4.6
  around ANTXR2's 1.6.

Checks (IBDverse, 437 old programs × 58 cells): new vs old per cell,
rho r = 0.80, excess r = 0.70 (median). Collagens: activity_sd 0.01–0.06 and
beta 0.00 in every non-fibroblast population, 1.18 in fibroblasts. E2F
targets: activity_sd 0.7–1.0 in dividing epithelium with ANTXR2 at the
2nd–5th percentile, as before.

### 2b. Fibroblasts under the new estimator, both datasets

IBDverse = mean over 3 states; GCA = normal. 597 library programs.
- **Interferon is the only thing shared at the top.** Programs at >= 90th
  percentile in both datasets (activity_sd >= 0.2 in both; 20 in IBDverse,
  38 in GCA): IRF2 regulon, CytoSig IFNG, CytoSig IL27, hallmark complement,
  hallmark IFN-α, hallmark IFN-γ — six, all interferon-type. GCA percentiles
  for the pre-named set: IFNG 97.5, IFN-γ 96, IFN-α 93, IRF1 82, IFN1 78,
  STAT1 64 (IBDverse 96–99). GCA excess is about a third of IBDverse's.
- **Nothing shared at the bottom** (<= 10th percentile in both: none).
- **Collagen / matrix:** NABA collagens 16th (IBDverse) and 23rd (GCA)
  percentile, ECM regulators 19th / 21st, TGF-β 14th / 24th: ANTXR2 couples
  a little less than matched genes, not an outlier. Core matrisome,
  glycoproteins, basement membrane disagree between datasets.
- Whole profile: Spearman 0.10 (excess) between datasets; IBDverse states
  agree with each other at only 0.38–0.46, so the ceiling is low.
- Program activity_mean agrees well across datasets (Spearman 0.76).

### 3. Epithelial co-expression, both datasets (reference-scaled runs)

IBDverse: 9 epithelial populations × 3 states = 26 cells. GCA: 7 populations
(ANTXR2 floor 0.05). Per program, the median anchor percentile over the
cells where it is active (activity_sd >= 0.2, in at least half the cells):
511 programs in both. No pre-registered test (testing strategy deferred).
- **The epithelial profile carries over:** Spearman 0.51 between datasets on
  median percentile (0.56 on mean z), against 0.10 for fibroblasts.
- **Low in both (median pct <= 20; 67 IBDverse, 13 GCA, 12 shared):** hallmark
  E2F targets (3rd / 13th), E2F3 regulon (5th / 15th), HMGA1, TCF7, and
  chromatin / repressor regulons (HDAC1, DNMT3A, KMT2B, HIC1, ZBTB17, SNAI2,
  PAX8, ZKSCAN7).
- **High in both (>= 80; 92 IBDverse, 30 GCA, 13 shared):** differentiation
  regulons CDX2, HNF1A, KLF4, KLF5, PPARA, NR2F2, plus RB1, IRF1, NFATC3,
  MAZ, DBP, FOXC2, LHX2.
- **Proliferation set by GCA population** (percentile): TA — E2F 1, DNA
  repair 2, MYC V1 3, MYC V2 4; stem — MYC V1 5, E2F 11; cycling goblet —
  MYC V1 4, E2F 14. In GCA's non-dividing populations (Enterocyte, Goblet,
  Paneth) the same programs are mid-range (35–83), whereas in IBDverse they
  are low in 60–90% of all epithelial cells. So: replicated where cells
  divide, not as an epithelium-wide tilt.
- **Not replicated:** CytoSig TGF-β3 (IBDverse 90th, GCA 63rd), PROGENy
  hypoxia (85th / 52nd), ZNF331 (86th / 37th).

- **Whole-profile scatter** (`scripts/program_landscape_ref_scatter.py`,
  `landscape_gca_ref/report/epithelial_pct_ibd_vs_gca.png`): 513 programs,
  Spearman 0.51, Pearson 0.55; 356 (69%) on the same side of the 50th
  percentile in both; slope 0.36 (GCA compressed toward the middle, SD 17
  vs 25). The agreement runs through the whole cloud, not the extremes.
  User's read (2026-11-10): favourable, justifies the cross-dataset design.

- **Same scatter, other compartments** (`report/scatter_summary.csv`):
  fibroblast 558 programs, Spearman 0.07 (1 GCA population); all
  lymphocytes 503, 0.13; T 507, -0.04 (GCA Trm_CD4, Trm_Th17 at ANTXR2
  0.05-0.07); B naive/memory 494, 0.30 (GCA ANTXR2 0.10-0.14).
  Caveat on comparing compartments: the epithelial 0.51 is a median over
  7 GCA populations; single matched epithelial pairs give 0.15-0.32
  (all single cross pairs: median 0.19). Fibroblast single pairs: -0.03,
  0.04, 0.21 against the three IBDverse states (which agree 0.38-0.46
  among themselves). So fibroblasts are weaker than epithelium, but less
  so than 0.07 vs 0.51 suggests; T cells show nothing.

- **Scope decision (user, 2026-11-10): lymphocytes are dropped going
  forward.** GCA Extended+ composition (level_1_annot): postnatal small
  intestine 385,271 cells, 58% epithelial, 22% T/NK, 13% B/plasma, 2.7%
  mesenchymal; whole atlas 1,596,200 cells, 41% epithelial.

- **Fetal fibroblasts added to the GCA side (user, 2026-11-10).** Extended+
  holds 59,192 fetal SI fibroblasts (26 donors; ANTXR2 0.23), in the same
  file as the postnatal ones (`program_extract.py --source gca_fetal_fib`;
  scored with `--pops fibroblast_fetal fibroblast_postnatal_fetal` on the
  existing GCA reference SDs, earlier rows unchanged). Cross-dataset
  Spearman vs IBDverse fibroblasts: postnatal 0.07, fetal 0.00,
  postnatal + fetal -0.04 (63,362 cells, 49 donors; 0.92 with fetal
  alone, since fetal is 93% of the pool). GCA fetal vs GCA postnatal
  profiles agree only 0.05, so the two are not the same profile with
  different noise. Fetal ANTXR2 sits above its anchors across the board
  (median pct 67 vs 46 postnatal), as in the pilot. Interferon programs
  stay high in the pool (IFN-g response 88, IFN-a 88, IFNG 80,
  complement 93) but IRF2 drops to 61. More cells did not rescue the
  fibroblast profile.

- **Donor split-half inside each dataset (user, 2026-11-11;
  `program_landscape_ref_donor_split.py`, `<run>/donor_split/`).** 200
  random halvings of the donors; Spearman between the halves' program
  percentile profiles (median, 10th-90th pct of splits), and the
  Spearman-Brown reliability of the full-donor profile:
  - GCA postnatal fibroblasts (24 donors, 4,225 cells): 0.04 (-0.04 to
    0.11), positive in 73% of splits; reliability 0.07. **Not
    reproducible against itself**, so it cannot replicate anything.
  - GCA fetal fibroblasts (25 donors, 59,137 cells): 0.20; reliability 0.34.
  - IBDverse fibroblasts: normal 0.29 (25 donors, 867 cells), CD
    uninflamed 0.46 (18, 1,139), CD inflamed 0.21 (13, 532); reliability
    0.35-0.63.
  - Context, epithelium: GCA Enterocyte 0.25, Goblet 0.14, TA 0.12;
    IBDverse top-villus enterocyte 0.73-0.89, MKI67 stem 0.32-0.70.
  So the fibroblast cross-dataset failure is a GCA-side noise problem in
  postnatal fibroblasts, not evidence of different biology; the fetal
  profile is modestly reproducible and still does not match IBDverse.

- **ECM programs in the most reliable fibroblast cell (IBDverse CD
  uninflamed; raw rho, anchors set aside at the user's request,
  2026-11-11).** 563 active programs, rho median 0.07, range -0.37 to
  0.46; the top is immune/interferon (NPM1, STAT5A, IRF2, RELB, IRF1,
  JAK-STAT, complement), not matrix. ECM programs, PC1 rho (rank of
  563): TGF-b signalling 0.26 (53), coagulation 0.24, ECM glycoproteins
  0.23 (82), proteoglycans 0.20; collagens -0.06 (436), ECM regulators
  -0.18, basement membranes -0.18, EMT -0.22 (549), core matrisome
  -0.22 (550). Signs hold in all three states. The matrix PC1s are
  contrast axes, not 'more matrix': core matrisome/EMT + = COL6A1,
  COL6A2, COL4A5/6, TNC, POSTN, LAMA5; - = DCN, LUM, FBLN1, TNXB; the
  glycoprotein PC1 is the same axis flipped. ANTXR2 is on the
  DCN/LUM/FBLN1 side in every case. Equal-weight rho is small and
  unstable across states (core matrisome 0.07-0.28, collagens
  0.00-0.15). Not yet checked at gene level.

- **Testing strategy, first version: donor bootstrap (user, 2026-11-11;
  `program_landscape_ref_test.py`, `<run>/test/`).** Compares two
  population x state cells of one atlas per program: d_beta = beta_A -
  beta_B and d_excess = (beta - mean anchor beta)_A - (...)_B. Donors
  are resampled (2,000 draws; union of the two cells' donors, shared
  donors kept paired) and the stored per-group sufficient statistics
  re-summed; no within-donor cell resampling (planned as a second layer
  in memento's GPU formulation once the host has a GPU). se = SD of
  draws, normal-approximation p, BH q over programs.
  - *Calibration* (20 random donor halvings of one cell tested against
    each other, truth = 0; share of active programs with p < 0.05,
    median over splits; nominal 0.05). IBDverse: top-villus enterocyte
    (145 donors) 0.04, LGR5 stem (46) 0.05, fibroblast CD uninflamed (18)
    0.06, CD inflamed (13) 0.04-0.06. GCA: Enterocyte (67) 0.04-0.05,
    stem (12) 0.08, fibroblast (24) 0.02-0.04. So per-test size is about
    right, mildly anti-conservative at 12 donors. But programs are
    correlated: a single null split reached 13-51% of programs at
    p < 0.05, so the COUNT of hits is not evidence. BH held better:
    d_excess q < 0.05 never exceeded 5% of programs in a null split;
    d_beta once hit 39% (GCA Enterocyte). Prefer d_excess.
  - *First comparison, top-villus enterocyte vs LGR5 stem, normal.*
    IBDverse (145 vs 46 donors, 37 shared): d_excess q < 0.05 for 120 of
    527 programs (108 higher in enterocyte: RXRB, ETV4, NR2F2, KLF6,
    CytoSig TNFA / IL1 / EGF / HGF; 12 lower: DNMT3A, SIRT1, SNAI2, ZEB1,
    TCF4, TGFB3). GCA Enterocyte vs Epithelial_stem (67 vs 12 donors):
    13 of 556. Across datasets d_excess correlates only 0.08 over 525
    programs; of the 120 IBDverse hits, 65% have the same sign in GCA
    and 21 are also p < 0.05 there. The contrast does not replicate as a
    profile; GCA stem (12 donors, ANTXR2 0.058) is the weak arm.

- **Testing moved into memento (user, 2026-11-14).** Branch
  `feature/program-correlation` of scrna-parameter-estimation, in the
  worktree `/home/ubuntu/Github/memento-program-correlation` (commit
  84dbf6b, local only): `memento.ht_program_correlation` (pooled
  within-group program coupling, reference scale, anchors as
  `references` -> `excess`, contrast between two conditions, NumPy and
  CUDA kernels). ANTXR2 adapter: `program_landscape_ref_memento.py`
  (`compare`, `calibrate`, `--mode donors|both|cells`), outputs in
  `<run>/test_memento/`; needs `PYTHONPATH=<worktree>`. Its observed
  values reproduce coupling.csv and its donor-mode results reproduce
  `program_landscape_ref_test.py` (294 / 117 programs at q < 0.05 for
  d_beta / d_excess, top-villus vs LGR5 stem). 49k cells, 600 programs,
  201 targets, 2,000 draws: 30 s on the T4, 108 s on CPU.
  - *Which layer to resample* (20 null donor-half splits each; share of
    programs with d_excess p < 0.05, SD of z; nominal 0.05 and 1):
    donors only 0.04-0.07, 0.93-1.01; cells only 0.05-0.06, 0.93-1.01;
    donors + cells 0.000-0.005, 0.58-0.67 (LGR5 stem 46 donors,
    fibroblast CD uninflamed 18, CD inflamed 13). The two-level
    bootstrap double counts cell sampling noise (SE ~1.5x too large);
    donors only is the default. Cells only is also calibrated on this
    null, i.e. donor-to-donor variation in coupling is small next to
    cell sampling noise here, but it assumes independent cells.
  - Correlated programs still make hit counts unreliable: one null
    split reached 17% of programs at q < 0.05 (donors only).
  - Host: the T4 needed the 580-server NVIDIA driver for kernel
    6.8.0-1066 (535 had modules for 6.5 only); torch 2.5.1+cu121 is in
    the venv.
- **Donors as strata (user, 2026-11-15).** User's view: when both cell
  types exist in each donor, resample cells within donor x cell type only
  (= `--mode cells`); `compare --shared` restricts to shared donors.
  Top-villus vs LGR5 stem, 37 shared donors, 17,384 of 49k cells:
  d_excess q < 0.05: cells 35, donors 3, both 0 (median se 0.022 / 0.024 /
  0.033); d_beta 250 / 212 / 102. Donor/cell SE ratio per program median
  1.10 (IQR 0.99-1.25), so donor x cell-type variation in coupling is
  small. Restricting to shared donors costs 26% in SE; all donors, cells
  mode: 134 (the 35 are a subset). Memento's `resample_rep` is the same
  two-stage scheme as `both`; toy simulation: calibrated when donor
  variance >= 10x cell noise, sqrt(2) SE inflation when ~0 (our case,
  ratio 1.37-1.39).
- **Mixed-model-type test (user, 2026-11-15; prototype
  `program_landscape_ref_mixed.py`, outputs `<run>/test_mixed/`).**
  Coupling written as a per-cell mean of z = S~ T~ / (sd_T sigma_ref);
  random donor intercept on z, moment estimate of between-donor variance
  (truncated at 0), GLS mean with weights 1 / (sb2 + se2 / n_i). sd_T
  fixed; donor effects of the two cell types treated as independent.
  - Top-villus vs LGR5 stem, d_excess q < 0.05 of 527: mixed 66 (median
    se 0.022) vs donor bootstrap 117 (0.018), cell bootstrap 134 (0.017),
    analytic cells-only 113 (0.019). d_beta: 250 vs 294 / 319 / 289.
  - Null donor halves (200 splits): mixed is conservative, 1.9-3.2% at
    p < 0.05, z SD 0.81-0.91 (top villus, stem, fibroblast CD uninflamed
    / inflamed); analytic cells-only 2.4-6.8%, z SD 0.84-1.05.
  - sb2 > 0 in 87% of top-villus programs yet cells-only is calibrated
    there, so the moment estimate is biased up (truncation; probably
    unequal within-donor variances).
  - Donor-specific within-donor variances (DerSimonian-Laird, `--het`):
    SE halves (0.010) but the estimate moves (r 0.60 with pooled for
    d_excess) because precision weights favour low-variance donors; null
    still conservative in median (2.4-3.6%) with unstable splits (max 57%
    of programs at p < 0.05). Not usable.
  - No power gain here: between-donor variance in coupling is ~0, so the
    model can only add a noisy non-negative term.
- **Decision (user, 2026-11-15): headline test = cells resampled within
  donor x subtype strata, all donors** (`--mode cells`, now the adapter
  default; files end `__cells`). Donor bootstrap is the check.
  - IBDverse top-villus vs LGR5 stem: d_excess q < 0.05 134 of 527
    (d_beta 319). GCA Enterocyte vs stem (119,093 cells): 11 of 556
    (d_beta 86), median se 0.029. Fibroblast CD uninflamed vs inflamed:
    0 (d_beta 1).
  - Across datasets (525 shared active programs): d_excess Spearman 0.08,
    3 programs significant in both (all same sign), 64% of IBDverse hits
    same sign in GCA; d_beta Spearman -0.02, 48 significant in both, only
    19 same sign. The enterocyte-vs-stem difference does not replicate,
    whatever the test.
- **Enterocyte vs fibroblast (user, 2026-11-15: "an easier one"), cells
  mode, normal tissue.** IBDverse top-villus vs pooled fibroblasts (867
  fibroblast cells in groups): d_excess q < 0.05 39 of 574 (d_beta 48),
  median se 0.067. GCA Enterocyte vs pooled fibroblasts (4,225 cells):
  2 of 590 (d_beta 7), se 0.051.
  - The contrast is the fibroblast side: corr(d_excess, -fibroblast
    excess) 0.97 in both datasets (SD of excess: enterocyte 0.03 / 0.02,
    fibroblast 0.13 / 0.06); 20 of the 39 IBDverse hits are programs
    active only in fibroblasts.
  - IBDverse hits: ANTXR2 couples in fibroblasts more than in enterocytes
    to coagulation, ECM glycoproteins, NR1H2, TBX2, SNAI2, SMAD7; less to
    (negative in fibroblasts) collagens, core matrisome, basement
    membranes, EMT, TGFb, GLI2/3, ZEB2, SMAD3 - the same fibroblast
    contrast axis as the ECM sanity check.
  - No replication: d_excess Spearman 0.07 across datasets (0.06 on
    programs active in both), 0 significant in both, 28 of 39 IBDverse
    hits same sign in GCA. GCA estimates are noise: variance of d_excess
    across programs / mean se^2 = 0.99 (IBDverse 2.3).
- **Skin fibroblasts (user, 2026-11-15: replicate where fibroblasts are
  plentiful).** Human dermal fibroblast atlas (CELLxGENE dataset a19d1667,
  version 8ae4471a, 25 studies, fibroblasts only; `/data/antxr2/raw/skin/`,
  1.8 GB). Healthy, 10x 3' v2 / v3: 78,643 cells, 105 donors (96 with a
  >= 20-cell group). `ANTXR2_LANDSCAPE=skin` (+ `_RUN=ref`), IBDverse
  weights, own reference SDs -> `landscape_skin_ref/`. Two populations:
  pooled within donor, and within donor x author cluster (`cluster_robust`).
  ANTXR2 raw mean 0.29 / cell (IBDverse fibroblasts 1.57, GCA 0.31).
  - Reliable: donor split-half 0.58 (reliability 0.73; by cluster 0.48 /
    0.65), vs IBDverse fibroblasts 0.21-0.46 and GCA 0.04. Excess q < 0.05
    (cells mode): 173 of 600 (by cluster 104); median se 0.011. Pooled vs
    by-cluster profiles agree 0.91.
  - Whole profile does not match gut: skin vs IBDverse fibroblasts
    (precision-weighted over the 3 states, 551 programs) Spearman -0.05;
    159 IBDverse hits, 85 same sign in skin; 61 significant in both, 38
    same sign. By collection: matrisome +0.50 (8), cytosig +0.31 (42),
    hallmark +0.12 (47), CollecTRI -0.11 (438; 24 of 44 same sign).
  - The ECM axis does replicate in sign: negative in both for collagens
    (skin -0.051 / IBDverse -0.138), core matrisome (-0.026 / -0.195),
    basement membranes (-0.044 / -0.222), EMT (-0.024 / -0.227),
    myogenesis (-0.037 / -0.221); positive in both for coagulation
    (0.043 / 0.359), proteoglycans (0.060 / 0.149), IFN-alpha (0.030 /
    0.173), IFN-gamma (0.023 / 0.176). All q < 0.05 in both. Skin effects
    are 4-10x smaller. Not matched: ECM glycoproteins, ECM regulators,
    ECM-affiliated, progeny TGFb (skin n.s.).
  - `program_landscape_ref_memento.py estimate --cell` gives per-cell
    coupling / excess with cell-bootstrap SE (`<run>/test_memento/
    *__estimate__cells.csv`).
- **ANTXR2 vs individual collagens (user, 2026-11-18).**
  `program_landscape_ref_memento.py genes --cell C` (matrisome collagens
  clearing the 0.07 floor; each gene its own one-gene "program", scale None
  -> noise-corrected correlation pooled within donor x subtype; excess over
  the cell's 200 matched anchors; cells mode, 2,000 draws, CPU because the
  GPU driver is not loaded on kernel 6.8.0-1069). Table:
  `landscape_skin_ref/test_memento/collagen_genes_skin_vs_gut.csv` (gut =
  IBDverse fibroblasts, precision-weighted over normal / CD uninflamed /
  CD inflamed).
  - Skin (25 genes, median se 0.013): raw correlations all small (-0.01 to
    0.15). Excess below anchors (q < 0.05) for COL1A1 -0.053, COL1A2
    -0.044, COL3A1 -0.052, COL5A1 -0.053, COL5A2, COL5A3, COL4A1, COL6A1
    -0.039, COL6A2 -0.044, COL18A1 -0.057; above for COL15A1 +0.103 (pct
    97). COL6A3 -0.015 n.s.
  - Gut (median se 0.045 combined, 0.07-0.10 per state): COL6A1 -0.178,
    COL6A2 -0.172, COL6A3 -0.121, COL1A1 -0.112, COL4A5 -0.191, COL12A1
    -0.157 (p < 0.05); COL14A1 +0.123, COL25A1 +0.099. Raw correlation
    with COL6A1/A2 negative in all three states (-0.11 to -0.24).
  - Skin vs gut over 22 shared genes: excess Spearman 0.26, 17 same sign;
    COL6A1/COL6A2 and COL1A1 negative in both. GCA: nothing (0 of 20).
  - Skin pooled within donor x author cluster (`cluster_robust`, 17
    clusters): COL6A1 / COL6A2 excess -0.001 / -0.002 (n.s.), COL1A1 /
    COL1A2 -0.027 / -0.024 (q 0.058), COL3A1 -0.027 (q 0.025); COL15A1
    +0.094 stays (q < 0.001). Only 2 of 25 q < 0.05 (pooled: 11); excess
    profiles agree 0.62. So the skin COL6 result is between-subtype
    composition, not within-subtype covariation.
  - Gut, pooled within donor x author subtype (LP ADAMDEC1, CPM, Crypt
    PI16, C3; 20-cell floor keeps 456 / 877 / 307 cells of 10 / 10 / 7
    donors; scratchpad driver over `cmd_genes`, outputs
    `landscape_ref/test_memento/by_author_subtype/`): COL6A1 excess -0.178
    -> -0.087 (p 0.11), COL6A2 -0.172 -> -0.069 (p 0.24), COL1A1 -0.112 ->
    -0.041, COL12A1 -0.157 -> -0.067; COL4A5 -0.124 stays (p 0.02).
    Composition explains about half of the gut COL6 effect; the remainder
    is not significant with ~40% fewer cells.
  - What the subtypes are (user asked, 2026-11-18). Skin `cluster_robust`
    = the atlas authors' consensus populations (Ascension & Izeta, bioRxiv
    2024.09.05.611379): 5 axes A-E (A ECM-secreting, A2 papillary; B
    immune; C hair-follicle, C1 sheath / C2 papilla; D, E1) plus T / U,
    which the authors call unreliable. They describe A as a continuum and
    A3, B3, C5, E1 as bridges. In our cells A4 is 92% one study with
    mito-RNA markers, U is mito / KRT14, and B1 (CXCL1/2/3, IL6, NFKBIA,
    MT2A) looks like an activation / dissociation state (Tabula Sapiens 92%
    B1, Boothby 80%). Between-cluster driver: ANTXR2-high E1 (1.70 per 10k),
    C2 (1.08), B4 (0.87) are COL6-low (20-35); A2 is ANTXR2 0.62, COL6 98.
  - Gut author subtypes are distinct populations: Fibroblast CPM (SOX6, F3,
    BMP4, POSTN, COL4A5; subepithelial-type) vs LP fibroblast ADAMDEC1
    (ADAMDEC1, APOE, CFD). Only these two reach 20 cells per donor (12 / 24
    donors; PI16 1, C3 0; 10 donors in >= 2 subtypes). Donor-median COL6
    55 vs 16 per 10k, ANTXR2 1.89 vs 2.49: COL6/ANTXR2 29 vs 6.5.
  - In skin the raw ANTXR2-COL6 correlation is slightly positive (0.04-0.05);
    the negative excess is relative to anchors (0.09), i.e. ANTXR2 tracks
    COL6 less than a matched gene does. In gut it is negative outright.
- **Single subtype per tissue (user, 2026-11-18: use the largest, most
  reliable subtype, donor x batch groups).** One sample per donor in both
  atlases, so groups = donor. Gut: IBDverse LP fibroblast ADAMDEC1, normal
  + CD uninflamed pooled (`fibroblast ADAMDEC1 uninflamed`; each state
  alone has 9-10 donors), 983 cells in 19 donor groups. Skin: dermal atlas
  A3 (`fibroblast A3`; B1 excluded as stress state), 16,931 cells, 70
  donors.
  - Donor split-half of the program profile: gut 0.04, skin 0.22
    (pooled skin fibroblasts were 0.58).
  - Programs, excess q < 0.05: 0 of 561 (gut), 0 of 600 (skin); p < 0.05
    41 / 90. Skin vs gut excess Spearman -0.02 (561); raw beta 0.67, i.e.
    shared with matched genes. Nominal in both: IFN-alpha (skin 0.056, gut
    0.122) and IFN-gamma (0.053 / 0.129), coagulation (0.039 / 0.123, gut
    n.s.).
  - Collagens: none q < 0.05 in either. COL6A1/A2/A3 excess skin +0.031 /
    +0.011 / +0.049, gut +0.034 / +0.055 / -0.052. COL1A1 skin -0.038 (q
    0.53), gut +0.009.
  - Reading: within one fibroblast subtype ANTXR2 has no detectable
    excess coupling to collagens or ECM programs in either tissue; the
    earlier negative ECM / COL6 signal was between-subtype composition. The
    COL6/ANTXR2 question is a subtype-level mean question.
- **Can the matched-gene bar be beaten? (user, 2026-11-18: "maybe the
  negative control is a bit too hard to beat").** Positive controls with
  fibroblast ECM roles, none in the collagen/core matrisome programs (MRC2,
  ITGA11, DDR2, FAP, PDGFRB, LRP1, ANTXR1), each against its own 200 anchors
  chosen as ANTXR2's are, same subtypes, cell resampling
  (`program_landscape_ref_controls.py`; `{cell}__controls.csv`).
  - Skin A3 (70 donors), programs at excess q < 0.05: PDGFRB 484, MRC2
    467, LRP1 464, ANTXR1 238, FAP 196, ITGA11 27, DDR2 14, ANTXR2 0.
    NABA_COLLAGENS anchor percentile: MRC2 100, PDGFRB 98, ITGA11 98,
    LRP1 96, ANTXR1 95, FAP 70, ANTXR2 52, DDR2 44.
  - Gut ADAMDEC1 (19 donors): PDGFRB 29, MRC2 11, DDR2 3, all others 0
    (FAP below floor). Collagens pct: LRP1 97, MRC2 97, ANTXR2 92, PDGFRB
    82, ANTXR1 76, DDR2 73, ITGA11 32.
  - The top passers (PDGFRB, MRC2, LRP1) beat their anchors on ~80% of all
    programs (median pct 93-94), i.e. they sit on the dominant fibroblast
    activation axis rather than having ECM-specific coupling.
  - IFN programs: in gut ANTXR2 is the highest of the 7 genes (pct 93-99 on
    IFN-alpha/gamma, IFN1, JAK-STAT, IRF1/4); in skin it is high (78-94)
    but LRP1/PDGFRB match it.
  - Reading: in skin the bar is beatable (5 of 7 controls clear it on
    collagens), so ANTXR2's null there is real: it couples to collagens like
    a typical fibroblast gene of its size. The gut single-subtype test is
    underpowered (even PDGFRB barely clears it); ANTXR2's gut point
    estimates sit with the controls.
- **COL6/ANTXR2 demoted to a side note (user, 2026-11-18); interferon is the
  lead.**
- **Interferon across cell types (user, 2026-11-18: "is that only in
  fibroblasts?").** Pre-specified set of 8 programs (hallmark IFN-alpha,
  IFN-gamma; CytoSig IFN1, IFNG; CollecTRI IRF1, IRF2, STAT1; PROGENy
  JAK-STAT), `estimate --units ... --tag ifn` (q over the 8), cells mode, every
  non-lymphocyte filled cell (`*__estimate__ifn__cells.csv`). Median excess
  (programs q < 0.05; none negative anywhere):
  - IBDverse, normal / CD uninflamed / CD inflamed: fibroblasts 0.12 / 0.18 /
    0.22 (8/8 each); top-villus enterocytes 0.045 / 0.077 / 0.137 (8/8 each);
    crypt OLFM4lo progenitors 0.024 / 0.048 / 0.127 (8/8 each); goblet
    MUC2hi 0.029 / 0.046 / 0.058 (7, 8, 6); MKI67 stem 0.012 / 0.032 /
    0.094 (1, 0, 8); LGR5 stem -0.011 / -0.008 / -0.028 (0). Rises with
    inflammation in every population.
  - GCA (normal): goblet 0.156 (8/8); fibroblasts postnatal 0.052 (0, se
    0.036), postnatal + fetal 0.033 (5); enterocyte 0.015 (0); TA -0.013,
    Paneth -0.037.
  - Skin: A3 0.046 (5/8), pooled 0.017 (4/8).
  - IBDverse Ileum goblet MUC2lo BCAS1int (3 states) hit the estimate's
    beta-vs-landscape assert; not resolved.
- **Stimulated airway epithelium (user, 2026-11-18: memento's HTECs).**
  `program_htec_ifn.py`, GEO GSE274751 (`GSE274751_HTEC_counts.h5ad`,
  `/data/antxr2/raw/htec/`): 69,958 cells, 2 donors, IFN-alpha / beta /
  gamma / lambda x 3-48 h + 2,360 unstimulated. Outputs
  `output/htec_ifn/`.
  - ANTXR2 is low: 0.137 counts / cell in ciliated, ~0.06 in basal and club
    (below the 0.07 floor). Mean vs control (UMI-pooled per donor, donor
    mean), ciliated: alpha 1.2-2.2x, beta 1.25-1.5x at every time point
    (pooled alpha 1.65x / 1.56x in the two donors); gamma, lambda ~1.
    Basal down (gamma 0.2-0.5x), club ~0.8-1.
  - Ciliated coupling (noise-corrected correlation, scale None, strata donor
    x time, 200 matched anchors, cell bootstrap 2,000): 8/8 IFN programs
    excess q < 0.05 under alpha (excess 0.13-0.18, pct 90-94) and beta
    (0.15-0.27, pct 91-97); lambda same direction, n.s. (0.05-0.15, q
    0.14-0.43); gamma ~0 (pct 60-86); unstimulated 0.01-0.13 n.s. (739
    cells, se 0.095). Anchors' mean correlation turns negative under
    stimulation (~-0.05): compositional, ISGs take read share.
  - All 600 programs (500 draws): IFN programs rank 3-68 under beta
    (JAK-STAT 3rd, IFN1 9th), 20-62 under alpha; 34 / 48 programs q < 0.05
    overall (29 / 36 positive), many IFN-adjacent (STAT4, IL27, IFNL,
    allograft rejection). ANTXR2's median excess over all programs is
    +0.02-0.03 under alpha / beta, so IFN (~0.2) stands well above its own
    baseline. Gamma and lambda: nothing.
  - Reading: in ciliated cells ANTXR2 is a weak type I ISG (mean and
    cell-to-cell); matches the tissue picture (coupling largest where IFN
    tone is highest, absent in stem / TA). Basal cells go the other way.
- **Type I IFN gut single-cell data (user, 2026-11-20).** None found with
  direct recombinant type I / III treatment. Closest: GSE171620 (Triana 2021,
  human ileal organoids, astrovirus mock / 4 h / 16 h, endogenous IFN
  response) and GSE156760 (colon / ileal organoids, SARS-CoV-2 mock / 12 /
  24 hpi; ISGs in bystanders). GSE142166 (IFN-beta / lambda organoids) is
  mouse bulk.
- **IFN-gamma colon organoids (user, 2026-11-20).**
  `program_organoid_ifng.py`, GEO GSE255366: mouse, one organoid line,
  untreated / IFN-gamma 24 h / 48 h, 60,134 cells (>= 1,000 UMI, mt <= 25%;
  median 1.3-1.9k UMI). Mouse -> human by MGI one-to-one homology
  (`annotation/HOM_MouseHumanSequence.rpt`). Outputs `output/organoid_ifng/`.
  - Strong response (Stat1 ~17x, Cxcl10 ~300x, Cd74 up); colonocytes mostly
    lost. Leiden clusters on non-induced genes still split by sample.
  - Antxr2 sits in stem / cycling clusters (1.3-2.0 counts / cell), ~0.02 in
    colonocytes, ~0.01 in the treated-only clusters. Per 10k: stem 0.79 /
    0.89 / 1.77, cycling 0.29 / 0.23 / 0.41 (untreated / 24 h / 48 h): no
    consistent induction.
  - Coupling, IFN set, strata sample x cluster (>= 20 cells): stem / cycling
    (6,367 cells, raw mean 0.92) excess -0.11 to +0.07, pct 10-78, nothing
    q < 0.05 (24 h leans negative on hallmark IFN / JAK-STAT, q 0.11); all
    cells excess -0.03 to +0.05, q >= 0.19.
  - Reading: no IFN-gamma coupling or induction, as in HTEC IFN-gamma. Caveats:
    mouse, one line, no replication; Antxr2 is in the stem compartment, where
    tissue coupling was also absent (IBDverse LGR5 stem).
- **Type I vs IFN-gamma split, and fibroblast type I (user, 2026-11-21).**
  `program_ifn_type_split.py` (outputs `output/ifn_type_split/`). Coupling here
  is the noise-corrected correlation (scale None), so tissue numbers are not on
  the beta scale of the earlier IFN-set table.
  - Sets from HTEC pseudobulk (6-48 h vs unstimulated, median log2FC over 6
    cell type x donor): type I-biased 49 genes (MX1/2, OAS1-3, IFIT1-3, RSAD2,
    ISG15, IFI44L, IRF7, USP18, CXCL10/11); gamma-biased 22 (CIITA, HLA-II,
    CD74, CXCL9, IRF1, GBP1/2/5, IDO1, WARS); shared 26 (STAT1, TAP1/2, PSMB8-10,
    HLA-E/F). Strict "only" sets are 2-15 genes: gamma induces ISGF3 genes ~2x
    in airway epithelium. Residual programs I|G and G|I (one set minus its
    within-stratum, noise-corrected regression on the other) stay linear.
  - **Within-cell coupling cannot name the cytokine.** HTEC ciliated (positive
    control): under pure IFN-beta the coupling sits on G|I (0.17, q 0.03), not
    I|G (0.07, n.s.); alpha splits weakly. Gamma-biased genes vary cell-to-cell
    with the type I response. In tissue both residuals are positive in places
    (I|G q < 0.05 in top-villus x3, crypt OLFM4lo uninflamed / inflamed; G|I in
    crypt normal / inflamed, top-villus inflamed): uninformative by design.
  - **Between samples, means do name it.** `donors`: per donor x state sample,
    log2 ANTXR2 per 10k ~ standardized type I + gamma set scores (+ state;
    donor-clustered SE). HTEC ciliated (control, samples donor x stim x time):
    type I +0.24 (p < 1e-4), gamma +0.06 (n.s.); basal gamma -0.26 (p 0.02).
    IBDverse: crypt OLFM4lo type I +0.14 (p < 0.001; inflamed +0.45, uninflamed
    +0.28, normal +0.07), gamma ~0 in every state; LGR5 stem type I +0.20 (p
    0.004), gamma -0.07; top-villus pooled both +0.10, inflamed type I +0.36 (p
    < 0.001) vs gamma +0.09; fibroblasts n.s. (13-25 donors per state); goblet
    nothing. Set scores correlate 0.6-0.9 across samples.
  - ANTXR2 per 10k by state (normal / uninflamed / inflamed): top-villus 0.86 /
    0.89 / 1.17, crypt OLFM4lo 0.82 / 0.83 / 1.07, fibroblasts 2.30 / 2.46 /
    2.35, goblet MUC2hi 1.44 / 1.37 / 1.26, LGR5 stem 0.53 / 0.50 / 0.56. Set
    activity rises with inflammation everywhere (fibroblast type I 3.3 / 4.0 /
    4.2, gamma 0.9 / 1.3 / 1.5).
  - Fibroblasts, IFN-beta (Hagai 2018, E-MTAB-7051; Smart-seq2, 5 donors, 681
    cells, 2 / 6 h, no unstimulated cells; reads, not UMIs): 6 h, 7 of 8 IFN
    programs q < 0.05 (excess 0.15-0.22, pct 99-100; IRF2 null), type I-biased
    0.20 (q 0.002), I|G 0.17 (q 0.001), G|I 0.05 (n.s.). All 605 programs: 37
    at q < 0.05, median excess 0.012; top are TNF / NF-kB (CytoSig TNFA 0.28,
    PROGENy NFkB 0.27), IFN programs rank 5-40. 2 h: nothing (1 program q <
    0.05). ANTXR2 per 10k 2 h -> 6 h not consistent across donors (0.36-0.82).
    No unstimulated cells in this experiment (checked SDRF / IDF); sibling
    Smart-seq2 controls (E-MTAB-7037, 306 cells, pooled wuye / wetu / vass,
    other plate; E-MTAB-5920 / 6066, ~30-60 cells) belong to dsRNA runs.
    The memento paper's Hagai2018 sets (memento repo, publication/.../
    bulk_comparison/hagai) are the Squair-compendium mononuclear phagocytes
    (mouse / pig / rabbit / rat, LPS / poly I:C), not these fibroblasts.
    E-MTAB-7167 (Hang et al., unstimulated HipSci dermal fibroblasts,
    Smart-seq2, 30 donors incl. joxm 83 / vass 66 / wetu 89 cells) is the only
    unstimulated single-cell match for IFN-beta donors: other study and library
    prep, FASTQ only (no matrix).
  - **Induction: bulk, paired (Hagai E-MTAB-5919, `bulk`).** Untreated vs
    IFN-beta 8 h, per donor, median-of-ratios. Columns are lane ids with no key
    in the SDRF; treatment called from MX1 / IFIT1 / RSAD2 (ISG log2FC 5-9),
    pairs by correlation on non-responsive genes (always the adjacent lane).
    Human ANTXR2 up in 6 / 6 donors, log2FC +0.53 (+0.39 to +0.66; 1.3-1.6x,
    joxm and pelm among them). Rhesus -0.11 (6 / 6 slightly down), mouse -0.19
    (3), rat +0.08 (3): the induction is human-specific here, consistent with
    the no-induction mouse organoid result.
  - **Immune contamination (user, 2026-11-21, light-touch; `contam`,
    `output/ifn_type_split/contam/`).** Fibroblasts, top-villus, crypt
    OLFM4lo x 3 states, 1,000 draws, IFN-set median excess (scale None):
    all cells / without >= 2 UMI of pan-immune genes (PTPRC, CD52, LYZ,
    TYROBP, ...; keeps 77-99%) / without each donor's top quarter by ambient
    plasma-cell Ig (IGKC, IGHA1, JCHAIN, ...; keeps 75-90%). Fibroblasts
    0.245 / 0.250 / 0.250 (normal), 0.242 / 0.184 / 0.256 (uninflamed),
    0.261 / 0.277 / 0.290 (inflamed); top-villus 0.090 / 0.082 / 0.076,
    0.135 / 0.118 / 0.128, 0.164 / 0.159 / 0.163; crypt 0.079 / 0.077 /
    0.079, 0.112 / 0.095 / 0.102, 0.224 / 0.193 / 0.228. All 8 programs stay
    q < 0.05 in epithelium (fibroblast counts drop to 3-5 / 8 where n falls
    to 650-880). ANTXR2's own coupling to the contamination programs is ~0:
    ambient Ig -0.02 to +0.02 in epithelium (fibroblast normal 0.25 +- 0.16,
    n.s.), immune-doublet -0.04 to +0.07. Crypt inflamed carries the most
    immune reads (4.1 UMI / cell, 0.1 after the filter) and keeps 0.19.
    Reading: not carried by immune-derived reads.
  - **Cross-species replication (Shaw et al. 2017 PLoS Biol, type I IFN,
    10 species; webserver isg.data.cvr.ac.uk, POST
    ISGwebServer/queryByGeneNameOrEnsemblId; saved in
    `output/ifn_type_split/shaw2017/`).** ANTXR2 log2FC (FDR): human +1.13
    (5e-18, called an ISG); rat +0.10 (0.36), cow -0.50 (0.03, down), sheep
    +0.02, pig +0.11, horse -0.33 (0.05), dog +0.04, flying fox -0.30, little
    brown bat -0.04, chicken +0.02 (all n.s. unless noted). Up only in human
    of 10; independent of Hagai (human 6 / 6 up, macaque / mouse / rat not).
    No macaque or mouse in Shaw. ANTXR1 (paralog): human -0.06, n.s.
    - *Ortholog flags checked 2026-11-21* (raw `q_ANTXR2.json`, fields 6-8 of
      each gene array): the 9 nulls are not detection failures. ANTXR2 is
      `present` AND `expressed` in all 10 species, one ortholog cluster
      (HS9951), protein identity to human 86-91% in the mammals (69% chicken).
      dN/dS 0.18-0.26 wherever reported (rat 0.206, cow 0.196, sheep 0.191,
      pig 0.181, horse 0.196, dog 0.248, little brown bat 0.259; None for
      human / flying fox / chicken): strong purifying selection on the
      protein alongside human-restricted induction. This closes the main hole
      in the cross-species claim and supplies the "conserved sequence,
      unconserved regulation" framing.
  - **Promoter (user, 2026-11-21, "internet search for #2").** No literature
    on ANTXR2 promoter IFN regulation found. ChIP-Atlas (hg38, +-5 kb of TSS;
    `annotation/chipatlas/`): ANTXR2 has weak STAT1 / STAT2 / IRF1 peaks in a
    minority of experiments (STAT1 20% of 103, avg score 55, 58th pct;
    IRF1 14% of 56; IRF9 1 of 56), the level of GAPDH, not of ISGs (MX1 /
    IFIT1 99th pct). Motif scan of -1 kb..+0.5 kb around the Ensembl gene
    start: human ANTXR2 carries a full canonical ISRE at +72
    (AGTTTCTTTTTCT, consensus AGTTTCNNTTTCN) and a core at -179; mouse
    Antxr2 has none in the same window. Weak evidence (one window around the
    annotated start; AT-rich human sequence, GC 0.39 vs 0.53; the same scan
    misses human MX1's ISREs, which sit outside this window).
  - Fibroblasts, dsRNA (E-MTAB-5988 / 5989; 10x, 3 donors pooled, not
    demultiplexed; 2,553 unstimulated, 1,944 poly I:C 6 h): ANTXR2 0.23 -> 0.35
    per 10k (1.5x, one pooled sample). 6 h: IFN programs excess 0.44-0.48 (pct
    98), but 285 / 605 programs q < 0.05, median excess 0.107, NF-kB / TNF /
    inflammatory as high: a responder-state effect, not IFN-specific.
    Unstimulated: nothing.
  - Reading: cytokine identity has to come from stimulation and from means.
    Both point to type I: alpha / beta induce and couple (HTEC), IFN-beta
    couples in dermal fibroblasts at 6 h, gamma does neither (HTEC, mouse colon
    organoids), and across tissue samples ANTXR2's mean follows type I, not
    gamma, activity in crypt progenitors, stem and inflamed enterocytes.
    Fibroblast tissue means are underpowered. The fibroblast IFN-beta coupling
    shares the top with NF-kB / TNF.

**Reading.** In both datasets ANTXR2 co-varies with the differentiation
regulons and against E2F / proliferation within epithelial populations —
the within-population counterpart of the step 1 mean result, probably the
same maturation gradient seen two ways rather than an independent finding.
