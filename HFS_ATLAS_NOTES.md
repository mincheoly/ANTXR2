# HFS cross-tissue atlas arm — methods, findings, and standing corrections

Companion to `ANALYSIS_SUMMARY.md`, which remains authoritative for the
co-expression / memento work. This file covers only the CELLxGENE cross-tissue
arm (`scripts/hfs_*.py`) and follows the same convention: **the corrections
section wins over anything stated earlier in this file or in session logs.**

## What this arm asks

Not "where is *ANTXR2* expressed" but "where is there the most collagen VI per
receptor to clear". Two readouts per tissue, both from library-size-normalised
pseudobulk CPM:

| metric | meaning |
|---|---|
| `col6_frac_of_collagen` | COL6 as a share of the tissue's total collagen output — composition, robust to absolute output |
| `col6_per_antxr2` | COL6 CPM / *ANTXR2* CPM — substrate:receptor ratio, the "load per receptor" |

Per the project's standing rule, these are reported side by side and **never
combined into a composite score.**

## Pipeline

```
scripts/hfs_config.py             paths, dataset registry, gene groups, thresholds
scripts/hfs_pseudobulk.py         streaming h5ad readers + extract() + pseudobulk()
scripts/hfs_download_atlases.py   fetch 11 CELLxGENE h5ads, verify identity
scripts/hfs_build_pseudobulk.py   per-atlas cell filters + grouping keys
scripts/hfs_tissue_specificity.py cross-tissue table (pooled + per-donor)
scripts/hfs_suspension_bias.py    nucleus/cell control + figure
scripts/hfs_gut_anomaly.py        gut compartment figure + rankings
```

Run order: download → build → specificity → (suspension_bias, gut_anomaly).

Two implementation notes worth keeping:

- **Blocks are sized by nonzero count, not row count.** That is what lets a
  10 GB atlas be aggregated on a 16 GB host without an out-of-core framework.
- **Dataset identity is verified against each file's own `uns["citation"]`**,
  not against the portal UI and not against transfer success. A completed
  download is not evidence of correct content; `hfs_download_atlases.py` aborts
  on a version mismatch.

## Findings

Numbers below are read back from the saved tables, not restated from memory.
23 whole-cell fibroblast populations + 2 single-nucleus.

**Load per receptor is highest in mucosal, cutaneous and synovial fibroblasts.**
Top five by `col6_per_antxr2`: gingiva 108, buccal mucosa 83, synovium (JIA) 80,
skin 52, hard palate 48. The next tissue down is trachea at 30, so there is a
visible gap between this group and everything else.

**Skin is the outlier on collagen composition**, at 0.65 COL6 share of total
collagen against a whole-cell median near 0.35. The composition metric and the
load metric therefore do *not* rank tissues the same way — gingiva leads on load
but sits low on share (0.15). Both are reported; neither is privileged.

**The gut is the inverse configuration, not a marginal anomaly.** Large
intestine is second-lowest of all 23 on load (8.5) while carrying high
*ANTXR2* (148 CPM, upper third of the panel); small intestine is similar
(16.2, 125 CPM). Whatever drives the intestinal phenotype, the collagen VI
clearance account does not explain it.

Within the gut the receptor is **stromal and glial, not epithelial** — every
epithelial label sits at the bottom of the per-cell-type ranking. This is in
tension with the mouse work placing the Wnt role in epithelium
cell-autonomously, and is the direct motivation for the inflamed-vs-healthy
epithelial arm: if the Wnt function is injury-conditional, a healthy atlas
cannot show it by construction.

## Standing corrections

**1. The organ-involvement enrichment does not survive a pre-specified test.**
Grouping tissues by HFS involvement class and testing core-affected against the
rest gives Mann–Whitney one-sided *p* = 0.11 on `col6_per_antxr2` (n = 5 vs 18)
and *p* = 0.51 on `col6_frac_of_collagen`. Against the "not reported" class
only: *p* = 0.13 and *p* = 0.67.

A much smaller *p* (3 × 10⁻⁵) is obtainable by grouping "mucocutaneous + joint"
against the rest — but that group is exactly the top five of the ranking being
tested, so it is selection on the outcome variable and **is not evidence.** Any
previously recorded significant enrichment for this table came from a post-hoc
grouping and should be treated as withdrawn.

The ranking itself is a real, reproducible ordering and is worth reporting as
**exploratory**. The claim "affected organs have higher collagen VI load per
receptor" is **not supported** at this sample size.

**2. The earlier claim that single-nucleus inflation is *ANTXR2*-specific was
wrong** and is withdrawn. It compared tendon nuclei to *other tissues'* cells,
confounding platform with tissue, and never asked whether *ANTXR2* moves more
than genes in general. The matched control (kidney, donor-restricted to donors
sampled both ways; heart, region-matched plus a single-donor pair; 9 strata)
shows the direction reproduces — nuclei recover relatively more *ANTXR2* and
relatively less collagen — but *ANTXR2* sits at the **53rd–87th percentile** of
the genome-wide shift, i.e. upper-middle, not the tail. The ratio is biased from
both ends by **2.2–15.7×**, dataset-dependent.

Consequence: the two tendon rows (Achilles 2.9, quadriceps 3.5) are
**unplaceable, not low.** Corrected for the observed range they span a band that
overlaps both the never-reported ceiling and the core-affected median, so tendon
cannot be ordered against the whole-cell tissues at all. They are excluded from
every cross-tissue comparison and retained only to make the exclusion explicit.

**3. Pooled tissue rows are cell-weighted, which deviates from the project's
donor-equal-weight convention.** `col6_*` on the pooled rows sums counts across
donors before normalising, so a high-cell-count donor contributes more. This was
chosen because several tissues have too few cells per donor to give stable
donor-level CPM. `hfs_col6_specificity_by_donor.csv` carries the per-donor rows
so the donor-mean-of-means version can be recomputed; the donor-level table is
the one to use for any claim where donor weighting could matter.

## Data gaps

Three organs with frequent HFS involvement have **no fibroblast single-cell data
in any atlas surveyed**: gingiva at 29 of 45 clinical reports is covered only by
the oral atlas (3′ v3 fibroblasts, three sites, small n), and **joint (28) and
bone (18) have none at all.** The synovium row is juvenile idiopathic arthritis,
not healthy tissue, and the HFS joint pathology is periarticular rather than
synovial — so it is an imperfect proxy and labelled as such in the figures (†).

Tendon is present but unplaceable (see correction 2). This is the standing
limitation of the arm: the tissues where the disease is most visible are the
ones the atlases cover worst.

## Provenance

All 11 datasets are CELLxGENE Discover, registered with dataset-version UUIDs in
`hfs_config.DATASETS` and verified against each downloaded file's `uns` citation.
Tabula Sapiens rows are pinned to 10x 3′ v3 so oral sites are comparable rather
than confounded with capture chemistry; the oral atlas is likewise 3′ v3 and
healthy-site only.
