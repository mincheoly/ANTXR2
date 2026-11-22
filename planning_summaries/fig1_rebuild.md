# Figure 1 rebuild and the intron/exon comparability problem (2026-11-22)

## State

Figure 1 (a-e) is rebuilt from atlases on disk; it replaces the old genome-wide
memento means (`/data/ANTXR2/celltype_expression`, lost).

| Step | Script | Output (`/data/antxr2/fig1/`) |
|---|---|---|
| per-cell panel extract (62 genes + library size), streamed from h5ad | `scripts/fig1_extract.py` | `cells/<atlas>.parquet`, `.meta.json` |
| donor x cell-type units | `scripts/fig1_aggregate.py` | `donor_units.parquet` |
| tables (platform test, 1a, 1b/c, 1d) | `scripts/fig1_tables.py` | `tables/*.csv` |
| figure | `scripts/fig1_panels.py` | `panels/fig1.png` |

Atlases (`/data/antxr2/raw/`): Tabula Sapiens v2 All Cells (CELLxGENE
version 1c88f927..., 45 GB, identity verified), cross-tissue immune Global (3 GB,
not used in any panel), GCA Healthy reference + Extended+, IBDverse, dermal
fibroblast atlas. Run with `/home/ubuntu/my_python/bin/python` from `scripts/`.
No `is_primary_data` filter: GCA Healthy reference is entirely `False` (its
primary copies are in "Extended - All genes"); atlases are never pooled.

## The comparability problem

Absolute ANTXR2 and every ratio with ANTXR2 in it differ by quantification mode.

- Intron-inclusive: Tabula Sapiens v2 (Cell Ranger 7.0.1, "default parameters"
  -> introns counted), IBDverse (Cell Ranger 7.2.0).
- Exon-only: GCA (Oliver et al. 2024 remapped everything with STARsolo and made
  both Gene and GeneFull, but the released objects are `gene_cellbender`, i.e.
  STARsolo Gene after CellBender -- confirmed in Teichlab/PanGIAtlas
  `notebooks/1_QC_and_integration/1_make_healthy_reference_nh3.ipynb`); the
  dermal fibroblast atlas (empirically, no long-gene excess).
- Measured (`tables/platform.csv`, fibroblasts vs GCA 10x 3' v3): log-ratio vs
  genomic span Spearman 0.52 (TS 3' v3), 0.51 (IBDverse), ~0 for GCA 3' v2 and
  skin; ANTXR2 (224 kb) at the 91st-95th percentile of the excess. TS 5' v2
  only 0.24 (5' reads land at the transcript start).
- Consequence seen in 1a: immune cells outrank stroma only in TS 3' v3; stroma is
  above immune in TS Smart-seq2, TS 5' v2 and GCA.
- No released GeneFull / spliced-unspliced matrices exist for GCA (checked the
  Pan-GI page objects and the 77 GB objects zip index). TS v2 BAMs and
  molecule_info.h5 are in `s3://czb-tabula-sapiens/TabulaSapiens_v2/<donor>/alignment_gencode41/`
  but downloads need TS data-access credentials (listing is public, GET is 403).

## Next: re-align a GCA subset to measure the factor

Goal: same reads counted exon-only (Gene) and intron-inclusive (GeneFull), so
the inflation factor is measured per gene x cell type, not modelled.

1. Pick adult GCA healthy-reference samples from Elmentaite 2021 (ArrayExpress
   E-MTAB-9543, adult GEX; related E-MTAB-9536/9532/9533/10386) that carry
   fibroblasts, epithelium and myeloid cells; prefer 10x 3' v2 (GCA's main
   chemistry) plus one 3' v3 and one 5' v2 if available. Match by `sampleID` in
   the GCA obs.
2. Reference: 10x `refdata-gex-GRCh38-2020-A` (what the GCA STARsolo run matched).
   STAR genome generate needs ~32 GB RAM (this is why the host is being resized).
3. STARsolo (>= 2.7.10a), Cell Ranger-matching settings:
   `--soloType CB_UMI_Simple --soloCBwhitelist <737K-august-2016 (v2) | 3M-february-2018 (v3)>
   --soloUMIlen <10 v2 | 12 v3> --soloFeatures Gene GeneFull Velocyto
   --soloCellFilter EmptyDrops_CR --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR
   --clipAdapterType CellRanger4 --outFilterScoreMin 30
   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outSAMtype None`
   (5' v2: `--soloStrand Reverse`, read 1 cDNA handling per chemistry).
4. Map barcodes to GCA cell annotations (obs index / sampleID), then per gene x
   cell type: factor = GeneFull / Gene (pooled UMIs, donor-equal weight).
5. Checks: (a) Gene counts reproduce the released GCA counts for the same cells
   (before CellBender, so expect close, not identical); (b) log(GeneFull/Gene)
   vs genomic span reproduces the TS-vs-GCA slope; (c) apply the factors to put
   GCA on the intron-inclusive scale (or divide TS/IBDverse by them) and re-run
   `fig1_tables.platform()`: the span correlation should go to ~0; then see
   whether immune still outranks stroma.

Also differs between atlases (legend, not figure): GCA and IBDverse counts are
CellBender-corrected; TS raw counts are not.
