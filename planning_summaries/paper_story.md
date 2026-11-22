# Paper story: ANTXR2 as a human-specific interferon-responsive gene

Synthesis written 2026-11-21 at the user's request, pulling together
`ANALYSIS_SUMMARY.md` (durable framing + partner-availability work) and
`program_coexpression.md` (the co-expression landscape through the interferon
work). Nothing here is a new analysis except the ortholog-flag check in §2.1.
Numbers are quoted from those two documents; their provenance is recorded there.

---

## 1. The claim

> **ANTXR2/CMG2 is a weak type I interferon-stimulated gene in human cells, and
> that responsiveness is not detectable in any other species tested — including
> mouse, in which essentially all of CMG2's mechanistic biology was
> established. In human tissue, ANTXR2's expression and its cell-to-cell
> covariation track type I (not γ) interferon tone, rising with inflammation
> across gut epithelial populations.**

Three qualifiers in that sentence are load-bearing and belong in the abstract:
**weak**, **type I**, and **human**.

**"Weak" is a measurement, not a hedge** — and saying so plainly is what makes
the rest of the evidence internally consistent. In Shaw's human fibroblasts,
the same experiment that puts ANTXR2 at +1.13 log2 (2.2×) puts **MX1 at +9.89
log2, ~950×** — a roughly 430-fold difference in response size. In the Hagai
bulk, ANTXR2 is +0.53 log2 while MX1 / IFIT1 / RSAD2 run log2 5–9. ANTXR2 is a
real but small responder, not a canonical ISG. Every other observation lines up
with that reading rather than against it, including the ChIP-Atlas result that
looked like a problem under a stronger claim (§6, item 3).

The coding sequence is conserved and the regulation is not: dN/dS ≈ 0.2 and
86–91% protein identity across mammals, against induction in human alone. That
asymmetry is the paper — a well-studied receptor carrying a human-specific
regulatory input that the standard model organism cannot report on.

### The three acts (figures in §3)

1. **ANTXR2 expression is ubiquitous and consistent with the known HFS
   mechanism — but that is not the full story.** The cross-tissue expression
   and ratio map sets up a partition that expression alone cannot complete.
2. **ANTXR2 participates in gene-regulatory programs in a context- and
   cell-type-specific, partially replicable way**, with collagen as the
   pertinent negative and interferon as the pertinent positive.
3. **That interferon positive is a weak, human- and type I-specific ISG
   signal.**

Two framings follow, and they are what make the claim worth publishing rather
than merely true:

1. **Model-organism blind spot.** Bracq 2025 (DSS regeneration), Bürgi 2017
   (collagen VI clearance in uteri) and the rest of the CMG2 literature are
   mouse. If interferon drives ANTXR2 in human and not in mouse, an entire
   input to human CMG2 biology is invisible by construction in those models.
2. **A state axis for a tissue-restricted disease.** HFS's germline LoF is in
   every cell but the phenotype is not. The project has been looking for the
   partition that explains the gap (`ANALYSIS_SUMMARY.md`, "Central question").
   Interferon supplies a *state* axis rather than another cell-type axis:
   ANTXR2 is highest where and when interferon tone is high.

Framing 1 is the stronger of the two and should lead. Framing 2 is a discussion
point, not a result — see §5.

### 1.1 The inferential chain, and why the order matters

The paper's logic runs **observational discovery → stimulation confirmation →
species scope**, and the first arrow is load-bearing:

1. **Discovery, observational.** A program-level screen in human gut tissue
   (~600 programs, noise-corrected coupling against 200 span-matched anchors)
   put interferon at the top of ANTXR2's fibroblast profile — the one thing
   that sits at the top in *both* gut atlases. Nothing told the screen to look
   at interferon.
2. **That nomination is what licensed the stimulation data.** IFN-β and the
   HTEC α/β/γ/λ panel were pulled in *because* the observational screen named
   type I interferon. They were not the starting point.
3. **Species scope** then came from asking whether the induction found in (2)
   generalizes — and it does not, which is the headline.

**Consequence for how the stimulation experiments are read.** They are
confirmatory tests of a pre-specified program set, not fresh screens. They do
not have to re-discover interferon as ANTXR2's top-ranked program, and they do
not have to show interferon outranking every co-activated pathway. The stimulus
is a known, pure cytokine; the set was fixed in advance. Demanding exclusivity
there would also demand the one thing the method provably cannot deliver
(§4.1): within-cell coupling cannot name the upstream signal, so a stimulation
experiment that co-activates NF-κB is not thereby uninformative about
interferon.

This is also the answer to "it's just a co-expression screen." The screen is
the hypothesis generator; the causal input comes from recombinant cytokine, and
the cross-species comparison comes from paired designs. Each step does the job
the previous step could not.

---

## 2. Evidence, strongest first

### 2.1 Cross-species induction (strongest; two independent datasets)

**Shaw et al. 2017** (PLoS Biol; `isg.data.cvr.ac.uk`), fibroblasts from 10
species treated with type I interferon. ANTXR2 log2FC (FDR):

| Species | log2FC | FDR | Call |
|---|---|---|---|
| **Human** | **+1.13** | **5.1e-18** | **ISG** |
| Rat | +0.10 | 0.36 | — |
| Cow | −0.50 | 0.034 | down |
| Sheep | +0.02 | 0.96 | — |
| Pig | +0.11 | 0.40 | — |
| Horse | −0.33 | 0.052 | — |
| Dog | +0.04 | 0.90 | — |
| Flying fox | −0.30 | 0.15 | — |
| Little brown bat | −0.04 | 0.91 | — |
| Chicken | +0.02 | 0.91 | — |

**The nine nulls are informative, not detection failures** (checked
2026-11-21 against the raw webserver JSON, `shaw2017/q_ANTXR2.json`): ANTXR2 is
flagged `present=True` **and** `expressed=True` in all 10 species, all in one
ortholog cluster (`HS9951`), protein identity to human 86–91% across the
mammals (69% chicken). This was the obvious hole in the claim and it closes.

**The same JSON carries the key supporting point.** dN/dS is 0.18–0.26 in every
species that reports it — strong purifying selection on the protein. So:
conserved coding sequence, human-restricted regulation. That is the sentence
the abstract is built around, and it also pre-empts the first reviewer
objection ("maybe it's just a poorly conserved gene").

**Hagai et al. 2018 bulk** (ArrayExpress E-MTAB-5919; our reanalysis, paired
untreated vs IFN-β 8 h per donor, median-of-ratios): human +0.53 log2, up in
**6/6 donors** (+0.39 to +0.66, i.e. 1.3–1.6×); rhesus −0.11 (6/6 slightly
down), mouse −0.19 (3 donors), rat +0.08 (3). Independent lab, independent
species panel, paired design — and it adds **macaque and mouse**, which Shaw
lacks. Mouse and rat agree between the two datasets.

Two datasets, 12 species between them, one species up.

### 2.2 Cytokine identity: type I, not γ (strong)

- **HTEC** (Kim 2024 / memento, GSE274751; primary human tracheal epithelium,
  2 donors, α/β/γ/λ × 3–48 h + unstimulated). Ciliated cells, mean ANTXR2 vs
  control: α 1.2–2.2×, β 1.25–1.5× at every time point; **γ and λ ~1**. Basal
  cells go the *other* way under γ (0.2–0.5×). Coupling: 8/8 IFN programs
  q<0.05 under α (excess 0.13–0.18) and β (0.15–0.27); γ ≈ 0; λ same direction,
  n.s.
- **Mouse colon organoids + IFN-γ** (GSE255366): no induction (per-10k
  0.79/0.89/1.77 stem, no consistent direction), no coupling (nothing q<0.05).
  Consistent with the γ-null, though it is mouse and therefore also consistent
  with §2.1 — it is not a clean mouse type-I test.
- **Human dermal fibroblasts + IFN-β** (Hagai E-MTAB-7051, Smart-seq2, 5
  donors, 681 cells, 2/6 h): at 6 h, 7/8 IFN programs q<0.05 (excess 0.15–0.22,
  anchor percentile 99–100); type I-biased set 0.20 (q 0.002), with the type I
  residual I|G 0.17 (q 0.001) against G|I 0.05 (n.s.). 2 h: nothing.
  TNF/NF-κB programs also score high here (CytoSig TNFA 0.28, PROGENy NFkB
  0.27). That does not weaken the test, for two reasons. The stimulus is a
  known, pure cytokine, so nothing about the inference depends on interferon
  outranking other programs — and the IFN set was pre-specified from the
  observational discovery, not read off this ranking. IFN-β engaging NF-κB arms
  is expected biology. The response is also **specific rather than a global
  shift**: 37 of 605 programs reach q<0.05, median excess 0.012 — against
  **285 of 605**, median 0.107, in the dsRNA experiment (poly I:C 6 h). That
  contrast is precisely why dsRNA was dropped as interferon evidence and this
  was kept.
- **Between-sample means name the cytokine** where within-cell coupling cannot
  (§4.1). Regressing log2 ANTXR2 per 10k on standardized type I and γ set
  scores, donor-clustered SEs: HTEC ciliated (control samples) type I **+0.24**
  (p<1e-4) vs γ +0.06 (n.s.); IBDverse crypt OLFM4lo type I +0.14 (p<0.001;
  inflamed +0.45, uninflamed +0.28, normal +0.07) with γ ≈ 0 in every state;
  LGR5 stem type I +0.20 (p 0.004), γ −0.07; inflamed top-villus type I +0.36
  vs γ +0.09.

So: α/β induce and couple; γ does neither, in two systems and two species.

### 2.3 Human tissue covariation (moderate — this is context, not the hook)

IBDverse ileum, ANTXR2's coupling to 8 pre-specified IFN programs (excess over
200 span-matched anchors, cells resampled within donor × subtype), by state
(normal / CD-uninflamed / CD-inflamed):

| Population | Normal | Uninflamed | Inflamed |
|---|---|---|---|
| Fibroblasts | 0.12 | 0.18 | 0.22 |
| Top-villus enterocytes | 0.045 | 0.077 | 0.137 |
| Crypt OLFM4lo progenitors | 0.024 | 0.048 | 0.127 |
| Goblet MUC2hi | 0.029 | 0.046 | 0.058 |
| MKI67 stem | 0.012 | 0.032 | 0.094 |
| LGR5 stem | −0.011 | −0.008 | −0.028 |

8/8 programs q<0.05 in fibroblasts, top-villus and crypt in all three states;
nothing negative anywhere; **it rises with inflammation in every population**.
Mean level moves much less (ANTXR2 per 10k, same ordering of states):
top-villus 0.86/0.89/1.17, crypt 0.82/0.83/1.07, fibroblasts 2.30/2.46/2.35,
goblet 1.44/1.37/1.26, LGR5 stem 0.53/0.50/0.56. Report both — the coupling
result is the larger effect and the mean is the more interpretable one.

**Controls that this survived:**
- *Not immune contamination* (2026-11-21). Removing cells with ≥2 pan-immune
  UMI, or each donor's top quarter by ambient plasma-cell Ig, leaves the
  coupling intact in all 9 population × state cells tested (e.g. fibroblast
  inflamed 0.261 → 0.277 / 0.290; crypt inflamed 0.224 → 0.193 / 0.228).
  ANTXR2's own coupling to the contamination programs is ~0 (−0.04 to +0.07).
  Crypt-inflamed cells carry the most immune reads (4.1 UMI/cell → 0.1 after
  filtering) and keep 0.19.
- *Not simply "any gene of this size."* Everything is excess over 200 anchors
  matched on genomic span, raw mean and detection. The span match matters and
  was a real correction: ANTXR2 spans 224 kb and intronic reads create a
  length axis that inflated its apparent coupling before the fix
  (`program_coexpression.md`, "Correction 2026-11-03").
- *Among ECM-role positive controls, ANTXR2 leads on IFN in gut.* Of 7 genes
  (MRC2, ITGA11, DDR2, FAP, PDGFRB, LRP1, ANTXR1), ANTXR2 is the highest on
  IFN programs in gut fibroblasts (pct 93–99).

Replication is partial and must be stated as such: GCA goblet 0.156 (8/8),
GCA fibroblasts weak (0.052, n.s.), skin fibroblast A3 0.046 (5/8). ANTXR2's
program profile *as a whole* does not transfer across datasets (fibroblast
Spearman 0.07–0.12) — **interferon is the one thing that sits at the top in
both gut atlases**, which is the honest way to put it.

### 2.4 Promoter (weakest; supporting only)

- Motif scan, −1 kb..+0.5 kb around the Ensembl gene start: human ANTXR2
  carries a full canonical ISRE at **+72** (`AGTTTCTTTTTCT`; consensus
  `AGTTTCNNTTTCN`) plus a core match at −179. Mouse *Antxr2* has none in the
  same window.
- ChIP-Atlas (hg38, ±5 kb of TSS): weak STAT1 (20% of 103 experiments, 58th
  percentile), IRF1 (14% of 56), IRF9 (1 of 56) — **GAPDH-level, not
  ISG-level** (MX1, IFIT1 sit at the 99th percentile).

Caveats are serious: one window around one annotated start, AT-rich human
sequence (GC 0.39 vs mouse 0.53) inflating match probability, and the same scan
*misses* MX1's known ISREs because they fall outside the window. No literature
on ANTXR2 promoter interferon regulation was found.

**Use this as a hypothesis generator and the motivation for the reporter
experiment in §4, not as evidence.** The ChIP-Atlas occupancy is *concordant*
with a ~2× responder rather than contradicting it (§6, item 3); what stays open
is whether the +72 element is the operative one.

### 2.5 What the co-expression landscape gave, and did not

Worth including as a figure because it establishes where the gene is and
constrains what the interferon result means — but the project should be candid
that most of it was negative.

**Held up (replicated in two gut atlases):** ANTXR2 dips in dividing
progenitors (excess −0.80 to −0.96 log2 vs stem, 0–2nd anchor percentile) and
rises in differentiated enterocytes (+0.63 to +1.53, 82–96th percentile), on
both the between-population ladder and within-population maturation quartiles.
Mean-level counterpart of a co-expression tilt away from E2F/MYC/mTORC1 that
survives dropping MKI67+ cells, i.e. it tracks a graded growth axis rather
than division.

**Did not hold up:**
- *Wnt.* H1 (crypt Wnt-output coupling higher in inflamed, from Bracq's model)
  not supported — point estimate the wrong sign in both references. H2
  untestable: the fetal-like/revival state Bracq's KO blocks is not sampled in
  these annotated crypt populations (3/9 YAP genes clear the floor).
- *ECM/collagen.* Within a single fibroblast subtype there is no detectable
  excess coupling to collagens or ECM programs in gut or skin. The earlier
  negative COL6 signal was between-subtype composition (gut COL6A1 −0.178 →
  −0.087 n.s. after conditioning; skin cluster-robust ≈ 0). This is a
  *consistent* result, not a disappointing one — Bürgi showed collagen VI mRNA
  does not change, so transcriptional coupling was never the expected readout
  (`ANALYSIS_SUMMARY.md` standing correction).
- *ANTXR2's profile is not distinctive overall.* It is no more structured than
  a span-matched gene's; its identity is *which* programs sit at the extremes.

---

## 3. Figure plan

Three figures (restructured 2026-11-21, user). Each act does what the previous
one could not, which is also the inferential chain of §1.1. Everything below
exists already; nothing in the plan requires new data.

### Figure 1 — Expression is ubiquitous and fits the known HFS mechanism, but does not explain it

*Theme: necessary, insufficient. This figure earns the pivot to co-expression.*

| Panel | Content |
|---|---|
| 1a | ANTXR2 across 238 cell types, 3 collections, donor-equal-weight means (`antxr2_celltype_table.csv`) |
| 1b | Cross-tissue fibroblast expression: gut vs skin vs other stroma |
| 1c | Ratios as raw ingredients side by side — COL6/ANTXR2, ANTXR1/ANTXR2, MRC2 |
| 1d | MRC2 at the detection floor in gut epithelium: the one partition that does work |
| 1e | **Within-tissue state structure.** ANTXR2 along the gut differentiation ladder: down in dividing progenitors (excess −0.80 to −0.96 log2 vs stem, 0–2nd anchor percentile), up in differentiated enterocytes (+0.63 to +1.53, 82–96th) — replicated in both atlases, on the between-population ladder *and* within-population maturation quartiles |

Claims: the stromal/mesenchymal enrichment replicates the classic CMG2
picture; terminally differentiated barrier epithelia (keratinocyte, corneal,
RPE) sit at the bottom; gut epithelium is the exception, upper-middle. The
receptor/substrate ratio differs sharply between gut fibroblast subtypes
(COL6/ANTXR2 ≈ 29 in Fibroblast CPM vs 6.5 in LP ADAMDEC1).

Caveats that belong **on the figure**, not only in the legend:
- **Donor-equal-weight aggregation and `n_donors`.** The literal #1 rank is 2
  donor-units / 1,415 cells, while the well-powered generic "enterocyte" label
  ranks 150/238 at ~63× lower. Thin labels must be marked as thin.
- **Exon-only vs intron-counted atlases.** ANTXR2 spans 224 kb and reads ~25×
  lower in an exon-only atlas; the per-gene IBDverse/GCA ratio tracks genomic
  span at Spearman 0.70–0.74. Every cross-platform comparison carries this.
- **No composite vulnerability score, ever** (standing project rule).
- **Co-presence is not redundancy.** Skin fibroblasts show robust ANTXR1
  co-presence *and* HFS's most visible pathology; the paralogs bind different
  collagen VI domains.

Why it is not the full story, in two steps — and 1e is what makes the second
step possible:

1. Expression is broad and fairly flat across tissues that do and do not
   develop HFS pathology, so a presence/absence map cannot partition the
   phenotype. MRC2 is the only clean absence, and it is in gut epithelium.
2. *Within* a tissue, ANTXR2 is not flat at all. It is graded along
   differentiation — ~1.6× between LGR5 stem (0.53 per 10k) and top-villus
   enterocytes (0.86) in normal tissue, reaching +0.63 to +1.53 log2 above
   matched genes on the ladder — variation that a cell-type-level average
   discards entirely.

So the gene's variation lives at the level of **cell state**, which is what
Act 2 is built to read. 1e is the last panel of Act 1 and the reason Act 2
exists.

### Figure 2 — Context- and cell-type-specific, partially replicable program participation

*Theme: what ANTXR2 travels with, where, and how reliably. Collagen is the
pertinent negative; interferon is the pertinent positive.*

| Panel | Content |
|---|---|
| 2a | **The same gradient, now in co-expression.** Within epithelial populations ANTXR2 couples *negatively* to E2F / MYC / G2M programs (ranks 1–2 of 51 matched genes), and the tilt survives removing MKI67+ cells (E2F −1.17 → −1.10, MYC V1 −1.02 → −0.95; 20 of 23 epithelial cells still negative, median z −1.08) — a graded growth axis, not a division artefact |
| 2b | Replicability, stated honestly: donor split-half reliability by population — skin fibroblasts 0.58, IBDverse fibroblasts 0.21–0.46, GCA postnatal fibroblasts 0.04; cross-dataset epithelial profile Spearman 0.51 vs fibroblast 0.07–0.12 |
| 2c | The landscape: ~600 programs × population × state (IBDverse, 59 of 81 cells filled), excess over 200 span-matched anchors |
| 2d | **Pertinent negative — collagen.** Within-subtype null in gut and skin; COL6A1 −0.178 → −0.087 (n.s.) after conditioning on subtype; skin cluster-robust ≈ 0 |
| 2e | The bar is beatable: 7 ECM-role control genes, 5 of 7 clear their anchors on collagens in skin — so ANTXR2's collagen null is a real null, not a power failure |
| 2f | **Pertinent positive — interferon.** IFN coupling by population × state, rising with inflammation in every population; LGR5 stem ≈ 0 |
| 2g | Contamination controls: immune-doublet and ambient-Ig filters leave it intact; ANTXR2's coupling to the contamination programs ≈ 0 |

**2a is 1e seen the other way, and the legend must say so.** The two are *not*
independent evidence — most likely one maturation gradient read in two
registers, a mean difference and a within-population covariation. Stating that
is what keeps the pair honest; splitting them across the figure break displays
the non-independence instead of hiding it inside a single panel. 2a also earns
Act 2 its credibility alongside 2b, which admits where the landscape is unreliable: it is
the cleanest replication in the paper, holding in both atlases and in both
registers.

Scope note: the same growth tilt appears in lymphocytes, more strongly than in
epithelium. Lymphocytes are out of scope by standing project decision, so the
MKI67-negative control above is quoted on epithelium only.

The two "pertinent" results are what make this a result rather than a survey.
The collagen negative is **predicted**: Bürgi showed collagen VI mRNA does not
change in Antxr2−/− uteri, so transcriptional coupling was never the expected
readout — and 2e demonstrates the assay could have found it had it been there.
The interferon positive is the only thing at the top of the fibroblast profile
in *both* gut atlases, and it is what nominates Figure 3.

### Figure 3 — A weak, human- and type I-specific ISG signal

| Panel | Content |
|---|---|
| 3a | HTEC dose/time: ANTXR2 mean under α / β / γ / λ × 3–48 h — α 1.2–2.2×, β 1.25–1.5×, γ and λ ≈ 1; basal cells *down* under γ |
| 3b | HTEC coupling: 8/8 IFN programs q<0.05 under α and β; γ ≈ 0; λ same direction, n.s. |
| 3c | Type I vs γ from between-sample means: HTEC ciliated +0.24 vs +0.06; IBDverse crypt +0.14 (inflamed +0.45), LGR5 stem +0.20, inflamed top-villus +0.36, γ ≈ 0 throughout — with the mouse colon organoid IFN-γ null alongside |
| 3d | Cross-species forest plot, Shaw 10 species, annotated with present/expressed, % identity and dN/dS; human alone up |
| 3e | Hagai paired bulk, 4 species: human +0.53 in 6/6 donors; rhesus, mouse, rat not |
| 3f | **Scale panel.** ANTXR2 against MX1 / IFIT1 / RSAD2 in the same experiments: Shaw human ANTXR2 +1.13 vs MX1 +9.89 log2; Hagai ANTXR2 +0.53 vs ISGs log2 5–9 |
| 3g | Promoter: human-vs-mouse alignment showing the +72 ISRE, beside the ChIP-Atlas occupancy panel (ANTXR2 at GAPDH level, MX1/IFIT1 at the 99th percentile) |

**Panel 3f is new in this plan and should not be cut.** It converts "weak" from
a hedge into a number, and it is what makes 3g coherent instead of
contradictory: a ~2× gene should not show ISG-level factor occupancy, and the
~950× gene should. Occupancy and fold change agree, so those two panels belong
side by side.

### Dropped from the earlier plan

The old Figure 5 (wet-lab validation: human vs mouse fibroblasts ± IFN-β, ISRE
reporter) is **not** in this plan. That is consistent: the hedged claim — a
weak, possible, human-specific ISG signal, concordant across two archives and
12 species — is a defensible computational result, and the experiment is what
would remove the word "possible" rather than what licenses publication
(§4, item 1).

Methods figure / supplement: the program-coupling estimator, anchor matching,
the gene-span correction, and the negative result in §4.1.

---

## 4. What is missing, ranked

**1. Human-vs-mouse induction, side by side.** Human and mouse primary
fibroblasts (or matched lines), ± IFN-β, time course, ANTXR2 by qPCR or
RNA-seq. Under the hedged claim this is **no longer make-or-break**: "a weak,
possible, human-specific ISG signal, concordant across two independent archives
spanning 12 species" is a defensible computational result, and the figure plan
has no experimental panel. What the experiment buys is deleting the word
*possible* — it converts agreement between two reanalyses into a direct
observation, and a ~2× effect is exactly the size a reviewer will want shown
rather than inferred. A macaque or chimp arm would settle human-alone vs
primate-wide at the same time (§6, item 5).

**2. The ISRE.** Reporter assay with and without the +72 element, and/or
STAT1/IRF9 ChIP-qPCR at the human ANTXR2 promoter. This is what turns a motif
scan into a mechanism, and it is the natural resolution of the ChIP-Atlas
tension (§6, item 3). A mouse-sequence reporter in human cells (and vice
versa) would be the cleanest demonstration that the *element*, not the cellular
context, carries the species difference.

**3. A human gut type I stimulation dataset.** The tissue arm (§2.3) is
observational and the stimulation arm (§2.2) is airway epithelium and dermal
fibroblasts. The gap is a gut cell type stimulated with type I interferon.
No public dataset with direct recombinant type I/III treatment of human gut
cells was found; the closest are **GSE171620** (Triana 2021, human ileal
organoids, astrovirus mock/4 h/16 h — endogenous type I/III response) and
GSE156760 (SARS-CoV-2, ISGs in bystander cells). GSE171620 was offered and
never run; it is the cheapest remaining item on this list.

**4. Functional consequence.** Nothing in this work shows that interferon-raised
ANTXR2 *does* anything — not collagen VI uptake, not Wnt transduction, not
anything. The paper can be honest that it establishes regulation and not
function, but a reviewer will ask, and one experiment (e.g. collagen VI
internalization ± IFN-β in human fibroblasts) would pre-empt it.

**5. Loose ends.** The IBDverse Ileum goblet MUC2lo BCAS1int AssertionError
(unresolved, 3 states); E-MTAB-7167 alignment for unstimulated dermal
fibroblast controls (**not recommended** — different study and library prep,
FASTQ only, so the comparison would be confounded by prep).

---

## 5. Stories considered and rejected

- **"Partner availability explains HFS tissue restriction"** — *not rejected,
  relocated.* The project's original frame, and under the three-act structure
  it becomes **Act 1**: the expression and ratio map is the setup, and MRC2's
  absence in gut epithelium is its one clean partition. What it cannot do is
  finish the argument, and that failure is the pivot the rest of the paper
  turns on. The specific mechanistic tests did fail — the Wnt arm failed its
  pre-registered test, and the regenerative state Bracq's KO concerns is not
  sampled in these annotations — so this frames and motivates rather than
  concludes.
- **A methods paper on program coupling.** There is a genuine methodological
  result here (§4.1 below, plus the span-matching correction, the resampling-layer
  calibration work, and the mixed-model comparison). But it belongs in the
  memento branch's own write-up or as this paper's methods supplement, not as
  the headline — the biology is the more interesting half.
- **"ANTXR2 is a stromal gene."** Ruled out early and standing policy:
  ANTXR2 is not stroma-specific, and the compartment ordering even flips
  between exon-only and intron-counted atlases.

### 4.1 One methodological result worth stating as a result

**Within-cell coupling cannot identify the upstream cytokine; between-sample
means can.** The positive control proves it: in HTEC ciliated cells under
*pure IFN-β*, the coupling sits on the γ-residual program (0.172, q 0.03), not
the type I residual (0.065, n.s.) — because γ-biased genes still vary
cell-to-cell with the type I response. The between-sample regression on the
same data recovers the right answer (type I +0.24, γ +0.06). Anyone using
single-cell co-expression to infer which signal is driving a gene needs this,
and it is a clean, falsifiable negative result with a worked positive control.

---

## 6. What a reviewer will attack

1. **"It's all reanalysis."** True as of today. §4 item 1 is the answer and
   there is no way around it.
2. **"Two donors in HTEC."** The cytokine-identity arm rests partly on a
   2-donor experiment. Mitigated by the dermal fibroblast IFN-β result and the
   tissue between-sample means pointing the same way, but it should be stated,
   not finessed.
3. **"Your own ChIP data argue against an ISRE mechanism."** ANTXR2's STAT1
   occupancy is GAPDH-level while MX1/IFIT1 are at the 99th percentile. This
   was previously flagged as the story's central tension; **under the "weak"
   framing it is not a tension but a consistency check that passes.** A gene
   induced ~2× should not show ISG-level factor occupancy, and MX1 — induced
   ~950× in that same Shaw experiment — should. Occupancy and fold change point
   the same way, which is why panels 3f and 3g belong on one figure. What stays
   genuinely open is whether the +72 element is the operative one: a reporter
   question (§4 item 2), not a contradiction.
4. **"Why does the gut tissue coupling rise with inflammation if the mean
   barely moves?"** Fair, and the honest answer is that the coupling and the
   mean are different statistics with different power here. Report both
   everywhere (standing project rule) and do not imply a large expression
   change in inflamed tissue.
5. **"Human-specific how — primate-specific, or human alone?"** Shaw has no
   macaque; Hagai's rhesus is flat (−0.11, 6/6 slightly down), which argues for
   human alone rather than primate. That rests on one dataset and 6 animals, so
   phrase it as "not seen in rhesus macaque" rather than "human-unique among
   primates." A chimp or macaque fibroblast experiment would settle it and is a
   cheap addition to §4 item 1.
6. **Intronic reads.** ANTXR2 is 224 kb and most of its IBDverse signal is
   pre-mRNA; the gene reads ~25× lower in an exon-only atlas. Span-matched
   anchors make the within-atlas comparisons fair, but any cross-platform
   statement needs this caveat attached.
7. **"NF-κB outranks interferon in your stimulated fibroblasts."** It does
   (0.27–0.28 vs 0.15–0.22), and the answer is that exclusivity was never the
   claim. The stimulus is recombinant IFN-β, the IFN set was pre-specified from
   the observational screen, and IFN-β engages NF-κB arms — so co-movement is
   expected, not confounding. The specificity evidence is the program count:
   37 of 605 at q<0.05 (median excess 0.012) under IFN-β versus 285 of 605
   (median 0.107) under dsRNA, which is why dsRNA was dropped and this was not.
   Pointing at the ranking also asks the within-cell estimator to identify the
   upstream cytokine, which §4.1 shows it cannot do for anyone.

---

## 7. Suggested abstract shape

Mapped to the three acts, with the hedges kept.

1. CMG2/ANTXR2 is a collagen VI receptor and Wnt co-receptor; loss causes HFS,
   whose pathology is tissue-restricted; its mechanistic biology is almost
   entirely mouse-derived.
2. *(Act 1)* Across 238 human cell types ANTXR2 is broadly expressed, with the
   stromal enrichment its collagen VI role predicts — but expression is too
   ubiquitous to partition the phenotype, and MRC2 is the only clean absence.
3. *(Act 2)* In a ~600-program co-expression screen across gut and skin,
   ANTXR2's program participation is context- and cell-type-specific and only
   partially replicable. It does **not** couple to collagen programs within
   fibroblast subtypes — as expected, since collagen VI mRNA is unchanged in
   Antxr2−/− tissue — while interferon programs are the one positive shared
   across gut atlases, rising with inflammation.
4. *(Act 3)* That signal is a weak type I interferon response: IFN-α/β, but not
   IFN-γ, induce ANTXR2 ~1.3–2.2× in primary human airway epithelium and dermal
   fibroblasts — some ~400× less than MX1 in the same experiments.
5. Across 12 species in two independent datasets the induction is restricted to
   human, while the protein is under strong purifying selection (dN/dS ≈ 0.2,
   86–91% identity): a human-specific regulatory input to a conserved receptor.
6. A canonical ISRE present in the human and absent in the mouse promoter
   nominates a mechanism, at occupancy levels consistent with a weak responder.
7. These results identify a regulatory layer of human CMG2 biology that the
   standard model organism cannot report on.
