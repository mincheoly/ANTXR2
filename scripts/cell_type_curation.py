"""Curated row lists and disease-tissue annotations for the ECM-clearance-pathway
heatmap (see plot_ecm_clearance_heatmap.py). Pure data, no I/O -- kept separate
from query/plotting code so the row selection (the single biggest judgment call
in that figure) is easy to review and edit on its own.

Disease context: ANTXR2 (CMG2) loss causes hyaline fibromatosis syndrome (HFS,
including the severe infantile-systemic-hyalinosis/ISH form); its paralog ANTXR1
(TEM8) loss causes GAPO syndrome -- a different clinical phenotype entirely.
Categories below are grouped by which disease affects that tissue and how
(nodular/fibrotic vs. diffuse/functional matters -- these are different kinds of
lesions, not the same pathology at different severity; see COVERAGE_GAP_NOTE and
the GI category's own note).
"""

# cell_type (whole-body atlas) -> category. A dict (not lists-of-lists) so every
# row has exactly one category -- e.g. `fibroblast` sits only in the HFS-nodular
# bucket even though it's also the literature-reported GAPO mechanistic proxy
# (GAPO fibroblasts show disrupted actin cytoskeleton / reduced ECM turnover);
# noted in the figure caption rather than duplicated as a second row.
WHOLE_BODY_CATEGORIES = {
    # HFS: discrete papular/nodular skin, gingival, and perianal lesions --
    # matrix-accumulation pathology, myofibroblasts are the classic fibrosis
    # effector cell.
    "fibroblast": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "fibroblast of gingiva": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "myofibroblast cell": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "adventitial cell": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "keratinocyte": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "melanocyte": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "sebocyte": "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    # HFS: chronic diarrhea + protein-losing enteropathy, with histologically
    # documented *diffuse* (non-nodular) hyaline deposition in the gut wall --
    # a functional/barrier phenotype, not a matrix-accumulation lesion the way
    # skin/gingiva are. Thin site-specific splits (duodenum/ileum/jejunum
    # enterocyte, 1-4 donors each) deliberately excluded; only well-powered
    # generic labels shown.
    "enterocyte": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "colonocyte": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "intestine goblet cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "intestinal crypt stem cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "paneth cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "intestinal tuft cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "M cell of gut": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    # DROPPED: "gastrointestinal tract (lamina propria) macrophage" (CL:0000865).
    # It is not a macrophage. All 42,302 cells carrying that label across the
    # entire Gut Cell Atlas are authored level_1=Mesenchymal, level_2=Fibroblast,
    # level_3=Lamina_propria_fibroblast_ADAMDEC1, and markers side decisively
    # with the authors: PTPRC 0.4% positive, CD68 2.5%, LYZ 0.5%, C1QA 0.5%,
    # against COL1A1 89%. The CELLxGene ontology mapping is wrong, systematically,
    # across all 10+ constituent studies. Removed rather than relabelled: the
    # generic `fibroblast` row already represents that compartment here, and a
    # second fibroblast row under a GI heading would imply a distinct GI-resident
    # population the label cannot support.
    "interstitial cell of Cajal": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "enteric neuron": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "enteroglial cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    # HFS (severe/ISH form specifically): additional hyaline deposition reported
    # in muscle, lymph node, spleen -- fibroblastic reticular cell / follicular
    # dendritic cell are lymphoid-organ stromal proxies (no dedicated "lymph
    # node"/"spleen" cell type exists in this atlas). Kept in the same category
    # as the GI rows above (rather than a 4th category) to stay within the
    # dataviz skill's validated 3-hue all-pairs-safe categorical set for a
    # row-annotation strip whose adjacency isn't preserved once rows are
    # reordered by clustering -- both are the "diffuse/systemic" side of HFS
    # pathology, as opposed to the nodular/fibrotic side above.
    "smooth muscle cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "skeletal muscle satellite stem cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "fibroblastic reticular cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "follicular dendritic cell": "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    # GAPO: progressive optic atrophy/glaucoma/keratopathy (eye proxies, not
    # optic-nerve-specific -- no such cell type in this atlas) and dilated
    # scalp veins (vascular).
    "corneal epithelial cell": "GAPO-relevant",
    "retinal pigment epithelial cell": "GAPO-relevant",
    "retinal blood vessel endothelial cell": "GAPO-relevant",
    "endothelial cell": "GAPO-relevant",
    # Comparison / neither disease -- baseline for reading the rest of the
    # figure. tissue-resident macrophage included deliberately: CTSB/CTSK/
    # LAMP1/MRC2 are canonically myeloid/lysosomal-high, so this row shows that
    # baseline explicitly rather than letting high clearance-gene signal
    # elsewhere be misread as disease-tissue-specific.
    "tissue-resident macrophage": "Comparison / neither disease",
    "hepatocyte": "Comparison / neither disease",
    "pancreatic acinar cell": "Comparison / neither disease",
    "naive B cell": "Comparison / neither disease",
    "neuron": "Comparison / neither disease",
}

WHOLE_BODY_CATEGORY_ORDER = [
    "HFS - nodular/fibrotic (skin/gingiva/perianal)",
    "HFS - diffuse/systemic (GI + muscle + lymphoid organs)",
    "GAPO-relevant",
    "Comparison / neither disease",
]

# Category -> (light hex, dark hex). Validated all-pairs-safe (both modes) via
# the dataviz skill's palette validator -- the skill's categorical slots 1-3
# (blue, orange, aqua) are the specific triple it documents as all-pairs-safe;
# reused here rather than an arbitrary 3-hue pick. "Comparison" is deliberately
# neutral gray, not a 4th hue -- it's the de-emphasized baseline, not a category
# competing for identity with the three disease-relevant hues.
WHOLE_BODY_CATEGORY_COLORS = {
    "HFS - nodular/fibrotic (skin/gingiva/perianal)": ("#2a78d6", "#3987e5"),
    "HFS - diffuse/systemic (GI + muscle + lymphoid organs)": ("#eb6834", "#d95926"),
    "GAPO-relevant": ("#1baf7a", "#199e70"),
    "Comparison / neither disease": ("#85837c", "#8f8d83"),
}

# All 13 fibroblast_subtype values in the skin atlas.
SKIN_SUBTYPES = [
    "F1: Superficial",
    "F2: Universal",
    "F2/3: Stroma_PPARG+",
    "F3: FRC-like",
    "F4: DP_HHIP+",
    "F4: DS_DPEP1+",
    "F4: TNN+COCH+",
    "F5: NGFR+",
    "F5: RAMP1+",
    "F6: Myofibroblast",
    "F6: Inflammatory myofibroblast",
    "F7: Fascia-like myofibroblast",
    "F_Fascia",
]

# Disease-emergent annotation, checked against the published paper (Steele et al.
# 2025) rather than left as a naming-pattern guess: the paper states three
# fibroblast populations have NO healthy-skin counterpart -- inflammatory
# myofibroblasts, (plain) myofibroblasts, and fascia-like myofibroblasts (the
# paper's own F6/F7/F8; our downloaded object uses a different F-number scheme
# for the same three groups, mapped below).
#
# `F_Fascia` is left "unresolved", not coerced to either side: the paper states
# healthy fascial fibroblasts were merged into F2: Universal rather than kept as
# their own group, yet our object has a distinct `F_Fascia` label with its own
# cells (2,978 cells, 3 donors) -- a real discrepancy between the object's
# labeling and the (secondary, not primary-supplementary-table) source used to
# check it. Don't silently resolve this either direction.
SKIN_DISEASE_STATUS = {
    "F1: Superficial": "healthy",
    "F2: Universal": "healthy",
    "F2/3: Stroma_PPARG+": "healthy",
    "F3: FRC-like": "healthy",
    "F4: DP_HHIP+": "healthy",
    "F4: DS_DPEP1+": "healthy",
    "F4: TNN+COCH+": "healthy",
    "F5: NGFR+": "healthy",
    "F5: RAMP1+": "healthy",
    "F6: Myofibroblast": "disease-emergent",
    "F6: Inflammatory myofibroblast": "disease-emergent",
    "F7: Fascia-like myofibroblast": "disease-emergent",
    "F_Fascia": "unresolved",
}

# status -> (light hex, dark hex). Only one real hue (violet, "disease-emergent")
# against two neutral grays -- an "emphasis" pattern (one accent + de-emphasis
# gray), not a 3-way categorical identity problem, so no all-pairs validation is
# needed here (trivial to separate a single hue from grays). "unresolved" gets
# its own lighter gray step, distinct from "healthy", so genuine uncertainty
# about F_Fascia isn't visually collapsed into a confident "healthy" reading.
SKIN_STATUS_COLORS = {
    "healthy": ("#85837c", "#8f8d83"),
    "disease-emergent": ("#4a3aa7", "#9085e9"),
    "unresolved": ("#c7c5be", "#4a4a46"),
}

COVERAGE_GAP_NOTE = (
    "Not represented in this atlas (no matching Cell Ontology label): bone/"
    "cartilage/joint/synovium (GAPO metaphyseal dysplasia; HFS joint "
    "contractures and osteopenia/osteolysis), teeth/odontogenic tissue (GAPO "
    "pseudoanodontia), thyroid and adrenal gland (HFS-ISH severe-form hyaline "
    "deposition sites), gonadal tissue (GAPO hypogonadism). Eye/vascular rows "
    "shown are proxies, not optic-nerve- or scalp-vein-specific."
)

SCORE_CAVEAT_NOTE = (
    "Descriptive juxtaposition of raw expression values and independently-"
    "sourced disease-tissue annotations; not a composite risk, redundancy, or "
    "predictive score."
)

SCALE_CAVEAT_NOTE = (
    "Color scales are independent per panel (own log floor, own colorbar); "
    "compare patterns within a panel, not color intensity across panels."
)

NORMALIZATION_NOTE = (
    "Cell color is min-max normalized per gene (0=that gene's lowest value "
    "shown, 1=its highest), so uniformly-high genes like LAMP1 don't dominate "
    "the color scale and drown out other genes' own cross-cell-type pattern. "
    "In-cell numbers (skin panel) and the CSV's mean_expression/"
    "log10_mean_expression columns are the real, unnormalized values -- only "
    "cell color is rescaled; see the CSV's normalized_0_1 column for the "
    "displayed color value itself."
)
