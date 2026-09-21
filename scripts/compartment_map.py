"""cell_type -> compartment for the non-Tabula-Sapiens fig5 atlases.

The TS objects are already per-compartment, so only oral / synovium /
tendon_ach / tendon_quad need a mapping. Curated from the full label
vocabulary of each file (every label present is listed; nothing falls through
silently). 'other' = mural/neural/adipose/muscle/skeletal lineages that belong
to none of the four compartments under study.

Notes on the judgement calls:
  Langerhans cell      -> immune  (epidermal dendritic cell, myeloid lineage)
  Merkel cell          -> epithelial (epidermal mechanoreceptor, keratin+)
  myoepithelial cell   -> epithelial
  osteoclast           -> other   (myeloid-derived but functionally skeletal;
                                   53 cells, below the group floor anyway)
  stromal cell         -> other   (unresolved label, not callable as fibroblast)
  mural cell/pericyte  -> other   (perivascular, not endothelium)
"""

EPITHELIAL = {
    "keratinocyte", "acinar cell of salivary gland", "salivary gland glandular cell",
    "myoepithelial cell", "Merkel cell",
}

ENDOTHELIAL = {
    "blood vessel endothelial cell", "endothelial cell of lymphatic vessel",
    "endothelial cell", "endothelial cell of venule", "endothelial cell of arteriole",
    "capillary endothelial cell", "endothelial cell of vascular tree",
}

IMMUNE = {
    "CD4-positive helper T cell", "CD8-positive, alpha-beta T cell", "plasma cell",
    "mast cell", "macrophage", "natural killer cell", "B cell", "dendritic cell",
    "neutrophil", "Langerhans cell",
    "CD8-positive, alpha-beta memory T cell", "CD4-positive, alpha-beta memory T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell", "classical monocyte",
    "naive B cell", "gamma-delta T cell", "CD1c-positive myeloid dendritic cell",
    "regulatory T cell", "mucosal-associated invariant T cell",
    "CD141-positive myeloid dendritic cell", "plasmacytoid dendritic cell",
    "precursor B cell", "T cell", "Peyer's patch macrophage",
    "CD34-negative, CD117-positive innate lymphoid cell, human",
    "group 2 innate lymphoid cell", "group 3 innate lymphoid cell",
    "granulocyte", "monocyte", "mononuclear phagocyte", "leukocyte",
}

FIBROBLAST = {"fibroblast"}

OTHER = {
    "mural cell", "melanocyte", "Schwann cell", "adult skeletal muscle myoblast",
    "slow muscle cell", "fast muscle cell", "blood vessel smooth muscle cell",
    "adipocyte", "skeletal muscle fiber", "chondrocyte",
    "skeletal muscle satellite cell", "pericyte", "neural cell", "stromal cell",
    "osteoblast", "osteoclast",
}

BY_COMPARTMENT = {"epithelial": EPITHELIAL, "endothelial": ENDOTHELIAL,
                  "immune": IMMUNE, "fibroblast": FIBROBLAST, "other": OTHER}

LOOKUP = {ct: comp for comp, s in BY_COMPARTMENT.items() for ct in s}


def compartment_of(labels):
    """Vectorised lookup; unseen labels raise so nothing is silently dropped."""
    missing = sorted({x for x in labels if x not in LOOKUP})
    if missing:
        raise KeyError(f"unmapped cell_type labels: {missing}")
    return [LOOKUP[x] for x in labels]
