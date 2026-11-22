"""Stream the pilot's cell populations x ALL genes out of a source h5ad into one
CSR .npz per source (program_config.cells_npz):

    data, indices, indptr, shape, genes, total_counts,
    donor, cell_type (population name), fine (subtype label), study, age,
    region, assay, disease

Each source has an adapter (ADAPTERS below) that maps its own obs schema onto
those columns and names the populations: every author subtype of each
mesenchymal lineage, plus enterocytes. Subtype identity is handled at scoring
time by donor x subtype grouping (Decision 4) and by program_config.CONTRASTS,
not by narrowing the extraction.

    gca       Gut Cell Atlas releases (healthy, extended): author annotations,
              Elmentaite2021 SI only. enterocyte = level_3 "Enterocyte"
              (BEST4+ cells are a separate label and are excluded).
    ibdverse  standard CELLxGENE schema, ileum only (both disease states;
              the disease slice is applied at scoring time). disease is
              "<disease>|<author_inflammation_status>", e.g.
              "Crohn disease|uninflamed". fibroblast =
              "fibroblast" + "subepithelial intestinal fibroblast", matching
              the discovery's all-subtype fibroblast population;
              enterocyte = "enterocyte of epithelium proper of ileum".
    ibdverse_crypt  same file, ileum crypt epithelium only: crypt = the
              author's ileal Stem category (Enterocyte stem OLFM4lo LGR5hi,
              Enterocyte stem MKI67, Enterocyte progenitor crypt OLFM4hi),
              fine = author_cell_type. For the inflamed-vs-uninflamed state
              contrast (planning note, "State contrast plan").
    ibdverse_landscape  same file, ileum, the landscape populations
              (program_config.LANDSCAPE_POPS by author_cell_type, fine =
              author_cell_type; plus pooled fibroblasts, fine = ontology
              class), one npz per population in LANDSCAPE_CELLS.
    gca_landscape  GCA Extended+, postnatal small intestine, all studies and
              disease labels (the slice is applied at scoring), the
              replication populations (GCA_LANDSCAPE_POPS by level_3_annot),
              one npz per population in GCA_LANDSCAPE_CELLS.

One superset per source is extracted and slices (age, disease, region) are
applied at scoring time, so trying a slice never re-streams the h5ad.

Before extracting, the file's uns["citation"] is checked against the expected
dataset version and the matrix's indptr[-1] against its stored nnz -- a
completed download is not evidence of correct content.

Usage:
    python program_extract.py --source {healthy,extended,ibdverse,ibdverse_crypt,
                                        ibdverse_landscape,gca_landscape} [--dry-run]
"""
import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fib_cell_extract as F                       # noqa: E402
from program_config import (                       # noqa: E402
    DONOR_COL, FINE_COL, GCA_FIB_ALL, GCA_FIB_FETAL, GCA_FIB_LEVEL2, GCA_LANDSCAPE_CELLS,
    GCA_LANDSCAPE_POPS, IBDVERSE_LANDSCAPE_CELLS,
    LANDSCAPE_FIB, LANDSCAPE_FIB_TYPES, LANDSCAPE_POPS, MIN_GROUP_CELLS, REGION, SLICES, SOURCES,
    SKIN_ASSAYS, SKIN_FIB_CLUSTER, SKIN_LANDSCAPE_CELLS,
    STUDY, cells_npz, pop_slug,
)


def gca(f):
    cols = ["study", DONOR_COL, "level_2_annot", FINE_COL, "age_unified",
            "organ_groups", "organ_unified", "assay"]
    o = pd.DataFrame({c: F.rc(f["obs"], c) for c in cols})
    pops = {"fibroblast": ("level_2_annot", "Fibroblast"),
            "smooth_muscle": ("level_2_annot", "Smooth_muscle"),
            "pericyte": ("level_2_annot", "Pericyte"),
            "myofibroblast": ("level_2_annot", "Myofibroblast"),
            "enterocyte": (FINE_COL, "Enterocyte")}
    ages = sorted({a for s in SLICES.values() if s["source"] in ("healthy", "extended")
                   for a in s.get("ages", [])})
    base = ((o.study == STUDY) & o.age_unified.isin(ages) & (o.organ_groups == REGION)).values
    ct = np.full(len(o), "", dtype=object)
    for name, (col, label) in pops.items():
        ct[base & (o[col].values == label)] = name
    return ct, pd.DataFrame({"donor": o[DONOR_COL], "fine": o[FINE_COL], "study": o.study,
                             "age": o.age_unified, "region": o.organ_unified,
                             "assay": o.assay, "disease": "normal"})


def ibdverse(f):
    cols = ["donor_id", "cell_type", "tissue", "disease", "development_stage", "assay",
            "author_inflammation_status"]
    o = pd.DataFrame({c: F.rc(f["obs"], c) for c in cols})
    pops = {"fibroblast": {"fibroblast", "subepithelial intestinal fibroblast"},
            "myofibroblast": {"myofibroblast cell"},
            "pericyte": {"pericyte"},
            "smooth_muscle": {"vascular associated smooth muscle cell"},
            "enterocyte": {"enterocyte of epithelium proper of ileum"}}
    base = (o.tissue == "ileum").values
    ct = np.full(len(o), "", dtype=object)
    for name, labels in pops.items():
        ct[base & o.cell_type.isin(labels).values] = name
    return ct, pd.DataFrame({"donor": o.donor_id, "fine": o.cell_type, "study": "IBDverse",
                             "age": o.development_stage, "region": o.tissue,
                             "assay": o.assay,
                             "disease": o.disease + "|" + o.author_inflammation_status})


CRYPT_TYPES = ["Enterocyte stem OLFM4lo LGR5hi", "Enterocyte stem MKI67",
               "Enterocyte progenitor crypt OLFM4hi"]


def ibdverse_crypt(f):
    _, obs = ibdverse(f)
    act = pd.Series(F.rc(f["obs"], "author_cell_type"))
    ct = np.where(((obs.region == "ileum") & act.isin(CRYPT_TYPES)).values, "crypt", "")
    return ct.astype(object), obs.assign(fine=act.values)


def ibdverse_landscape(f):
    _, obs = ibdverse(f)
    act = pd.Series(F.rc(f["obs"], "author_cell_type")).values
    ileum = (obs.region == "ileum").values
    ct = np.full(len(obs), "", dtype=object)
    fib = ileum & obs.fine.isin(LANDSCAPE_FIB_TYPES).values
    ct[fib] = LANDSCAPE_FIB
    auth = ileum & pd.Series(act).isin(LANDSCAPE_POPS).values
    ct[auth] = act[auth]
    return ct, obs.assign(fine=np.where(auth, act, obs.fine.values))


def gca_landscape(f):
    cols = ["study", DONOR_COL, "level_2_annot", FINE_COL, "age_unified", "organ_groups",
            "organ_unified", "assay", "disease"]
    o = pd.DataFrame({c: F.rc(f["obs"], c) for c in cols})
    base = ((o.organ_groups == REGION)
            & ~o.age_unified.astype(str).str.contains("trimester")).values
    ct = np.where(base & o[FINE_COL].isin(GCA_LANDSCAPE_POPS).values,
                  o[FINE_COL].values, "").astype(object)
    ct[base & (o.level_2_annot == GCA_FIB_LEVEL2).values] = LANDSCAPE_FIB
    return ct, pd.DataFrame({"donor": o[DONOR_COL], "fine": o[FINE_COL], "study": o.study,
                             "age": o.age_unified, "region": o.organ_unified,
                             "assay": o.assay, "disease": o.disease})


def gca_fetal_fib(f):
    """Fetal small-intestine fibroblasts of Extended+ (same file and genes as
    gca_landscape); main() then writes the postnatal + fetal union."""
    cols = [DONOR_COL, "level_2_annot", FINE_COL, "age_unified", "organ_groups", "organ_unified",
            "study", "assay", "disease"]
    o = pd.DataFrame({c: F.rc(f["obs"], c) for c in cols})
    fetal = ((o.organ_groups == REGION) & (o.level_2_annot == GCA_FIB_LEVEL2)
             & o.age_unified.astype(str).str.contains("trimester")).values
    ct = np.where(fetal, GCA_FIB_FETAL, "").astype(object)
    return ct, pd.DataFrame({"donor": o[DONOR_COL], "fine": o[FINE_COL], "study": o.study,
                             "age": o.age_unified, "region": o.organ_unified,
                             "assay": o.assay, "disease": o.disease})


def skin_fib(f):
    """Human dermal fibroblast atlas: healthy skin, 10x 3' v2 / v3; every cell
    is a fibroblast. fine = one class (pooled within donor); the author
    cluster is kept in `cluster` for the by-cluster population."""
    cols = ["donor_id", "disease", "assay", "tissue", "Author", "Age", "cluster_robust"]
    o = pd.DataFrame({c: F.rc(f["obs"], c) for c in cols})
    keep = ((o.disease == "normal") & o.assay.isin(SKIN_ASSAYS)).values
    ct = np.where(keep, LANDSCAPE_FIB, "").astype(object)
    return ct, pd.DataFrame({"donor": o.donor_id, "fine": "skin fibroblast", "study": o.Author,
                             "age": o.Age, "region": o.tissue, "assay": o.assay,
                             "disease": o.disease, "cluster": o.cluster_robust})


def refine_npz(src, out, cell_type, fine_key):
    """Copy a population npz with `fine` replaced by another stored column."""
    z = dict(np.load(src, allow_pickle=False))
    z["fine"] = z[fine_key]
    z["cell_type"] = np.full(len(z["fine"]), cell_type)
    np.savez_compressed(out, **z)
    print(f"wrote {out}: fine = {fine_key}")


def merge_npz(paths, out, cell_type):
    """Row-concatenate per-population npz files written by F.extract."""
    import scipy.sparse as sp
    zs = [np.load(p, allow_pickle=False) for p in paths]
    assert all((z["genes"] == zs[0]["genes"]).all() for z in zs)
    X = sp.vstack([sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
                   for z in zs]).tocsr()
    rest = {k: np.concatenate([z[k].astype(str) if z[k].dtype.kind == "U" else z[k] for z in zs])
            for k in zs[0].files if k not in ("data", "indices", "indptr", "shape", "genes")}
    rest["cell_type"] = np.full(X.shape[0], cell_type)
    np.savez(out, data=X.data, indices=X.indices.astype(np.int32), indptr=X.indptr.astype(np.int64),
             shape=np.array(X.shape, dtype=np.int64), genes=zs[0]["genes"], **rest)
    print(f"wrote {out}: {X.shape[0]} cells, {X.nnz} nnz")


ADAPTERS = {"healthy": gca, "extended": gca, "ibdverse": ibdverse,
            "ibdverse_crypt": ibdverse_crypt, "ibdverse_landscape": ibdverse_landscape,
            "gca_landscape": gca_landscape, "gca_fetal_fib": gca_fetal_fib, "skin_fib": skin_fib}
# sources written as one npz per population instead of one per source
PER_POPULATION = {"ibdverse_landscape": IBDVERSE_LANDSCAPE_CELLS,
                  "gca_landscape": GCA_LANDSCAPE_CELLS, "gca_fetal_fib": GCA_LANDSCAPE_CELLS,
                  "skin_fib": SKIN_LANDSCAPE_CELLS}


def verify(f, version_id):
    cit = f["uns/citation"][()]
    cit = cit.decode() if isinstance(cit, bytes) else str(cit)
    assert version_id in cit, f"citation does not name {version_id}: {cit}"
    X = f["raw/X"] if "raw" in f and "X" in f["raw"] else f["X"]
    assert X["indptr"][-1] == X["data"].shape[0], "indptr[-1] != stored nnz (truncated?)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, choices=sorted(SOURCES))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    src = SOURCES[a.source]
    assert os.path.getsize(src["h5ad"]) == src["size"], "file size differs from the API's"

    with h5py.File(src["h5ad"], "r") as f:
        verify(f, src["version_id"])
        ct, obs = ADAPTERS[a.source](f)
    keep = ct != ""
    sel = obs[keep].assign(cell_type=ct[keep])
    g = sel.groupby(["disease", "cell_type", "fine", "donor"]).size()
    g = g[g >= MIN_GROUP_CELLS]
    print(f"{a.source} superset: {len(sel)} cells, {sel.donor.nunique()} donors; donors "
          f"with a >= {MIN_GROUP_CELLS}-cell group per subtype:")
    print(g.groupby(level=[0, 1, 2]).size().to_string())
    if a.dry_run:
        return

    obs_cols = {"cell_type": ct, **{k: obs[k].astype(str).values for k in obs.columns}}
    if a.source in PER_POPULATION:
        files = np.array([pop_slug(c) for c in ct], dtype=object)
        out_dir = PER_POPULATION[a.source]
    else:
        files = np.full(len(obs), a.source)
        out_dir = cells_npz(a.source).parent
    written = F.extract(str(src["h5ad"]), keep, files, obs_cols, "cells", out_dir=str(out_dir))
    for c, n, nnz, path in written:
        print(f"wrote {path}: {n} cells, {nnz} nnz")
    if a.source == "skin_fib":
        refine_npz(out_dir / f"cells_{pop_slug(LANDSCAPE_FIB)}.npz",
                   out_dir / f"cells_{pop_slug(SKIN_FIB_CLUSTER)}.npz", SKIN_FIB_CLUSTER, "cluster")
    if a.source == "gca_fetal_fib":
        merge_npz([out_dir / f"cells_{pop_slug(p)}.npz" for p in (LANDSCAPE_FIB, GCA_FIB_FETAL)],
                  out_dir / f"cells_{pop_slug(GCA_FIB_ALL)}.npz", GCA_FIB_ALL)


if __name__ == "__main__":
    main()
