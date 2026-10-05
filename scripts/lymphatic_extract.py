"""Pool lymphatic endothelial cells (LEC) from every local atlas into one
genome-wide population for the within-LEC ANTXR2 co-expression scan.

Sources and memento grouping label (cell_type column; groups are donor x
cell_type, so subtype / tissue / chemistry is fixed within every group --
pooling across subtypes inside a group is what produced the endothelial ECM
artefact):

  gut   gut.h5ad (Elmentaite 2021 Gut Cell Atlas), author_cell_type LEC1-LEC6,
        label = "<subtype>|<assay>" (the atlas mixes 3' v2 and 5' v2 within
        donors)
  ts    ts_endothelium.h5ad, cell_type == lymphatic, 10x 3' v3 only,
        label = tissue_in_publication; Stomach excluded (the whole stomach
        endothelium is labelled lymphatic but only 42% is PROX1/CCL21+)
  site  oral / synovium / tendon LECs sliced from comp_cell_extract.py output,
        label = site

Donor ids are prefixed with the source. Gene spaces differ between atlases
(32k / 61k / 35k features), so the pooled matrix is restricted to the gene
symbols present in every source.

Usage: python lymphatic_extract.py <workspace dir>
Writes <workspace>/lymphatic/lymphatic_pooled.npz and lymphatic_sources.csv.
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fib_cell_extract as F                       # noqa: E402

LEC = "endothelial cell of lymphatic vessel"
SITES = ["oral_buccal_mucosa", "oral_gingiva", "oral_hard_palate", "syn_synovium_JIA",
         "tach_tendon_Achilles", "tquad_tendon_quadriceps"]


def load_npz(path):
    z = np.load(path, allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return X, z["genes"].astype(str), z["total_counts"], z["donor"].astype(str), z["cell_type"].astype(str)


def extract_h5ad(path, mask, donor, label, tag, out_dir):
    F.extract(path, mask, np.array([tag] * len(mask)),
              {"donor": donor, "cell_type": label}, tag, out_dir=out_dir)
    return load_npz(os.path.join(out_dir, f"{tag}_{tag}.npz"))


def main(ws):
    out = os.path.join(ws, "lymphatic")
    os.makedirs(out, exist_ok=True)
    parts = {}

    gut = os.path.join(ws, "data", "gut.h5ad")
    with h5py.File(gut, "r") as f:
        ob = {c: F.rc(f["obs"], c) for c in ["author_cell_type", "donor_id", "assay"]}
    m = np.char.startswith(ob["author_cell_type"].astype(str), "LEC")
    lab = np.array([f"{s}|{a}" for s, a in zip(ob["author_cell_type"], ob["assay"])])
    parts["gut"] = extract_h5ad(gut, m, ob["donor_id"], lab, "gut", out)

    ts = os.path.join(ws, "data", "ts_endothelium.h5ad")
    with h5py.File(ts, "r") as f:
        ob = {c: F.rc(f["obs"], c) for c in ["cell_type", "assay", "tissue_in_publication", "donor_id"]}
    m = ((ob["cell_type"] == LEC) & (ob["assay"] == "10x 3' v3")
         & (ob["tissue_in_publication"] != "Stomach"))
    parts["ts"] = extract_h5ad(ts, m, ob["donor_id"], ob["tissue_in_publication"], "ts", out)

    # each site extract carries its own gene space; the intersection below
    # handles any differences between them
    for site in SITES:
        X, g, tot, d, ct = load_npz(os.path.join(ws, "comp_cells", f"endothelial_{site}.npz"))
        k = ct == LEC
        parts[f"site:{site}"] = (X[k], g, tot[k], d[k], np.array([site] * int(k.sum())))

    common = None
    for X, g, *_ in parts.values():
        u = pd.Index(g).drop_duplicates(keep=False)
        common = u if common is None else common.intersection(u)
    common = common.sort_values()
    mats, rows = [], []
    for src, (X, g, tot, d, ct) in parts.items():
        uniq = ~pd.Index(g).duplicated(keep=False)   # duplicated symbols are not in `common`
        pos = np.flatnonzero(uniq)[pd.Index(g[uniq]).get_indexer(common)]
        assert len(pos) == len(common)
        mats.append(X[:, pos])
        sname = src.split(":")[0]
        rows.append(pd.DataFrame({"source": sname, "donor": [f"{sname}:{x}" for x in d],
                                  "cell_type": ct, "total_counts": tot}))
    X = sp.vstack(mats).tocsr()
    obs = pd.concat(rows, ignore_index=True)
    np.savez(os.path.join(out, "lymphatic_pooled.npz"), data=X.data.astype(np.float32),
             indices=X.indices.astype(np.int32), indptr=X.indptr, shape=np.array(X.shape),
             genes=np.asarray(common, dtype=str), total_counts=obs.total_counts.values,
             donor=np.asarray(obs.donor, dtype=str), cell_type=np.asarray(obs.cell_type, dtype=str))
    s = obs.groupby(["source", "donor", "cell_type"]).size().rename("n_cells").reset_index()
    s.to_csv(os.path.join(out, "lymphatic_sources.csv"), index=False)
    print(f"{X.shape[0]} cells x {X.shape[1]} shared genes; "
          f"{(s.n_cells >= 100).sum()} donor x label groups with >=100 cells")
    print(s[s.n_cells >= 100].groupby("source").size().to_string())


if __name__ == "__main__":
    main(sys.argv[1])
