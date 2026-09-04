"""Carve a three-group Elmentaite2021 working set out of the Gut Cell Atlas and
attach two independently-derived capture-rate columns, as preparation for
single-cell co-expression analysis.

This is preparation only -- no correlation is computed here. Each run writes
/data/ANTXR2/coexpression/README_<variant>.md with the full rationale; see
config.py's COEXPR_* / Q_CHEM_* constants for the parameters.

Two variants, differing only in which annotation column defines the groups:
    celltype -- CELLxGene-harmonised `cell_type` (pools author subtypes)
    level3   -- the authors' own `level_3_annot` (narrower, cleaner populations)
See config.COEXPR_VARIANTS. Both write their own h5ad; neither overwrites the
other, and the variant is recorded in uns['provenance'].

Reading follows the established pattern from compute_means.py: raw h5py +
anndata.io.sparse_dataset()/read_elem(), never anndata.read_h5ad(backed="r").
The source file's obsp holds a 1.6M x 1.6M neighbour graph -- the exact
structure that OOM-killed the Tabula Sapiens run -- and backed mode does not
keep it lazy. We touch only obs/var/raw.X/obsm.

Usage:
    python subset_coexpression_data.py [--variant {celltype,level3}] [--dry-run]
"""
import argparse
import os
import time

import anndata
import anndata.io as anndata_io
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from capture_rate import (
    chem_capture_rate, gate_masks, median_umi_by_gate, pbmc_scaled_capture_rate,
)
from compute_means import available_mem_gb, peak_rss_gb
from config import (
    ASSUMED_SATURATION, COEXPR_DIR, COEXPR_IMMUNE_LEVEL1, COEXPR_STUDY,
    COEXPR_VARIANTS, DONOR_UNIFIED_COL, GUT_ATLAS_H5AD, PBMC_GATES,
    PBMC_REFERENCE_DIR, PBMC_REFERENCES, Q_CHEM_BY_ASSAY, Q_CHEM_DEFAULT,
)

RAW_SLOT = "raw/X"
GENE_NAME_COL = "feature_name"

# Markers used only for the post-hoc sanity report -- these are what caught the
# "gastrointestinal tract (lamina propria) macrophage" curation error, so the
# check earns its place rather than being decorative.
SANITY_MARKERS = ["ANTXR2", "ANTXR1", "MRC2", "COL1A1", "EPCAM", "PTPRC", "C1QA", "CD68"]


def read_markers_chunked(sparse_ds, positions, gene_idx, chunk=20000):
    """Read only total-UMI and a handful of marker columns for `positions`.

    The immune compartment used for capture-rate estimation is ~84k cells, and
    we need just 5-6 gene columns plus row sums from it -- materialising the
    full rows would cost several GB for no reason. CSR row-slicing returns whole
    rows regardless, so we chunk and discard.

    Returns (total_umi, {gene: counts}) as dense 1-D arrays over `positions`.
    """
    total = np.zeros(len(positions), dtype=np.float64)
    cols = {g: np.zeros(len(positions), dtype=np.float64) for g in gene_idx}
    for start in range(0, len(positions), chunk):
        sl = slice(start, min(start + chunk, len(positions)))
        block = sparse_ds[positions[sl]]
        block = block.tocsr() if sp.issparse(block) else sp.csr_matrix(block)
        total[sl] = np.asarray(block.sum(axis=1)).ravel()
        for g, i in gene_idx.items():
            cols[g][sl] = np.asarray(block[:, i].todense()).ravel()
        del block
    return total, cols


def gate_masks_from_columns(marker_cols, n_cells, gates):
    """Same gating logic as capture_rate.gate_masks, but over pre-extracted
    marker columns rather than a full matrix. Kept in step with that function
    deliberately -- both express `all positive markers > 0 and all negative
    markers == 0`."""
    out = {}
    for name, spec in gates.items():
        needed = list(spec["pos"]) + list(spec["neg"])
        if any(g not in marker_cols for g in needed):
            continue
        mask = np.ones(n_cells, dtype=bool)
        for g in spec["pos"]:
            mask &= marker_cols[g] > 0
        for g in spec["neg"]:
            mask &= marker_cols[g] == 0
        out[name] = mask
    return out


def load_pbmc_reference(assay):
    """Read a 10x public PBMC reference into (counts_csr, var_names, total_umi)."""
    spec = PBMC_REFERENCES[assay]
    path = os.path.join(PBMC_REFERENCE_DIR, spec["filename"])
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} missing -- run download_pbmc_reference.py first")

    if spec["kind"] == "h5":
        ad = sc.read_10x_h5(path)
    elif spec["kind"] == "mtx_tar":
        import tarfile
        extract_dir = os.path.join(PBMC_REFERENCE_DIR, spec["name"] + "_extracted")
        if not os.path.isdir(extract_dir):
            with tarfile.open(path) as tf:
                tf.extractall(extract_dir)
        # the tarball nests as filtered_gene_bc_matrices/<genome>/
        mtx_dir = None
        for root, _dirs, files in os.walk(extract_dir):
            if any(f.startswith("matrix.mtx") for f in files):
                mtx_dir = root
                break
        if mtx_dir is None:
            raise RuntimeError(f"no matrix.mtx found under {extract_dir}")
        ad = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols")
    else:
        raise ValueError(f"unknown reference kind {spec['kind']!r}")

    ad.var_names_make_unique()
    counts = ad.X.tocsr() if sp.issparse(ad.X) else sp.csr_matrix(ad.X)
    total = np.asarray(counts.sum(axis=1)).ravel()
    return counts, list(ad.var_names), total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="select cells and report, but don't write the h5ad")
    ap.add_argument("--variant", default="celltype", choices=sorted(COEXPR_VARIANTS),
                    help="which annotation column defines the groups (see config.COEXPR_VARIANTS)")
    args = ap.parse_args()

    variant = COEXPR_VARIANTS[args.variant]
    label_col = variant["label_col"]
    labels = variant["labels"]
    out_h5ad = variant["h5ad"]

    t0 = time.time()
    print(f"=== variant: {args.variant} (group by {label_col!r}) -> {out_h5ad}")
    print(f"=== source: {GUT_ATLAS_H5AD} ===")
    f = h5py.File(GUT_ATLAS_H5AD, "r")
    obs = anndata_io.read_elem(f["obs"])
    var = anndata_io.read_elem(f["var"])
    raw_X = anndata_io.sparse_dataset(f[RAW_SLOT])
    print(f"  full shape: {raw_X.shape}, obs cols: {len(obs.columns)}")

    study = obs["study"].astype(str).values
    group_label = obs[label_col].astype(str).values
    mask = (study == COEXPR_STUDY) & np.isin(group_label, labels)
    positions = np.where(mask)[0]
    print(f"  selected {len(positions):,} cells "
          f"(study={COEXPR_STUDY}, {label_col} in {labels})")
    for ct in labels:
        n = int(((study == COEXPR_STUDY) & (group_label == ct)).sum())
        print(f"    {ct:12s} {n:7,}")
    if len(positions) == 0:
        raise SystemExit(f"no cells selected -- check COEXPR_STUDY / {label_col} labels")

    sub_obs = obs.iloc[positions].copy()
    # Drop unused categories on every categorical. Not cosmetic: unused
    # categories inherited from the full 1.6M-cell file leak into memento's
    # create_groups as phantom empty groups (documented in compute_skin_means.py).
    n_dropped = 0
    for col in sub_obs.columns:
        if isinstance(sub_obs[col].dtype, pd.CategoricalDtype):
            before = len(sub_obs[col].cat.categories)
            sub_obs[col] = sub_obs[col].cat.remove_unused_categories()
            n_dropped += before - len(sub_obs[col].cat.categories)
    print(f"  dropped {n_dropped} unused categories across categorical columns")

    print(f"  reading raw counts... (peak_rss={peak_rss_gb():.1f}GB, "
          f"avail={available_mem_gb():.1f}GB)")
    counts = raw_X[positions]
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    print(f"  counts: {counts.shape}, nnz={counts.nnz:,}, "
          f"peak_rss={peak_rss_gb():.1f}GB")

    obsm = {}
    for k in f["obsm"].keys():
        arr = anndata_io.read_elem(f["obsm"][k])
        obsm[k] = np.asarray(arr)[positions]
        print(f"  obsm/{k}: {obsm[k].shape}")

    var_names = list(var[GENE_NAME_COL].astype(str).values)
    # NB: the file stays open past this point -- the PBMC-scaled capture rate
    # below reads marker columns for this study's immune compartment, which is
    # a different (larger) set of rows than the working set. Closed at the end.
    raw_X_open = raw_X

    # --- capture rate column 1: chemistry x sequencing saturation -------------
    print(f"\n=== capture rate 1: chemistry x saturation "
          f"(ASSUMED_SATURATION={ASSUMED_SATURATION}) ===")
    assays = sub_obs["assay"].astype(str).values
    chem_q = np.zeros(len(sub_obs), dtype=float)
    chem_detail = {}
    for a in sorted(set(assays)):
        q, lam, det = chem_capture_rate(a, ASSUMED_SATURATION, Q_CHEM_BY_ASSAY, Q_CHEM_DEFAULT)
        chem_q[assays == a] = q
        chem_detail[a] = {"q": q, "lambda": lam, "detected_fraction": det,
                          "q_chem": Q_CHEM_BY_ASSAY.get(a, Q_CHEM_DEFAULT)}
        print(f"  {a:12s} q_chem={chem_detail[a]['q_chem']:.4f} "
              f"lambda={lam:.3f} detected={det:.4f} -> q={q:.5f}  "
              f"(naive product would be {chem_detail[a]['q_chem']*ASSUMED_SATURATION:.5f})")

    # --- capture rate column 2: PBMC-scaled -----------------------------------
    #
    # Estimated from this study's IMMUNE compartment, not from the working set:
    # the working set is fibroblast/enterocyte/macrophage and holds essentially
    # no lymphocytes, so T/B gates applied to it would match only doublets and
    # ambient RNA. Capture rate is a property of the assay, so we measure it
    # where matched populations exist and apply it to the cells we keep.
    print("\n=== capture rate 2: PBMC-scaled (matched immune gates) ===")
    total_umi = np.asarray(counts.sum(axis=1)).ravel()

    gate_genes = sorted({g for spec in PBMC_GATES.values()
                         for g in list(spec["pos"]) + list(spec["neg"])})
    name_to_idx = {g: i for i, g in enumerate(var_names)}
    missing = [g for g in gate_genes if g not in name_to_idx]
    if missing:
        raise SystemExit(f"gate markers missing from source var: {missing}")
    gene_idx = {g: name_to_idx[g] for g in gate_genes}

    level1 = obs["level_1_annot"].astype(str).values
    imm_pos = np.where((study == COEXPR_STUDY) & np.isin(level1, COEXPR_IMMUNE_LEVEL1))[0]
    print(f"  immune reference compartment: {len(imm_pos):,} cells "
          f"({', '.join(COEXPR_IMMUNE_LEVEL1)})")
    imm_total, imm_cols = read_markers_chunked(raw_X_open, imm_pos, gene_idx)
    imm_assays = obs["assay"].astype(str).values[imm_pos]
    imm_masks = gate_masks_from_columns(imm_cols, len(imm_pos), PBMC_GATES)
    print(f"  gate sizes (whole immune compartment): "
          + ", ".join(f"{g}={int(m.sum()):,}" for g, m in sorted(imm_masks.items())))

    pbmc_q = np.zeros(len(sub_obs), dtype=float)
    pbmc_detail = {}
    for a in sorted(set(assays)):
        if a not in PBMC_REFERENCES:
            print(f"  {a}: no chemistry-matched PBMC reference -- skipping")
            continue
        ref_counts, ref_var, ref_total = load_pbmc_reference(a)
        ref_masks = gate_masks(ref_counts, ref_var, PBMC_GATES)
        ref_umi = median_umi_by_gate(ref_total, ref_masks)

        in_assay_imm = imm_assays == a
        ours_masks_a = {g: (m & in_assay_imm) for g, m in imm_masks.items()}
        our_umi = median_umi_by_gate(imm_total, ours_masks_a)
        gate_n = {g: int((m & in_assay_imm).sum()) for g, m in imm_masks.items()}

        q_ref = Q_CHEM_BY_ASSAY.get(a, Q_CHEM_DEFAULT)
        q, ratios = pbmc_scaled_capture_rate(our_umi, ref_umi, q_ref)
        pbmc_q[assays == a] = q
        pbmc_detail[a] = {"q": q, "q_ref": q_ref, "ratios": ratios,
                          "our_umi": our_umi, "ref_umi": ref_umi,
                          "gate_n": gate_n,
                          "reference": PBMC_REFERENCES[a]["name"]}
        print(f"  {a:12s} vs {PBMC_REFERENCES[a]['name']}")
        for g in sorted(ratios):
            print(f"      {g:9s} n={gate_n[g]:7,}  ours={our_umi[g]:8.0f}  "
                  f"ref={ref_umi[g]:8.0f}  ratio={ratios[g]:.3f}")
        print(f"      -> median ratio {np.median(list(ratios.values())):.3f} "
              f"x q_ref {q_ref:.4f} = q {q:.5f}")

    f.close()

    if (pbmc_q == 0).any():
        # Any assay without a matched reference would leave zeros, which memento
        # cannot use. Fail loudly rather than writing an unusable column.
        missing = sorted(set(assays[pbmc_q == 0]))
        raise SystemExit(f"capture_rate_pbmc unset for assays {missing} -- "
                         "add a chemistry-matched reference or handle explicitly")

    # Canonical donor key. `donor_id` is passed through unmodified from the
    # source (so nothing about the raw file is silently rewritten), but it is
    # NOT safe to group on: it merges two people under `A34 (417C)` and splits
    # three others across several strings. `donor` carries the resolved
    # identity under a name that does not invite reaching for it by reflex --
    # use this one, and `donorID_unified` remains available as the source.
    #
    # Note the Phase 1 parquet resolves this differently: there `donor_id`
    # itself holds the resolved identity, because it is a pipeline output
    # rather than a passthrough subset. Same name, different convention in the
    # two artifacts -- hence this explicitly-named column here.
    sub_obs["donor"] = sub_obs[DONOR_UNIFIED_COL].astype(str)

    sub_obs["capture_rate_chem"] = chem_q
    sub_obs["capture_rate_pbmc"] = pbmc_q
    for name, arr in [("capture_rate_chem", chem_q), ("capture_rate_pbmc", pbmc_q)]:
        if not ((arr > 0).all() and (arr < 1).all()):
            raise SystemExit(f"{name} outside (0,1) -- memento asserts max<1")

    # --- sanity report --------------------------------------------------------
    # obs['n_counts'] is carried from the source but is NOT exactly the row sum
    # of the matrix we ship: it runs slightly higher (median ~2 counts, ~0.04%),
    # because it was computed before this release dropped genes -- the release is
    # titled "18485 genes" but var holds 18,370, and feature_is_filtered is
    # all-False so it isn't a masking artifact. Immaterial in size, but reported
    # so nobody later assumes the two are interchangeable. Everything in this
    # script that needs depth uses matrix row sums, never n_counts.
    src_n_counts = sub_obs["n_counts"].values.astype(float)
    d = src_n_counts - total_umi
    print(f"\n  note: obs['n_counts'] exceeds matrix row sums by median "
          f"{np.median(d):.1f} counts ({np.median(d/src_n_counts)*100:.3f}%), "
          f"max {d.max():.0f} ({(d/src_n_counts).max()*100:.1f}%) -- see README")

    print("\n=== marker sanity by cell type (raw mean / % positive) ===")
    gi = {g: i for i, g in enumerate(var_names)}
    ct_sub = sub_obs[label_col].astype(str).values
    rows = []
    for ct in labels:
        m = ct_sub == ct
        s = counts[m]
        rec = {label_col: ct, "n_cells": int(m.sum()),
               "medUMI": float(np.median(total_umi[m]))}
        for g in SANITY_MARKERS:
            if g not in gi:
                continue
            col = np.asarray(s[:, gi[g]].todense()).ravel()
            rec[f"{g}_mean"] = round(float(col.mean()), 4)
            rec[f"{g}_pct"] = round(float((col > 0).mean() * 100), 1)
        rows.append(rec)
    rep = pd.DataFrame(rows).set_index(label_col)
    pd.set_option("display.width", 250)
    print(rep[["n_cells", "medUMI"] + [f"{g}_mean" for g in SANITY_MARKERS if g in gi]].to_string())
    print()
    print(rep[[f"{g}_pct" for g in SANITY_MARKERS if g in gi]].to_string())

    # Donor identity uses donorID_unified, NOT donor_id. In this atlas `donor_id`
    # is wrong both ways: `A34 (417C)` covers two different people (D11, 55-74y;
    # D12, 18-34y), and three people appear under several ids (D2 -> T036 /
    # T036NEG / T036POS). Counting by donor_id would both merge donors and
    # pseudo-replicate them -- and memento's 2D path bootstraps over donors.
    print("\n=== donor overlap (obs['donor'] = donorID_unified) ===")
    donor_sets = {ct: set(sub_obs.loc[ct_sub == ct, "donor"].astype(str))
                  for ct in labels}
    raw_sets = {ct: set(sub_obs.loc[ct_sub == ct, "donor_id"].astype(str))
                for ct in labels}
    for ct in labels:
        print(f"  {ct:24s} {len(donor_sets[ct]):3d} donors "
              f"(raw donor_id would say {len(raw_sets[ct])})")
    shared = set.intersection(*donor_sets.values())
    shared_raw = set.intersection(*raw_sets.values())
    print(f"  shared across all three: {len(shared)} "
          f"(raw donor_id would say {len(shared_raw)})")
    n_u, n_r = sub_obs["donor"].nunique(), sub_obs["donor_id"].nunique()
    print(f"  total: {n_u} donors ({n_r} distinct donor_id strings)")

    print("\n=== chemistry split ===")
    print(pd.crosstab(sub_obs[label_col].astype(str), sub_obs["assay"].astype(str)).to_string())

    if args.dry_run:
        print(f"\n[dry run] not writing. elapsed {time.time()-t0:.0f}s")
        return

    # --- write ----------------------------------------------------------------
    os.makedirs(COEXPR_DIR, exist_ok=True)
    adata = anndata.AnnData(X=counts, obs=sub_obs, var=var.copy(), obsm=obsm)
    adata.uns["provenance"] = {
        "source_h5ad": GUT_ATLAS_H5AD,
        "source_slot": RAW_SLOT,
        "study": COEXPR_STUDY,
        "variant": args.variant,
        "label_col": label_col,
        "labels": list(labels),
        "assumed_saturation": ASSUMED_SATURATION,
        "capture_rate_chem_detail": chem_detail,
        "capture_rate_pbmc_detail": {
            a: {"q": d["q"], "q_ref": d["q_ref"], "reference": d["reference"],
                "ratios": d["ratios"]}
            for a, d in pbmc_detail.items()
        },
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    print(f"\nwriting {out_h5ad} ...")
    adata.write_h5ad(out_h5ad, compression="gzip", compression_opts=4)
    size = os.path.getsize(out_h5ad)
    print(f"  wrote {size/1e9:.2f} GB, peak_rss={peak_rss_gb():.1f}GB")

    write_readme(rep, chem_detail, pbmc_detail, sub_obs, shared, size,
                 args.variant, label_col, labels, out_h5ad)
    print(f"\ndone in {time.time()-t0:.0f}s")


def write_readme(rep, chem_detail, pbmc_detail, sub_obs, shared_donors, size_bytes,
                 variant, label_col, labels, out_h5ad):
    path = os.path.join(COEXPR_DIR, f"README_{variant}.md")
    ct_counts = sub_obs[label_col].astype(str).value_counts()
    chem_lines = "\n".join(
        f"| `{a}` | {d['q_chem']:.4f} | {d['lambda']:.3f} | {d['detected_fraction']:.4f} | **{d['q']:.5f}** |"
        for a, d in sorted(chem_detail.items()))
    pbmc_lines = "\n".join(
        f"| `{a}` | {d['reference']} | "
        + ", ".join(f"{g} {r:.2f}" for g, r in sorted(d["ratios"].items()))
        + f" | **{d['q']:.5f}** |"
        for a, d in sorted(pbmc_detail.items()))

    with open(path, "w") as fh:
        fh.write(f"""# Co-expression working set: Elmentaite2021 three-cell-type subset

Generated by `scripts/subset_coexpression_data.py`. **Preparation only** -- no
correlation has been computed. See
`/home/ubuntu/.claude/plans/clean-up-the-tissue-enchanted-axolotl.md` for the
plan this implements and `ANALYSIS_SUMMARY.md` for how we got here.

## What this is

`elmentaite2021_trio.h5ad` ({size_bytes/1e9:.2f} GB) -- {len(sub_obs):,} cells x
18,370 genes, carved from the Gut Cell Atlas "Extended+" release (1,596,200
cells), filtered to `study == "Elmentaite2021"` and three cell types.

| cell type | cells | role |
|---|---:|---|
| fibroblast | {ct_counts.get('fibroblast', 0):,} | canonical case: ANTXR1 + ANTXR2 + MRC2 all present |
| enterocyte | {ct_counts.get('enterocyte', 0):,} | the unexplained case: ANTXR2 present, ANTXR1/MRC2 absent |
| macrophage | {ct_counts.get('macrophage', 0):,} | myeloid calibration control |

**{len(shared_donors)} donors are shared across all three groups** (counted by
`donorID_unified`, not `donor_id` -- see "Donor identity" below) -- one lab,
one protocol, 100% `Non_pathological` samples. That is the reason for this
particular slice: it removes donor, cohort, and protocol confounds from the
three-way comparison simultaneously.

## Curation error corrected here

`gastrointestinal tract (lamina propria) macrophage` (CL:0000865) is **not a
macrophage** and is deliberately excluded. All 42,302 cells carrying that label
across the entire Gut Cell Atlas are authored as
`level_3_annot=Lamina_propria_fibroblast_ADAMDEC1` (Mesenchymal / Fibroblast),
and markers agree: PTPRC 0.4% positive, CD68 2.5%, C1QA 0.5%, COL1A1 **89%**
positive. The myeloid member of this trio is the plain `macrophage` label.

This error also affects earlier project outputs (the ECM-clearance heatmap and
the ANTXR1 ranking in `ANALYSIS_SUMMARY.md`, where it is described as a
macrophage) -- tracked as an open item there, not fixed by this script.

## Donor identity -- use `donorID_unified`, not `donor_id`

`donor_id` is unreliable in this atlas, in both directions:

- **collision**: `A34 (417C)` maps to two different people -- `D11` (55-74y,
  31,370 cells) and `D12` (18-34y, 5,814 cells). Different age brackets;
  unambiguously two donors sharing one id string.
- **aliases**: three people appear under several ids -- `D12` -> `A32 (411C)` /
  `A34 (417C)`, `D2` -> `T036` / `T036NEG` / `T036POS`, `D5` -> `T110NEG` /
  `T110POS` (the NEG/POS pairs are sorted fractions of one donor).

So `donor_id` merges two donors *and* splits three others: 41 id strings for 38
actual people in this study. This matters beyond bookkeeping -- memento's 2D
path bootstraps over donors, so keying on `donor_id` would pseudo-replicate.

**Use `obs['donor']`.** It is added by the pipeline and holds the resolved
identity (a copy of `donorID_unified`). `donor_id` and `donorID_unified` are
both passed through unmodified from the source, so nothing is rewritten -- but
`donor_id` is not safe to group on and `donor` exists so that the safe choice
is the one with the obvious name.

Note the Phase 1 parquet (`/data/ANTXR2/celltype_expression/`) resolves this
*differently*: there, `donor_id` itself holds the resolved identity, because it
is a pipeline output rather than a passthrough subset. Same column name, two
conventions across the two artifacts -- which is exactly why this file carries
an explicitly-named `donor` column.

(Separately: 4 donors elsewhere in the Gut Cell Atlas span multiple `sourceID`s
-- `390C`, `HT-228`, `HT-234`, `HT-236` -- but none are in Elmentaite2021, so
they do not affect this working set.)

## Contents

- `X` -- **raw integer counts**, copied from the source's `raw/X`. The source's
  `X` is log-normalized and is deliberately not carried: memento models the
  count-sampling process and requires raw counts.
- `var` -- unmodified from source (18,370 genes; index is Ensembl ID,
  `feature_name` carries the symbol -- note this differs from the skin atlas,
  where the index *is* the symbol). Genes are not subset: the co-expression
  scope is still open.
- `obs` -- all source columns for these cells, with unused categoricals dropped
  (they otherwise leak into memento's `create_groups` as phantom empty groups),
  plus the two capture-rate columns below.
- `obsm` -- `X_umap`, `X_scANVI` (row-sliced).
- **Excluded**: `obsp` (the 1.6M x 1.6M neighbour graph -- the structure that
  OOM-killed the Tabula Sapiens run), `uns` colour arrays (their lengths are
  tied to full-category counts and would be wrong after subsetting).

## Capture rate columns

Two independent estimates, kept separate rather than averaged. They rest on
different assumptions, so disagreement between them is informative.

### `capture_rate_chem` -- chemistry x sequencing saturation

| assay | q_chem | lambda | detected fraction | **q** |
|---|---:|---:|---:|---:|
{chem_lines}

`q_chem` is memento's own broad droplet rate (0.07), from their CELLxGene
pipeline (`publication/cellxgene/make_cube.py`) and their PBMC co-expression
analysis. Both chemistries get the same value because **no published
per-chemistry q exists** -- rather than invent one, the 5'-vs-3' difference is
measured empirically by the other column.

The saturation term is *not* a naive product. 10x "sequencing saturation" is a
read-duplicate fraction, not a molecule-loss fraction, so it is inverted through
the Poisson read->molecule relation (`saturation = 1 - (1-e^-L)/L`, then
`detected = 1-e^-L`). At high saturation this is nearly inert, which is the
honest result: the naive product would have understated q by ~15%.

**`ASSUMED_SATURATION = {ASSUMED_SATURATION}` is an assumption, not a
measurement.** The paper never states saturation, and
ENA has `read_count=0` for all 89 runs of the adult GEX accession (E-MTAB-9543),
so it cannot be recovered for these samples. Only E-MTAB-8901 (developing gut,
HiSeq 4000, ~44K reads/cell) supports the ~80-85% figure it is based on. Change
`ASSUMED_SATURATION` in `config.py` and re-run to test sensitivity.

### `capture_rate_pbmc` -- transferred from a 10x public PBMC reference

| assay | reference | UMI ratios by gate (ours/ref) | **q** |
|---|---|---|---:|
{pbmc_lines}

For a given cell type, detected UMIs scale with capture efficiency, so the ratio
of UMI depth between our data and a reference of assumed q transfers that q. We
use the median ratio across matched marker-gated immune populations (T, B, NK,
Mono/Mac) so one atypical population cannot dominate, with the gates applied
identically to both datasets.

**Main limitation**: this assumes a given immune cell type carries the same
absolute mRNA content in gut tissue as in blood. Tissue-resident lymphocytes are
plausibly more activated and RNA-richer than circulating ones, which would
inflate the apparent q. This is precisely why it is a second, independent column
and not a replacement for the chemistry-based one.

## `n_counts` is not the row sum -- don't treat them as interchangeable

`obs['n_counts']` is carried through from the source, but it runs slightly
*higher* than the row sums of the matrix shipped here -- median ~2 counts
(~0.04%), p99 ~1.7%, max ~14% on one cell. It is not a masking artifact
(`feature_is_filtered` is all-False); the release is titled "18485 genes" while
`var` holds 18,370, so `n_counts` evidently predates that gene subsetting.

The size is immaterial, but the distinction matters for anything depth-related:
every calculation in `subset_coexpression_data.py` that needs UMI depth (both
capture-rate columns, all the gate medians) uses **matrix row sums**, never
`n_counts`.

## Marker sanity check

These are what caught the curation error above, so they are recomputed on every
run rather than trusted from documentation.

```
{rep.to_string()}
```

## Not in scope

No correlation / 2D moments; no Kong2023 or disease contrast (every
Elmentaite2021 sample is `Non_pathological`, so its Crohn cells are
non-inflamed tissue from Crohn donors); no colonocyte (11 donors, 92% 5' v2 vs
enterocyte's 61%, so that contrast would be partly a chemistry contrast); no
FASTQ re-alignment; no change to the Phase 1 `DEFAULT_CAPTURE_RATE`.
""")
    print(f"  wrote {path}")


if __name__ == "__main__":
    main()
