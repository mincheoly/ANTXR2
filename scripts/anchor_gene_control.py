"""Task 2 of prompts/partner_availability.md: anchor-gene specificity control.

Establishes the nonspecific-covariation floor that gates every positive
correlation claim in the project (the fibroblast ECM result and the
enterocyte digestion/absorption result specifically). For each cell type,
~20 anchor genes are matched to ANTXR2 on both mean expression and detection
rate, excluding ECM/ribosomal/mitochondrial genes (so the null isn't
contaminated by the modules under test), then run through the IDENTICAL
pipeline ANTXR2 went through:
  - `coexpression_full_ht.py`'s one-sample bootstrap (reusing
    `coexpression_discovery_replication.run_discovery_ht`, generalized here
    to accept an arbitrary `target_symbol`), same num_boot per cell type as
    ANTXR2's own run (fibroblast 5000, enterocyte 3000 -- see
    ANALYSIS_SUMMARY.md), same adaptive min_perc_group retry.
  - `coexpression_gsea.py --metric zscore`'s prerank GSEA, same gene-set
    libraries (GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022),
    same FDR_THRESH=0.25.

Macrophage is excluded (standing correction: exploratory-only everywhere).

Per-cell-type invocation, per-anchor checkpointing to disk (each anchor's HT
CSV and GSEA CSV are written as soon as computed, and skipped on a re-run if
already present) -- this project's own scripts note that a multi-minute-per-
gene job run inside one long-lived process is less reliable in this tooling
than one invocation per cell type; the per-anchor checkpoint additionally
means an interrupted run loses at most the one anchor in flight, not the
whole batch.

Run in the `antxr2` conda env:
`conda run -n antxr2 python scripts/anchor_gene_control.py --cell-types fibroblast`
(one cell type at a time; run again with `--cell-types enterocyte`, then with
`--assemble-only` across both to build the aggregate report once every anchor
for both cell types has a cached HT+GSEA result).
"""
import argparse
import os
import pickle
import re

import anndata
import gseapy as gp
import numpy as np
import pandas as pd

from config import COEXPR_FIGURES_DIR, COEXPR_VARIANTS, COEXPR_TARGET_GENE, GENE_PANEL_ARMS
from coexpression_discovery_replication import GENE_NAME_COL, LEVEL3_LABEL, run_discovery_ht
from coexpression_full_ht import all_usable_donors
from coexpression_gsea import GENE_SETS, FDR_THRESH

OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "anchor_gene_control")
HT_DIR = os.path.join(OUT_DIR, "full_ht")
GSEA_DIR = os.path.join(OUT_DIR, "gsea_zscore")
ANALYSIS_SUMMARY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ANALYSIS_SUMMARY.md")

INPUT_H5AD = COEXPR_VARIANTS["level3"]["h5ad"]
DONOR_COL = "donor"
CELL_TYPES = ["fibroblast", "enterocyte"]  # macrophage excluded -- standing correction, exploratory-only
N_ANCHORS = 20
# Matches ANTXR2's OWN full_dataset_ht run per cell type (ANALYSIS_SUMMARY.md,
# "Full-dataset one-sample test" section) -- "same num_boot" per the prompt.
NUM_BOOT_BY_CT = {"fibroblast": 5000, "enterocyte": 3000}
NUM_CPUS = 12
RANDOM_STATE_BASE = 1000  # offset from ANTXR2's own random_state=42 -- distinct bootstrap draws per anchor
MEAN_TOL_LOG10_INIT = 0.25   # ~1.8-fold band, initial
DETECT_TOL_INIT = 0.05       # +/-5 percentage points, initial
BAND_WIDEN_FACTOR = 1.5
MAX_WIDEN_STEPS = 6
ANCHOR_SEED = 20260906

# Exclusion regexes -- ECM/matrisome-core, cytoplasmic ribosomal, mitochondrial
# (mtDNA-encoded). Matches the two blocks the project's own fibroblast pairwise
# heatmap visually identified (ECM/collagen module; ribosomal-protein block;
# mitochondrial-gene block) plus a standard matrisome-core keyword list, so the
# anchor null isn't drawn from the modules under test.
ECM_REGEX = re.compile(
    r"^(COL\d|LAMA\d|LAMB\d|LAMC\d|FN1$|ELN$|DCN$|BGN$|LUM$|VCAN$|SPARC$|POSTN$|"
    r"FBN\d|TNC$|THBS\d|MMP\d|TIMP\d|ADAMTS|ADAM\d|ITGA\d|ITGB\d|FAP$|PLOD\d|LOX$|LOXL\d)"
)
RIBO_REGEX = re.compile(r"^(RPS|RPL)\d")
MITO_REGEX = re.compile(r"^MT-")
# Also exclude genes that are themselves part of the project's own hypothesis
# panels -- an anchor should be an arbitrary expression-matched gene, not one
# of the partners/paralogs already under test.
EXTRA_EXCLUDE = set(g for arm in GENE_PANEL_ARMS.values() for g in arm) | {COEXPR_TARGET_GENE}

# Key terms checked for recovery, taken verbatim from ANTXR2's own significant
# z-score GSEA output (gsea_zscore/{ct}_gsea_significant.csv) -- these are the
# terms the pre-registration and the project's log specifically call out.
KEY_TERMS = {
    "fibroblast": [
        ("GO_Biological_Process_2023", "Extracellular Matrix Organization (GO:0030198)"),
        ("KEGG_2021_Human", "ECM-receptor interaction"),
        ("Reactome_2022", "Collagen Formation R-HSA-1474290"),
        ("KEGG_2021_Human", "Focal adhesion"),
    ],
    "enterocyte": [
        ("KEGG_2021_Human", "Protein digestion and absorption"),
        ("KEGG_2021_Human", "Fat digestion and absorption"),
    ],
}


def is_excluded(symbol):
    return bool(ECM_REGEX.match(symbol) or RIBO_REGEX.match(symbol) or MITO_REGEX.match(symbol)) \
        or symbol in EXTRA_EXCLUDE


# ---------------------------------------------------------------------------
# Step 1: anchor selection
# ---------------------------------------------------------------------------

def select_anchors(ct, n_anchors=N_ANCHORS, seed=ANCHOR_SEED):
    adata = anndata.read_h5ad(INPUT_H5AD)
    label = LEVEL3_LABEL[ct]
    sub = adata[adata.obs["level_3_annot"].astype(str) == label]
    sym = adata.var[GENE_NAME_COL].astype(str).values
    X = sub.X
    mean_expr = np.asarray(X.mean(axis=0)).ravel()
    detect_rate = np.asarray((X > 0).mean(axis=0)).ravel()
    df_all = pd.DataFrame({
        "gene_id": adata.var.index, "gene_symbol": sym, "mean_expr": mean_expr, "detect_rate": detect_rate,
    })

    # ANTXR2's own stats must come from the FULL per-gene table, not the
    # tested-universe-filtered one below -- ANTXR2 is excluded from its own
    # "tested against every other gene" CSV (coexpression_full_ht.py removes
    # the target from test_genes), so intersecting first would silently drop it.
    target_row = df_all[df_all["gene_symbol"] == COEXPR_TARGET_GENE]
    if target_row.empty:
        raise SystemExit(f"{COEXPR_TARGET_GENE} not found in {ct}'s var at all")
    target_log_mean = np.log10(target_row["mean_expr"].iloc[0] + 1e-9)
    target_detect = target_row["detect_rate"].iloc[0]
    print(f"  {ct}: {COEXPR_TARGET_GENE} mean_expr={target_row['mean_expr'].iloc[0]:.5f} "
          f"(log10={target_log_mean:.3f}), detect_rate={target_detect:.4f}")

    tested = pd.read_csv(os.path.join(COEXPR_FIGURES_DIR, "full_dataset_ht", f"{ct}_full_dataset_ht.csv"),
                          comment="#")
    tested_ids = set(tested["gene_id"])
    df = df_all[df_all["gene_id"].isin(tested_ids)].reset_index(drop=True)
    print(f"  {ct}: {len(df)} genes in the already-tested universe available as anchor candidates")

    df["log_mean"] = np.log10(df["mean_expr"] + 1e-9)
    df["mean_diff"] = (df["log_mean"] - target_log_mean).abs()
    df["detect_diff"] = (df["detect_rate"] - target_detect).abs()

    excluded_mask = df["gene_symbol"].apply(is_excluded)
    n_excluded = int(excluded_mask.sum())
    candidates = df[~excluded_mask].reset_index(drop=True)
    print(f"  {ct}: {n_excluded} candidate genes excluded (ECM/ribosomal/mitochondrial/panel-overlap), "
          f"{len(candidates)} remain as eligible anchor candidates")

    mean_tol, detect_tol = MEAN_TOL_LOG10_INIT, DETECT_TOL_INIT
    for step in range(MAX_WIDEN_STEPS + 1):
        band = candidates[(candidates["mean_diff"] <= mean_tol) & (candidates["detect_diff"] <= detect_tol)]
        print(f"  {ct}: band attempt {step} (mean_tol=log10±{mean_tol:.3f}, detect_tol=±{detect_tol:.3f}): "
              f"{len(band)} candidates in band")
        if len(band) >= n_anchors or step == MAX_WIDEN_STEPS:
            break
        mean_tol *= BAND_WIDEN_FACTOR
        detect_tol *= BAND_WIDEN_FACTOR

    if len(band) < n_anchors:
        print(f"  {ct}: WARNING -- only {len(band)} candidates available even after widening to "
              f"mean_tol=log10±{mean_tol:.3f}, detect_tol=±{detect_tol:.3f}; using all of them")
        selected = band.copy()
    else:
        rng = np.random.default_rng(seed)
        selected = band.sample(n=n_anchors, random_state=rng.integers(1e9)).reset_index(drop=True)

    selected = selected.sort_values("gene_symbol").reset_index(drop=True)
    band_record = {"mean_tol_log10": mean_tol, "detect_tol": detect_tol, "n_candidates_in_band": len(band),
                    "n_selected": len(selected)}
    return selected, band_record


# ---------------------------------------------------------------------------
# Step 2: run each anchor through the identical full-HT pipeline
# ---------------------------------------------------------------------------

def run_anchor_ht(full, ct, anchor_symbol, anchor_idx, num_boot, num_cpus):
    out_path = os.path.join(HT_DIR, f"{ct}_{anchor_symbol}_full_dataset_ht.csv")
    if os.path.exists(out_path):
        print(f"    [{ct}/{anchor_symbol}] cached HT result found, skipping computation")
        return pd.read_csv(out_path, comment="#")

    label = LEVEL3_LABEL[ct]
    sub_full = full[full.obs["level_3_annot"].astype(str) == label]
    donors = all_usable_donors(sub_full, DONOR_COL)
    random_state = RANDOM_STATE_BASE + anchor_idx

    ht_df, mpg_used = run_discovery_ht(full, ct, donors, num_boot, num_cpus, random_state,
                                        target_symbol=anchor_symbol)
    ht_df["z"] = ht_df["coef"] / ht_df["se"]
    sym_map = dict(zip(full.var.index, full.var[GENE_NAME_COL].astype(str)))
    ht_df.insert(1, "gene_symbol", ht_df["gene_id"].map(sym_map))

    header = (f"# anchor-gene full-dataset one-sample bootstrap test, gene={anchor_symbol}, "
              f"cell_type={ct}, num_boot={num_boot}, min_perc_group_final={mpg_used:.2f}, "
              f"random_state={random_state}\n")
    with open(out_path, "w") as f:
        f.write(header)
        ht_df.sort_values("pval").to_csv(f, index=False)
    print(f"    [{ct}/{anchor_symbol}] {len(donors)} donors, min_perc_group={mpg_used:.2f}, "
          f"{len(ht_df)} genes tested -- wrote {out_path}")
    return ht_df


# ---------------------------------------------------------------------------
# Step 3: GSEA (zscore metric) on each anchor's ranking
# ---------------------------------------------------------------------------

def run_anchor_gsea(ct, anchor_symbol):
    out_path = os.path.join(GSEA_DIR, f"{ct}_{anchor_symbol}_gsea_all.csv")
    if os.path.exists(out_path):
        print(f"    [{ct}/{anchor_symbol}] cached GSEA result found, skipping computation")
        return pd.read_csv(out_path)

    ht_path = os.path.join(HT_DIR, f"{ct}_{anchor_symbol}_full_dataset_ht.csv")
    df = pd.read_csv(ht_path, comment="#")
    df = df.dropna(subset=["z"]).rename(columns={"z": "score"})
    df["abs_score"] = df["score"].abs()
    df = df.sort_values("abs_score", ascending=False).drop_duplicates("gene_symbol", keep="first")
    rnk = df.set_index("gene_symbol")["score"].sort_values(ascending=False)

    pre = gp.prerank(rnk=rnk, gene_sets=GENE_SETS, organism="human", min_size=5, max_size=1000,
                      permutation_num=1000, outdir=None, seed=0, threads=4, no_plot=True)
    res = pre.res2d.copy()
    split = res["Term"].str.split("__", n=1, expand=True)
    res["Gene_set"] = split[0]
    res["Term"] = split[1]
    res["FDR q-val"] = pd.to_numeric(res["FDR q-val"], errors="coerce")
    res["NES"] = pd.to_numeric(res["NES"], errors="coerce")
    res.to_csv(out_path, index=False)
    print(f"    [{ct}/{anchor_symbol}] GSEA done, {len(res)} terms tested -- wrote {out_path}")
    return res


# ---------------------------------------------------------------------------
# main per-cell-type run
# ---------------------------------------------------------------------------

def run_cell_type(ct):
    os.makedirs(HT_DIR, exist_ok=True)
    os.makedirs(GSEA_DIR, exist_ok=True)

    anchors, band_record = select_anchors(ct)
    anchor_csv = os.path.join(OUT_DIR, f"{ct}_selected_anchors.csv")
    anchors.to_csv(anchor_csv, index=False)
    band_json = os.path.join(OUT_DIR, f"{ct}_anchor_band.pkl")
    with open(band_json, "wb") as f:
        pickle.dump(band_record, f)
    print(f"  wrote {anchor_csv} ({band_record})")

    print(f"\nloading {INPUT_H5AD} for {ct} HT runs...")
    full = anndata.read_h5ad(INPUT_H5AD)
    num_boot = NUM_BOOT_BY_CT[ct]

    for i, row in anchors.iterrows():
        sym = row["gene_symbol"]
        print(f"\n  -- anchor {i+1}/{len(anchors)}: {sym} --")
        run_anchor_ht(full, ct, sym, i, num_boot, NUM_CPUS)
        run_anchor_gsea(ct, sym)
    del full
    print(f"\n{ct}: all {len(anchors)} anchors done (HT + GSEA).")


# ---------------------------------------------------------------------------
# assemble: aggregate report across cell types
# ---------------------------------------------------------------------------

def antxr2_gsea_lookup(ct, gene_set, term):
    """ANTXR2's own {ct}_gsea_all.csv (coexpression_gsea.py) never got the
    Gene_set/Term split applied to the _all.csv (only _significant.csv gets
    it, inside summarize()) -- Term there is still "<library>__<term>"."""
    path = os.path.join(COEXPR_FIGURES_DIR, "gsea_zscore", f"{ct}_gsea_all.csv")
    df = pd.read_csv(path)
    if "Gene_set" not in df.columns:
        split = df["Term"].str.split("__", n=1, expand=True)
        df["Gene_set"], df["Term"] = split[0], split[1]
    row = df[(df["Gene_set"] == gene_set) & (df["Term"] == term)]
    if row.empty:
        return None, None
    return float(row["NES"].iloc[0]), float(row["FDR q-val"].iloc[0])


def assemble():
    lines = []
    lines.append("\n### Task 2 results: anchor-gene specificity control\n\n")
    lines.append(
        "**Confirmatory** (pre-registered interpretation rule, see pre-registration above) -- this "
        "gates the fibroblast ECM-correlation and enterocyte digestion/absorption-correlation "
        "exploratory GSEA findings reported earlier in this log; those remain exploratory GSEA "
        "characterizations, and this section is what determines whether they carry any ANTXR2-"
        "specific evidentiary weight.\n\n"
    )

    for ct in CELL_TYPES:
        anchor_csv = os.path.join(OUT_DIR, f"{ct}_selected_anchors.csv")
        band_pkl = os.path.join(OUT_DIR, f"{ct}_anchor_band.pkl")
        if not (os.path.exists(anchor_csv) and os.path.exists(band_pkl)):
            lines.append(f"#### {ct}: NOT YET RUN\n\n")
            continue
        anchors = pd.read_csv(anchor_csv)
        with open(band_pkl, "rb") as f:
            band = pickle.load(f)

        lines.append(f"#### {ct}\n\n")
        lines.append(
            f"**Anchor selection**: {band['n_selected']} genes randomly drawn (seed={ANCHOR_SEED}) from "
            f"{band['n_candidates_in_band']} candidates matched to ANTXR2 within log10±"
            f"{band['mean_tol_log10']:.3f} mean expression and ±{band['detect_tol']:.3f} detection rate "
            f"(ECM/ribosomal/mitochondrial genes and this project's own panel genes excluded from the "
            f"candidate pool). Selected: {', '.join(sorted(anchors['gene_symbol']))}. Full candidate "
            f"table: `{ct}_selected_anchors.csv`.\n\n"
        )

        n_anchors_run = sum(
            os.path.exists(os.path.join(GSEA_DIR, f"{ct}_{sym}_gsea_all.csv")) for sym in anchors["gene_symbol"]
        )
        lines.append(f"{n_anchors_run}/{len(anchors)} anchors have a completed HT+GSEA run.\n\n")
        if n_anchors_run == 0:
            continue

        for gene_set, term in KEY_TERMS[ct]:
            antxr2_nes, antxr2_fdr = antxr2_gsea_lookup(ct, gene_set, term)
            anchor_nes = []
            n_recovered = 0
            for sym in anchors["gene_symbol"]:
                gsea_path = os.path.join(GSEA_DIR, f"{ct}_{sym}_gsea_all.csv")
                if not os.path.exists(gsea_path):
                    continue
                df = pd.read_csv(gsea_path)
                row = df[(df["Gene_set"] == gene_set) & (df["Term"] == term)]
                if row.empty:
                    continue
                nes, fdr = float(row["NES"].iloc[0]), float(row["FDR q-val"].iloc[0])
                anchor_nes.append(nes)
                if fdr < FDR_THRESH and nes > 0:
                    n_recovered += 1

            if not anchor_nes:
                lines.append(f"- **{term}** ({gene_set}): no anchor GSEA results available yet.\n")
                continue
            anchor_nes_arr = np.array(anchor_nes)
            if antxr2_nes is not None:
                percentile = float((anchor_nes_arr < antxr2_nes).mean() * 100)
            else:
                percentile = float("nan")
            lines.append(
                f"- **{term}** ({gene_set}): ANTXR2 NES={antxr2_nes:.2f} (FDR={antxr2_fdr:.3g}) vs. "
                f"**{n_recovered}/{len(anchor_nes)} anchors also recover this term** at the same FDR<"
                f"{FDR_THRESH} threshold (anchor NES range [{anchor_nes_arr.min():.2f}, "
                f"{anchor_nes_arr.max():.2f}], mean {anchor_nes_arr.mean():.2f}). ANTXR2's NES sits at "
                f"the **{percentile:.0f}th percentile** of the anchor distribution.\n"
            )

        # interpretation call, applying the pre-registered rule mechanically
        majority_recover = []
        for gene_set, term in KEY_TERMS[ct]:
            anchor_nes = []
            n_recovered = 0
            for sym in anchors["gene_symbol"]:
                gsea_path = os.path.join(GSEA_DIR, f"{ct}_{sym}_gsea_all.csv")
                if not os.path.exists(gsea_path):
                    continue
                df = pd.read_csv(gsea_path)
                row = df[(df["Gene_set"] == gene_set) & (df["Term"] == term)]
                if row.empty:
                    continue
                nes, fdr = float(row["NES"].iloc[0]), float(row["FDR q-val"].iloc[0])
                anchor_nes.append(nes)
                if fdr < FDR_THRESH and nes > 0:
                    n_recovered += 1
            if anchor_nes:
                majority_recover.append(n_recovered > len(anchor_nes) / 2)
        if majority_recover:
            if any(majority_recover):
                lines.append(
                    f"\n**Interpretation (pre-registered rule applied): at least one key term is "
                    f"recovered by a MAJORITY of anchors in {ct}** -- for that term/those terms, the "
                    f"positive ANTXR2 correlation is a property of {ct}'s transcriptional program, not "
                    f"specific evidence about ANTXR2. This is a real, useful negative for those terms, "
                    f"reported as such rather than suppressed.\n\n"
                )
            else:
                lines.append(
                    f"\n**Interpretation (pre-registered rule applied): no key term is recovered by a "
                    f"majority of anchors in {ct}** -- the ANTXR2 correlation with this program is not "
                    f"simply a generic property of {ct}'s cells at ANTXR2's expression level, supporting "
                    f"(but not proving causally) a specific relationship.\n\n"
                )

    lines.append(check_z_vs_mean())

    with open(ANALYSIS_SUMMARY_PATH, "a") as f:
        f.write("".join(lines))
    print(f"\nappended Task 2 results to {ANALYSIS_SUMMARY_PATH}")


def check_z_vs_mean():
    """Task 2 step 6: does |z| correlate with per-gene mean expression across
    the tested universe? Uses ANTXR2's own already-computed full_dataset_ht.csv
    (no new computation needed) plus mean expression from the level3 h5ad."""
    lines = ["#### |z| vs. mean-expression check (Task 2, step 6)\n\n"]
    adata = anndata.read_h5ad(INPUT_H5AD)
    sym = adata.var[GENE_NAME_COL].astype(str).values
    for ct in CELL_TYPES:
        label = LEVEL3_LABEL[ct]
        sub = adata[adata.obs["level_3_annot"].astype(str) == label]
        mean_expr = np.asarray(sub.X.mean(axis=0)).ravel()
        mean_df = pd.DataFrame({"gene_id": adata.var.index, "mean_expr": mean_expr})

        ht = pd.read_csv(os.path.join(COEXPR_FIGURES_DIR, "full_dataset_ht", f"{ct}_full_dataset_ht.csv"),
                          comment="#")
        merged = ht.merge(mean_df, on="gene_id", how="left").dropna(subset=["z", "mean_expr"])
        r = float(np.corrcoef(np.log10(merged["mean_expr"] + 1e-9), merged["z"].abs())[0, 1])
        r_signed = float(np.corrcoef(np.log10(merged["mean_expr"] + 1e-9), merged["z"])[0, 1])
        lines.append(
            f"- **{ct}**: corr(log10 mean expression, |z|) = {r:.3f} across {len(merged)} tested genes "
            f"(corr with signed z = {r_signed:.3f}). "
        )
        if abs(r) > 0.2:
            lines.append(
                "Non-trivial positive relationship -- consistent with memento's bootstrap se shrinking "
                "with cell count/detection rate, so highly-expressed genes can get inflated |z| in both "
                "tails. Treat this as a caveat on z-ranked results generally (this project's GSEA "
                "z-score ranking included); the anchor-matching above at least holds mean expression "
                "roughly fixed between ANTXR2 and its anchors, which limits (but does not eliminate) "
                "this confound's effect on the recovery-count comparison specifically.\n"
            )
        else:
            lines.append("No strong relationship -- the |z| vs. expression-level confound does not "
                         "appear to be a major factor for this cell type's ranking.\n")
    return "".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-types", default=None, help="comma-separated (default: fibroblast,enterocyte)")
    ap.add_argument("--assemble-only", action="store_true")
    args = ap.parse_args()

    if args.assemble_only:
        assemble()
        return

    cell_types = args.cell_types.split(",") if args.cell_types else CELL_TYPES
    for ct in cell_types:
        print(f"\n{'='*70}\n{ct}\n{'='*70}")
        run_cell_type(ct)


if __name__ == "__main__":
    main()
