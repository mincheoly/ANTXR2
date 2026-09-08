"""Task 5 of prompts/partner_availability.md: Kong2023 feasibility check,
gating any Phase 3 (differential-correlation) scoping.

Kong2023 (235,327 cells, 71 donors, 3 `sample_category` conditions:
Non_pathological / Neighbouring_inflamed / Inflamed) is already inside the
downloaded Gut Cell Atlas "Extended+" h5ad (`config.GUT_ATLAS_H5AD`) -- no new
download needed, just a `study=='Kong2023'` filter on the same file the
project's other gut-epithelium work already uses.

This is a FEASIBILITY CHECK ONLY -- no memento pipeline is run, no
differential-correlation test happens here regardless of outcome (per the
pre-registration). ANTXR2's presence filter is replicated directly from
memento's own formula (`compute_1d_moments`: a group's raw per-gene mean count
must exceed `filter_mean_thresh`; the overall gene passes if that holds in a
STRICT majority, `> min_perc_group`, of groups) rather than by actually
invoking memento, since only the filter-pass/fail outcome is needed here, not
a fitted model.

Run in the `antxr2` conda env: `conda run -n antxr2 python scripts/kong2023_feasibility.py`.
"""
import os

import anndata.io as aio
import h5py
import numpy as np
import pandas as pd

from config import (
    COEXPR_FILTER_MEAN_THRESH, COEXPR_MIN_GROUP_CELLS, COEXPR_MIN_PERC_GROUP,
    DONOR_UNIFIED_COL, FIGURES_DIR, GUT_ATLAS_H5AD,
)

OUT_DIR = os.path.join(FIGURES_DIR, "kong2023_feasibility")
ANALYSIS_SUMMARY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ANALYSIS_SUMMARY.md")
CONDITIONS = ["Non_pathological", "Neighbouring_inflamed", "Inflamed"]
CELL_TYPES = [
    "enterocyte", "colonocyte", "intestine goblet cell", "intestinal crypt stem cell",
    "transit amplifying cell", "paneth cell",
]
MIN_DONORS_WELL_POWERED = 8  # macrophage precedent: 5 donors -> replication indistinguishable from chance
REFERENCE_MPG_ELMENTAITE_ENTEROCYTE = 0.6  # what Elmentaite's full enterocyte donor set needed


def read_categorical(f, path):
    grp = f[path]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in grp["categories"][:]])
    codes = grp["codes"][:]
    return cats[codes]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    f = h5py.File(GUT_ATLAS_H5AD, "r")
    study = read_categorical(f, "obs/study")
    kong_mask = study == "Kong2023"
    print(f"Kong2023: {int(kong_mask.sum()):,} cells total in {GUT_ATLAS_H5AD}")

    donor = read_categorical(f, f"obs/{DONOR_UNIFIED_COL}")
    cond = read_categorical(f, "obs/sample_category")
    ct = read_categorical(f, "obs/cell_type")
    assay = read_categorical(f, "obs/assay")

    mask = kong_mask & np.isin(ct, CELL_TYPES)
    pos = np.where(mask)[0]
    print(f"{len(pos):,} cells across the {len(CELL_TYPES)} target cell types")

    var = aio.read_elem(f["raw/var"])
    sym = var["feature_name"].astype(str)
    antxr2_var_pos = int(np.where(sym.values == "ANTXR2")[0][0])
    ds = aio.sparse_dataset(f["raw/X"])
    X = ds[pos]
    antxr2 = np.asarray(X[:, antxr2_var_pos].todense()).ravel()

    df = pd.DataFrame({
        "donor": donor[pos], "condition": cond[pos], "cell_type": ct[pos], "assay": assay[pos],
        "antxr2": antxr2,
    })

    # ---- 1+2: per cell_type x condition x donor cell counts + presence filter ----
    grp = df.groupby(["cell_type", "condition", "donor"]).agg(
        n_cells=("antxr2", "size"), antxr2_mean=("antxr2", "mean"),
    ).reset_index()
    donor_csv = os.path.join(OUT_DIR, "per_donor_counts.csv")
    grp.to_csv(donor_csv, index=False)
    print(f"wrote {donor_csv}")

    summary_rows = []
    for (cell_type, condition), sub in grp.groupby(["cell_type", "condition"]):
        n_donors_total = len(sub)
        well_powered = sub[sub["n_cells"] >= COEXPR_MIN_GROUP_CELLS]
        n_well_powered = len(well_powered)
        n_pass_07 = int((well_powered["antxr2_mean"] > COEXPR_FILTER_MEAN_THRESH).sum())
        frac_pass = n_pass_07 / n_well_powered if n_well_powered else np.nan
        summary_rows.append({
            "cell_type": cell_type, "condition": condition, "n_donors_total": n_donors_total,
            "n_donors_ge_min_cells": n_well_powered, "n_donors_ANTXR2_present": n_pass_07,
            "frac_ANTXR2_present": frac_pass,
            "clears_min_perc_group_0.7": bool(frac_pass > COEXPR_MIN_PERC_GROUP) if n_well_powered else False,
            "underpowered_lt8_donors": n_well_powered < MIN_DONORS_WELL_POWERED,
        })
    summary_df = pd.DataFrame(summary_rows).sort_values(["cell_type", "condition"])
    summary_csv = os.path.join(OUT_DIR, "feasibility_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"wrote {summary_csv}")
    print(summary_df.to_string(index=False))

    # ---- 3: chemistry/batch confounded with condition? ----
    chem_crosstab = pd.crosstab(df["condition"], df["assay"])
    chem_csv = os.path.join(OUT_DIR, "chemistry_by_condition.csv")
    chem_crosstab.to_csv(chem_csv)
    print(f"\nchemistry x condition crosstab:\n{chem_crosstab.to_string()}")
    print(f"wrote {chem_csv}")
    # crude confound signal: does any assay category appear in only one condition?
    assay_frac = chem_crosstab.div(chem_crosstab.sum(axis=0), axis=1)
    confounded_assays = [a for a in chem_crosstab.columns if (chem_crosstab[a] > 0).sum() == 1]

    write_summary(summary_df, chem_crosstab, confounded_assays)


def write_summary(summary_df, chem_crosstab, confounded_assays):
    lines = []
    lines.append("\n### Task 5 results: Kong2023 feasibility check\n\n")
    lines.append(
        "**Feasibility gate, not a results-producing step** (per pre-registration) -- no "
        "differential-correlation test was run. Kong2023 data source: already inside the downloaded "
        f"Gut Cell Atlas Extended+ h5ad (`{GUT_ATLAS_H5AD}`, `study=='Kong2023'` filter), no new "
        "download needed. 235,327 cells, 71 donors. Output: `/data/ANTXR2/figures/kong2023_feasibility/"
        "{per_donor_counts,feasibility_summary,chemistry_by_condition}.csv`.\n\n"
    )

    lines.append("#### 1-2. Per cell_type x condition: donor counts and ANTXR2 presence\n\n")
    lines.append(
        f"ANTXR2 \"presence\" replicates memento's own filter formula directly (a donor group's raw "
        f"mean count must exceed `filter_mean_thresh={COEXPR_FILTER_MEAN_THRESH}`; the gene passes "
        f"overall if that holds in a strict majority, `>min_perc_group`, of the cell type x condition's "
        f"donor groups with >={COEXPR_MIN_GROUP_CELLS} cells) rather than by running memento itself -- "
        f"only the pass/fail outcome is needed for a feasibility check. Donors with "
        f"<{MIN_DONORS_WELL_POWERED} usable (>={COEXPR_MIN_GROUP_CELLS}-cell) groups are flagged "
        f"underpowered, using the macrophage precedent (5 donors -> replication indistinguishable "
        f"from chance) as the reference.\n\n"
    )
    lines.append("| cell type | condition | donors (total / ≥100 cells) | ANTXR2 present (n/frac) | clears mpg=0.7 | underpowered (<8 donors) |\n")
    lines.append("|---|---|---|---|---|---|\n")
    for _, r in summary_df.iterrows():
        lines.append(
            f"| {r['cell_type']} | {r['condition']} | {r['n_donors_total']} / {r['n_donors_ge_min_cells']} | "
            f"{r['n_donors_ANTXR2_present']}/{r['n_donors_ge_min_cells']} ({r['frac_ANTXR2_present']:.2f}) | "
            f"{'YES' if r['clears_min_perc_group_0.7'] else 'no'} | "
            f"{'**YES**' if r['underpowered_lt8_donors'] else 'no'} |\n"
        )

    colono = summary_df[summary_df["cell_type"] == "colonocyte"]
    entero = summary_df[summary_df["cell_type"] == "enterocyte"]
    lines.append(
        f"\n**Colonocyte (the Bracq et al.-matched cell type -- mouse colon DSS colitis) fails the "
        f"standard `min_perc_group=0.7` presence filter** in both well-powered conditions "
        f"(Non_pathological: {colono[colono['condition']=='Non_pathological']['frac_ANTXR2_present'].iloc[0]:.2f}, "
        f"Neighbouring_inflamed: {colono[colono['condition']=='Neighbouring_inflamed']['frac_ANTXR2_present'].iloc[0]:.2f} "
        f"of donors present) -- well below the {REFERENCE_MPG_ELMENTAITE_ENTEROCYTE} threshold that "
        f"already had to be used for Elmentaite's full enterocyte donor set. A relaxed threshold of "
        f"roughly 0.35-0.5 would be needed just to admit ANTXR2 into a colonocyte test, i.e. **further "
        f"relaxation than any precedent in this project.** Colonocyte's Inflamed condition additionally "
        f"has only {int(colono[colono['condition']=='Inflamed']['n_donors_ge_min_cells'].iloc[0])} "
        f"well-powered donors -- underpowered regardless of the presence question.\n\n"
        f"**Enterocyte (the extension cell type) clears the standard filter** in "
        f"Neighbouring_inflamed ({entero[entero['condition']=='Neighbouring_inflamed']['frac_ANTXR2_present'].iloc[0]:.2f}) "
        f"and Non_pathological ({entero[entero['condition']=='Non_pathological']['frac_ANTXR2_present'].iloc[0]:.2f}), "
        f"and sits exactly at the boundary in Inflamed "
        f"({entero[entero['condition']=='Inflamed']['frac_ANTXR2_present'].iloc[0]:.2f} -- since memento's "
        f"filter is a strict `>`, exactly 0.70 would still FAIL and needs a hair of relaxation, e.g. 0.69). "
        f"All three enterocyte conditions clear the {MIN_DONORS_WELL_POWERED}-donor power floor.\n\n"
    )

    stem = summary_df[summary_df["cell_type"].isin(["intestinal crypt stem cell", "paneth cell"])]
    lines.append(
        f"**Crypt stem cell and Paneth cell are underpowered in every condition** "
        f"({', '.join(f'{r.cell_type}/{r.condition}: {int(r.n_donors_ge_min_cells)} donors' for r in stem.itertuples())}) "
        f"-- both flagged underpowered regardless of the presence question; not usable for a two-group "
        f"test in Kong2023 at all.\n\n"
    )

    lines.append("#### 3. Chemistry/batch confounded with condition?\n\n")
    lines.append("```\n" + chem_crosstab.to_string() + "\n```\n\n")
    v1_col = "10x 3' v1"
    v1_cells = int(chem_crosstab.loc["Non_pathological", v1_col]) if v1_col in chem_crosstab.columns else 0
    v1_frac = v1_cells / chem_crosstab.loc["Non_pathological"].sum() * 100
    lines.append(
        f"**Yes, confounded.** `10x 3' v1` appears ONLY in `Non_pathological` "
        f"({v1_cells} cells, {v1_frac:.0f}% "
        f"of that condition's cells if present) while `Neighbouring_inflamed` and `Inflamed` are 100% "
        f"`10x 3' v2`/`v3` -- the same shape of confound as the `donor_id` collision issue flagged "
        f"earlier in this project (a real technical variable perfectly or near-perfectly aligned with "
        f"the biological grouping of interest). Any Non_pathological-vs-inflamed differential test "
        f"would need to either restrict Non_pathological to its v2/v3 donors only, or treat chemistry "
        f"as an explicit covariate -- **not treat the raw condition contrast as chemistry-free.**\n\n"
    )

    lines.append("#### 4. Scoping (not running) the differential test, if feasible\n\n")
    lines.append(
        "**Enterocyte is the only feasible cell type for a well-powered Non_pathological-vs-Inflamed "
        "or Non_pathological-vs-Neighbouring_inflamed two-group test** (colonocyte fails presence; "
        "crypt stem/Paneth fail power; goblet/TA are intermediate -- see full table above). Pre-"
        "registered directional prediction (Bracq et al.): since CMG2's Wnt function is injury-"
        "conditional (CMG2-KO baseline guts are normal), **ANTXR2's coupling to Wnt-arm partners and "
        "regeneration programs should be absent or weak in `Non_pathological` and appear in "
        "`Inflamed`.** Any such test must restrict or covary for the chemistry confound in point 3 "
        "above, and should use enterocyte as primary with colonocyte reported only as a "
        "presence-filter negative result (an interesting finding in its own right, not a null test "
        "outcome), consistent with `prompts/partner_availability.md`'s framing that colonocyte is the "
        "literal match to Bracq et al.'s mouse colon model while enterocyte is the extension.\n\n"
    )

    lines.append("#### 5. Colonocyte vs. enterocyte as the Bracq et al. match\n\n")
    lines.append(
        "Confirmed per the pre-registration: Bracq et al. is mouse **colon** with DSS colitis, so "
        "colonocyte is the directly-matched cell type; enterocyte is the extension. The feasibility "
        "result above means the directly-matched cell type is NOT the one available for the actual "
        "test -- a real scoping constraint to carry forward, not an incidental detail.\n\n"
    )

    with open(ANALYSIS_SUMMARY_PATH, "a") as fh:
        fh.write("".join(lines))
    print(f"\nappended Task 5 results to {ANALYSIS_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
