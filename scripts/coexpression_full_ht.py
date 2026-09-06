"""Full-dataset one-sample significance test for ANTXR2-vs-gene correlations,
per user request this session: the discovery/replication split in
coexpression_discovery_replication.py deliberately used only HALF the donors
per cell type (to keep gene selection non-circular); this script runs the same
memento one-sample bootstrap test (`ht_2d_moments` with a constant treatment
column -- see that script's docstring for why this is a real one-sample test,
not a hack) on ALL usable donors per cell type, for maximum power. This is the
authoritative, best-powered significance test in this project -- downstream
work (pairwise heatmaps, GSEA z-score ranking, cross-cell-type gene overlap)
should build on ITS output, not the discovery-half-only results, which exist
only to validate the method itself.

Output per cell type: a CSV of EVERY gene with a valid one-sample test (coef,
se, z=coef/se, pval, fdr) -- not just FDR-significant ones, since GSEA and any
future re-analysis need the full ranked universe, not a pre-filtered subset
(this project's established GSEA convention, see coexpression_gsea.py).

Run in the `antxr2` conda env: `conda run -n antxr2 python
scripts/coexpression_full_ht.py --cell-types fibroblast` (one at a time --
each is a multi-minute job; see coexpression_discovery_replication.py's
docstring for why this project's tooling runs these one cell type per
invocation rather than as one long-lived job).
"""
import argparse
import os

import anndata
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from config import COEXPR_MIN_GROUP_CELLS, COEXPR_FIGURES_DIR
from coexpression_discovery_replication import (
    INPUT_H5AD, GENE_NAME_COL, DONOR_COL, LEVEL3_LABEL, run_discovery_ht,
)

OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "full_dataset_ht")

CSV_HEADER = """\
# ANTXR2-vs-gene one-sample bootstrap significance test (memento ht_2d_moments,
# constant-treatment trick -- see coexpression_discovery_replication.py's docstring),
# run on ALL usable donors for this cell type (not a discovery/replication half --
# that split exists only in coexpression_discovery_replication.py, to validate this
# method; this file is the best-powered, authoritative result for downstream use).
# Every gene with a valid test is included, not just FDR-significant ones -- GSEA
# and any other full-ranking analysis should use this file directly.
# z = coef / se (a signed significance-weighted score -- large |coef| with large se,
# i.e. a noisy/imprecise estimate, gets pulled toward 0 relative to a gene with the
# same coef but tighter se).
"""


def all_usable_donors(sub_full, donor_col):
    n_cells_by_donor = sub_full.obs[donor_col].value_counts()
    usable = n_cells_by_donor[n_cells_by_donor >= COEXPR_MIN_GROUP_CELLS].index.tolist()
    return n_cells_by_donor.reindex(usable).sort_values(ascending=False).index.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-types", default="fibroblast,enterocyte,macrophage")
    ap.add_argument("--num-boot", type=int, default=5000)
    ap.add_argument("--num-cpus", type=int, default=12)
    ap.add_argument("--random-state", type=int, default=42)
    args = ap.parse_args()
    cell_types = args.cell_types.split(",")

    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"loading {INPUT_H5AD} ...")
    full = anndata.read_h5ad(INPUT_H5AD)

    for ct in cell_types:
        print(f"\n{'='*70}\n{ct} (full dataset, all usable donors)\n{'='*70}")
        label = LEVEL3_LABEL[ct]
        sub_full = full[full.obs["level_3_annot"].astype(str) == label]
        donors = all_usable_donors(sub_full, DONOR_COL)
        print(f"{len(donors)} usable donors: {donors}")

        ht_df, mpg_used = run_discovery_ht(full, ct, donors, args.num_boot, args.num_cpus, args.random_state)
        ht_df["z"] = ht_df["coef"] / ht_df["se"]
        sym_map = dict(zip(full.var.index, full.var[GENE_NAME_COL].astype(str)))
        ht_df.insert(1, "gene_symbol", ht_df["gene_id"].map(sym_map))

        n_sig_01 = int((ht_df["fdr"] < 0.1).sum())
        n_sig_05 = int((ht_df["fdr"] < 0.05).sum())
        print(f"{len(ht_df)} genes with a valid one-sample test (min_perc_group={mpg_used:.2f}, "
              f"{len(donors)} donors); {n_sig_01} FDR<0.1, {n_sig_05} FDR<0.05")

        out_path = os.path.join(OUT_DIR, f"{ct}_full_dataset_ht.csv")
        with open(out_path, "w") as f:
            f.write(CSV_HEADER)
            ht_df.sort_values("pval").to_csv(f, index=False)
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
