"""Downstream products of the full-dataset one-sample significance test
(coexpression_full_ht.py): cross-cell-type gene overlap and a pairwise
correlation heatmap over the union of significant genes -- replacing the
original pipeline's top-50-by-|mean_corr| gene panel (coexpression_pipeline.py
select_top_genes) with a significance-based panel, per user request.

Reuses coexpression_pipeline.py's compute_pairwise_correlations /
average_pairwise_across_donors / plot_heatmaps unchanged -- only the gene
panel construction changes (FDR<thresh from the full-dataset one-sample test,
not |mean_corr| magnitude ranking).

Run in the `antxr2` conda env: `conda run -n antxr2 python
scripts/coexpression_full_ht_downstream.py`.
"""
import os

import anndata

from config import COEXPR_CELL_TYPES, COEXPR_FIGURES_DIR, COEXPR_TARGET_GENE, COEXPR_VARIANTS
from coexpression_discovery_replication import GENE_NAME_COL
from coexpression_pipeline import (
    average_pairwise_across_donors, compute_pairwise_correlations, plot_heatmaps,
)
from coexpression_full_ht import OUT_DIR as HT_DIR

INPUT_H5AD = COEXPR_VARIANTS["level3"]["h5ad"]
FDR_THRESH = 0.1
MIN_PERC_GROUP, SHRINKAGE, TRIM_PERCENT = 0.7, 0.0, 0.5
# The union of ALL FDR<0.1 genes (542 across the 3 cell types) is the right universe for
# overlap counts and GSEA (no size constraint there), but far too large to RENDER as a
# legible heatmap (542x542 blew past image size limits and would be illegible even if it
# hadn't). Cap the heatmap panel to each cell type's top-N by p-value (not magnitude) among
# its own significant genes -- same "top-N-per-cell-type, then union" structure as the
# original pipeline's select_top_genes, just with significance instead of |mean_corr| as
# the ranking criterion. Purely a rendering concession; overlap/GSEA are unaffected.
HEATMAP_TOP_N_PER_CT = 50


def load_sig_sets():
    import pandas as pd
    sig = {}
    tables = {}
    for ct in COEXPR_CELL_TYPES:
        df = pd.read_csv(os.path.join(HT_DIR, f"{ct}_full_dataset_ht.csv"), comment="#")
        tables[ct] = df
        sig[ct] = set(df.loc[df["fdr"] < FDR_THRESH, "gene_id"])
        print(f"  {ct}: {len(df)} genes tested, {len(sig[ct])} FDR<{FDR_THRESH}")
    return sig, tables


def report_overlap(sig, sym_map):
    import itertools
    import pandas as pd
    print(f"\nPairwise overlap between cell types' FDR<{FDR_THRESH} significant gene sets:")
    rows = []
    for a, b in itertools.combinations(COEXPR_CELL_TYPES, 2):
        inter = sig[a] & sig[b]
        frac = len(inter) / max(1, min(len(sig[a]), len(sig[b])))
        print(f"  {a} ({len(sig[a])}) vs {b} ({len(sig[b])}): {len(inter)} shared ({frac*100:.0f}% of the smaller set)")
        rows.append({"cell_type_a": a, "cell_type_b": b, "n_a": len(sig[a]), "n_b": len(sig[b]),
                      "n_shared": len(inter), "frac_of_smaller": frac,
                      "shared_gene_symbols": ";".join(sorted(sym_map.get(g, g) for g in inter))})
    all3 = sig["fibroblast"] & sig["enterocyte"] & sig["macrophage"]
    print(f"  significant in all three cell types: {len(all3)}")
    out_path = os.path.join(HT_DIR, "significant_gene_overlap.csv")
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"  wrote {out_path}")


def main():
    print("Loading full-dataset one-sample test results...")
    sig, tables = load_sig_sets()
    sym_map_from_tables = {}
    for df in tables.values():
        sym_map_from_tables.update(dict(zip(df["gene_id"], df["gene_symbol"])))
    report_overlap(sig, sym_map_from_tables)

    print(f"\nNOTE: macrophage's full-dataset significant genes (n={len(sig['macrophage'])}) use all "
          f"5 donors for BOTH selection and reporting -- unlike fibroblast/enterocyte, this has NOT "
          f"been independently validated by the discovery/replication split (which found macrophage's "
          f"agreement indistinguishable from chance, p=0.25). Read macrophage's presence in the panel "
          f"below as unvalidated/exploratory, not confirmed.")

    print(f"\nloading {INPUT_H5AD} ...")
    raw_adata = anndata.read_h5ad(INPUT_H5AD)
    sym_to_id = dict(zip(raw_adata.var[GENE_NAME_COL].astype(str), raw_adata.var.index))
    target_id = sym_to_id[COEXPR_TARGET_GENE]
    sym_map = dict(zip(raw_adata.var.index, raw_adata.var[GENE_NAME_COL].astype(str)))

    heatmap_sets = {}
    for ct in COEXPR_CELL_TYPES:
        top = tables[ct][tables[ct]["fdr"] < FDR_THRESH].sort_values("pval").head(HEATMAP_TOP_N_PER_CT)
        heatmap_sets[ct] = set(top["gene_id"])
        print(f"  heatmap panel from {ct}: {len(heatmap_sets[ct])} of {len(sig[ct])} significant genes "
              f"(top {HEATMAP_TOP_N_PER_CT} by p-value)")

    union_ids = sorted(set().union(*heatmap_sets.values()) | {target_id})
    print(f"\nheatmap union panel: {len(union_ids)} genes ({len(union_ids)-1} significant + ANTXR2) "
          f"-- capped for legibility; overlap counts above and GSEA use the FULL significant sets")

    pairwise_results = compute_pairwise_correlations(raw_adata, union_ids, MIN_PERC_GROUP, SHRINKAGE, TRIM_PERCENT)
    donor_avg_2d = average_pairwise_across_donors(pairwise_results)

    os.makedirs(HT_DIR, exist_ok=True)
    plot_heatmaps(donor_avg_2d, union_ids, sym_map, target_id, HT_DIR)
    print(f"\nwrote heatmaps to {HT_DIR}/ ({{cell_type}}_pairwise_correlation.png)")


if __name__ == "__main__":
    main()
