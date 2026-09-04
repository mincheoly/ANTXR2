"""Gene set enrichment (Enrichr over-representation, via gseapy) on each cell
type's top-N ANTXR2-correlated gene list from coexpression_pipeline.py, plus
a detailed overlap report between cell types' top-gene lists (not just the
counts already stored in memento_correlations['top_gene_overlap']).

This is over-representation analysis (ORA / Enrichr's hypergeometric test on
a fixed gene set against the whole-genome background), not ranked GSEA
(Subramanian et al.) -- appropriate here because the input is a fixed top-N
gene list, not a fully ranked gene list with a running-sum statistic. Uses
Enrichr's default whole-genome background, not the memento-filtered ~4,009-
gene universe those genes were selected from -- a known limitation, noted in
the output.

Requires internet access to Enrichr (maayanlab.cloud). Run in the `antxr2`
conda env: `conda run -n antxr2 python scripts/coexpression_enrichment.py`.
"""
import itertools
import os

import anndata
import gseapy as gp
import pandas as pd

from config import COEXPR_CELL_TYPES, COEXPR_FIGURES_DIR, COEXPR_OUTPUT_H5AD

GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
ADJ_P_THRESH = 0.05
OUT_DIR = os.path.join(COEXPR_FIGURES_DIR, "enrichment")


def load_top_genes():
    a = anndata.read_h5ad(COEXPR_OUTPUT_H5AD)
    mc = a.uns["memento_correlations"]
    return {ct: mc["donor_averaged"][ct]["top_genes"].copy()
            for ct in COEXPR_CELL_TYPES if ct in mc["donor_averaged"]}


def run_enrichment(genes, name):
    print(f"\n  querying Enrichr ({', '.join(GENE_SETS)}) for {name} ({len(genes)} genes)...")
    try:
        enr = gp.enrichr(gene_list=genes, gene_sets=GENE_SETS, organism="human", outdir=None)
    except Exception as e:
        print(f"  ENRICHMENT FAILED for {name}: {e}")
        return pd.DataFrame()
    return enr.results


def summarize(df, name, top_n_terms=15):
    if df.empty:
        return df
    sig = df[df["Adjusted P-value"] < ADJ_P_THRESH].sort_values("Adjusted P-value").reset_index(drop=True)
    print(f"\n=== {name}: {len(sig)} significant terms (adj p < {ADJ_P_THRESH}) "
          f"across {df['Gene_set'].nunique()} libraries ({len(df)} terms tested) ===")
    for gs in GENE_SETS:
        sub = sig[sig["Gene_set"] == gs].head(top_n_terms)
        print(f"\n  -- {gs} (top {len(sub)} of {len(sig[sig['Gene_set']==gs])} significant) --")
        for _, row in sub.iterrows():
            print(f"    {row['Term']:65s} adj_p={row['Adjusted P-value']:.2e}  "
                  f"overlap={row['Overlap']:6s}  genes={row['Genes']}")
    return sig


def overlap_report(gene_sets):
    cts = list(gene_sets)
    print("\n\n=== Overlap between top-gene lists (gene symbols, not just counts) ===")
    for a, b in itertools.combinations(cts, 2):
        shared = gene_sets[a] & gene_sets[b]
        print(f"\n  {a} vs {b}: {len(shared)} shared / {len(gene_sets[a])} and {len(gene_sets[b])}")
        if shared:
            print(f"    {sorted(shared)}")
    if len(cts) >= 3:
        for combo in itertools.combinations(cts, 3):
            triple = set.intersection(*[gene_sets[c] for c in combo])
            print(f"\n  {' & '.join(combo)} (all three): {len(triple)} shared")
            if triple:
                print(f"    {sorted(triple)}")
    union = set.union(*gene_sets.values())
    unique_per_ct = {ct: gene_sets[ct] - set.union(*(gene_sets[o] for o in cts if o != ct)) for ct in cts}
    print(f"\n  union across all cell types: {len(union)} genes")
    for ct in cts:
        print(f"  unique to {ct} only: {len(unique_per_ct[ct])} / {len(gene_sets[ct])}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    top_genes = load_top_genes()

    all_sig = {}
    for ct, df in top_genes.items():
        genes = df["gene_symbol"].tolist()
        print(f"\n\n######## {ct} ({len(genes)} genes) ########")
        res = run_enrichment(genes, ct)
        if res.empty:
            continue
        res.to_csv(os.path.join(OUT_DIR, f"{ct}_enrichr_all.csv"), index=False)
        sig = summarize(res, ct)
        sig.to_csv(os.path.join(OUT_DIR, f"{ct}_enrichr_significant.csv"), index=False)
        all_sig[ct] = sig

    gene_sets = {ct: set(df["gene_symbol"]) for ct, df in top_genes.items()}
    overlap_report(gene_sets)

    print(f"\n\nwrote per-cell-type Enrichr results to {OUT_DIR}/")


if __name__ == "__main__":
    main()
