"""Proliferation tilt, steps 1 and 3: is ANTXR2 lower in cycling cells
(--axis cycle) or in high-growth cells (--axis growth) at the mean level?
(planning_summaries/program_coexpression.md, "Proliferation tilt").

The landscape's clearest lineage-wide signal is that ANTXR2 couples less than
matched genes to E2F / MYC / G2M / DNA-repair programs within populations.
This asks the same thing on means, where it is directly interpretable.

A. Within population x state, per donor (>= MIN_GROUP_CELLS cells; the
   fibroblast floor for pooled fibroblasts), cells are split into quartiles of
   a cell-cycle score (Tirosh et al. 2016 S + G2M genes; mean per-gene z of
   log-normalised expression within the cell). For ANTXR2 and each of the
   population's 50 landscape anchors (span + raw mean + detection matched):
   lfc = log2 mean(X/sf) in Q4 over Q1, averaged over donors. --axis growth
   uses the hallmark MTORC1_SIGNALING + MYC_TARGETS_V1 genes minus the cycle
   genes instead; --no-mki67 first removes cells with any MKI67 count, so
   the axis is read among non-dividing cells. --all-pops (exploratory) also
   takes population x state cells the landscape left unfilled because ANTXR2
   is below the raw-mean floor, with anchors matched here (span window + raw
   mean + detection, candidates above LOW_FLOOR) instead of the landscape's.
   excess = lfc_ANTXR2 - mean lfc_anchor, z = excess / anchor SD; donors =
   share of donors where ANTXR2's lfc is below the median anchor's.
B. (IBDverse, cycle axis) Cycling vs quiescent stem cells as author types:
   "Enterocyte stem MKI67" vs "Enterocyte stem OLFM4lo LGR5hi", per donor
   with >= MIN_GROUP_CELLS of each, size factors computed jointly, anchors
   matched across both.

Writes LANDSCAPE_OUT/cycling/quartiles_<axis>[_nomki67][_allpops].csv and
stem_pair.csv.
ANTXR2_LANDSCAPE=gca runs it on the replication atlas.

Usage:
    python program_landscape_cycling.py [--axis {cycle,growth}] [--no-mki67] [--all-pops]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import (                       # noqa: E402
    LANDSCAPE_ATLAS, LANDSCAPE_CELLS, LANDSCAPE_FIB, LANDSCAPE_FIB_FLOOR, LANDSCAPE_LINEAGES,
    LANDSCAPE_OUT, LANDSCAPE_STATES, MIN_GROUP_CELLS, N_ANCHORS, TARGET, gene_spans, pop_slug,
)
from program_scores import gene_stats, load_cells, load_library, size_factors  # noqa: E402

OUT = LANDSCAPE_OUT / "cycling"
OUT.mkdir(parents=True, exist_ok=True)
STATE_OF = {v: k for k, v in LANDSCAPE_STATES.items()}

# Tirosh et al. 2016 (Seurat cc.genes.updated.2019 symbols)
S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
           "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "RPA2", "NASP",
           "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
           "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
           "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
G2M_GENES = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
             "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2",
             "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
             "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
             "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
             "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]
STEM_PAIR = ("Enterocyte stem MKI67", "Enterocyte stem OLFM4lo LGR5hi")
GROWTH_PROGRAMS = ["hallmark:MTORC1_SIGNALING", "hallmark:MYC_TARGETS_V1"]
LOW_FLOOR = 0.01              # --all-pops anchor candidates: raw mean floor


def axis_genes(axis):
    if axis == "cycle":
        return S_GENES + G2M_GENES
    lib = load_library()
    g = set(lib[lib.program_id.isin(GROWTH_PROGRAMS)].gene)
    return sorted(g - set(S_GENES + G2M_GENES))


def normalise(X, genes):
    sf = size_factors(X, genes)
    return (sp.diags(1 / sf) @ X).tocsr(), np.median(sf)


def axis_score(Y, scale, genes, rows, gene_list):
    """Mean per-gene z of log1p(Y * scale) over the axis genes detected in
    >= 1% of `rows`."""
    idx = [genes.get_loc(g) for g in gene_list if g in genes]
    L = Y[rows][:, idx].toarray() * scale
    keep = (L > 0).mean(0) >= 0.01
    L = np.log1p(L[:, keep])
    return ((L - L.mean(0)) / (L.std(0) + 1e-9)).mean(1)


def lfc(Y, rows_hi, rows_lo, T):
    hi = np.asarray(Y[rows_hi][:, T].mean(0)).ravel()
    lo = np.asarray(Y[rows_lo][:, T].mean(0)).ravel()
    c = 0.05 * (hi + lo) / 2 + 1e-12                     # relative pseudocount
    return np.log2((hi + c) / (lo + c))


def summarise(L, **extra):
    """L: donors x targets log2 ratios, column 0 = ANTXR2."""
    a, anc = L[:, 0], L[:, 1:]
    m = anc.mean(0)
    excess = a.mean() - m.mean()
    return dict(extra, n_donors=len(L), antxr2_lfc=a.mean(), anchor_lfc=m.mean(),
                anchor_sd=m.std(ddof=1), excess=excess, excess_z=excess / m.std(ddof=1),
                antxr2_pct=(m < a.mean()).mean() * 100,
                donors_below=(a < np.median(anc, 1)).mean())


def matched_here(X, genes, rows):
    from program_landscape import span_matched   # binds LANDSCAPE_OUT; import late
    raw_mean, det = gene_stats(X, rows)
    f = pd.DataFrame({"log10_mean": np.log10(raw_mean + 1e-6), "det": det}, index=genes)
    anchors, _ = span_matched({"c": f}, [set(genes[raw_mean >= LOW_FLOOR])],
                              np.log10(gene_spans()), N_ANCHORS)
    return anchors


def quartiles(axis, no_mki67, all_pops):
    gene_list = axis_genes(axis)
    meta = json.loads((LANDSCAPE_OUT / "score_meta.json").read_text())
    cells = pd.read_csv(LANDSCAPE_OUT / "cells.csv")
    lineage_of = {p: lin for lin, ps in LANDSCAPE_LINEAGES.items() for p in ps}
    rows = []
    for pop in sorted(lineage_of if all_pops else meta):
        X, obs, genes = load_cells(LANDSCAPE_CELLS / f"cells_{pop_slug(pop)}.npz")
        Y, scale = normalise(X, genes)
        mki = genes.get_loc("MKI67")
        floor = LANDSCAPE_FIB_FLOOR if pop == LANDSCAPE_FIB else MIN_GROUP_CELLS
        state = obs.disease.map(STATE_OF)
        pc = cells[cells.population == pop]
        for st, filled in zip(pc.state, pc.filled):
            if not (filled or all_pops):
                continue
            anchors = (meta[pop]["anchors"] if filled
                       else matched_here(X, genes, np.where(state == st)[0]))
            T = np.array([genes.get_loc(g) for g in [TARGET] + anchors])
            L, mk = [], []
            ok = state == st
            if no_mki67:
                ok &= X[:, mki].toarray().ravel() == 0
            for d, ix in obs[ok].groupby("donor").groups.items():
                ix = np.asarray(ix)
                if len(ix) < floor:
                    continue
                s = axis_score(Y, scale, genes, ix, gene_list)
                q = pd.qcut(s, 4, labels=False, duplicates="drop")
                hi, lo = ix[q == q.max()], ix[q == 0]
                L.append(lfc(Y, hi, lo, T))
                mk.append([(X[hi, mki] > 0).mean(), (X[lo, mki] > 0).mean()])
            if len(L) < 3:
                continue
            mk = np.array(mk).mean(0)
            rows.append(summarise(np.array(L), population=pop, lineage=lineage_of[pop],
                                  state=st, mki67_det_q4=mk[0], mki67_det_q1=mk[1],
                                  anchors="landscape" if filled else "matched_here",
                                  antxr2_raw_mean=pc.antxr2_raw_mean[pc.state == st].iloc[0]))
        print(f"{pop}: done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / f"quartiles_{axis}{'_nomki67' if no_mki67 else ''}"
                     f"{'_allpops' if all_pops else ''}.csv", index=False)
    return out


def stem_pair():
    from program_landscape import span_matched   # binds LANDSCAPE_OUT; import late
    cyc, qui = STEM_PAIR
    parts = [load_cells(LANDSCAPE_CELLS / f"cells_{pop_slug(p)}.npz") for p in STEM_PAIR]
    assert parts[0][2].equals(parts[1][2])
    genes = parts[0][2]
    X = sp.vstack([p[0] for p in parts]).tocsr()
    obs = pd.concat([p[1] for p in parts], ignore_index=True)
    obs["state"] = obs.disease.map(STATE_OF)
    Y, _ = normalise(X, genes)
    span = np.log10(gene_spans())
    rows = []
    for st, o in obs[obs.state.notna()].groupby("state"):
        n = o.groupby(["donor", "cell_type"]).size().unstack(fill_value=0)
        donors = n.index[(n[cyc] >= MIN_GROUP_CELLS) & (n[qui] >= MIN_GROUP_CELLS)]
        if len(donors) < 3:
            continue
        oo = o[o.donor.isin(donors)]
        feats, passing = {}, []
        for p in STEM_PAIR:                 # anchors matched in both types
            r = oo.index[oo.cell_type == p].values
            raw_mean, det = gene_stats(X, r)
            feats[p] = pd.DataFrame({"log10_mean": np.log10(raw_mean + 1e-6), "det": det},
                                    index=genes)
            passing.append(set(genes[raw_mean >= 0.07]))
        anchors, _ = span_matched(feats, passing, span, N_ANCHORS)
        T = np.array([genes.get_loc(g) for g in [TARGET] + anchors])
        L = np.array([lfc(Y, oo.index[(oo.donor == d) & (oo.cell_type == cyc)].values,
                          oo.index[(oo.donor == d) & (oo.cell_type == qui)].values, T)
                      for d in donors])
        rows.append(summarise(L, state=st, comparison=f"{cyc} / {qui}"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "stem_pair.csv", index=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", choices=["cycle", "growth"], default="cycle")
    ap.add_argument("--no-mki67", action="store_true")
    ap.add_argument("--all-pops", action="store_true")
    a = ap.parse_args()
    pd.set_option("display.width", 250)
    print(quartiles(a.axis, a.no_mki67, a.all_pops).round(3).to_string())
    if a.axis == "cycle" and not (a.no_mki67 or a.all_pops) and LANDSCAPE_ATLAS == "ibdverse":
        print(stem_pair().round(3).to_string())


if __name__ == "__main__":
    main()
