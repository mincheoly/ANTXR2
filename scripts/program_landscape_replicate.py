"""Proliferation tilt, step 4: the pre-registered replication tests in the
GCA landscape (planning_summaries/program_coexpression.md, "Step 4 plan").

Reads the IBDverse landscape report (OUT_DIR/landscape/report) and the GCA
one (OUT_DIR/landscape_gca/report, after ANTXR2_LANDSCAPE=gca
program_landscape.py all + program_landscape_report.py, and
program_landscape_cycling.py for P3).

    P1  per lineage, each named program: GCA consensus z < 0 and ANTXR2 in
        the bottom P1_RANK of the lineage's genes (ANTXR2 + 50 anchors)
    P2  per lineage, the named IBDverse proliferation set (PLAN; scored in
        GCA): mean GCA consensus z vs N_DRAWS random same-size sets of
        GCA-scored library programs, one-sided "less"
    P3  cells whose top cycle quartile is >= 10 points richer in MKI67+ cells:
        ANTXR2 excess < 0 count, two-sided sign test

Epithelium is primary; lymphocytes and the positive side are secondary.
Writes OUT_DIR/landscape_gca/replication.json and prints the tables.

Usage:
    python program_landscape_replicate.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import OUT_DIR, SEED          # noqa: E402

IBD, GCA = OUT_DIR / "landscape", OUT_DIR / "landscape_gca"
N_DRAWS = 10_000
P1_RANK = 5
PLAN = {
    "epithelial": dict(p1=["hallmark:E2F_TARGETS", "hallmark:MYC_TARGETS_V1"],
                       set=["collectri:E2F2", "collectri:E2F3", "collectri:E2F4",
                            "hallmark:E2F_TARGETS", "hallmark:MYC_TARGETS_V1",
                            "hallmark:DNA_REPAIR", "hallmark:G2M_CHECKPOINT"]),
    "lymphocyte": dict(p1=["hallmark:MTORC1_SIGNALING", "hallmark:MYC_TARGETS_V1"],
                       set=["hallmark:MTORC1_SIGNALING", "hallmark:E2F_TARGETS",
                            "hallmark:MYC_TARGETS_V1", "hallmark:MYC_TARGETS_V2"]),
}
POSITIVE = ["cytosig:TGFB3", "progeny:Hypoxia", "collectri:ZNF331"]


def consensus(root, lin):
    f = root / "report" / f"consensus_{lin}.csv"
    if not f.exists():                       # < 3 filled cells: no consensus
        return None
    c = pd.read_csv(f, index_col=0)
    return c[~c.index.str.startswith(("discovery:", "bracq:"))]


def lineage(lin, plan):
    ibd, gca = consensus(IBD, lin), consensus(GCA, lin)
    if gca is None:
        return dict(verdict="untestable", note="fewer than 3 filled GCA cells"), None
    missing = [u for u in plan["p1"] + plan["set"] if u not in gca.index]
    if set(plan["p1"]) <= set(missing) and set(plan["set"]) <= set(missing):
        return dict(verdict="untestable", note="no named program scored in GCA",
                    n_programs_gca=len(gca), programs_gca=list(gca.index)), None
    p1 = {u: dict(consensus_z=float(gca.consensus_z.get(u, np.nan)),
                  rank=int(gca.rank_among_genes[u]) if u in gca.index else None,
                  n_genes=int(gca.n_genes.iloc[0]),
                  passed=bool(u in gca.index and gca.consensus_z[u] < 0
                              and gca.rank_among_genes[u] <= P1_RANK))
          for u in plan["p1"]}
    S_ibd = plan["set"]
    S = [u for u in S_ibd if u in gca.index]
    z = gca.consensus_z.dropna()
    obs = z[S].mean()
    rng = np.random.default_rng(SEED)
    null = np.array([rng.choice(z.values, len(S), replace=False).mean() for _ in range(N_DRAWS)])
    p2 = dict(set=S, missing_in_gca=[u for u in S_ibd if u not in gca.index], mean_z=float(obs),
              null_mean=float(null.mean()), p=float((1 + (null <= obs).sum()) / (N_DRAWS + 1)))
    p2["passed"] = p2["p"] < 0.05
    pos = {u: float(gca.consensus_z.get(u, np.nan)) for u in POSITIVE}
    both = ibd.join(gca, lsuffix="_ibd", rsuffix="_gca", how="inner")
    rho = both.consensus_z_ibd.corr(both.consensus_z_gca, method="spearman")
    if all(r["passed"] for r in p1.values()) and p2["passed"]:
        verdict = "replicated"
    elif (all(r["consensus_z"] < 0 for r in p1.values()) and obs < 0) \
            or any(r["passed"] for r in p1.values()):
        verdict = "partial"
    else:
        verdict = "not replicated"
    table = both.loc[[u for u in S_ibd + POSITIVE if u in both.index],
                     ["consensus_z_ibd", "rank_among_genes_ibd", "consensus_z_gca",
                      "rank_among_genes_gca"]]
    return dict(p1=p1, p2=p2, positive=pos, profile_spearman=float(rho),
                n_programs_both=len(both), verdict=verdict), table


def p3():
    q = pd.read_csv(GCA / "cycling" / "quartiles_cycle.csv")
    r = q[q.mki67_det_q4 - q.mki67_det_q1 >= 0.10]
    n_neg = int((r.excess < 0).sum())
    return dict(n_cells=len(r), n_negative=n_neg, median_z=float(r.excess_z.median()),
                p=float(binomtest(n_neg, len(r)).pvalue) if len(r) else None,
                cells=r[["population", "state", "excess", "excess_z", "donors_below"]]
                .round(3).to_dict("records"))


def main():
    pd.set_option("display.width", 250)
    res = {}
    for lin, plan in PLAN.items():
        res[lin], table = lineage(lin, plan)
        print(f"== {lin}: {res[lin]['verdict']}")
        print(json.dumps({k: v for k, v in res[lin].items() if k != "p2"}, indent=1))
        if "p2" in res[lin]:
            print(json.dumps(res[lin]["p2"], indent=1))
        if table is not None:
            print(table.round(2).to_string())
    if (GCA / "cycling" / "quartiles_cycle.csv").exists():
        res["p3"] = p3()
        print("P3:", json.dumps({k: v for k, v in res["p3"].items() if k != "cells"}))
    (GCA / "replication.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
