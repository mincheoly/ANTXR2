"""Gene-level ANTXR2 co-expression table within lymphatic endothelium.

Pathway enrichment is not used here: at ~23 donor x subtype groups the
question is which individual genes track ANTXR2, with each gene's own
specificity control.

Per partner gene g (pooled within-LEC scan, coexpr_point_scan.py output with
ANTXR2 + 30 matched anchors):

  antxr2_r      ANTXR2-g correlation, residualized on g's expression exactly
                as in coexpr_pathway_gsea.py (running median over genes sorted
                by raw mean, computed per target)
  anchor_mean/sd  the same residualized correlation of g with each anchor
  z             (antxr2_r - anchor_mean) / anchor_sd: does ANTXR2 track g
                more than an expression-matched arbitrary gene does
  r_gut/r_ts/r_site  ANTXR2's residualized correlation in each source alone.
                Only the gut groups hold lymphatic subtype fixed (donor x
                subtype); ts and site groups are donor x tissue and can mix
                capillary / collecting / valve LECs, so a gene supported only
                there may be composition, not co-regulation.
  lec6_log2fc   donor-paired LEC6 vs other-LEC fold change (gut): a large
                value means g marks the ANTXR2-high subtype, the composition
                route by which it could correlate with ANTXR2.

Annotation (gene name, summary) from mygene.info, cached.

Usage: python lymphatic_gene_table.py <lymphatic dir> [scan_dir] [by_source_dir] [out_suffix]
"""
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from endothelial_source_check import residualize   # noqa: E402

SOURCES = ["gut", "ts", "site"]


def wide(scan_csv):
    d = residualize(pd.read_csv(scan_csv))
    return d.pivot_table(index="gene", columns="target", values="resid")


def annotate(genes, cache):
    ann = json.load(open(cache)) if os.path.exists(cache) else {}
    todo = [g for g in genes if g not in ann]
    for i in range(0, len(todo), 200):
        body = ("q=" + ",".join(todo[i:i + 200]) +
                "&scopes=symbol&fields=name,summary&species=human").encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                     headers={"Content-Type": "application/x-www-form-urlencoded"})
        for hit in json.load(urllib.request.urlopen(req, timeout=60)):
            q = hit["query"]
            if q not in ann or "notfound" in ann[q]:
                ann[q] = {"name": hit.get("name", ""), "summary": hit.get("summary", "")} \
                    if not hit.get("notfound") else {"notfound": True}
    json.dump(ann, open(cache, "w"))
    return ann


def main(lym, scan_dir="scan", src_dir="scan_by_source", suffix=""):
    w = wide(os.path.join(lym, scan_dir, "lymphatic_pooled.csv"))
    anchors = [c for c in w.columns if c != "ANTXR2"]
    t = pd.DataFrame({"antxr2_r": w["ANTXR2"], "anchor_mean": w[anchors].mean(1),
                      "anchor_sd": w[anchors].std(1), "n_anchors": w[anchors].notna().sum(1)})
    t["z"] = (t.antxr2_r - t.anchor_mean) / t.anchor_sd
    for s in SOURCES:
        ws = wide(os.path.join(lym, src_dir, f"lymphatic_{s}.csv"))
        t[f"r_{s}"] = ws["ANTXR2"].reindex(t.index)
    t["n_sources_same_sign"] = (np.sign(t[[f"r_{s}" for s in SOURCES]]).eq(np.sign(t.antxr2_r), axis=0)).sum(1)
    lec6 = pd.read_csv(os.path.join(lym, "LEC6_vs_otherLEC_paired.csv")).set_index("gene")
    t["lec6_log2fc"] = lec6.log2fc.reindex(t.index)
    t = t.dropna(subset=["z"]).sort_values("z", ascending=False)
    t.to_csv(os.path.join(lym, f"antxr2_lymphatic_gene_table{suffix}.csv"))

    # candidates: specific vs anchors (|z| >= 3), same sign in all three
    # sources, and same sign within the subtype-fixed gut groups
    cand = t[(t.z.abs() >= 3) & (t.n_sources_same_sign == 3)]
    ann = annotate(list(cand.index), os.path.join(lym, "mygene_cache.json"))
    cand = cand.assign(name=[ann.get(g, {}).get("name", "") for g in cand.index],
                       summary=[(ann.get(g, {}).get("summary", "") or "")[:160] for g in cand.index])
    cand.to_csv(os.path.join(lym, f"antxr2_lymphatic_gene_candidates{suffix}.csv"))
    pd.set_option("display.width", 260)
    pd.set_option("display.max_colwidth", 60)
    cols = ["z", "antxr2_r", "anchor_mean", "r_gut", "r_ts", "r_site", "lec6_log2fc", "name"]
    print(f"{len(t)} genes; |z|>=3 & same sign in all 3 sources: {len(cand)} "
          f"({(cand.z > 0).sum()} positive / {(cand.z < 0).sum()} negative)")
    print(f"null expectation for |z|>=3 alone: ~{0.0027 * len(t):.0f} genes "
          f"(observed |z|>=3 before replication filter: {(t.z.abs() >= 3).sum()})")
    print(cand[cand.z > 0][cols].head(40).round(3).to_string())
    print(cand[cand.z < 0].sort_values("z")[cols].head(25).round(3).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:])
