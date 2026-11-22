"""Figure 1 (Act 1) tables from the donor units (fig1_aggregate.py).

Comparability rule (measured 2026-11-22, see fig1_platform.csv): Tabula Sapiens
10x 3' v3 and IBDverse read long genes high relative to the Gut Cell Atlas
(log-ratio vs genomic span, Spearman ~0.5 over the panel; ANTXR2, 224 kb, at the
91st-95th percentile), consistent with intronic reads being counted; GCA and the
dermal fibroblast atlas do not (rho ~0). Absolute ANTXR2 and any ratio with
ANTXR2 in it are therefore compared only within a class:
  intron-inclusive: tabula_sapiens (10x 3' v3), ibdverse
  exon-only:        gca_healthy, skin_fibroblast

Donor-equal weight throughout: a donor's value pools all of its cells in the
group (units of >= MIN_CELLS cells); the group value is the mean over donors.
Ratios are per donor (pooled UMI sums) and summarised as the geometric mean.

    -> /data/antxr2/fig1/tables/{platform, 1a_celltypes, 1bc_fibroblasts, 1d_gut_mrc2}.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from program_config import gene_spans

ROOT = Path("/data/antxr2/fig1")
TAB = ROOT / "tables"
TAB.mkdir(parents=True, exist_ok=True)
MIN_CELLS = 20
THIN_DONORS = 3
ADULT = {"18-34", "35-54", "47-80", "55-74", "75+"}
COL6 = ["COL6A1", "COL6A2", "COL6A3"]
# HFS involvement fixed from the clinical literature before looking at any value:
# skin (nodules, thickening), gingiva (hypertrophy), intestine (protein-losing
# enteropathy). Joints are involved but not sampled here.
HFS_TISSUES = {"Skin", "Skin (dermal atlas)", "Oral mucosa", "Small intestine", "Large intestine"}


def load():
    u = pd.read_parquet(ROOT / "donor_units.parquet")
    u = u[u.n_cells >= MIN_CELLS].copy()
    ex = u.extra.str.split("|", expand=True)
    u["x0"], u["x1"], u["x2"], u["x3"] = (ex[i] if i in ex else None for i in range(4))
    return u


def summarise(d, by):
    """donor value = pooled over the donor's units in the group; group = donor mean."""
    sums = [c for c in d.columns if c.startswith("sum_")]
    d = d.assign(**{f"w_{c[5:]}": d[c] * d.n_cells for c in d.columns if c.startswith("mean_")})
    w = [c for c in d.columns if c.startswith("w_")]
    don = d.groupby(by + ["donor"], observed=True)[sums + w + ["n_cells", "total"]].sum().reset_index()
    for c in w:
        don[f"cp10k_{c[2:]}"] = don[c] / don.n_cells             # mean per-cell CP10k
        don[f"raw_{c[2:]}"] = don[f"sum_{c[2:]}"] / don.n_cells   # mean UMI per cell
        don[f"det_{c[2:]}"] = np.nan
    don["col6_sum"] = don[[f"sum_{g}" for g in COL6]].sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        don["lr_col6_antxr2"] = np.log2(don.col6_sum / don.sum_ANTXR2)
        don["lr_antxr1_antxr2"] = np.log2(don.sum_ANTXR1 / don.sum_ANTXR2)
    don = don.replace([np.inf, -np.inf], np.nan)
    keep = [c for c in don.columns if c.startswith(("cp10k_", "raw_"))] + ["lr_col6_antxr2", "lr_antxr1_antxr2"]
    g = don.groupby(by, observed=True)
    out = g[keep].mean()
    for c in ["lr_col6_antxr2", "lr_antxr1_antxr2"]:
        out[c.replace("lr_", "gm_")] = 2 ** out[c]
        out[c + "_sd"] = g[c].std()
    out["cp10k_ANTXR2_sd"] = g["cp10k_ANTXR2"].std()
    out["n_donors"] = g.size()
    out["n_cells"] = g.n_cells.sum()
    return out.reset_index(), don


def platform(u):
    span = gene_spans()
    genes = [c[4:] for c in u.columns if c.startswith("sum_")]
    fib = u[u.cell_type.str.contains("fibroblast")]

    def pooled(a, assay):
        d = fib[(fib.atlas == a) & (fib.assay == assay)]
        return pd.Series(d[[f"sum_{g}" for g in genes]].sum().values / d.total.sum() * 1e4, index=genes)

    ref = pooled("gca_healthy", "10x 3' v3")
    rows = []
    for a, assay in [("tabula_sapiens", "10x 3' v3"), ("ibdverse", "10x 3' v3"), ("tabula_sapiens", "10x 5' v2"),
                     ("skin_fibroblast", "10x 3' v3"), ("skin_fibroblast", "10x 3' v2"), ("gca_healthy", "10x 3' v2"),
                     ("cross_tissue_immune", "10x 3' v3")]:
        x = pooled(a, assay)
        ok = (x > 0.05) & (ref > 0.05) & x.index.isin(span.index)
        if ok.sum() < 10:
            continue
        lr = np.log2(x[ok] / ref[ok])
        r, p = spearmanr(lr, np.log10(span[lr.index]))
        rows.append(dict(atlas=a, assay=assay, vs="gca_healthy 10x 3' v3 fibroblasts", n_genes=int(ok.sum()),
                         rho_span=r, p=p, antxr2_log2_ratio=lr.get("ANTXR2"),
                         antxr2_pct=(lr < lr["ANTXR2"]).mean() * 100))
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "platform.csv", index=False)
    print(t.round(3).to_string(index=False))


def fig1a(u):
    ts = u[(u.atlas == "tabula_sapiens") & (u.assay == "10x 3' v3")]
    t, _ = summarise(ts, ["cell_type"])
    comp = ts.groupby(["cell_type", "x0"], observed=True).n_cells.sum().reset_index() \
        .sort_values("n_cells").drop_duplicates("cell_type", keep="last").set_index("cell_type").x0
    t["compartment"] = t.cell_type.map(comp)
    t["thin"] = t.n_donors < THIN_DONORS
    t = t.sort_values("cp10k_ANTXR2", ascending=False).reset_index(drop=True)
    t.insert(0, "rank", t.index + 1)
    # replication: same ontology labels in GCA healthy (adult, exon-only) -- ranks only
    g = u[(u.atlas == "gca_healthy") & u.x1.isin(ADULT)]
    gt, _ = summarise(g, ["cell_type"])
    m = t.merge(gt[["cell_type", "cp10k_ANTXR2", "n_donors"]], on="cell_type", suffixes=("", "_gca"))
    m = m[(m.n_donors >= THIN_DONORS) & (m.n_donors_gca >= THIN_DONORS)]
    r, p = spearmanr(m.cp10k_ANTXR2, m.cp10k_ANTXR2_gca)
    t["cp10k_ANTXR2_gca"] = t.cell_type.map(gt.set_index("cell_type").cp10k_ANTXR2)
    t.to_csv(TAB / "1a_celltypes.csv", index=False)
    print(f"1a: {len(t)} TS cell types ({t.thin.sum()} thin, <{THIN_DONORS} donors); "
          f"rank agreement with GCA on {len(m)} shared labels: Spearman {r:.2f} (p={p:.1e})")
    print(t.groupby("compartment").cp10k_ANTXR2.median().round(3).to_string())
    print(t[["rank", "cell_type", "compartment", "cp10k_ANTXR2", "n_donors", "n_cells"]].head(15).to_string(index=False))


def fib_rows(u):
    """Fibroblast populations by tissue, tagged with comparability class."""
    ts = u[(u.atlas == "tabula_sapiens") & (u.assay == "10x 3' v3") & u.cell_type.str.contains("fibroblast")
           & (u.x0 == "Stromal")].assign(group=lambda d: d.x1.str.replace("_", " ").str.capitalize(), cls="intron")
    gca = u[(u.atlas == "gca_healthy") & u.x1.isin(ADULT) & (u.x3 == "Fibroblast")] \
        .assign(group=lambda d: d.tissue.str.replace("_", " ").str.capitalize(), cls="exon")
    skin = u[(u.atlas == "skin_fibroblast") & (u.disease == "normal") & u.assay.str.startswith("10x 3'")
             & ~u.x1.str.contains("trimester|post-fertilization|fetal|unknown", na=True)] \
        .assign(group="Skin (dermal atlas)", cls="exon")
    ibd = u[(u.atlas == "ibdverse") & (u.disease == "normal") & (u.x1 == "Mesenchymal")
            & u.annot.str.contains("ibroblast")].assign(group=lambda d: "IBDverse " + d.tissue, cls="intron")
    return pd.concat([ts, ibd, gca, skin])


def fig1bc(u):
    d = fib_rows(u)
    t, don = summarise(d, ["cls", "atlas", "group"])
    t["hfs"] = t.group.str.replace("IBDverse ileum", "Small intestine").str.replace("IBDverse rectum", "Large intestine") \
        .isin(HFS_TISSUES)
    t = t[t.n_donors >= 2].sort_values(["cls", "atlas", "cp10k_ANTXR2"], ascending=[False, False, False])
    t.to_csv(TAB / "1bc_fibroblasts.csv", index=False)
    don.to_csv(TAB / "1bc_fibroblasts_donors.csv", index=False)
    cols = ["cls", "atlas", "group", "hfs", "n_donors", "n_cells", "cp10k_ANTXR2", "gm_col6_antxr2",
            "gm_antxr1_antxr2", "cp10k_MRC2"]
    print(t[cols].round(2).to_string(index=False))


def fig1d(u):
    ibd = u[(u.atlas == "ibdverse") & (u.disease == "normal") & (u.tissue != "blood")
            & u.x1.isin(["Stem", "Enterocyte", "Colonocyte", "Secretory", "Mesenchymal", "Myeloid"])] \
        .assign(group=lambda d: d.x1, cls="intron")
    gca = u[(u.atlas == "gca_healthy") & u.x1.isin(ADULT) & u.tissue.isin(["Small_intestine", "Large_intestine"])
            & u.x3.isin(["Epithelial_stem", "Transit_amplifying", "Absorptive", "Secretory", "Enteroendocrine",
                         "Fibroblast", "Macrophage"])].assign(group=lambda d: d.x3, cls="exon")
    t, don = summarise(pd.concat([ibd, gca]), ["cls", "atlas", "group"])
    # detection: fraction of cells with >=1 UMI, cell-weighted within donor, donor mean
    d = pd.concat([ibd, gca])
    for gname in ["ANTXR2", "MRC2"]:
        d[f"wd_{gname}"] = d[f"det_{gname}"] * d.n_cells
    dd = d.groupby(["atlas", "group", "donor"])[["wd_ANTXR2", "wd_MRC2", "n_cells"]].sum()
    for gname in ["ANTXR2", "MRC2"]:
        dd[f"det_{gname}"] = dd[f"wd_{gname}"] / dd.n_cells
    det = dd.groupby(["atlas", "group"])[["det_ANTXR2", "det_MRC2"]].mean().reset_index()
    t = t.merge(det, on=["atlas", "group"])
    t.to_csv(TAB / "1d_gut_mrc2.csv", index=False)
    don.to_csv(TAB / "1d_gut_mrc2_donors.csv", index=False)
    print(t[["atlas", "group", "n_donors", "raw_ANTXR2", "raw_MRC2", "det_ANTXR2", "det_MRC2",
             "cp10k_ANTXR2", "cp10k_MRC2"]].round(3).to_string(index=False))


if __name__ == "__main__":
    u = load()
    platform(u)
    fig1a(u)
    fig1bc(u)
    fig1d(u)
