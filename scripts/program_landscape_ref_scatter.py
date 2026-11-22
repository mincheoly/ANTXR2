"""Cross-dataset scatter per compartment: ANTXR2 median anchor percentile per
program, IBDverse (landscape_ref) vs Gut Cell Atlas (landscape_gca_ref), PC1
scheme, programs active (activity_sd >= 0.2) in at least half the compartment's
population x state cells in each dataset.
Writes landscape_gca_ref/report/{compartment}_pct_ibd_vs_gca.{png,csv} and
scatter_summary.csv.

Usage:
    python program_landscape_ref_scatter.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import OUT_DIR  # noqa: E402

IBD, GCA = OUT_DIR / "landscape_ref", OUT_DIR / "landscape_gca_ref"
REPORT = GCA / "report"
REPORT.mkdir(parents=True, exist_ok=True)
MIN_ACTIVITY_SD = 0.2

IFN = ["hallmark:INTERFERON_GAMMA_RESPONSE", "hallmark:INTERFERON_ALPHA_RESPONSE", "cytosig:IFNG"]
FIB = lambda p: p == "fibroblast (pooled)"   # noqa: E731
# compartment -> (IBDverse population test, GCA populations, programs to label)
# Lymphocytes dropped from the comparison (user, 2026-11-10).
COMPARTMENTS = {
    "epithelial": (lambda p: p.startswith(("Enterocyte", "Ileum goblet")),
                   ["Enterocyte", "TA", "Goblet", "Goblet_progenitor", "Goblet_cycling", "Paneth",
                    "Epithelial_stem"],
                   ["hallmark:E2F_TARGETS", "hallmark:MYC_TARGETS_V1", "collectri:E2F3", "collectri:HDAC1",
                    "collectri:ZNF331", "cytosig:TGFB3", "progeny:Hypoxia"]),
    "fibroblast": (FIB, ["fibroblast (pooled)"], IFN),                              # GCA postnatal
    "fibroblast_fetal": (FIB, ["fibroblast (fetal)"], IFN),
    "fibroblast_postnatal_fetal": (FIB, ["fibroblast (postnatal + fetal)"], IFN),
}
TITLES = {"fibroblast": "fibroblasts (atlas: postnatal)", "fibroblast_fetal": "fibroblasts (atlas: fetal)",
          "fibroblast_postnatal_fetal": "fibroblasts (atlas: postnatal + fetal)", "epithelial": "epithelium"}


def median_pct(d, keep):
    c = pd.read_csv(d / "coupling.csv")
    c = c[(c.scheme == "pc1") & c.population.map(keep) & c.pct.notna()]
    n = c.cell.nunique()
    g = c[c.activity_sd >= MIN_ACTIVITY_SD].groupby("unit").pct.agg(["median", "size"])
    return g.loc[g["size"] >= n / 2, "median"], n


SHORT = {"hallmark:INTERFERON_GAMMA_RESPONSE": "IFN-γ and IFN-α response",
         "hallmark:INTERFERON_ALPHA_RESPONSE": "", "cytosig:IFNG": "IFNG (CytoSig)"}


def short(u):
    return SHORT.get(u, u.split(":")[-1].replace("HALLMARK_", "").replace("_", " "))


def scatter(name, keep_ibd, pops_gca, labels):
    x, nx = median_pct(IBD, keep_ibd)
    y, ny = median_pct(GCA, lambda p: p in pops_gca)
    df = pd.concat([x.rename("ibd"), y.rename("gca")], axis=1).dropna()
    df.to_csv(REPORT / f"{name}_pct_ibd_vs_gca.csv")
    rs, rp = spearmanr(df.ibd, df.gca)[0], pearsonr(df.ibd, df.gca)[0]
    slope, icpt = np.polyfit(df.ibd, df.gca, 1)
    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=160)
    ax.axhline(50, color="#d0d0d0", lw=0.8, zorder=0)
    ax.axvline(50, color="#d0d0d0", lw=0.8, zorder=0)
    ax.plot([0, 100], [0, 100], color="#b0b0b0", lw=0.8, ls=(0, (4, 4)), zorder=0)
    ax.scatter(df.ibd, df.gca, s=14, color="#3b6fb6", alpha=0.55, lw=0)
    ax.plot([0, 100], [icpt, icpt + 100 * slope], color="#1f3f6e", lw=2)
    for u in [u for u in labels if u in df.index]:
        ax.scatter(df.ibd[u], df.gca[u], s=26, color="#1f3f6e", edgecolor="white", lw=1, zorder=3)
        right = df.ibd[u] > 75
        ax.annotate(short(u), (df.ibd[u], df.gca[u]), xytext=(-5 if right else 5, 5 if "and" in short(u) else -9 if right else 4),
                    ha="right" if right else "left", textcoords="offset points", fontsize=7.5, color="#222")
    ax.set(xlim=(0, 100), ylim=(0, 100), aspect="equal",
           xlabel=f"IBDverse: median percentile over {nx} population × state cells",
           ylabel=f"Gut Cell Atlas: median percentile over {ny} population(s)")
    ax.set_title(f"ANTXR2 coupling percentile per program, {TITLES.get(name, name)}\n"
                 f"{len(df)} programs · Spearman {rs:.2f} · Pearson {rp:.2f}", fontsize=10, loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=8.5)
    ax.xaxis.label.set_size(9)
    ax.yaxis.label.set_size(9)
    fig.tight_layout()
    fig.savefig(REPORT / f"{name}_pct_ibd_vs_gca.png")
    plt.close(fig)
    same = ((df.ibd > 50) == (df.gca > 50)).mean()
    return dict(compartment=name, ibd_cells=nx, gca_cells=ny, programs=len(df), spearman=rs, pearson=rp,
                slope=slope, same_side=same, sd_ibd=df.ibd.std(), sd_gca=df.gca.std())


if __name__ == "__main__":
    out = pd.DataFrame([scatter(k, *v) for k, v in COMPARTMENTS.items()])
    out.to_csv(REPORT / "scatter_summary.csv", index=False)
    print(out.round(3).to_string(index=False))
