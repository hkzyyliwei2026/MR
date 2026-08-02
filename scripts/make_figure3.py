# -*- coding: utf-8 -*-
"""
Figure 3 - instrument-threshold instability (two panels).

A: rank-rank scatter of P values under the two thresholds (conduction disorders), marking the 5 traits that reach FDR significance under the genome-wide threshold
B: distribution of instrument counts (nIV) per trait under both thresholds, annotating the nIV of the genome-wide-threshold FDR signals

Reads only from derived/ and results/tables/; writes to results/figures/.
"""
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
FIG = PROJ / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

OUTCOME = "Cardiac conduction disorders"

# --- Panel A data: ranks of the traits analysable under both thresholds ---
inst = pd.read_csv(PROJ / "derived" / "threshold_instability.csv")
inst = inst[inst["outcome"] == OUTCOME].copy()
inst["rank_main"] = inst["p_main"].rank(method="average")
inst["rank_gws"] = inst["p_gws"].rank(method="average")
rho = inst["rank_main"].corr(inst["rank_gws"], method="pearson")  # Pearson correlation of ranks = Spearman

cross = pd.read_csv(PROJ / "derived" / "table1_threshold_crosswalk.csv")
cross = cross[cross["outcome"] == OUTCOME]
sig_traits = set(cross["trait"])
inst["sig"] = inst["trait"].isin(sig_traits)

# --- Panel B data: nIV distribution under the genome-wide threshold (full GWS table, n=607, matching Section 3.4) ---
gws = pd.read_csv(PROJ / "results" / "tables" / "MR_immune_GWS_CONDUCTIO.csv")
n_gws = len(gws)
# Primary-threshold instrument counts, plotted alongside so the contrast is visible in the figure
main = pd.read_csv(PROJ / "results" / "tables" / "MR_immune_CONDUCTIO.csv")
n_main = len(main)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))
plt.rcParams.update({
    "font.family": "Arial",
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
})

def lowess_line(x, y, frac=0.28):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    n = len(xs)
    k = max(20, int(np.ceil(frac * n)))
    fitted = np.empty(n)
    for i, x0 in enumerate(xs):
        dist = np.abs(xs - x0)
        h = np.partition(dist, min(k, n - 1))[min(k, n - 1)]
        if h == 0:
            fitted[i] = np.mean(ys[dist == 0])
            continue
        w = (1 - (dist / h) ** 3) ** 3
        w[dist > h] = 0
        X = np.column_stack([np.ones(n), xs - x0])
        XtW = X.T * w
        try:
            beta = np.linalg.pinv(XtW @ X) @ (XtW @ ys)
            fitted[i] = beta[0]
        except np.linalg.LinAlgError:
            fitted[i] = np.average(ys, weights=w)
    return xs, fitted

# ===== Panel A =====
axA.plot([1, len(inst)], [1, len(inst)], ls="--", c="#999999", lw=0.8, zorder=1)
bg = inst[~inst["sig"]]
axA.scatter(bg["rank_main"], bg["rank_gws"], s=8, c="#bbbbbb", alpha=0.5,
            linewidths=0, label="Other phenotypes")
hi = inst[inst["sig"]]
xs, ys = lowess_line(inst["rank_main"], inst["rank_gws"])
axA.plot(xs, ys, color="#0072B2", lw=1.5, alpha=0.9, label="LOESS trend")
axA.scatter(hi["rank_main"], hi["rank_gws"], s=60, c="#E69F00", edgecolors="black",
            linewidths=0.6, zorder=5, label=r"FDR-significant ($P<5\times10^{-8}$)")
axA.set_xlabel(r"Primary-screen $P$-value rank ($P<1\times10^{-5}$)")
axA.set_ylabel(r"Genome-wide-threshold $P$-value rank ($P<5\times10^{-8}$)")
axA.set_title("A  Rank–rank concordance of the two thresholds",
              loc="left", fontsize=11, fontweight="bold")
axA.text(0.03, 0.97, r"Spearman $\rho$ = %.2f" % rho, transform=axA.transAxes,
         va="top", fontsize=10, bbox=dict(boxstyle="round", fc="white", ec="#888"))
axA.legend(loc="lower right", fontsize=8, framealpha=0.9)

# ===== Panel B =====
cats = ["1", "2", "3", r"$\geq$4"]
def _bucket(series):
    return [int((series == 1).sum()), int((series == 2).sum()),
            int((series == 3).sum()), int((series >= 4).sum())]
c_gws, c_main = _bucket(gws["nIV"]), _bucket(main["nIV"])
# The two thresholds have different denominators (607 vs 731 analysable phenotypes),
# so the comparison is plotted as a percentage of phenotypes rather than as raw counts.
p_gws = [100*c/n_gws for c in c_gws]
p_main = [100*c/n_main for c in c_main]
x = np.arange(len(cats)); w = 0.38
bars_g = axB.bar(x - w/2, p_gws, w, color="#E69F00", edgecolor="black", linewidth=0.5,
                 label=r"Genome-wide ($P<5\times10^{-8}$), n=%d" % n_gws)
bars_m = axB.bar(x + w/2, p_main, w, color="#56B4E9", edgecolor="black", linewidth=0.5,
                 label=r"Primary ($P<1\times10^{-5}$), n=%d" % n_main)
for bars, pcts, counts in ((bars_g, p_gws, c_gws), (bars_m, p_main, c_main)):
    for b, pc, c in zip(bars, pcts, counts):
        label = f"{pc:.1f}%" if 0 < pc < 99.95 else f"{pc:.0f}%"
        axB.text(b.get_x() + b.get_width()/2, pc + 1.5,
                 "%s\n(%d)" % (label, c), ha="center", va="bottom", fontsize=7)
axB.set_xticks(x); axB.set_xticklabels(cats)
axB.set_ylim(0, 118)
axB.set_xlabel("Instrumental variants per phenotype")
axB.set_ylabel("Phenotypes (%)")
axB.set_title("B  Instrument sparsity under the two thresholds",
              loc="left", fontsize=11, fontweight="bold")
axB.legend(loc="upper left", fontsize=7.5, framealpha=0.95)
sig_nIV = ", ".join(str(int(v)) for v in sorted(cross["gws_nIV"]))
axB.text(0.97, 0.62, "Genome-wide-threshold FDR signals\nused " + sig_nIV + " variants",
         transform=axB.transAxes, ha="right", va="top", fontsize=7.5,
         bbox=dict(boxstyle="round", fc="#FFF7E6", ec="#E69F00"))

plt.tight_layout()
out = FIG / "Figure_3_ThresholdInstability.png"
plt.savefig(out, dpi=600, bbox_inches="tight")          # 600 dpi, matching the other figures
outpdf = FIG / "Figure_3_ThresholdInstability.pdf"
plt.savefig(outpdf, bbox_inches="tight")                # vector version (sharpest for the text annotations)
print(f"wrote {out}")
print(f"wrote {outpdf}")
print(f"Spearman rho = {rho:.3f}; genome-wide nIV buckets {c_gws} of {n_gws}; "
      f"primary nIV buckets {c_main} of {n_main}; sig nIV = {sig_nIV}")
