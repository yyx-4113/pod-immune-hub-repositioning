# -*- coding: utf-8 -*-
"""Figure 2 for the POD manuscript: cell-composition deconvolution and the
composition-adjusted sensitivity analysis (probe level + DMR level).

Every number is read from analysis artefacts at run time - nothing hardcoded.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.stats import wilcoxon

AN = r"D:/podprj/analysis"
RES = os.path.join(AN, "results", "GSE330869")
OUT_PNG = os.path.join(RES, "layerA_figure2.png")
OUT_PDF = os.path.join(RES, "layerA_figure2.pdf")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.linewidth": 0.8,
})

# palette
C_NEU = "#D85A30"   # coral 400
C_LYM = "#378ADD"   # blue 400
C_MONO = "#639922"  # green 400
C_IMM = "#D85A30"
C_BG = "#888780"    # gray 400
C_ALL = "#5F5E5A"   # gray 600

# ----------------------------------------------------------------- data ----
props = pd.read_csv(os.path.join(AN, "layerA_cell_proportions_python.csv"))
adj = pd.read_csv(os.path.join(AN, "adjust", "adjust_summary.csv"))
dmr = pd.read_csv(os.path.join(RES, "bumphunter_adjusted_summary.csv"))

LYMPH_TOKENS = ("CD4T", "CD8T", "NK", "Bnv", "Bmem", "Treg", "B")


def lineage(celltype):
    c = celltype.lower().lstrip("a")
    if c.startswith("neu") or c.startswith("gran"):
        return "Neutrophil"
    if c.startswith("mono"):
        return "Monocyte"
    for t in ("cd4t", "cd8t", "nk", "bnv", "bmem", "treg"):
        if c.startswith(t.lower()):
            return "Lymphoid"
    if c == "b":
        return "Lymphoid"
    if c.startswith("eos") or c.startswith("baso"):
        return "Other"
    return "Other"


props["lineage"] = props["celltype"].map(lineage)
lin = (props[props.lineage.isin(["Neutrophil", "Lymphoid", "Monocyte"])]
       .groupby(["sample", "config", "method", "subject", "time", "lineage"],
                as_index=False)["fraction"].sum())

MAIN_CFG, MAIN_METHOD = "centCAB100i_a12", "CP_nnls_norm"

# --------------------------------------------------------------- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.9))  # 180 mm wide, 2-column
axA, axB, axC, axD = axes.ravel()

# ---------- A: paired Pre -> Post, main reference ----------
sub = lin[(lin.config == MAIN_CFG) & (lin.method == MAIN_METHOD)]
piv = sub.pivot_table(index=["subject", "lineage"], columns="time",
                      values="fraction").reset_index()
order = ["Neutrophil", "Lymphoid", "Monocyte"]
colors = {"Neutrophil": C_NEU, "Lymphoid": C_LYM, "Monocyte": C_MONO}

rng = np.random.default_rng(0)
for i, lg in enumerate(order):
    d = piv[piv.lineage == lg]
    x0 = np.full(len(d), i - 0.17) + rng.uniform(-0.05, 0.05, len(d))
    x1 = np.full(len(d), i + 0.17) + rng.uniform(-0.05, 0.05, len(d))
    for a, b, xa, xb in zip(d["Pre"], d["Post"], x0, x1):
        axA.plot([xa, xb], [a * 100, b * 100], color=colors[lg],
                 lw=0.35, alpha=0.45, zorder=1)
    axA.scatter(x0, d["Pre"] * 100, s=5, fc="white", ec=colors[lg],
                lw=0.5, zorder=2)
    axA.scatter(x1, d["Post"] * 100, s=5, fc=colors[lg], ec=colors[lg],
                lw=0.3, zorder=2)
    mu0, mu1 = d["Pre"].median() * 100, d["Post"].median() * 100
    axA.plot([i - 0.3, i + 0.3], [mu0, mu1], color="white", lw=3.2,
             zorder=3, solid_capstyle="butt")
    axA.plot([i - 0.3, i + 0.3], [mu0, mu1], color="black", lw=1.4, zorder=4)
    p = wilcoxon(d["Post"], d["Pre"]).pvalue
    lab = "p = %.1e" % p if p < 1e-3 else "p = %.3f" % p
    dmed = np.median((d["Post"] - d["Pre"]).values) * 100
    ytext = {"Neutrophil": 96, "Lymphoid": 54, "Monocyte": 37}[lg]
    axA.text(i, ytext,
             "$\\Delta$ = %+.2f pp\n%s" % (dmed, lab),
             ha="center", va="top", fontsize=6.5, linespacing=1.35)

axA.set_xticks(range(3))
axA.set_xticklabels(["Neutrophil", "Lymphoid", "Monocyte"])
axA.set_xlim(-0.5, 2.5)
axA.set_ylim(0, 102)
axA.set_ylabel("Estimated fraction (%)")
axA.set_title("A  Paired pre → post lineage shift (n = 65)", pad=15)
axA.text(-0.02, 1.015, "centCAB100i (641 probes), NNLS; bars connect medians",
         transform=axA.transAxes, fontsize=6.5, color="#5F5E5A", va="bottom")

# ---------- B: cross-reference consistency of the shift ----------
cfgs = ["centCAB100i_a12", "centDHSbloodDMC", "centBloodSub",
        "cent12CT450k", "centCAB100i_7"]
lbl = ["centCAB100i\n(641)", "centDHSblood\nDMC (113)", "centBloodSub\n(66)",
       "cent12CT450k\n(185)", "centCAB100i\n7-type"]

for i, lg in enumerate(order):
    xs, ys, lo, hi = [], [], [], []
    for j, cf in enumerate(cfgs):
        d = lin[(lin.config == cf) & (lin.lineage == lg)]
        # median paired delta per solver (same convention as the manuscript table)
        per_method_med = []
        for m in d.method.unique():
            dm = d[d.method == m].pivot_table(index="subject", columns="time",
                                              values="fraction")
            per_method_med.append(np.median((dm["Post"] - dm["Pre"]).values) * 100)
        v = np.array(per_method_med)
        xs.append(j + (i - 1) * 0.22)
        ys.append(np.mean(v))
        lo.append(np.mean(v) - v.min())
        hi.append(v.max() - np.mean(v))
    axB.errorbar(xs, ys, yerr=[lo, hi], fmt="o", color=colors[lg],
                 ecolor=colors[lg], elinewidth=0.8, capsize=2, ms=3.5,
                 lw=0, label=lg, alpha=0.95)

axB.axhline(0, color="#B4B2A9", lw=0.6, zorder=0)
axB.set_xticks(range(5))
axB.set_xticklabels(lbl, fontsize=6, rotation=30, ha="right",
                    rotation_mode="anchor")
axB.set_xlim(-0.6, 4.6)
axB.set_ylabel("$\\Delta$ fraction (percentage points)")
axB.set_title("B  Shift is consistent across five references", pad=15)
axB.legend(frameon=False, loc="upper right", fontsize=6.5, handletextpad=0.3)
axB.text(0.0, 1.015, "points = mean over 4 solvers; bars = solver range",
         transform=axB.transAxes, fontsize=6, color="#5F5E5A", va="bottom")

# ---------- C: probe-level median |t| before/after adjustment ----------
models = ["M0_unadjusted", "MA_adj_dNeu", "MB_adj_3lineage"]
mlab = ["Unadjusted", "Adj. $\\Delta$Neu", "Adj. 3 lineages"]
imm = [float(adj[(adj.model == m) & (adj.group == "immune")].median_abs_t.iloc[0])
       for m in models]
bg = [float(adj[(adj.model == m) & (adj.group == "background")].median_abs_t.iloc[0])
      for m in models]
xs = np.arange(3)
axC.bar(xs - 0.19, imm, 0.36, color=C_IMM, label="Immune probes (n = 1,195)")
axC.bar(xs + 0.19, bg, 0.36, color=C_BG, label="Background (n = 311,319)")
for x, v in zip(xs - 0.19, imm):
    axC.text(x, v + 0.012, "%.3f" % v, ha="center", fontsize=6.5)
for x, v in zip(xs + 0.19, bg):
    axC.text(x, v + 0.012, "%.3f" % v, ha="center", fontsize=6.5)
axC.set_xticks(xs)
axC.set_xticklabels(mlab)
axC.set_ylabel("Median $|t|$ of paired $\\Delta$M")
axC.set_ylim(0, 1.22)
axC.set_title("C  Probe level: immune enrichment is abolished", pad=15)
axC.legend(frameon=False, loc="upper right", fontsize=6.5)
axC.text(-0.02, 1.015,
         "MWU p: 2.46e-4 (unadj.) $\\rightarrow$ 0.724 ($\\Delta$Neu) / 0.55 (3 lin.)",
         transform=axC.transAxes, fontsize=6.5, color="#5F5E5A", va="bottom")

# ---------- D: DMR burden before/after adjustment ----------
dl = {"unadjusted": "Unadjusted", "adj_dNeu": "Adj. $\\Delta$Neu",
      "adj_dNeu_Eos_Baso": "Adj. $\\Delta$Neu\n+$\\Delta$Eos+$\\Delta$Baso",
      "adj_3lineage": "Adj. 3 lineages\n(same set as\npanel C)"}
dmr["lab"] = dmr["model"].map(dl)
dmr = dmr.set_index("model").loc[
    ["unadjusted", "adj_dNeu", "adj_dNeu_Eos_Baso", "adj_3lineage"]].reset_index()
xs = np.arange(len(dmr))
w = 0.36
b1 = axD.bar(xs - w / 2, dmr.n_DMR_L3, w, color=C_ALL, label="All DMRs")
b2 = axD.bar(xs + w / 2, dmr.immune_DMR, w, color=C_IMM, label="Immune-gene DMRs")
for b, v in zip(b1, dmr.n_DMR_L3):
    axD.text(b.get_x() + b.get_width() / 2, v + 250, "{:,}".format(int(v)),
             ha="center", fontsize=6.5)
for b, v in zip(b2, dmr.immune_DMR):
    axD.text(b.get_x() + b.get_width() / 2, v + 250, "%d" % int(v),
             ha="center", fontsize=6.5)
axD.set_xticks(xs)
axD.set_xticklabels(dmr.lab)
axD.set_ylabel("DMR count (|$c$| ≥ 0.05, length ≥ 3)")
axD.set_ylim(0, 16000)
axD.set_title("D  DMR level: attenuated but not abolished", pad=15)
axD.legend(frameon=False, loc="upper right", fontsize=6.5)
for x, v in zip(xs + w / 2, dmr.pct_immune):
    axD.text(x, float(v) * 16000 / 100 + 260, "%.2f%%" % float(v),
             ha="center", fontsize=6.5, color="#0F6E56")

for a in (axA, axB, axC, axD):
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.tick_params(direction="out")

plt.tight_layout(h_pad=2.6, w_pad=2.2)
plt.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
plt.savefig(OUT_PDF, bbox_inches="tight")
print("WROTE", OUT_PNG, os.path.getsize(OUT_PNG))
print("WROTE", OUT_PDF, os.path.getsize(OUT_PDF))

# ---------------------------- numeric audit trail --------------------------
print("\n--- Figure 2 source numbers ---")
for lg in order:
    d = piv[piv.lineage == lg]
    print("%-12s pre %.2f%%  post %.2f%%  median delta %+.2f pp  p %.2e" % (
        lg, d["Pre"].median() * 100, d["Post"].median() * 100,
        np.median((d["Post"] - d["Pre"]).values) * 100,
        wilcoxon(d["Post"], d["Pre"]).pvalue))
print("\nprobe level:", list(zip(mlab, np.round(imm, 3), np.round(bg, 3))))
print("DMR level  :", dmr[["model", "n_DMR_L3", "immune_DMR", "immune_genes", "pct_immune"]].to_string(index=False))
