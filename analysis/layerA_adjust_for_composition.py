# -*- coding: utf-8 -*-
"""
Sensitivity analysis: does the pre->post immune methylation signature survive
adjustment for the estimated cell-composition shift?

For every probe p and subject j we have paired deltas  D_pj = Post - Pre
(in M-values for testing, in beta for effect size).

  Model 0 (unadjusted):  D_pj = a_p           + e_pj
  Model A (primary):     D_pj = a_p + b_p*dNeu_j + e_pj
  Model B (full lineage):D_pj = a_p + b1*dNeu_j + b2*dLymph_j + b3*dMono_j + e_pj

`a_p` is the composition-ADJUSTED mean change. We compare immune-gene probes vs
all other probes before and after adjustment (competitive gene-set test).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from scipy.stats import mannwhitneyu, wilcoxon

BASE = r"D:/podprj/analysis"
GSE = rf"{BASE}/results/GSE330869"
DEL = rf"{BASE}/delta_ready"
OUT = rf"{BASE}/adjust"

os.makedirs(OUT, exist_ok=True)
lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    lines.append(s)


def bh_fdr(p):
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(q, 0, 1)
    return out


# ------------------------------------------------------------------ inputs
probes = [l.strip() for l in open(rf"{DEL}/probes.txt", encoding="utf-8") if l.strip()]
subjects = [l.strip() for l in open(rf"{DEL}/subjects.txt", encoding="utf-8") if l.strip()]
NP, NS = len(probes), len(subjects)
say(f"probes={NP}  subjects={NS}")

def load_bin(name):
    a = np.fromfile(rf"{DEL}/{name}", dtype=np.float32)
    assert a.size == NP * NS, f"{name}: {a.size} != {NP*NS}"
    return a.reshape((NP, NS), order="F")

DM_float = load_bin("delta_M.bin")
DB = load_bin("delta_beta.bin")
DM = np.ascontiguousarray(DM_float.T, dtype=np.float64)   # NS x NP
Ybeta = np.ascontiguousarray(DB.T, dtype=np.float32)
del DM_float, DB

# ------------------------------------------------- cell-composition deltas
frac = pd.read_csv(rf"{BASE}/layerA_cell_proportions_python.csv")
sel = frac[(frac["config"] == "centCAB100i_a12") & (frac["method"] == "CP_penalty_K1000")].copy()
w = sel.pivot_table(index="subject", columns="celltype", values="fraction")
met = sel[["subject", "time"]].drop_duplicates().set_index("subject")["time"]
LY = ["aNK", "aCD4Tnv", "aCD4Tmem", "aBnv", "aBmem", "aTreg", "aCD8Tmem", "aCD8Tnv"]
w["lymphoid"] = w[LY].sum(axis=1)
pre_sub = met[met == "Pre"].index
dNeu_sub = (w.loc[[s for s in subjects if s in w.index], "aNeu"] -
            w.loc[[s for s in subjects if s in w.index], "aNeu"]).copy()
dcomp = pd.DataFrame(index=subjects, dtype=float)
for col in ["aNeu", "lymphoid", "aMono", "aEos"]:
    pre_v = w.loc[[s for s in subjects if s in w.index], col]
    post_v = pre_v.copy()
    dcomp[col] = np.nan
# recompute properly: need per-subject pre and post values
w2 = sel.pivot_table(index="subject", columns="celltype", values="fraction")
w2["lymphoid"] = w2[LY].sum(axis=1)
fr = sel.set_index(["subject", "time"]).pivot_table(index="subject", columns="time", values="fraction")
LYc = [c for c in LY if c in sel["celltype"].unique()]
lymph_wide = {}
for s in subjects:
    d_post = {}
    d_pre = {}
    for ct in sel["celltype"].unique():
        v_post = sel[(sel["subject"] == s) & (sel["time"] == "Post") & (sel["celltype"] == ct)]["fraction"]
        v_pre = sel[(sel["subject"] == s) & (sel["time"] == "Pre") & (sel["celltype"] == ct)]["fraction"]
        if len(v_post) and len(v_pre):
            d_post[ct] = float(v_post.iloc[0]); d_pre[ct] = float(v_pre.iloc[0])
    if not d_post:
        continue
    lymph_wide[s] = dict(
        dNeu=d_post.get("aNeu", np.nan) - d_pre.get("aNeu", np.nan),
        dLymph=sum(d_post.get(c, 0) for c in LYc) - sum(d_pre.get(c, 0) for c in LYc),
        dMono=d_post.get("aMono", np.nan) - d_pre.get("aMono", np.nan),
        dEos=d_post.get("aEos", np.nan) - d_pre.get("aEos", np.nan),
    )
dcomp = pd.DataFrame(lymph_wide).T.loc[subjects]
assert dcomp.shape[0] == NS
say("\nper-subject composition deltas (mean +- SD):")
for c in dcomp.columns:
    say(f"  {c:<8} mean={dcomp[c].mean()*100:+7.3f}pp  sd={dcomp[c].std(ddof=1)*100:6.3f}pp")

# ------------------------------------------------- probe -> gene -> immune
ann = pd.read_csv(rf"{GSE}/EPICv2_anno_small.csv", dtype=str,
                  usecols=["ID", "chr", "UCSC_RefGene_Name"])
ann["clean"] = ann["ID"].str.replace(r"^([^_]+)_.*$", r"\1", regex=True)
ann_map = {}
for pid, gid in zip(ann["clean"], ann["UCSC_RefGene_Name"]):
    if not isinstance(gid, str) or gid == "":
        continue
    s = ann_map.get(pid)
    if s is None:
        s = ann_map[pid] = set()
    for g in gid.split(";"):
        if g:
            s.add(g)
say(f"\nannotation rows={len(ann)}  probes with >=1 gene={len(ann_map)}")

ig = [l.strip() for l in open(rf"{BASE}/signatures/immune_genes.txt", encoding="utf-8") if l.strip()]
ig = sorted({g.split("#")[0].strip() for g in ig if g.strip() and not g.startswith("#")})
ig = [g for g in ig if g]
imm_gene_set = set(ig)
say(f"immune gene set size={len(imm_gene_set)}")

is_imm = np.zeros(NP, dtype=bool)
imm_gene_of = np.array([""] * NP, dtype=object)
for i, p in enumerate(probes):
    gs = ann_map.get(p)
    if not gs:
        continue
    hit = gs & imm_gene_set
    if hit:
        is_imm[i] = True
        imm_gene_of[i] = ";".join(sorted(hit))
say(f"immune probes (map to >=1 immune gene): {is_imm.sum()} / {NP} "
    f"({100*is_imm.mean():.2f}%)")

# ------------------------------------------------- fitted models
def fit(X, Ymat):
    """OLS across all probes. Returns (intercept t-stats, slopesymean effect (None))."""
    XtXinv = np.linalg.pinv(X.T @ X)
    P = XtXinv @ X.T                    # (k+1) x NS
    B = P @ Ymat                        # (k+1) x NP
    resid = Ymat - X @ B
    df = X.shape[0] - X.shape[1]
    sigma2 = (resid ** 2).sum(axis=0) / df
    var_int = np.maximum(sigma2 * XtXinv[0, 0], 1e-300)
    t = B[0] / np.sqrt(var_int)
    slope_n = B[1] / np.sqrt(np.maximum(sigma2 * XtXinv[1, 1], 1e-300))
    return t, B[0], slope_n, df

ones = np.ones(NS)
results = {}

# Model 0: intercept only
mean0 = DM.mean(axis=0)
sd0 = DM.std(axis=0, ddof=1)
se0 = sd0 / np.sqrt(NS)
t0 = np.divide(mean0, se0, out=np.full(NP, np.nan), where=se0 > 0)
p0 = 2 * tdist.sf(np.abs(np.nan_to_num(t0, nan=0.0)), NS - 1)
p0 = np.where(np.isnan(t0), np.nan, p0)
results["M0_unadjusted"] = dict(t=t0, mean=mean0, p=p0, q=bh_fdr(np.nan_to_num(p0, nan=1.0)))

# Model A: + dNeu
XA = np.column_stack([ones, dcomp["dNeu"].values])
tA, bA, sA, dfA = fit(XA, DM)
pA = 2 * tdist.sf(np.abs(np.nan_to_num(tA, nan=0.0)), dfA)
pA = np.where(np.isnan(tA), np.nan, pA)
results["MA_adj_dNeu"] = dict(t=tA, mean=bA, p=pA, q=bh_fdr(np.nan_to_num(pA, nan=1.0)), slope=sA)

# Model B: + dNeu + dLymph + dMono
XB = np.column_stack([ones, dcomp[["dNeu", "dLymph", "dMono"]].values])
tB, bB, sB, dfB = fit(XB, DM)
pB = 2 * tdist.sf(np.abs(np.nan_to_num(tB, nan=0.0)), dfB)
pB = np.where(np.isnan(tB), np.nan, pB)
results["MB_adj_3lineage"] = dict(t=tB, mean=bB, p=pB, q=bh_fdr(np.nan_to_num(pB, nan=1.0)), slope=sB)

# ------------------------------------------------- report
say("\n\n########## 1. HOW MANY PROBES SURVIVE FDR<0.05 ##########")
say(f"{'model':<20}{'all probes':>14}{'immune probes':>16}{'% immune':>11}{'sig in immune':>16}")
for nm, r in results.items():
    q = r["q"]; p = r["p"]
    sig_all = int(np.nansum((q < 0.05) & ~np.isnan(p)))
    sig_imm = int(np.nansum((q < 0.05) & ~np.isnan(p) & is_imm))
    n_imm = int(is_imm.sum())
    say(f"{nm:<20}{sig_all:>14,}{n_imm:>16,}{100*n_imm/NP:>10.2f}%{sig_imm:>16,}")

say("\n\n########## 2. COMPETITIVE GENE-SET TEST (immune vs background) ##########")
say("  Mann-Whitney on |t|; effect = median |t| immune - median |t| background")
for nm, r in results.items():
    at = np.abs(r["t"])
    ok = ~np.isnan(at)
    a = at[ok & is_imm]; b = at[ok & ~is_imm]
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    say(f"{nm:<20} median|t| immune={np.median(a):7.3f}  background={np.median(b):7.3f}  "
        f"ratio={np.median(a)/np.median(b):5.3f}  MWU p={p:.4g}  (n imm={len(a)})")

say("\n\n########## 3. ATTENUATION OF THE IMMUNE SIGNAL ##########")
base = np.abs(results["M0_unadjusted"]["t"])
for nm in ["MA_adj_dNeu", "MB_adj_3lineage"]:
    at = np.abs(results[nm]["t"])
    rat_im = np.nanmedian(at[is_imm] / np.where(base[is_imm] > 0, base[is_imm], np.nan))
    rat_bg = np.nanmedian(at[~is_imm] / np.where(base[~is_imm] > 0, base[~is_imm], np.nan))
    say(f"{nm:<20} retained fraction of |t|:  immune={rat_im:6.3f}   background={rat_bg:6.3f}")

say("\n\n########## 4. TOP IMMUNE GENES FROM THE DMR ANALYSIS ##########")
cnt = pd.read_csv(rf"{GSE}/layerA_DMR_immune_gene_counts.csv")
top_genes = cnt.sort_values("n_DMR", ascending=False)["gene"].head(15).tolist()
say("genes: " + ", ".join(top_genes))
rows = []
for g in top_genes:
    idx = np.where(is_imm)[0]
    gi = np.array([i for i in idx if g in imm_gene_of[i].split(";")])
    if gi.size == 0:
        continue
    rr = dict(gene=g, n_probes=int(gi.size))
    for nm, key in [("M0", "M0_unadjusted"), ("MA", "MA_adj_dNeu"), ("MB", "MB_adj_3lineage")]:
        r = results[key]
        rr[f"n_q05_{nm}"] = int(np.nansum(r["q"][gi] < 0.05))
        rr[f"med_abst_{nm}"] = float(np.nanmedian(np.abs(r["t"][gi])))
    rr["mean_dBeta"] = float(np.nanmean(np.nanmean(Ybeta[:, gi], axis=0)))
    rr["mean_dBeta_adj"] = float(np.nanmean(results["MA_adj_dNeu"]["mean"][gi]))
    rows.append(rr)
topdf = pd.DataFrame(rows)
say("\n" + topdf.to_string(index=False))

# save per-immune-probe audit table
idx = np.where(is_imm)[0]
audit = pd.DataFrame(dict(
    probe=[probes[i] for i in idx],
    gene=imm_gene_of[idx],
    mean_dBeta=np.nanmean(Ybeta[:, idx], axis=0),
    t_unadj=results["M0_unadjusted"]["t"][idx],
    q_unadj=results["M0_unadjusted"]["q"][idx],
    t_adj_dNeu=results["MA_adj_dNeu"]["t"][idx],
    q_adj_dNeu=results["MA_adj_dNeu"]["q"][idx],
    t_adj_3lin=results["MB_adj_3lineage"]["t"][idx],
    q_adj_3lin=results["MB_adj_3lineage"]["q"][idx],
))
audit.to_csv(rf"{OUT}/adjust_immune_probes.csv", index=False)
say(f"\nwrote {OUT}/adjust_immune_probes.csv ({len(audit)} rows)")

sumrows = []
for nm, r in results.items():
    for grp, mask in [("immune", is_imm), ("background", ~is_imm)]:
        sumrows.append(dict(model=nm, group=grp, n_probes=int(mask.sum()),
                            n_FDR05=int(np.nansum((r["q"] < 0.05) & ~np.isnan(r["p"]) & mask)),
                            median_abs_t=float(np.nanmedian(np.abs(r["t"])[mask]))))
pd.DataFrame(sumrows).to_csv(rf"{OUT}/adjust_summary.csv", index=False)
say(f"wrote {OUT}/adjust_summary.csv")

with open(rf"{OUT}/adjust_report.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("\n=== DONE ===")
