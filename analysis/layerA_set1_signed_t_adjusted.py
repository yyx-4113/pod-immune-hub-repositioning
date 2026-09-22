# -*- coding: utf-8 -*-
"""
T1-07 补分析: Set-1 signed-t (directional probe test) 在组成校正后的 MWU p.

主脚本 layerA_adjust_for_composition.py 第 195-200 行只对 |t| 做 MWU (Set-2)。
稿件 Results/Discussion 报 Set-1 signed-t 未校正 p=0.041, 但校正后未报。
本补丁补算 Set-1 signed-t (不取 abs) 在 Model 0 / A(dNeu) / B(3-lineage) 的 MWU p,
并写入 set1_signed_t_adjusted.csv。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from scipy.stats import mannwhitneyu

BASE = r"D:\2026.9\极速交付9月会员日优惠套路\05_多组学+虚拟敲除药物发现\术后谵妄免疫枢纽\pod-epigenetic-immune-repositioning\analysis"
DEL = os.path.join(BASE, "delta_ready")
GSE = os.path.join(BASE, "results", "GSE330869")
OUT = os.path.join(BASE, "adjust")
os.makedirs(OUT, exist_ok=True)

# -------------------------------------------------- inputs
probes = [l.strip() for l in open(os.path.join(DEL, "probes.txt"), encoding="utf-8") if l.strip()]
subjects = [l.strip() for l in open(os.path.join(DEL, "subjects.txt"), encoding="utf-8") if l.strip()]
NP, NS = len(probes), len(subjects)
print(f"probes={NP}  subjects={NS}")


def load_bin(name):
    a = np.fromfile(os.path.join(DEL, name), dtype=np.float32)
    assert a.size == NP * NS, f"{name}: {a.size} != {NP*NS}"
    return a.reshape((NP, NS), order="F")


DM_float = load_bin("delta_M.bin")
DM = np.ascontiguousarray(DM_float.T, dtype=np.float64)   # NS x NP
del DM_float
print(f"DM shape (NS x NP) = {DM.shape}")

# -------------------------------------------------- composition deltas
frac = pd.read_csv(os.path.join(BASE, "layerA_cell_proportions_python.csv"))
sel = frac[(frac["config"] == "centCAB100i_a12") & (frac["method"] == "CP_penalty_K1000")].copy()
LY = ["aNK", "aCD4Tnv", "aCD4Tmem", "aBnv", "aBmem", "aTreg", "aCD8Tmem", "aCD8Tnv"]
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
    )
dcomp = pd.DataFrame(lymph_wide).T.loc[subjects]
assert dcomp.shape[0] == NS
print("per-subject composition deltas (mean +- SD):")
for c in dcomp.columns:
    print(f"  {c:<8} mean={dcomp[c].mean()*100:+7.3f}pp  sd={dcomp[c].std(ddof=1)*100:6.3f}pp")

# -------------------------------------------------- probe -> gene -> immune
ann = pd.read_csv(os.path.join(GSE, "EPICv2_anno_small.csv"), dtype=str,
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

ig = [l.strip() for l in open(os.path.join(BASE, "signatures", "immune_genes.txt"), encoding="utf-8") if l.strip()]
ig = sorted({g.split("#")[0].strip() for g in ig if g.strip() and not g.startswith("#")})
ig = [g for g in ig if g]
imm_gene_set = set(ig)
print(f"immune gene set size={len(imm_gene_set)}")

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
print(f"immune probes (map to >=1 immune gene): {is_imm.sum()} / {NP} ({100*is_imm.mean():.2f}%)")

# -------------------------------------------------- fitted models
def fit(X, Ymat):
    XtXinv = np.linalg.pinv(X.T @ X)
    P = XtXinv @ X.T
    B = P @ Ymat
    resid = Ymat - X @ B
    df = X.shape[0] - X.shape[1]
    sigma2 = (resid ** 2).sum(axis=0) / df
    var_int = np.maximum(sigma2 * XtXinv[0, 0], 1e-300)
    t = B[0] / np.sqrt(var_int)
    return t, B[0], df


ones = np.ones(NS)

# Model 0: intercept only
mean0 = DM.mean(axis=0)
sd0 = DM.std(axis=0, ddof=1)
se0 = sd0 / np.sqrt(NS)
t0 = np.divide(mean0, se0, out=np.full(NP, np.nan), where=se0 > 0)

# Model A: + dNeu
XA = np.column_stack([ones, dcomp["dNeu"].values])
tA, bA, dfA = fit(XA, DM)

# Model B: + dNeu + dLymph + dMono
XB = np.column_stack([ones, dcomp[["dNeu", "dLymph", "dMono"]].values])
tB, bB, dfB = fit(XB, DM)

models = {
    "M0_unadjusted": t0,
    "MA_adj_dNeu": tA,
    "MB_adj_3lineage": tB,
}

# -------------------------------------------------- Set-1 signed-t MWU
print("\n########## Set-1 SIGNED-T MWU (immune vs background, signed t, two-sided) ##########")
rows = []
for nm, tvec in models.items():
    ok = ~np.isnan(tvec)
    a = tvec[ok & is_imm]
    b = tvec[ok & ~is_imm]
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    med_imm = float(np.median(a))
    med_bg = float(np.median(b))
    mean_imm = float(np.mean(a))
    mean_bg = float(np.mean(b))
    print(f"{nm:<20} n_imm={len(a):5}  n_bg={len(b):6}  "
          f"median_t imm={med_imm:+7.4f}  bg={med_bg:+7.4f}  "
          f"mean_t imm={mean_imm:+7.4f}  bg={mean_bg:+7.4f}  "
          f"MWU p={p:.4g}")
    rows.append(dict(
        model=nm,
        n_immune=len(a),
        n_background=len(b),
        median_t_immune=med_imm,
        median_t_background=med_bg,
        mean_t_immune=mean_imm,
        mean_t_background=mean_bg,
        mwu_p_signed=p,
    ))

df_out = pd.DataFrame(rows)
out_path = os.path.join(OUT, "set1_signed_t_adjusted.csv")
df_out.to_csv(out_path, index=False)
print(f"\nwrote {out_path}")
print("\n=== DONE ===")
