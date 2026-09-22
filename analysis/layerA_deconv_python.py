# -*- coding: utf-8 -*-
"""
Layer A cell-type deconvolution -- pure-Python re-implementation using EpiDISH
built-in reference centroids (NO external download, NO manifest needed).

Design: reference-based linear mixture model on beta values
        y_p (observed whole-blood beta) = sum_k R_pk * w_k + e_p
        with  w_k >= 0,  sum_k w_k = 1

References exported locally from the EpiDISH package (already installed):
  centBloodSub.m    188x7   (B, NK, CD4T, CD8T, Mono, Neutro, Eosino)
  centDHSbloodDMC.m 333x7   (same 7 types)
  cent12CT450k.m    600x12  (CD4Tnv, Baso, CD4Tmem, Bmem, Bnv, Treg,
                             CD8Tmem, CD8Tnv, Eos, NK, Neu, Mono)
  centCAB100i.m    1906x19  -> split into primary 7 (Gran lumped) and 'a*' 12

NOTE on probe coverage: the GSE330869 beta matrix we hold contains ~312k rows,
i.e. only ~33% of EPIC v2 probes (the authors' filtered deposition). Therefore
every reference retains only ~30-36% of its probes after intersection. This is
a row-coverage limitation of the DATA, not an address-space mismatch.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import wilcoxon, pearsonr
from sklearn.svm import NuSVR

OUT = r"D:/podprj/analysis/deconv_ref_ready"
DST = r"D:/podprj/analysis"

CONFIGS = [
    ("centBloodSub",    "ref_centBloodSub.csv",    "beta_centBloodSub.csv",    "none"),
    ("centDHSbloodDMC", "ref_centDHSbloodDMC.csv", "beta_centDHSbloodDMC.csv", "none"),
    ("cent12CT450k",    "ref_cent12CT450k.csv",    "beta_cent12CT450k.csv",    "none"),
    ("cent12CT",        "ref_cent12CT.csv",        "beta_cent12CT.csv",        "none"),
    ("centCAB100i_a12", "ref_centCAB100i.csv",     "beta_centCAB100i.csv",     "a_only"),
    ("centCAB100i_7",   "ref_centCAB100i.csv",     "beta_centCAB100i.csv",     "a_drop"),
]


def load_pair(refname, betaname, mode):
    ref = pd.read_csv(os.path.join(OUT, refname))
    bet = pd.read_csv(os.path.join(OUT, betaname))
    ref = ref.rename(columns={ref.columns[0]: "probe"})
    bet = bet.rename(columns={bet.columns[0]: "probe"})
    common = [p for p in ref["probe"] if p in set(bet["probe"])]
    ref = ref[ref["probe"].isin(common)].set_index("probe").loc[common]
    bet = bet[bet["probe"].isin(common)].set_index("probe").loc[common]
    cols = list(ref.columns)
    if mode == "a_only":
        cols = [c for c in cols if c.startswith("a")]
    elif mode == "a_drop":
        cols = [c for c in cols if not c.startswith("a")]
    ref = ref[cols]
    return ref, bet


# ---------------------------------------------------------------- estimators
def cp_nnls_norm(R, y):
    """NNLS then normalise to sum 1."""
    w, _ = nnls(R, y)
    s = w.sum()
    return w / s if s > 0 else np.full(R.shape[1], np.nan)


def cp_penalty(R, y, kappa=1000.0):
    """Constrained projection: NNLS with explicit sum-to-one penalty row."""
    K = R.shape[1]
    lam = kappa * np.sqrt(np.mean(R ** 2))
    Ra = np.vstack([R, lam * np.ones((1, K))])
    ya = np.concatenate([y, [lam]])
    w, _ = nnls(Ra, ya)
    return w


def cbs_svr(R, y, nu=0.5, C=1.0, drop_q=None):
    """CIBERSORT-style linear nu-SVR; clip negatives, normalise."""
    try:
        svr = NuSVR(kernel="linear", C=C, nu=nu, shrinking=False, max_iter=200000)
        svr.fit(R, y)
        coef = np.asarray(svr.coef_).ravel()
    except Exception:
        return np.full(R.shape[1], np.nan), False
    if coef.size != R.shape[1]:
        return np.full(R.shape[1], np.nan), False
    if drop_q is not None and drop_q > 0:
        thr = np.quantile(np.abs(coef), drop_q)
        coef[np.abs(coef) <= thr] = 0.0
    w = np.clip(coef, 0.0, None)
    s = w.sum()
    return (w / s if s > 0 else np.full(R.shape[1], np.nan)), True


def rpc_irls(R, y, n_iter=30, c=4.685, kappa=1000.0):
    """RPC-style robust regression: Tukey-biweight IRLS with sum-to-one."""
    K = R.shape[1]
    lam = kappa * np.sqrt(np.mean(R ** 2))
    w = cp_penalty(R, y, kappa=kappa)
    for _ in range(n_iter):
        resid = y - R @ w
        mad = np.median(np.abs(resid - np.median(resid)))
        s = 1.4826 * mad
        if s <= 0 or np.isnan(s):
            s = 1e-4
        s = max(s, 1e-4)
        u = resid / (c * s)
        wt = np.where(np.abs(u) < 1.0, (1.0 - u ** 2) ** 2, 0.0)
        sw = np.sqrt(wt)
        Ra = np.vstack([R * sw[:, None], lam * np.ones((1, K))])
        ya = np.concatenate([y * sw, [lam]])
        wn, _ = nnls(Ra, ya)
        if np.linalg.norm(wn - w) < 1e-7:
            w = wn
            break
        w = wn
    return w


METHODS = {
    "CP_nnls_norm": cp_nnls_norm,
    "CP_penalty_K1000": lambda R, y: cp_penalty(R, y, 1000.0),
    "CP_penalty_K100": lambda R, y: cp_penalty(R, y, 100.0),
    "CBS_svr_nu0.5": lambda R, y: cbs_svr(R, y, 0.5, 1.0)[0],
    "CBS_svr_nu0.25": lambda R, y: cbs_svr(R, y, 0.25, 1.0)[0],
    "RPC_irls": rpc_irls,
}


def fit_config(ref, bet, method_fn):
    R = ref.values.astype(float)
    Y = bet.values.astype(float)
    cols = list(ref.columns)
    samples = list(bet.columns)
    rows = []
    qual = []
    for j, smp in enumerate(samples):
        y = Y[:, j]
        w = method_fn(R, y)
        fit = R @ w
        rmse = float(np.sqrt(np.mean((y - fit) ** 2)))
        rmse0 = float(np.sqrt(np.mean((y - y.mean()) ** 2)))  # intercept-only baseline
        try:
            corr = float(pearsonr(y, fit)[0])
        except Exception:
            corr = np.nan
        qual.append(dict(sample=smp, rmse=rmse, rmse_baseline=rmse0,
                         r2_like=1 - rmse ** 2 / max(rmse0 ** 2, 1e-12),
                         corr=corr, sum_w=float(np.sum(w)),
                         min_w=float(np.min(w)), max_w=float(np.max(w))))
        for k, ct in enumerate(cols):
            rows.append(dict(sample=smp, celltype=ct, fraction=float(w[k])))
    return pd.DataFrame(rows), pd.DataFrame(qual)


def main():
    meta = pd.read_csv(os.path.join(OUT, "sample_meta.csv"))
    meta = meta.rename(columns={"Unnamed: 0": "drop"}) if "Unnamed: 0" in meta.columns else meta
    mmap = meta.set_index("sample")[["subject", "time"]].to_dict("index")

    all_long, all_qual, rep = [], [], []

    for cname, rf, bf, mode in CONFIGS:
        try:
            ref, bet = load_pair(rf, bf, mode)
        except Exception as e:
            rep.append(f"[{cname}] LOAD FAIL: {e}")
            continue
        rep.append(f"\n===== {cname} :: {ref.shape[0]} probes x {ref.shape[1]} types =====")
        rep.append("types: " + ", ".join(ref.columns))
        for mname, mfn in METHODS.items():
            try:
               flong, fq = fit_config(ref, bet, mfn)
            except Exception as e:
                rep.append(f"  {mname}: FIT FAIL {e}")
                continue
            flong["config"] = cname
            flong["method"] = mname
            fq["config"] = cname
            fq["method"] = mname
            all_long.append(flong)
            all_qual.append(fq)
            rep.append(f"  {mname:<18} mean sum(w)={fq['sum_w'].mean():.4f} "
                       f"rmse={fq['rmse'].mean():.4f} (baseline {fq['rmse_baseline'].mean():.4f}) "
                       f"r2={fq['r2_like'].mean():.3f} corr={fq['corr'].mean():.3f} "
                       f"min_w={fq['min_w'].mean():.4f}")

    long = pd.concat(all_long, ignore_index=True)
    qual = pd.concat(all_qual, ignore_index=True)
    long["subject"] = long["sample"].map(lambda s: mmap.get(s, {}).get("subject"))
    long["time"] = long["sample"].map(lambda s: mmap.get(s, {}).get("time"))
    qual["subject"] = qual["sample"].map(lambda s: mmap.get(s, {}).get("subject"))
    qual["time"] = qual["sample"].map(lambda s: mmap.get(s, {}).get("time"))

    long.to_csv(os.path.join(DST, "layerA_cell_proportions_python.csv"), index=False)
    qual.to_csv(os.path.join(DST, "layerA_deconv_quality.csv"), index=False)

    # ---------------- per cell type mean proportions and paired tests -------
    rep.append("\n\n===== ABSOLUTE FRACTIONS (mean over 65 Pre / 65 Post) =====")
    summ = []
    for (cname, mname), g in long.groupby(["config", "method"]):
        for ct, gc in g.groupby("celltype"):
            pre = gc[gc["time"] == "Pre"].set_index("subject")["fraction"]
            post = gc[gc["time"] == "Post"].set_index("subject")["fraction"]
            common = [s for s in pre.index if s in post.index]
            pre = pre.loc[common].values
            post = post.loc[common].values
            try:
                W, p = wilcoxon(post, pre, zero_method="wilcox")
            except Exception:
                W, p = np.nan, np.nan
            summ.append(dict(config=cname, method=mname, celltype=ct, n=len(common),
                             mean_pre=float(np.mean(pre)), mean_post=float(np.mean(post)),
                             sd_pre=float(np.std(pre, ddof=1)), sd_post=float(np.std(post, ddof=1)),
                             median_delta=float(np.median(post - pre)),
                             wilcoxon_p=float(p)))
    summdf = pd.DataFrame(summ)
    summdf.to_csv(os.path.join(DST, "layerA_deconv_summary.csv"), index=False)

    for (cname, mname), g in summdf.groupby(["config", "method"]):
        rep.append(f"\n--- {cname} / {mname} ---")
        gg = g.sort_values("mean_pre", ascending=False)
        for _, r in gg.iterrows():
            rep.append(f"  {r['celltype']:<10} Pre={r['mean_pre']*100:6.2f}%  Post={r['mean_post']*100:6.2f}%  "
                       f"delta={r['median_delta']*100:+6.2f}pp  p={r['wilcoxon_p']:.4g}")

    txt = "\n".join(rep)
    with open(os.path.join(DST, "layerA_deconv_report.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
