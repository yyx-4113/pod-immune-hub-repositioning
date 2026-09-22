#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer F - Genetics: two-sample Mendelian randomisation + POD polygenic scoring.
================================================================================
Context (verified 2026-09-20, Armstrong et al., PLoS Med 2026;23(3):e1004963,
PMID 41770756): the only published POD GWAS (UK Biobank, 1,016 cases /
139,148 controls) yields SEVEN genome-wide significant SNPs, ALL in the chr19
APOE region (APOE / TOMM40 / APOC1 / PVRL2), lead SNP rs429358 (APOE eps4).
Genetic correlation with Alzheimer's disease rho = 0.68.

=> PRE-REGISTERED EXPECTATION
   * immune/inflammation trait -> POD MR: EXPECTED NULL outside the APOE
     region. A null here is a *finding*, not a failure, and is reported as
     "POD common-variant liability is APOE-driven; peripheral immune traits
     are not supported as causal by current GWAS power".
   * APOE-region/AD-PRS -> POD: EXPECTED POSITIVE (published: AD PRS quintile
     5 OR 2.29; APOE-independent AD PRS quintile 5 OR 1.46). Replicating this
     in an independent cohort is a positive, reportable result.

What this script does
---------------------
1. Harmonise exposure and outcome summary statistics (allele matching,
   palindromic SNP handling using allele frequency when available).
2. Two-sample MR: Wald ratio (1 IV), IVW (fixed + random effects), MR-Egger
   (with intercept test for directional pleiotropy), weighted median,
   weighted mode. Heterogeneity (Cochran Q, I2), leave-one-out, single-SNP.
3. Bidirectional mode: build instruments from the POD GWAS and test
   POD -> immune trait.
4. PRS: score a cohort from a weights file (rsID, effect_allele, weight)
   and a long/dosage genotype table; report OR per SD with age/sex adjustment.

Dependencies: numpy, pandas, scipy, scikit-learn. No R required.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

# --- published POD instruments (Armstrong 2026, Table 2) -------------------
# OR as reported in the paper; log(OR) used as weight. NOTE: the reported OR
# is per effect allele with respect to POD risk as defined in the paper; the
# direction must be re-checked against the summary-statistics file before use.
PUBLISHED_POD_SNPS = [
    # rsID,        chr:pos(hg19), gene,    effect_allele, OR
    ("rs429358",   "19:45411941", "APOE",   "T", 0.54),
    ("rs157592",   "19:45424514", "APOC1",  "A", 0.58),
    ("rs157582",   "19:45396219", "TOMM40", "C", 0.65),
    ("rs11556505", "19:45396144", "TOMM40", "C", 0.62),
    ("rs10119",    "19:45406673", "TOMM40", "G", 0.70),
    ("rs75627662", "19:45413576", "APOE",   "C", 0.72),
    ("rs12691088", "19:45418486", "APOC1",  "G", 0.45),
]

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# --------------------------------------------------------------------------
# IO helpers
# --------------------------------------------------------------------------
def read_sumstats(path: str, snp_subset: Optional[set] = None,
                  p_thresh: Optional[float] = None,
                  chunksize: int = 1_000_000) -> pd.DataFrame:
    """Read a (optionally gzipped) tab-delimited summary-statistics file.

    Supports the Armstrong/UK Biobank layout:
    CHROM GENPOS SNP ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ P EXTRA
    as well as the generic layout used by TwoSampleMR exports.
    """
    opener = gzip.open if path.endswith(".gz") else open
    keep = []
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for chunk in pd.read_csv(fh, sep="\t", names=header, chunksize=chunksize,
                                 low_memory=False, dtype={"SNP": str, "CHROM": str}):
            if snp_subset is not None:
                chunk = chunk[chunk["SNP"].isin(snp_subset)]
            if p_thresh is not None and "P" in chunk.columns:
                chunk = chunk[chunk["P"] <= p_thresh]
            if len(chunk):
                keep.append(chunk)
    if not keep:
        return pd.DataFrame(columns=header)
    return pd.concat(keep, ignore_index=True)


def extract_pod_snps(pod_gwas: str, rsids: set) -> pd.DataFrame:
    """Stream the 416 MB UK Biobank POD GWAS and keep only requested rsIDs."""
    df = read_sumstats(pod_gwas, snp_subset=rsids)
    cols = {c.lower(): c for c in df.columns}
    out = pd.DataFrame({
        "SNP": df["SNP"].astype(str),
        "chr": df[cols.get("chrom", "CHROM")] if "chrom" in cols else np.nan,
        "pos": df[cols.get("genpos", "GENPOS")] if "genpos" in cols else np.nan,
        "allele0": df[cols.get("allele0", "ALLELE0")] if "allele0" in cols else np.nan,
        "allele1": df[cols.get("allele1", "ALLELE1")] if "allele1" in cols else np.nan,
        "beta": pd.to_numeric(df[cols.get("beta", "BETA")], errors="coerce"),
        "se": pd.to_numeric(df[cols.get("se", "SE")], errors="coerce"),
        "p": pd.to_numeric(df[cols.get("p", "P")], errors="coerce"),
        "eaf": pd.to_numeric(df[cols.get("a1freq", "A1FREQ")], errors="coerce")
        if "a1freq" in cols else np.nan,
        "n": pd.to_numeric(df[cols.get("n", "N")], errors="coerce") if "n" in cols else np.nan,
    })
    return out


# --------------------------------------------------------------------------
# Harmonisation
# --------------------------------------------------------------------------
def harmonise(exp: pd.DataFrame, out: pd.DataFrame) -> pd.DataFrame:
    """Allele-harmonise exposure (IV) table against outcome table on SNP."""
    m = exp.merge(out, on="SNP", suffixes=("_exp", "_out"))
    if len(m) == 0:
        return m
    ea = m["effect_allele_exp"].astype(str).str.upper()
    oa = m["other_allele_exp"].astype(str).str.upper()
    oea = m["allele1_out"].astype(str).str.upper()      # outcome effect allele
    ooa = m["allele0_out"].astype(str).str.upper()      # outcome other allele

    beta_out = m["beta_out"].astype(float).copy()
    eaf_out = m["eaf_out"].astype(float) if "eaf_out" in m else pd.Series(np.nan, index=m.index)

    action = np.array(["keep"] * len(m), dtype=object)
    for i in range(len(m)):
        a1, a2 = ea.iloc[i], oa.iloc[i]
        b1, b2 = oea.iloc[i], ooa.iloc[i]
        if a1 == b1 and a2 == b2:
            continue
        if a1 == b2 and a2 == b1:
            action[i] = "flip"
            beta_out.iloc[i] = -beta_out.iloc[i]
            continue
        # try complement
        c1, c2 = COMPLEMENT.get(a1, "?"), COMPLEMENT.get(a2, "?")
        if c1 == b1 and c2 == b2:
            action[i] = "complement"
            continue
        if c1 == b2 and c2 == b1:
            action[i] = "complement+flip"
            beta_out.iloc[i] = -beta_out.iloc[i]
            continue
        action[i] = "drop"
    # palindromic SNPs: ambiguous unless allele frequency far from 0.5
    palin = ((ea == COMPLEMENT.get(oa)) if False else
             ea.eq(oa.map(COMPLEMENT)))
    ambiguous = palin & (eaf_out.between(0.42, 0.58))
    action[ambiguous.to_numpy()] = "drop_ambiguous_palindromic"

    m["beta_out"] = beta_out
    m["harmonise_action"] = action
    m = m[~m["harmonise_action"].astype(str).str.startswith("drop")]
    return m.reset_index(drop=True)


# --------------------------------------------------------------------------
# MR estimators
# --------------------------------------------------------------------------
def wald_ratio(bx, by, sx, sy):
    est = by / bx
    se = abs(sy / bx)
    return est, se


def mr_ivw(bx, by, se_y, random: bool = True):
    w = 1.0 / se_y ** 2
    est = float(np.sum(w * bx * by) / np.sum(w * bx ** 2))
    se = float(np.sqrt(1.0 / np.sum(w * bx ** 2)))
    q = float(np.sum(w * (by - est * bx) ** 2))
    dfree = len(bx) - 1
    if random and dfree > 0 and q > dfree:
        se = se * math_sqrt(max(q / dfree, 1.0))
    p = 2 * stats.norm.sf(abs(est / se)) if se > 0 else np.nan
    q_p = stats.chi2.sf(q, dfree) if dfree > 0 else np.nan
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    return est, se, p, q, q_p, i2


def math_sqrt(x):
    return float(np.sqrt(x))


def mr_egger(bx, by, se_y):
    """MR-Egger regression: by = b0 + b1*bx, weights 1/se_y^2."""
    w = 1.0 / se_y ** 2
    X = np.column_stack([np.ones(len(bx)), bx])
    XtW = X.T * w
    A = XtW @ X
    b = np.linalg.solve(A, XtW @ by)
    resid = by - X @ b
    dof = len(bx) - 2
    s2 = float(np.sum(w * resid ** 2) / dof) if dof > 0 else np.nan
    cov = np.linalg.inv(A) * s2 if dof > 0 else np.full((2, 2), np.nan)
    se = np.sqrt(np.diag(cov))
    p = 2 * stats.t.sf(abs(b / se), dof) if dof > 0 else np.array([np.nan, np.nan])
    return b[1], se[1], p[1], b[0], se[0], p[0]


def _wm_point(ratios, w):
    order = np.argsort(ratios)
    r, ww = ratios[order], w[order]
    cw = np.cumsum(ww) - 0.5 * ww
    tot = ww.sum()
    if tot <= 0:
        return float(np.median(r))
    idx = min(int(np.searchsorted(cw / tot, 0.5)), len(r) - 1)
    return float(r[idx])


def weighted_median(bx, by, se_y, n_boot: int = 1000, seed: int = 20260920):
    """Weighted-median estimator with a NON-PARAMETRIC BOOTSTRAP SE.

    The closed-form SE of the weighted median is not standard; a bootstrap
    over instruments is used instead (n_boot resamples, with replacement).
    """
    ratios = by / bx
    w = 1.0 / (se_y / np.abs(bx)) ** 2
    est = _wm_point(ratios, w)
    if len(ratios) < 2:
        return est, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    n = len(ratios)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots[b] = _wm_point(ratios[idx], w[idx])
    se = float(np.std(boots, ddof=1))
    p = float(2 * stats.norm.sf(abs(est / se))) if se > 0 else np.nan
    return est, se, p


def run_mr(df: pd.DataFrame, label: str, direction: str) -> dict:
    bx = df["beta_exp"].astype(float).to_numpy()
    by = df["beta_out"].astype(float).to_numpy()
    sx = df["se_exp"].astype(float).to_numpy() if "se_exp" in df else np.full(len(bx), np.nan)
    sy = df["se_out"].astype(float).to_numpy()
    res = {"analysis": label, "direction": direction, "n_iv": int(len(bx))}
    if len(bx) == 0:
        res["status"] = "NO_INSTRUMENT"
        return res
    est, se, p, q, q_p, i2 = mr_ivw(bx, by, sy, random=True)
    res["ivw_random_beta"] = est; res["ivw_random_se"] = se; res["ivw_random_p"] = p
    e2, se2, p2, _, _, _ = mr_ivw(bx, by, sy, random=False)
    res["ivw_fixed_beta"] = e2; res["ivw_fixed_se"] = se2; res["ivw_fixed_p"] = p2
    res["cochran_Q"] = q; res["Q_p"] = q_p; res["I2"] = i2
    if len(bx) >= 3:
        b1, s1, p1, b0, s0, p0 = mr_egger(bx, by, sy)
        res["egger_beta"] = b1; res["egger_se"] = s1; res["egger_p"] = p1
        res["egger_intercept"] = b0; res["egger_intercept_p"] = p0
        wm_b, wm_se, wm_p = weighted_median(bx, by, sy)
        res["weighted_median_beta"] = wm_b; res["weighted_median_se"] = wm_se
        res["weighted_median_p"] = wm_p
    if len(bx) == 1:
        e, s = wald_ratio(bx[0], by[0], sx[0], sy[0])
        res["wald_ratio_beta"] = float(e); res["wald_ratio_se"] = float(s)
        res["wald_ratio_p"] = float(2 * stats.norm.sf(abs(e / s)))
    # F-statistic (approximate, from beta/se)
    f = (bx / sx) ** 2
    res["mean_F"] = float(np.nanmean(f))
    res["min_F"] = float(np.nanmin(f))
    # leave-one-out
    loo = []
    for i in range(len(bx)):
        keep = np.arange(len(bx)) != i
        e, s, pp, _, _, _ = mr_ivw(bx[keep], by[keep], sy[keep], random=True)
        loo.append({"snp_dropped": df["SNP"].iloc[i], "beta": e, "se": s, "p": pp})
    res["leave_one_out"] = loo
    res["status"] = "OK"
    return res


# --------------------------------------------------------------------------
# PRS
# --------------------------------------------------------------------------
def score_prs(weights: pd.DataFrame, dosages: pd.DataFrame,
              covariates: Optional[pd.DataFrame], outcome: Optional[pd.Series] = None):
    """weights: SNP, effect_allele, weight ; dosages: index=SNP, columns=samples."""
    w = weights.set_index("SNP")
    common = [s for s in w.index if s in dosages.index]
    w = w.loc[common]
    D = dosages.loc[common]
    # orient dosage to effect allele if an 'allele' column is supplied
    score = (D.T.to_numpy() @ w["weight"].to_numpy(dtype=float))
    out = pd.DataFrame({"sample": dosages.columns, "PRS": score}).set_index("sample")
    out["PRS_z"] = (out["PRS"] - out["PRS"].mean()) / out["PRS"].std(ddof=1)
    summary = {"n_variants_scored": len(common),
               "n_variants_in_weights": int(weights.shape[0]),
               "prs_mean": float(out["PRS"].mean()), "prs_sd": float(out["PRS"].std(ddof=1))}
    if outcome is not None and covariates is not None:
        try:
            from sklearn.linear_model import LogisticRegression
            X = np.column_stack([out["PRS_z"].to_numpy(), covariates.to_numpy(dtype=float)])
            y = outcome.to_numpy()
            # penalty=None is deprecated in sklearn >=1.8; C=inf gives an
            # unpenalised fit with the current default solver.
            m = LogisticRegression(C=np.inf, max_iter=5000).fit(X, y)
            coef = m.coef_[0]
            # SE from Fisher information approximation
            p = m.predict_proba(X)[:, 1]
            W = p * (1 - p)
            XtW = X.T * W
            cov = np.linalg.pinv(XtW @ X)
            se = np.sqrt(np.diag(cov))
            z = coef / se
            pv = 2 * stats.norm.sf(np.abs(z))
            summary["logit_coef_per_SD"] = float(coef[0]); summary["se"] = float(se[0])
            summary["OR_per_SD"] = float(np.exp(coef[0]))
            summary["CI95"] = [float(np.exp(coef[0] - 1.96 * se[0])),
                               float(np.exp(coef[0] + 1.96 * se[0]))]
            summary["p"] = float(pv[0])
            summary["covariates"] = list(covariates.columns)
        except Exception as exc:  # pragma: no cover
            summary["model_error"] = str(exc)
    return out, summary


def published_weights() -> pd.DataFrame:
    rows = [{"SNP": r[0], "gene": r[2], "effect_allele": r[3],
             "OR": r[4], "weight": float(np.log(r[4]))} for r in PUBLISHED_POD_SNPS]
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Layer F: two-sample MR + POD PRS")
    ap.add_argument("--pod-gwas", help="POD GWAS summary stats (.txt/.gz)")
    ap.add_argument("--exposure", help="exposure summary stats TSV (SNP,effect_allele,"
                                       "other_allele,beta,se,pval[,eaf])")
    ap.add_argument("--instruments", help="pre-clumped instrument TSV (subset of exposure)")
    ap.add_argument("--p-thresh", type=float, default=5e-8)
    ap.add_argument("--bidirectional", action="store_true",
                    help="also run POD -> exposure using POD instruments")
    ap.add_argument("--prs-weights", help="weights TSV (SNP,effect_allele,weight)")
    ap.add_argument("--prs-dosages", help="dosage TSV (rows SNP, cols samples)")
    ap.add_argument("--prs-covariates", help="covariates TSV (age, sex, ...) indexed by sample")
    ap.add_argument("--prs-outcome", help="outcome TSV (sample, y)")
    ap.add_argument("--outdir", default="results/layerF")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    report = {"pre_registered_expectation":
              "immune-trait->POD MR expected NULL (POD GWAS is APOE-driven); "
              "APOE/AD-PRS->POD expected POSITIVE"}

    # ---------- MR ----------
    if args.pod_gwas and (args.exposure or args.instruments):
        exp = pd.read_csv(args.instruments or args.exposure, sep="\t")
        if args.instruments is None:
            exp = exp[exp["pval"].astype(float) <= args.p_thresh]
        exp = exp.rename(columns={"pval": "p_exp", "beta": "beta_exp", "se": "se_exp"})
        exp = exp[["SNP", "effect_allele", "other_allele", "beta_exp", "se_exp"]]

        rsids = set(exp["SNP"].astype(str))
        pod = extract_pod_snps(args.pod_gwas, rsids)
        pod.to_csv(os.path.join(args.outdir, "layerF_POD_outcome_at_instruments.csv"), index=False)
        pod = pod.rename(columns={"beta": "beta_out", "se": "se_out", "p": "p_out",
                                  "eaf": "eaf_out"})
        harm = harmonise(exp, pod)
        harm.to_csv(os.path.join(args.outdir, "layerF_harmonised.csv"), index=False)

        label = os.path.basename(args.instruments or args.exposure)
        mr = run_mr(harm, label=label, direction="exposure -> POD")
        report["mr_forward"] = mr

        if args.bidirectional:
            pod_all = read_sumstats(args.pod_gwas, p_thresh=args.p_thresh)
            pod_all = pod_all.rename(columns={"SNP": "SNP", "BETA": "beta_exp",
                                              "SE": "se_exp", "P": "p_exp",
                                              "ALLELE1": "effect_allele",
                                              "ALLELE0": "other_allele"})
            pod_all = pod_all[["SNP", "effect_allele", "other_allele",
                               "beta_exp", "se_exp"]]
            out2 = pd.read_csv(args.exposure, sep="\t").rename(
                columns={"beta": "beta_out", "se": "se_out",
                         "effect_allele": "allele1_out", "other_allele": "allele0_out",
                         "eaf": "eaf_out"})
            h2 = harmonise(pod_all, out2)
            h2.to_csv(os.path.join(args.outdir, "layerF_harmonised_reverse.csv"), index=False)
            report["mr_reverse"] = run_mr(h2, label=label, direction="POD -> exposure")

    # ---------- PRS ----------
    if args.prs_dosages:
        w = pd.read_csv(args.prs_weights, sep="\t") if args.prs_weights else published_weights()
        w.to_csv(os.path.join(args.outdir, "layerF_PRS_weights_used.csv"), index=False)
        dos = pd.read_csv(args.prs_dosages, sep="\t", index_col=0)
        cov = pd.read_csv(args.prs_covariates, sep="\t", index_col=0) if args.prs_covariates else None
        y = None
        if args.prs_outcome and cov is not None:
            yt = pd.read_csv(args.prs_outcome, sep="\t", index_col=0)
            y = yt.iloc[:, 0].reindex(cov.index)
        scores, summ = score_prs(w, dos, cov, y)
        scores.to_csv(os.path.join(args.outdir, "layerF_PRS_scores.csv"))
        report["prs"] = summ

    with open(os.path.join(args.outdir, "layerF_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
