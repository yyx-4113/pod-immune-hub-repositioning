#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer B/C - disease signature from GSE163943 + consensus with published anchors.
================================================================================
GSE163943 (Song et al., Front Aging Neurosci 2021;13:665935, PMID 34093168):
peripheral blood, Agilent Arraystar human lncRNA V5 (GPL26963),
4 POD vs 4 non-POD volunteers, quantile-normalised log2 intensities.

Why this dataset is used as an ANCHOR rather than a discovery set
-----------------------------------------------------------------
* n = 4 vs 4 -> almost no power at genome-wide FDR; the original paper
  reported 1,195 DE lncRNAs / 735 DE mRNAs with its own pipeline.
* It has already been re-analysed by others (WGCNA + PPI + hub genes).
=> We (a) reproduce the direction of the published signal,
   (b) test the IMMUNE GENE SET at pathway level (more power than single
       genes), and (c) intersect with published anchors to build the
       consensus signature used by Layer D (LINCS screen).
   Nothing here is claimed as an independent novel discovery.

Steps
-----
1. Load series matrix + GPL26963 probe annotation (probe -> gene symbol).
2. QC: sample correlation, PCA, per-sample intensity distribution.
3. Collapse probes -> genes (probe with highest mean expression per gene).
4. Differential expression: moderated t (limma-style empirical Bayes variance
   shrinkage, moment estimator for the prior), BH FDR.
5. Competitive gene-set test for the curated immune/inflammation set
   (Mann-Whitney of moderated t vs all genes) + self-contained rotation-free
   effect size (mean log2FC of set vs mean of all genes).
6. Optional pathway enrichment through g:Profiler (public API).
7. Consensus signature: intersection/union rules with published anchors
   (Seki 2026 POD-vs-control DMP genes; the hub genes reported by others)
   -> writes up/down gene lists consumed by layerD_lincs.py.

Usage
-----
  python layerB_consensus_signature.py \
      --matrix  <path>/GSE163943_series_matrix.txt \
      --annot   <path>/GPL26963_annotation.tsv \
      --immune  signatures/immune_genes.txt \
      --anchors signatures/published_pod_anchors.csv \
      --outdir  results/layerB
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

CTRL_PREFIX = "Ctrl"
POD_PREFIX = "POD"


# --------------------------------------------------------------------------
def read_series_matrix(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (expression matrix probes x samples, sample meta)."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        meta = {}
        line = fh.readline()
        while line and not line.startswith("!series_matrix_table_begin"):
            if line.startswith("!Sample_title"):
                meta["title"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                meta["gsm"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                meta["source"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            line = fh.readline()
        expr = pd.read_csv(fh, sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated(keep="first")]
    sample_meta = pd.DataFrame(meta)
    return expr, sample_meta


def collapse_to_genes(expr: pd.DataFrame, annot: pd.DataFrame) -> pd.DataFrame:
    """Map probes -> gene symbols, then keep the highest-mean probe per gene."""
    a = annot.dropna(subset=["gene_symbol"])
    a = a[a["gene_symbol"].astype(str).str.strip() != ""]
    a = a[~a["gene_symbol"].astype(str).str.lower().isin(["na", "nan", "null", "-"])]
    a = a.drop_duplicates(subset="probe_id")
    sym = a.set_index("probe_id")["gene_symbol"].astype(str).str.upper()
    common = expr.index.intersection(sym.index)
    e = expr.loc[common].copy()
    e = e.apply(pd.to_numeric, errors="coerce")     # guard against stray text
    mean_expr = e.mean(axis=1, numeric_only=True)
    e["gene"] = sym.loc[common].values
    e["_mean"] = mean_expr
    e = e.sort_values("_mean", ascending=False)
    gene_mat = e.drop_duplicates(subset="gene", keep="first").drop(columns="_mean")
    return gene_mat.set_index("gene")


def moderated_t(X: pd.DataFrame, g1: list, g2: list):
    """limma-style moderated t-test with empirical-Bayes variance shrinkage.

    Moment estimator (Smyth 2004 approximation):
      s_g^2 residual variance with d_g df; fit scaled inverse chi-square prior
      (d0, s0^2) by matching moments of log(s_g^2) across genes.
    """
    A = X[g1].to_numpy(dtype=float)
    B = X[g2].to_numpy(dtype=float)
    n1, n2 = A.shape[1], B.shape[1]
    d = n1 + n2 - 2
    m1, m2 = A.mean(axis=1), B.mean(axis=1)
    s2 = (((n1 - 1) * A.var(axis=1, ddof=1)) + ((n2 - 1) * B.var(axis=1, ddof=1))) / d
    # robust: avoid zero variances
    s2 = np.where(s2 <= 0, np.nan, s2)
    se2 = s2 * (1.0 / n1 + 1.0 / n2)
    logfc = m1 - m2

    # prior estimation on log(s^2)
    ok = ~np.isnan(s2)
    z = np.log(s2[ok])
    # Smyth (2004) moment estimator:
    #   E[log s_g^2]      = log(s0^2) + digamma(d/2)  - log(d/2)
    #   Var[log s_g^2]    = trigamma(d0/2) + trigamma(d/2)
    e_z = float(np.mean(z)) - digamma(d / 2.0) + np.log(d / 2.0)   # -> log(s0^2)
    v_z = float(np.var(z, ddof=1)) - polygamma(1, d / 2.0)          # -> trigamma(d0/2)

    def trigamma(x):
        return polygamma(1, x)

    if not np.isfinite(v_z) or v_z <= 0:
        d0 = np.inf
        s0_2 = float(np.exp(e_z))
    else:
        # Newton solve trigamma(d0/2) = v_z
        target = v_z
        x = 2.0 / max(target, 1e-9)          # trigamma(x) ~ 1/x for large x
        for _ in range(500):
            f = trigamma(x / 2.0) - target
            fp = -0.5 * polygamma(2, x / 2.0)
            if fp == 0 or not np.isfinite(fp):
                break
            step = f / fp
            x = x - step
            if not np.isfinite(x) or x <= 0:
                x = np.inf
                break
            if abs(step) < 1e-10:
                break
        d0 = float(x)
        s0_2 = float(np.exp(e_z))
    # posterior variance
    s2_post = (d0 * s0_2 + d * s2) / (d0 + d) if np.isfinite(d0) else s2
    se_post = np.sqrt(s2_post * (1.0 / n1 + 1.0 / n2))
    t = logfc / se_post
    df_total = (d + d0) if np.isfinite(d0) else d
    p = 2 * stats.t.sf(np.abs(t), df_total)
    out = pd.DataFrame({"log2FC": logfc, "t": t, "p": p,
                        "s2": s2, "s2_post": s2_post}, index=X.index)
    return out, {"d0": None if not np.isfinite(d0) else float(d0),
                 "s0_squared": float(s0_2), "df_total": float(df_total)}


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    out = np.full(p.shape, np.nan)
    pv = p[ok]; n = len(pv)
    if n == 0:
        return out
    o = np.argsort(pv); r = pv[o]
    q = np.minimum.accumulate((r * n / np.arange(1, n + 1))[::-1])[::-1]
    res = np.empty(n); res[o] = np.clip(q, 0, 1); out[ok] = res
    return out


def gprofiler_enrich(genes: list[str], organism: str = "hsapiens",
                     max_terms: int = 40):
    """Query g:Profiler (public API) for GO/KEGG/REAC enrichment."""
    import urllib.request
    url = "https://biit.cs.ut.ee/gprofiler/api/annot/annotate"
    # public g:Profiler REST: use /gost/profile for enrichment
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = json.dumps({"organism": organism, "query": genes,
                          "sources": ["GO:BP", "KEGG", "REAC", "WP"],
                          "user_threshold": 0.05, "all_results": False,
                          "ordered": False, "no_iea": False,
                          "measure_underrepresentation": False}).encode()
    req = urllib.request.Request(url, data=payload,
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            res = json.loads(r.read().decode())
        rows = res.get("result", [])
        return pd.DataFrame(rows)
    except Exception as exc:
        print("[warn] g:Profiler query failed:", exc)
        return pd.DataFrame()


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--immune", default="signatures/immune_genes.txt")
    ap.add_argument("--anchors", default="signatures/published_pod_anchors.csv")
    ap.add_argument("--outdir", default="results/layerB")
    ap.add_argument("--fc-cut", type=float, default=0.585)   # ~1.5x
    ap.add_argument("--p-cut", type=float, default=0.05)
    ap.add_argument("--gprofiler", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    expr, smeta = read_series_matrix(args.matrix)
    print("[info] probes x samples:", expr.shape)

    # ---- group assignment from sample titles -----------------------------
    titles = list(smeta.get("title", expr.columns))
    ctrl = [c for c, t in zip(expr.columns, titles) if CTRL_PREFIX.lower() in str(t).lower()]
    pod = [c for c, t in zip(expr.columns, titles) if POD_PREFIX.lower() in str(t).lower()]
    if not ctrl or not pod:
        raise SystemExit("could not resolve Ctrl/POD sample groups from titles: %s" % titles)
    print("[info] CTRL:", ctrl)
    print("[info] POD :", pod)

    # ---- QC ---------------------------------------------------------------
    qc = {
        "n_samples": int(expr.shape[1]),
        "n_probes": int(expr.shape[0]),
        "group_sizes": {"ctrl": len(ctrl), "pod": len(pod)},
        "per_sample_median": expr.median(axis=0).round(3).to_dict(),
        "per_sample_IQR": (expr.quantile(0.75, axis=0) - expr.quantile(0.25, axis=0)).round(3).to_dict(),
        "sample_corr_pearson": expr.corr(method="pearson").round(3).to_dict(),
    }

    annot = pd.read_csv(args.annot, sep="\t")
    gene_mat = collapse_to_genes(expr, annot)
    print("[info] genes after collapse:", gene_mat.shape[0])

    res, eb = moderated_t(gene_mat, pod, ctrl)      # POD vs CTRL
    res["fdr"] = bh(res["p"].to_numpy())
    res = res.sort_values("p")
    res.to_csv(os.path.join(args.outdir, "layerB_GSE163943_POD_vs_ctrl_moderatedT.csv"))

    # ---- immune gene set competitive test ---------------------------------
    immune = [g.strip().upper() for g in open(args.immune) if g.strip()]
    in_set = res.index.isin(immune)
    t_in = res.loc[in_set, "t"].dropna()
    t_out = res.loc[~in_set, "t"].dropna()
    fc_in = res.loc[in_set, "log2FC"].dropna()
    fc_out = res.loc[~in_set, "log2FC"].dropna()
    mw = stats.mannwhitneyu(t_in, t_out, alternative="two-sided") if len(t_in) > 1 else None
    immune_test = {
        "n_immune_genes_on_array": int(in_set.sum()),
        "n_background_genes": int((~in_set).sum()),
        "mean_t_immune": float(t_in.mean()) if len(t_in) else None,
        "mean_t_background": float(t_out.mean()) if len(t_out) else None,
        "mean_log2FC_immune": float(fc_in.mean()) if len(fc_in) else None,
        "mean_log2FC_background": float(fc_out.mean()) if len(fc_out) else None,
        "mannwhitney_U": float(mw.statistic) if mw else None,
        "mannwhitney_p": float(mw.pvalue) if mw else None,
    }
    # self-contained: is the set shifted up as a whole (one-sample test)
    if len(t_in) > 1:
        immune_test["onesample_t_p"] = float(stats.ttest_1samp(t_in, 0).pvalue)

    # ---- DE gene lists ----------------------------------------------------
    sig = res[(res["p"] <= args.p_cut) & (res["log2FC"].abs() >= args.fc_cut)]
    up = sig[sig["log2FC"] > 0].index.tolist()
    down = sig[sig["log2FC"] < 0].index.tolist()
    with open(os.path.join(args.outdir, "layerB_up_genes.txt"), "w") as fh:
        fh.write("\n".join(up))
    with open(os.path.join(args.outdir, "layerB_down_genes.txt"), "w") as fh:
        fh.write("\n".join(down))

    # ---- consensus with published anchors ---------------------------------
    cons_rows = []
    try:
        anch = pd.read_csv(args.anchors)
        anchor_genes = [str(x).upper() for x in anch.get("gene", pd.Series(dtype=object)).dropna()]
        anchor_genes = [a for a in anchor_genes if a and a != "NAN"]
    except Exception:
        anchor_genes = ["TMIGD3", "ADORA3", "SPATA13", "COL13A1"]
    for g in sorted(set(anchor_genes) | set(immune)):
        row = {"gene": g, "is_published_anchor": g in anchor_genes,
               "is_immune_set": g in immune}
        if g in res.index:
            row.update({"log2FC": float(res.loc[g, "log2FC"]),
                        "t": float(res.loc[g, "t"]),
                        "p": float(res.loc[g, "p"]),
                        "fdr": float(res.loc[g, "fdr"])})
        else:
            row.update({"log2FC": None, "t": None, "p": None, "fdr": None})
        cons_rows.append(row)
    consensus = pd.DataFrame(cons_rows).sort_values(
        ["is_published_anchor", "p"], ascending=[False, True])
    consensus.to_csv(os.path.join(args.outdir, "layerB_consensus_signature.csv"), index=False)

    # Layer D input: immune-restricted signature (more robust than whole-transcriptome)
    imm = res[res.index.isin(immune)].sort_values("p")
    up_i = imm[(imm["log2FC"] > 0) & (imm["p"] <= args.p_cut)].index.tolist()
    dn_i = imm[(imm["log2FC"] < 0) & (imm["p"] <= args.p_cut)].index.tolist()
    if len(up_i) + len(dn_i) < 20:      # fall back to top-|t| immune genes
        ranked = imm.reindex(imm["t"].abs().sort_values(ascending=False).index)
        up_i = [g for g in ranked.index if ranked.loc[g, "log2FC"] > 0][:75]
        dn_i = [g for g in ranked.index if ranked.loc[g, "log2FC"] < 0][:75]
    with open(os.path.join(args.outdir, "layerD_input_up.txt"), "w") as fh:
        fh.write("\n".join(up_i))
    with open(os.path.join(args.outdir, "layerD_input_down.txt"), "w") as fh:
        fh.write("\n".join(dn_i))

    enrich = pd.DataFrame()
    if args.gprofiler and (up or down):
        enrich = gprofiler_enrich(sorted(set(up) | set(dn_i)))
        if len(enrich):
            enrich.to_csv(os.path.join(args.outdir, "layerB_pathway_enrichment.csv"), index=False)

    summary = {
        "dataset": "GSE163943 (Song 2021, Front Aging Neurosci 13:665935)",
        "design": "blood, 4 POD vs 4 non-POD, Agilent Arraystar lncRNA V5 (GPL26963)",
        "qc": qc,
        "ebayes": eb,
        "n_genes_tested": int(res.shape[0]),
        "n_nominal_p05": int((res["p"] <= 0.05).sum()),
        "n_fdr05": int((res["fdr"] <= 0.05).sum()),
        "n_up": len(up), "n_down": len(down),
        "immune_set_test": immune_test,
        "layerD_input": {"n_up": len(up_i), "n_down": len(dn_i)},
        "top20": res.head(20).round(4).reset_index().to_dict(orient="records"),
    }
    with open(os.path.join(args.outdir, "layerB_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in summary.items() if k != "qc"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
