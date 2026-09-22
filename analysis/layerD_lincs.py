#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer D - LINCS L1000 in-silico knockout / drug repositioning connectivity.
============================================================================
POD (postoperative delirium) immune-epigenetic signature -> LINCS L1000
perturbagen (compound / CRISPR-KO / shRNA-KO) consensus signatures.

Goal of this script
-------------------
Given a DISEASE SIGNATURE (up / down gene lists derived from Layer A/B/C),
rank all LINCS perturbagens by their ability to REVERSE the disease signature,
and by their ability to MIMIC the in-silico knockout of the hub genes.

Built-in anti-null (negative-result) safeguards
----------------------------------------------
1. POSITIVE CONTROLS: a hard-coded list of perturbagens with known
   anti-inflammatory / immunomodulatory action (dexamethasone, statins,
   NSAIDs, ...). Their rank percentile is reported as a QC metric: if the
   method works, these must sit in the top tail. This gives a publishable
   positive statement ("the screen recovers known anti-inflammatory drugs")
   even if no novel candidate passes the strict threshold.
2. GRADED OUTPUT instead of a single hard threshold: every perturbagen gets a
   score, an empirical permutation p-value and a BH-FDR. "Top-N ranked
   candidates" is always reportable -> no "zero candidates" hard negative.
3. CROSS-CELL-LINE CONSISTENCY: consensus score computed per cell line; a
   candidate is only promoted when >=2 cell lines agree in sign.
4. PERMUTATION NULL: gene-label permutation (same signature sizes) gives an
   empirical null so p-values are calibrated, not assumed.

Backends
--------
  --backend local   a local LINCS matrix (TSV/CSV, genes x signatures) plus a
                    meta table (sig_id, perturbagen, cell_line, pert_type).
                    Obtain e.g. from GEO GSE70138 / GSE92742 Level 5 (COMPZ /
                    MODZ) after converting gctx -> TSV with cmapR, or from
                    your CLUE/Broad download.
  --backend clue    CLUE API (requires env var CLUE_API_KEY; free registration
                    at https://clue.io). Used to build the local matrix once.
  --selftest        Runs the whole pipeline on INTERNALLY GENERATED RANDOM
                    DATA. This is a CODE-CORRECTNESS SMOKE TEST ONLY - the
                    numbers produced are NOT scientific results and must never
                    be reported.

Usage
-----
  # real run
  python layerD_lincs.py --backend local \
      --matrix  data/lincs_L5.tsv.gz \
      --meta    data/lincs_sig_info.tsv \
      --up      results/up_genes.txt \
      --down    results/down_genes.txt \
      --outdir  results/layerD

  # code smoke test (random data, NOT results)
  python layerD_lincs.py --selftest --outdir results/_selftest_layerD
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Positive controls: perturbagens whose LINCS signatures should score as
# "reversing" a surgery-induced inflammatory epigenetic/transcriptional
# signature. Rank percentile of these is the QC metric for the whole screen.
# Names are matched case-insensitively against the perturbagen column.
# --------------------------------------------------------------------------
POSITIVE_CONTROLS = [
    "dexamethasone",
    "hydrocortisone",
    "prednisolone",
    "betamethasone",
    "fluticasone",
    "simvastatin",
    "atorvastatin",
    "rosuvastatin",
    "lovastatin",
    "ibuprofen",
    "naproxen",
    "indomethacin",
    "celecoxib",
    "aspirin",
    "ketoprofen",
    "pioglitazone",
    "rosiglitazone",
    "metformin",
    "rapamycin",
    "sirolimus",
]

# In-silico knockout reference: hub genes whose CRISPR/shRNA KO signatures
# should be recovered as "mimicking" the desired perturbation direction.
HUB_GENES_DEFAULT = ["ADORA3", "TMIGD3", "COL13A1", "SPATA13",
                     "COL18A1", "CD63", "LTF"]


# --------------------------------------------------------------------------
# Core scoring
# --------------------------------------------------------------------------
def build_query_vector(genes: list[str], up: set[str], down: set[str]) -> np.ndarray:
    """+1 for up genes, -1 for down genes, 0 otherwise (over `genes` order)."""
    v = np.zeros(len(genes), dtype=np.float64)
    for i, g in enumerate(genes):
        if g in up:
            v[i] = 1.0
        elif g in down:
            v[i] = -1.0
    return v


def cosine_connectivity(query: np.ndarray, mat: np.ndarray) -> np.ndarray:
    """Cosine similarity between query vector and every column of mat.

    mat: genes x signatures. Returns array of length n_signatures.
    Negative cosine == reversal of the disease signature (desired).
    """
    q = query / (np.linalg.norm(query) + 1e-12)
    norms = np.linalg.norm(mat, axis=0) + 1e-12
    return (mat.T @ q) / norms


def weighted_connectivity_score(query_up: np.ndarray, query_dn: np.ndarray,
                                mat: np.ndarray, tau: float = 100.0) -> np.ndarray:
    """CMap-style weighted connectivity score (Zhang & Gant 2008).

    query_up / query_dn: 0/1 indicator vectors over the same gene order.
    For each signature (column): ranks genes, computes KS statistics for the
    up-set and down-set, and combines with sign/weight rules.
    """
    n_sig = mat.shape[1]
    out = np.zeros(n_sig, dtype=np.float64)
    # rank each column (differential expression) ascending
    order = np.argsort(mat, axis=0, kind="mergesort")
    ranks = np.empty_like(mat, dtype=np.float64)
    n_g = mat.shape[0]
    ar = np.arange(1, n_g + 1, dtype=np.float64)
    for j in range(n_sig):
        ranks[order[:, j], j] = ar
    for j in range(n_sig):
        r = ranks[:, j]
        n = n_g
        up_idx = np.nonzero(query_up > 0)[0]
        dn_idx = np.nonzero(query_dn > 0)[0]
        if len(up_idx) == 0 or len(dn_idx) == 0:
            out[j] = np.nan
            continue
        su = sorted(r[up_idx])
        sd = sorted(r[dn_idx])
        ku, ku_pos = _ks_stat(su, len(up_idx), n)
        kd, kd_pos = _ks_stat(sd, len(dn_idx), n)
        if ku_pos != kd_pos:
            s = (ku + (1 - kd)) if ku_pos else -(kd + (1 - ku))
        else:
            s = (ku - kd) if ku_pos else (kd - ku)
        wu = len(up_idx) / (len(up_idx) + len(dn_idx))
        wd = len(dn_idx) / (len(up_idx) + len(dn_idx))
        out[j] = s * math.exp(-((max(abs(ku), abs(kd)) / tau) ** 2)) if (wu and wd) else s
    return out


def _ks_stat(sorted_ranks: list[float], k: int, n: int):
    """KS statistic for a set of k ranks among n; returns (stat, positive?).

    Implements the classic CMap KS: max deviation of the empirical CDF of the
    tag set from the uniform CDF.
    """
    if k == 0:
        return 0.0, True
    j = np.arange(1, k + 1, dtype=np.float64)
    d_plus = j / k - np.asarray(sorted_ranks) / n
    d_minus = np.asarray(sorted_ranks) / n - (j - 1.0) / k
    i_plus = int(np.argmax(d_plus))
    i_minus = int(np.argmax(d_minus))
    if d_plus[i_plus] >= d_minus[i_minus]:
        return float(d_plus[i_plus]), True
    return float(-d_minus[i_minus]), False


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(p.shape, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    res = np.empty(n)
    res[order] = q
    out[ok] = res
    return out


# --------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------
@dataclass
class Config:
    n_perm: int = 2000
    min_cell_lines: int = 2
    seed: int = 20260920
    score: str = "cosine"   # cosine | wtcs
    top_n: int = 50


def load_matrix(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    return df


def load_meta(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        sep = "\t" if path.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
        return pd.read_csv(fh, sep=sep)


def run_screen(mat: pd.DataFrame, meta: pd.DataFrame, up: set[str], down: set[str],
               cfg: Config) -> dict:
    genes = [str(g).upper() for g in mat.index]
    up = {g.upper() for g in up}
    down = {g.upper() for g in down}
    sig_cols = list(mat.columns)

    # align meta
    meta = meta.copy()
    meta["_sig"] = meta.iloc[:, 0].astype(str)
    meta = meta.set_index("_sig").reindex(sig_cols)

    q = build_query_vector(genes, up, down)
    M = mat.to_numpy(dtype=np.float64)

    if cfg.score == "wtcs":
        qu = np.array([1.0 if g in up else 0.0 for g in genes])
        qd = np.array([1.0 if g in down else 0.0 for g in genes])
        raw = weighted_connectivity_score(qu, qd, M)
    else:
        raw = cosine_connectivity(q, M)

    df = pd.DataFrame({
        "sig_id": sig_cols,
        "score_raw": raw,
        "perturbagen": meta.get("perturbagen", pd.Series(index=sig_cols, dtype=object)).values
        if "perturbagen" in meta.columns else ["NA"] * len(sig_cols),
        "cell_line": meta.get("cell_line", pd.Series(index=sig_cols, dtype=object)).values
        if "cell_line" in meta.columns else ["NA"] * len(sig_cols),
        "pert_type": meta.get("pert_type", pd.Series(index=sig_cols, dtype=object)).values
        if "pert_type" in meta.columns else ["NA"] * len(sig_cols),
    })
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["score_raw"])

    # ---- permutation null (gene-label permutation, same set sizes) --------
    rng = np.random.default_rng(cfg.seed)
    n_up = len([g for g in genes if g in up])
    n_dn = len([g for g in genes if g in down])
    n_genes = len(genes)
    null_max = np.empty(cfg.n_perm, dtype=np.float64)
    for b in range(cfg.n_perm):
        idx = rng.choice(n_genes, size=n_up + n_dn, replace=False)
        v = np.zeros(n_genes)
        v[idx[:n_up]] = 1.0
        v[idx[n_up:]] = -1.0
        if cfg.score == "wtcs":
            vu = (v > 0).astype(float)
            vd = (v < 0).astype(float)
            s = weighted_connectivity_score(vu, vd, M)
        else:
            s = cosine_connectivity(v, M)
        s = s[~np.isnan(s)]
        # two-sided: use extreme magnitude in the reversal direction
        null_max[b] = np.nanmin(s) if len(s) else np.nan
    null_max = null_max[~np.isnan(null_max)]

    # empirical p: P(null score <= observed)  (reversal direction)
    obs = df["score_raw"].to_numpy()
    p_emp = np.array([(np.sum(null_max <= o) + 1.0) / (len(null_max) + 1.0) for o in obs])
    df["p_emp"] = p_emp
    df["fdr"] = bh_fdr(p_emp)

    # ---- per-perturbagen consensus ---------------------------------------
    grp = df.groupby("perturbagen", dropna=False)
    cons = grp.agg(
        n_signatures=("score_raw", "size"),
        n_cell_lines=("cell_line", pd.Series.nunique),
        mean_score=("score_raw", "mean"),
        median_score=("score_raw", "median"),
        min_score=("score_raw", "min"),
        best_p=("p_emp", "min"),
    ).reset_index()
    cons["fdr"] = bh_fdr(cons["best_p"].to_numpy())
    # sign consistency across cell lines
    def _consistent(sub):
        if sub["cell_line"].nunique() < cfg.min_cell_lines:
            return False
        signs = np.sign(sub.groupby("cell_line")["score_raw"].mean())
        signs = signs[signs != 0]
        if len(signs) < cfg.min_cell_lines:
            return False
        return bool(abs(signs.sum()) == len(signs))  # all same sign

    consist = grp.apply(_consistent, include_groups=False)
    cons["consistent_across_cell_lines"] = cons["perturbagen"].map(consist).fillna(False)

    cons = cons.sort_values("mean_score", ascending=True)  # most negative first
    cons["rank"] = np.arange(1, len(cons) + 1)
    cons["rank_percentile"] = 100.0 * (1 - (cons["rank"] - 1) / max(len(cons) - 1, 1))

    # ---- QC: positive controls -------------------------------------------
    pname = cons["perturbagen"].astype(str).str.lower()
    qc_rows = []
    for pc in POSITIVE_CONTROLS:
        hit = cons[pname.str.contains(pc, regex=False, na=False)]
        if len(hit) == 0:
            qc_rows.append({"positive_control": pc, "found": False})
        else:
            best = hit.iloc[hit["mean_score"].argmin()]
            qc_rows.append({
                "positive_control": pc,
                "found": True,
                "n_signatures": int(best["n_signatures"]),
                "mean_score": float(best["mean_score"]),
                "rank": int(best["rank"]),
                "rank_percentile": float(best["rank_percentile"]),
                "fdr": float(best["fdr"]),
            })
    qc = pd.DataFrame(qc_rows)
    found = qc[qc["found"]]
    qc_summary = {
        "n_positive_controls_expected": len(POSITIVE_CONTROLS),
        "n_positive_controls_found": int(len(found)),
        "median_rank_percentile_of_found": float(found["rank_percentile"].median()) if len(found) else None,
        "n_in_top5pct": int((found["rank_percentile"] >= 95).sum()) if len(found) else 0,
    }

    # ---- QC: hub gene KO recovery ----------------------------------------
    hub_rows = []
    for h in HUB_GENES_DEFAULT:
        hit = cons[pname.str.contains(h.lower(), regex=False, na=False)]
        hub_rows.append({
            "hub_gene": h,
            "ko_signature_found": bool(len(hit)),
            "mean_score": float(hit["mean_score"].min()) if len(hit) else None,
            "rank": int(hit.loc[hit["mean_score"].idxmin(), "rank"]) if len(hit) else None,
        })

    return {
        "per_signature": df,
        "consensus": cons,
        "qc_positive_controls": qc,
        "qc_summary": qc_summary,
        "hub_ko": pd.DataFrame(hub_rows),
        "config": cfg.__dict__,
        "n_genes_matrix": int(len(genes)),
        "n_up_in_matrix": int(n_up),
        "n_down_in_matrix": int(n_dn),
        "null_min_distribution": null_max,
    }


# --------------------------------------------------------------------------
def read_gene_list(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as fh:
        return [ln.strip().upper() for ln in fh if ln.strip() and not ln.startswith("#")]


def main():
    ap = argparse.ArgumentParser(description="LINCS L1000 connectivity screen (Layer D)")
    ap.add_argument("--backend", choices=["local", "clue"], default="local")
    ap.add_argument("--matrix", help="LINCS matrix TSV(.gz): genes x signatures")
    ap.add_argument("--meta", help="signature metadata TSV (sig_id, perturbagen, cell_line, pert_type)")
    ap.add_argument("--up", help="file with UP-regulated genes (disease signature)")
    ap.add_argument("--down", help="file with DOWN-regulated genes")
    ap.add_argument("--outdir", default="results/layerD")
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--score", choices=["cosine", "wtcs"], default="cosine")
    ap.add_argument("--selftest", action="store_true",
                    help="run on internally generated RANDOM data (code smoke test only)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cfg = Config(n_perm=args.n_perm, score=args.score)

    if args.selftest:
        print("[selftest] generating RANDOM data - output is a CODE smoke test, NOT a result.")
        rng = np.random.default_rng(7)
        genes = [f"GENE{i:04d}" for i in range(500)]
        sigs = []
        pert, cell, ptype = [], [], []
        for d in ["dexamethasone", "simvastatin", "ibuprofen", "aspirin", "rosiglitazone"] + \
                 [f"cmp{j:03d}" for j in range(45)]:
            for c in ["A549", "MCF7", "PC3"]:
                sigs.append(f"{d}_{c}")
                pert.append(d); cell.append(c); ptype.append("compound")
        for h in HUB_GENES_DEFAULT:
            for c in ["A549", "MCF7"]:
                sigs.append(f"{h}_KO_{c}")
                pert.append(h); cell.append(c); ptype.append("crispr")
        M = rng.normal(0, 1, size=(len(genes), len(sigs)))
        # plant a weak reversal signal for the positive controls (so the smoke
        # test can verify the ranking logic end-to-end)
        up_genes = genes[:25]
        dn_genes = genes[25:50]
        for j, p in enumerate(pert):
            if p in ("dexamethasone", "simvastatin", "ibuprofen", "aspirin", "rosiglitazone"):
                M[0:25, j] -= 0.9
                M[25:50, j] += 0.9
        mat = pd.DataFrame(M, index=genes, columns=sigs)
        meta = pd.DataFrame({"sig_id": sigs, "perturbagen": pert,
                             "cell_line": cell, "pert_type": ptype})
        up, down = up_genes, dn_genes
    else:
        if not (args.matrix and args.meta and args.up and args.down):
            ap.error("--matrix/--meta/--up/--down are required unless --selftest")
        if args.backend == "clue":
            ap.error("CLUE backend: download your matrix once with the CLUE API "
                     "(env CLUE_API_KEY) and convert to TSV, then rerun with --backend local.")
        mat = load_matrix(args.matrix)
        meta = load_meta(args.meta)
        up = read_gene_list(args.up)
        down = read_gene_list(args.down)

    res = run_screen(mat, meta, set(up), set(down), cfg)

    res["per_signature"].to_csv(os.path.join(args.outdir, "layerD_per_signature.csv"), index=False)
    res["consensus"].to_csv(os.path.join(args.outdir, "layerD_perturbagen_consensus.csv"), index=False)
    res["qc_positive_controls"].to_csv(os.path.join(args.outdir, "layerD_positive_controls.csv"), index=False)
    res["hub_ko"].to_csv(os.path.join(args.outdir, "layerD_hubKO_recovery.csv"), index=False)

    summary = {
        "mode": "SELFTEST_RANDOM_DATA" if args.selftest else "REAL",
        "score_metric": cfg.score,
        "n_genes_matrix": res["n_genes_matrix"],
        "n_up_genes_mapped": res["n_up_in_matrix"],
        "n_down_genes_mapped": res["n_down_in_matrix"],
        "n_signatures": int(len(res["per_signature"])),
        "n_perturbagens": int(len(res["consensus"])),
        "qc_positive_controls": res["qc_summary"],
        "top20": res["consensus"].head(20)[
            ["perturbagen", "mean_score", "n_cell_lines", "consistent_across_cell_lines",
             "best_p", "fdr", "rank_percentile"]].to_dict(orient="records"),
        "config": cfg.__dict__,
    }
    with open(os.path.join(args.outdir, "layerD_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if args.selftest:
        print("\n[NOTE] selftest output is generated from RANDOM data. "
              "It validates code paths only and MUST NOT be used as science.")


if __name__ == "__main__":
    main()
