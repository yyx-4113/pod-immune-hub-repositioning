#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer C2 - PBMC single-cell localisation of the POD immune hubs (GSE252572).
============================================================================
Pure-Python (numpy / scipy / pandas / sklearn) implementation - no scanpy,
no Seurat, no numba. Written this way because scanpy silently crashes in
sandboxed environments (numba cache) and R/Seurat may not be available.

Input: CellRanger output per sample, i.e. directories containing
       matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz
       (GSE252572_RAW.tar, 339 MB, contains 8 such samples).

Design decisions (pre-registered)
---------------------------------
* Cell-type calls are made by CANONICAL MARKER SCORES per cluster, and the
  script always prints the top marker genes per cluster so the annotation can
  be AUDITED rather than trusted.
* **NO cell-level differential statistics.** All group comparisons are done at
  the SAMPLE level (pseudobulk per sample, or per-sample cell-type
  proportions) to avoid pseudoreplication / double-dipping. With n=8 samples
  this is the only defensible test.
* If the deposited sample table does not carry POD status, the script stops
  and says so instead of guessing.

Usage
-----
  python layerC_scrna.py --data-dir data/GSE252572 --outdir results/layerC \
      [--sample-meta samples.tsv] [--hubs ADORA3,TMIGD3,COL18A1,CD63,LTF]
"""
from __future__ import annotations

import argparse
import gzip
import os
from glob import glob

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse as sp
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD

# Canonical PBMC marker sets (score = mean scaled expression within cluster)
MARKERS = {
    "CD4_T": ["IL7R", "CCR7", "LDHB", "CD3D", "CD3E", "LTB", "TRAC"],
    "CD8_T": ["CD8A", "CD8B", "GZMK", "CCL5", "NKG7", "CD3D"],
    "NK": ["GNLY", "NKG7", "GZMB", "KLRD1", "PRF1", "FGFBP2"],
    "B": ["MS4A1", "CD79A", "CD79B", "IGHM", "TCL1A", "BANK1"],
    "Mono_classical": ["LYZ", "LST1", "S100A8", "S100A9", "FCN1", "VCAN", "CD14"],
    "Mono_nonclassical": ["FCGR3A", "MS4A7", "CST3", "LST1", "CDKN1C"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "CD1C", "AIF1"],
    "Platelet": ["PPBP", "PF4", "GP9", "ITGA2B"],
}


def read_10x(sample_dir: str):
    """Return (genes x cells) sparse matrix + gene symbols + barcodes."""
    mtx = glob(os.path.join(sample_dir, "**", "matrix.mtx.gz"), recursive=True)
    if not mtx:
        mtx = glob(os.path.join(sample_dir, "**", "matrix.mtx"), recursive=True)
    if not mtx:
        raise FileNotFoundError("no matrix.mtx(.gz) under %s" % sample_dir)
    mtx_path = mtx[0]
    base = os.path.dirname(mtx_path)

    def opener(p):
        return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

    feats = glob(os.path.join(base, "features.tsv.gz")) + glob(os.path.join(base, "genes.tsv.gz"))
    feats += glob(os.path.join(base, "features.tsv")) + glob(os.path.join(base, "genes.tsv"))
    if not feats:
        raise FileNotFoundError("no features.tsv under %s" % base)
    ft = pd.read_csv(opener(feats[0]), sep="\t", header=None)
    sym = ft.iloc[:, 1].astype(str).str.upper().values

    bc = glob(os.path.join(base, "barcodes.tsv.gz")) + glob(os.path.join(base, "barcodes.tsv"))
    barcodes = None
    if bc:
        barcodes = pd.read_csv(opener(bc[0]), sep="\t", header=None)[0].astype(str).values

    with (gzip.open(mtx_path, "rb") if mtx_path.endswith(".gz") else open(mtx_path, "rb")) as fh:
        M = sio.mmread(fh).tocsr()
    return M, sym, barcodes


def qc_filter(M, sym, min_genes=200, max_mito=0.15):
    """Cell QC: gene count and mitochondrial fraction."""
    mito = np.array([g.startswith("MT-") for g in sym])
    n_genes = np.asarray((M > 0).sum(axis=0)).ravel()
    tot = np.asarray(M.sum(axis=0)).ravel()
    mito_frac = np.asarray(M[mito].sum(axis=0)).ravel() / np.maximum(tot, 1)
    keep = (n_genes >= min_genes) & (mito_frac <= max_mito)
    info = {"n_cells_raw": int(M.shape[1]), "n_cells_kept": int(keep.sum()),
            "median_genes_per_cell": float(np.median(n_genes)) if len(n_genes) else 0.0,
            "median_mito_frac": float(np.median(mito_frac)) if len(mito_frac) else 0.0}
    # Robustness: some datasets (e.g. shallow or targeted assays) would lose
    # every cell to the default threshold. Fall back to a data-driven cut-off
    # and record that it happened, rather than crashing downstream.
    if keep.sum() < 50 and len(n_genes):
        alt = max(50, int(np.median(n_genes) * 0.5))
        keep = (n_genes >= alt) & (mito_frac <= max_mito)
        info["qc_fallback_min_genes"] = int(alt)
        info["n_cells_kept"] = int(keep.sum())
    if keep.sum() == 0:
        raise ValueError("QC removed every cell for this sample - check the matrix")
    return keep, info


def normalise_and_cluster(M, n_hvg=2000, n_pcs=30, n_clusters=12, seed=12345):
    """CPM normalise -> log1p -> HVG -> SVD -> KMeans."""
    Mc = M.tocsc().astype(np.float64)
    tot = np.asarray(Mc.sum(axis=0)).ravel()
    tot[tot == 0] = 1.0
    X = Mc.multiply(1e4 / tot).tocsc()
    X = X.tocsr()
    X.data = np.log1p(X.data)

    mean = np.asarray(X.mean(axis=1)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=1)).ravel()
    var = sq - mean ** 2
    idx = np.argsort(-var)[:n_hvg]

    Xi = np.asarray(X[idx].todense()).T          # cells x hvg
    Xi = (Xi - Xi.mean(0)) / (Xi.std(0) + 1e-9)
    Xi = np.clip(Xi, -10, 10)
    svd = TruncatedSVD(n_components=min(n_pcs, Xi.shape[1] - 1), random_state=seed)
    pcs = svd.fit_transform(Xi)
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=seed).fit(pcs)
    return km.labels_, pcs, idx


def score_clusters(X_log, sym, labels, top_n=15):
    """Mean expression per cluster for canonical markers + top markers per cluster."""
    labs = np.unique(labels)
    rows = []
    sym_idx = {g: i for i, g in enumerate(sym)}
    for cl in labs:
        sel = np.nonzero(labels == cl)[0]
        sub = X_log[:, sel]
        m = np.asarray(sub.mean(axis=1)).ravel()
        row = {"cluster": int(cl), "n_cells": int(len(sel))}
        for ct, genes in MARKERS.items():
            gidx = [sym_idx[g] for g in genes if g in sym_idx]
            row["score_" + ct] = float(np.mean(m[gidx])) if gidx else np.nan
        top = np.argsort(-m)[:top_n]
        row["top_markers"] = ";".join(sym[top])
        rows.append(row)
    df = pd.DataFrame(rows)
    marker_cols = [c for c in df.columns if c.startswith("score_")]
    df["called_type"] = df[marker_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    return df


def align_genes(mats, syms):
    """Restrict every sample to the Intersection of gene symbols, in one order.

    Returns (list_of_genes_x_cells matrices, common_symbol_list).
    """
    common = set(syms[0])
    for s in syms[1:]:
        common &= set(s)
    order = [g for g in syms[0] if g in common]
    pos = {g: i for i, g in enumerate(order)}
    out = []
    for M, sym in zip(mats, syms):
        rows = np.array([pos[g] for g in sym if g in pos], dtype=int)
        cols = np.arange(len(rows), dtype=int)
        P = sp.csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)),
                          shape=(len(order), len(sym)))
        out.append((P @ M).tocsr())
    return out, order


def global_cluster(mats_genes_x_cells, n_hvg=2000, n_pcs=30, n_clusters=12, seed=12345):
    """POOLED clustering across all samples.

    WHY THIS MATTERS: an earlier version of this script clustered each sample
    separately. Cluster IDs produced that way are arbitrary per sample, so
    comparing "cluster k proportion" between groups compares different cell
    populations and is meaningless. Pooling guarantees that cluster k means the
    same population in every sample. Fixed 2026-09-20.
    """
    Xg = sp.hstack(mats_genes_x_cells).tocsc()            # genes x all cells
    Xc = Xg.T.tocsr()                                      # cells x genes
    mean = np.asarray(Xc.mean(axis=0)).ravel()
    sq = np.asarray(Xc.multiply(Xc).mean(axis=0)).ravel()
    var = sq - mean ** 2
    idx = np.argsort(-var)[:n_hvg]
    Xi = Xc[:, idx]
    svd = TruncatedSVD(n_components=min(n_pcs, Xi.shape[1] - 1), random_state=seed)
    pcs = svd.fit_transform(Xi)
    pcs = (pcs - pcs.mean(0)) / (pcs.std(0) + 1e-9)
    pcs = np.clip(pcs, -10, 10)
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=seed).fit(pcs)
    return km.labels_, idx, Xg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--outdir", default="results/layerC")
    ap.add_argument("--sample-meta", help="TSV: sample,group(POD/CTRL)")
    ap.add_argument("--hubs", default="ADORA3,TMIGD3,COL18A1,CD63,LTF")
    ap.add_argument("--n-clusters", type=int, default=12)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    samples = sorted([d for d in glob(os.path.join(args.data_dir, "*"))
                      if os.path.isdir(d)])
    if not samples:
        samples = [args.data_dir]
    print("[info] samples:", [os.path.basename(s) for s in samples])

    # ---- pass 1: read + QC + per-sample CPM/log1p -------------------------
    mats, syms, names, qc_rows = [], [], [], []
    for s in samples:
        name = os.path.basename(s.rstrip("/\\"))
        M, sym, bc = read_10x(s)
        keep, q = qc_filter(M, sym)
        q["sample"] = name
        qc_rows.append(q)
        Mf = M[:, keep].tocsc().astype(np.float64)
        tot = np.asarray(Mf.sum(axis=0)).ravel()
        tot[tot == 0] = 1.0
        Xl = Mf.multiply(1e4 / tot).tocsr()
        Xl.data = np.log1p(Xl.data)
        mats.append(Xl)
        syms.append(list(sym))
        names.append(name)

    # ---- pass 2: POOLED clustering (cluster IDs comparable across samples) --
    mats, sym_common = align_genes(mats, syms)
    sym_arr = np.array(sym_common)
    n_cells_pooled = sum(m.shape[1] for m in mats)
    print(f"[info] pooled: {n_cells_pooled} cells x {len(sym_common)} genes")
    labels, hvg, Xpool = global_cluster(mats, n_clusters=args.n_clusters)

    ct = score_clusters(Xpool, sym_arr, labels)
    ct.to_csv(os.path.join(args.outdir, "layerC_cluster_markers.csv"), index=False)
    print("[audit] global cluster annotation (MUST be checked against top_markers):")
    for _, r in ct.iterrows():
        print(f"   cluster {int(r['cluster']):>2}  n={int(r['n_cells']):>6}  "
              f"{r['called_type']:<18} {r['top_markers'][:90]}")
    type_of = dict(zip(ct["cluster"], ct["called_type"]))

    # ---- pass 3: per-sample composition, hubs, batch diagnostic -------------
    per_sample_props, per_sample_hub, comp_rows = [], [], []
    sym_idx = {g: i for i, g in enumerate(sym_common)}
    off = 0
    for name, M in zip(names, mats):
        n = M.shape[1]
        lab = labels[off:off + n]
        off += n
        cnt = pd.Series(lab).value_counts().reindex(sorted(ct["cluster"]), fill_value=0)
        props = cnt.astype(float) / max(int(cnt.sum()), 1)
        props.name = name
        per_sample_props.append(props)
        for c in sorted(ct["cluster"]):
            comp_rows.append({"cluster": int(c), "sample": name,
                              "n_cells": int((lab == c).sum()),
                              "frac_of_cluster": float((lab == c).sum() /
                                                       max(int((labels == c).sum()), 1))})
        hub_row = {"sample": name}
        for h in [x.strip().upper() for x in args.hubs.split(",") if x.strip()]:
            if h in sym_idx:
                hub_row[h] = float(np.asarray(M[sym_idx[h]].todense()).ravel().mean())
            else:
                hub_row[h] = np.nan
        per_sample_hub.append(hub_row)

    pd.DataFrame(qc_rows).to_csv(os.path.join(args.outdir, "layerC_qc_per_sample.csv"), index=False)
    props = pd.concat(per_sample_props, axis=1).T
    props.to_csv(os.path.join(args.outdir, "layerC_cluster_proportions.csv"))
    pd.DataFrame(comp_rows).to_csv(os.path.join(args.outdir, "layerC_cluster_batch_composition.csv"),
                                   index=False)
    # collapse clusters that received the same marker-based call
    tp = pd.DataFrame({name: {t: sum(float(p[c]) for c in p.index if type_of.get(c) == t)
                              for t in sorted(set(type_of.values()))}
                       for name, p in zip(names, per_sample_props)}).T
    tp.to_csv(os.path.join(args.outdir, "layerC_celltype_proportions.csv"))
    hubs = pd.DataFrame(per_sample_hub).set_index("sample")
    hubs.to_csv(os.path.join(args.outdir, "layerC_hub_expression_per_sample.csv"))

    summary = {"n_samples": len(samples), "qc": qc_rows,
               "note": "group comparison requires a sample table with POD status"}

    meta = None
    if args.sample_meta:
        meta = pd.read_csv(args.sample_meta, sep="\t").set_index("sample")
    if meta is not None and "group" in meta.columns:
        g = meta["group"].astype(str).str.upper().reindex(tp.index)
        out_rows = []
        # SAMPLE-LEVEL test of cell-type proportions (avoids pseudoreplication)
        for c in tp.columns:
            a = tp.loc[g == "POD", c].dropna()
            b = tp.loc[g != "POD", c].dropna()
            if len(a) >= 2 and len(b) >= 2:
                u = stats.mannwhitneyu(a, b, alternative="two-sided")
                out_rows.append({"cell_type": str(c),
                                 "n_POD_samples": int(len(a)), "n_CTRL_samples": int(len(b)),
                                 "mean_POD": float(a.mean()), "mean_CTRL": float(b.mean()),
                                 "log2_ratio": float(np.log2((a.mean() + 1e-6) / (b.mean() + 1e-6))),
                                 "mannwhitney_p": float(u.pvalue)})
        res = pd.DataFrame(out_rows).sort_values("mannwhitney_p")
        res.to_csv(os.path.join(args.outdir, "layerC_proportion_sample_level_test.csv"), index=False)
        summary["proportion_test"] = res.to_dict(orient="records")
        summary["annotation"] = ct[["cluster", "n_cells", "called_type", "top_markers"]].to_dict(orient="records")
        summary["warning"] = (
            "DESCRIPTIVE ONLY. Clustered POOLED across samples so that cluster IDs are "
            "comparable; proportions tested at SAMPLE level (no cell-level pseudoreplication), "
            "but with 4 POD vs 3 non-POD libraries derived from only 2+2 patients, Mann-Whitney "
            "has essentially no power. Report exact p; do NOT claim significance.")
    else:
        summary["status"] = "NO_GROUP_LABELS: sample-level comparison skipped"

    with open(os.path.join(args.outdir, "layerC_summary.json"), "w", encoding="utf-8") as fh:
        import json
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print("[done] wrote results to", args.outdir)
    print("[note] cluster annotation must be AUDITED against the top_markers column.")


if __name__ == "__main__":
    main()
