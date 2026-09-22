#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer F (real-data arm): gene-set / locus-level characterisation of the
PUBLIC postoperative-delirium (POD) GWAS (Armstrong et al., PLOS Medicine 2026;
PMID 41770756; DOI 10.5523/bris.1m83zai2e26yq2lro3tixz9kqq).

WHAT THIS SCRIPT ACTUALLY DOES (no individual-level data required):
  F1  QC of the summary statistics: variant count, genomic inflation factor,
      genome-wide-significant loci (p < 5e-8), lead SNP effect sizes.
  F2  Gene-level association scores: for every autosomal gene, the minimum
      GWAS p-value over all SNPs inside the gene body +/- 10 kb (hg19/UCSC
      refGene coordinates -- the GWAS is build GRCh37/hg19).
  F3  Competitive gene-set test for a curated peripheral immune/inflammation
      gene set, using a SNP-count-stratified permutation null (this removes
      the gene-length / SNP-density confound that a naive Mann-Whitney test
      would carry).
  F4  The same test restricted to autosomes EXCLUDING chromosome 19, i.e.
      "does any immune-gene common variation contribute to POD once APOE is
      removed?"  This is the pre-registered honest test of the peripheral
      immune hypothesis at the germline level.
  F5  Suggestive (5e-8 <= p < 1e-5) NON-chr19 loci, distance-clumped, with
      nearest gene -- reported explicitly as hypothesis-generating only.

Everything written to <outdir>. No number in the output is hard-coded; all
values are read from the downloaded files at run time.

Usage:
  python layerF_real_gwas.py \
      --gwas  <ukb_gwas_delirium_results.txt.gz> \
      --genes <refGene_hg19.txt.gz> \
      --immune signatures/immune_genes.txt \
      --outdir results/layerF
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

GW_SIG = 5e-8
SUGGESTIVE = 1e-5
PAD = 10_000          # bp padding around gene body for the gene-level score
MAX_GENE_SPAN = 3_000_000   # QC: max plausible gene span (bp)


# ----------------------------------------------------------------------------
# F1  load summary statistics / QC
# ----------------------------------------------------------------------------
def load_gwas(path: str, chunksize: int = 2_000_000) -> pd.DataFrame:
    """Stream the GWAS, keeping only the numeric columns we need."""
    chroms, poss, ps, chis = [], [], [], []
    n = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline().split()
        idx = {h: i for i, h in enumerate(header)}
        ci, pi, ppi, hi = idx["CHROM"], idx["GENPOS"], idx["P"], idx["CHISQ"]
        while True:
            c, po, p, ch = [], [], [], []
            for _ in range(chunksize):
                line = fh.readline()
                if not line:
                    break
                f = line.split()
                try:
                    c.append(f[ci]); po.append(int(f[pi]))
                    p.append(float(f[ppi])); ch.append(float(f[hi]))
                except (ValueError, IndexError):
                    continue
            if not c:
                break
            n += len(c)
            chroms.append(pd.Series(c)); poss.append(pd.Series(po, dtype="int64"))
            ps.append(pd.Series(p, dtype="float64")); chis.append(pd.Series(ch, dtype="float64"))
    df = pd.DataFrame({
        "CHROM": pd.concat(chroms, ignore_index=True),
        "GENPOS": pd.concat(poss, ignore_index=True),
        "P": pd.concat(ps, ignore_index=True),
        "CHISQ": pd.concat(chis, ignore_index=True),
    })
    df.attrs["n_raw"] = n
    return df


def top_loci_table(path: str, threshold: float = GW_SIG, keep: int = 500) -> pd.DataFrame:
    """Second light pass: keep the rsIDs / alleles of the top associations."""
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().split()
        for line in fh:
            f = line.split()
            try:
                if float(f[header.index("P")]) < threshold:
                    rows.append(f)
            except (ValueError, IndexError):
                continue
    out = pd.DataFrame(rows, columns=header)
    out["P"] = out["P"].astype(float)
    out["BETA"] = out["BETA"].astype(float)
    out["SE"] = out["SE"].astype(float)
    out["A1FREQ"] = out["A1FREQ"].astype(float)
    out["N"] = out["N"].astype(int)
    out = out.sort_values("P").head(keep).reset_index(drop=True)
    return out


# ----------------------------------------------------------------------------
# F2  gene windows (hg19)
# ----------------------------------------------------------------------------
def load_refgene(path: str) -> pd.DataFrame:
    cols = ["bin", "name", "chrom", "strand", "txStart", "txEnd",
            "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
            "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols,
                     dtype={"chrom": str, "name2": str}, low_memory=False)
    df = df[~df["chrom"].str.contains("_")]          # drop haplotype/random contigs
    df = df[df["name2"].notna() & (df["name2"] != "")]
    df["chrom"] = df["chrom"].str.replace(r"^chr", "", regex=True)
    g = (df.groupby(["chrom", "name2"], as_index=False)
           .agg(start=("txStart", "min"), end=("txEnd", "max")))
    g["chrom_code"] = pd.to_numeric(g["chrom"], errors="coerce")
    g = g[g["chrom_code"].notna()].copy()            # autosomes only (1..22)
    g["chrom_code"] = g["chrom_code"].astype(int)
    # QC: drop implausible gene spans. UCSC refGene contains a few records whose
    # txStart/txEnd span tens of Mb (annotation artefacts); those windows would
    # scoop up hundreds of thousands of SNPs and dominate any gene-level score.
    span = g["end"] - g["start"]
    n_drop = int((span > MAX_GENE_SPAN).sum())
    if n_drop:
        print(f"     dropped {n_drop} gene records with span > {MAX_GENE_SPAN/1e6:.1f} Mb "
              f"(annotation artefacts)")
    g = g[span <= MAX_GENE_SPAN]
    return g[["chrom_code", "name2", "start", "end"]].reset_index(drop=True)


def gene_level_scores(snps: pd.DataFrame, genes: pd.DataFrame, pad: int = PAD):
    """min GWAS p (and SNP count) inside gene body +/- pad, per gene."""
    snps = snps.sort_values(["chrom_code", "GENPOS"]).reset_index(drop=True)
    recs = np.full(len(genes), np.nan)
    counts = np.zeros(len(genes), dtype=int)
    for cc, sub in snps.groupby("chrom_code", sort=False):
        pos = sub["GENPOS"].to_numpy()
        pv = sub["P"].to_numpy()
        cump = np.minimum.accumulate(pv)          # NOT used directly (needs window)
        # prefix-min is wrong for arbitrary windows; use sparse approach below
        gsel = genes.index[genes["chrom_code"] == cc].to_numpy()
        if len(gsel) == 0:
            continue
        gs = genes.loc[gsel, "start"].to_numpy() - pad
        ge = genes.loc[gsel, "end"].to_numpy() + pad
        order = np.argsort(gs)
        gs_s, ge_s, gsel_s = gs[order], ge[order], gsel[order]
        # sweep: for each gene window take min p over SNPs inside it
        lo = np.searchsorted(pos, gs_s, side="left")
        hi = np.searchsorted(pos, ge_s, side="right")
        for k in range(len(gsel_s)):
            a, b = lo[k], hi[k]
            counts[gsel_s[k]] = b - a
            if b > a:
                recs[gsel_s[k]] = pv[a:b].min()
    out = pd.DataFrame({"gene": genes["name2"], "chrom": genes["chrom_code"],
                        "n_snps": counts, "min_p": recs})
    # collapse duplicate symbols (paralogous/aliased rows) -> keep strongest
    out = out.sort_values("min_p", na_position="last")
    out = out.groupby("gene", as_index=False).first()
    out["score"] = -np.log10(out["min_p"].fillna(1.0))
    return out


# ----------------------------------------------------------------------------
# F3/F4  stratified competitive gene-set test
# ----------------------------------------------------------------------------
def stratified_set_test(gene_scores: pd.DataFrame, gene_set: set,
                        n_perm: int = 5000, n_strata: int = 10, seed: int = 20260920):
    """Competitive test of `gene_set` against background genes matched on
    SNP-count decile (removes gene-length / marker-density confounding)."""
    gs = gene_scores.copy()
    gs = gs[gs["n_snps"] > 0].reset_index(drop=True)
    try:
        gs["stratum"] = pd.qcut(gs["n_snps"], n_strata, labels=False, duplicates="drop")
    except ValueError:
        gs["stratum"] = 0
    present = sorted(set(gs["gene"]) & gene_set)
    if not present:
        return {"n_set_genes_present": 0, "error": "no set genes present in gene table"}
    obs = gs.loc[gs["gene"].isin(present), "score"].mean()
    back = gs[~gs["gene"].isin(present)].reset_index(drop=True)
    by_stratum = back.groupby("stratum")["score"].apply(lambda s: s.to_numpy())
    set_strata = gs.loc[gs["gene"].isin(present), "stratum"].to_numpy()
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        draws = []
        for st in set_strata:
            arr = by_stratum.get(st)
            if arr is None or len(arr) == 0:
                draws.append(0.0)
            else:
                draws.append(arr[rng.integers(0, len(arr))])
        null[i] = np.mean(draws)
    p_one = (np.sum(null >= obs) + 1) / (n_perm + 1)
    # naive (unstratified) Mann-Whitney, reported only for transparency
    from scipy import stats as _st
    in_set = gs.loc[gs["gene"].isin(present), "score"].to_numpy()
    out_set = back["score"].to_numpy()
    try:
        mw = _st.mannwhitneyu(in_set, out_set, alternative="greater")
        mw_p = float(mw.pvalue)
    except Exception:
        mw_p = float("nan")
    return {
        "n_set_genes_present": len(present),
        "mean_score_set": float(obs),
        "mean_score_background": float(out_set.mean()),
        "n_perm": n_perm,
        "perm_p_one_sided": float(p_one),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
        "z_vs_null": float((obs - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) > 0 else float("nan"),
        "naive_mannwhitney_p": mw_p,
        "set_genes": present,
    }


# ----------------------------------------------------------------------------
# F5  suggestive non-chr19 loci (hypothesis generating only)
# ----------------------------------------------------------------------------
def suggestive_loci(path: str, genes: pd.DataFrame, pad: int = 250_000,
                    lo: float = GW_SIG, hi: float = SUGGESTIVE):
    """Independent loci with 5e-8 <= p < 1e-5 outside chromosome 19.

    lo/hi are NUMERIC bounds: lo = 5e-8 (smaller), hi = 1e-5 (larger).
    An earlier version had these swapped, which made the condition `1e-5 <= p <
    5e-8` unsatisfiable and silently returned zero loci -- fixed 2026-09-20.
    """
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().split()
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.split()
            try:
                if f[ix["CHROM"]] == "19":
                    continue
                p = float(f[ix["P"]])
            except (ValueError, IndexError):
                continue
            if lo <= p < hi:
                rows.append((f[ix["CHROM"]], int(f[ix["GENPOS"]]), f[ix["SNP"]],
                             float(f[ix["BETA"]]), float(f[ix["SE"]]), p))
    if not rows:
        return pd.DataFrame(columns=["chrom", "pos", "snp", "beta", "se", "p", "nearest_gene"])
    d = pd.DataFrame(rows, columns=["chrom", "pos", "snp", "beta", "se", "p"])
    d = d.sort_values("p").reset_index(drop=True)
    # greedy distance clumping: keep the strongest SNP within +/- pad
    keep = []
    for _, r in d.iterrows():
        if all(abs(int(r["chrom"]) != int(k["chrom"])) or abs(int(r["pos"]) - int(k["pos"])) > pad
               for k in keep):
            keep.append(r)
    d = pd.DataFrame(keep).reset_index(drop=True)
    # nearest gene
    gn = []
    for _, r in d.iterrows():
        cc = pd.to_numeric(r["chrom"], errors="coerce")
        sub = genes[genes["chrom_code"] == cc]
        if sub.empty:
            gn.append("NA"); continue
        mid = (sub["start"] + sub["end"]) / 2.0
        j = int((np.abs(mid.to_numpy() - int(r["pos"]))).argmin())
        gn.append(sub["name2"].iloc[j])
    d["nearest_gene"] = gn
    return d


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--immune", required=True)
    ap.add_argument("--outdir", default="results/layerF")
    ap.add_argument("--n-perm", type=int, default=5000)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("[F1] loading GWAS summary statistics ...")
    snps = load_gwas(args.gwas)
    snps["chrom_code"] = pd.to_numeric(snps["CHROM"], errors="coerce")
    n_var = len(snps)
    lam = float(np.median(snps["CHISQ"]) / 0.4549)
    n_sig = int((snps["P"] < GW_SIG).sum())
    chrom_sig = (snps.loc[snps["P"] < GW_SIG, "CHROM"].value_counts().to_dict())
    print(f"     variants={n_var}  lambda={lam:.4f}  p<5e-8={n_sig}  chroms={chrom_sig}")

    print("[F1] top loci table ...")
    top = top_loci_table(args.gwas)
    top.to_csv(outdir / "layerF_gwas_top_loci.csv", index=False)

    print("[F2] gene windows ...")
    genes = load_refgene(args.genes)
    print(f"     autosomal genes: {len(genes)}")
    gscore = gene_level_scores(snps[snps["chrom_code"].notna()].copy(), genes)
    gscore.to_csv(outdir / "layerF_gene_level_scores.csv", index=False)
    print(f"     gene scores computed: {len(gscore)}  (min p = {gscore['min_p'].min():.3e})")

    immune = {l.strip() for l in open(args.immune, encoding="utf-8") if l.strip()}

    print("[F3] stratified gene-set test (all autosomes) ...")
    full = stratified_set_test(gscore, immune, n_perm=args.n_perm)

    gs_no19 = gscore[gscore["chrom"] != 19].reset_index(drop=True)
    print("[F4] stratified gene-set test (chr19 / APOE excluded) ...")
    no19 = stratified_set_test(gs_no19, immune, n_perm=args.n_perm)

    print("[F5] suggestive non-chr19 loci ...")
    sug = suggestive_loci(args.gwas, genes)
    sug.to_csv(outdir / "layerF_suggestive_non_APOE_loci.csv", index=False)
    print(f"     suggestive independent non-chr19 loci: {len(sug)}")

    summary = {
        "gwas_file": os.path.basename(args.gwas),
        "n_variants": n_var,
        "genomic_inflation_lambda": lam,
        "n_genome_wide_significant": n_sig,
        "chromosome_distribution_of_hits": {str(k): int(v) for k, v in chrom_sig.items()},
        "lead_snp": {
            "rsid": str(top.iloc[0]["SNP"]), "chrom": str(top.iloc[0]["CHROM"]),
            "pos": int(top.iloc[0]["GENPOS"]), "beta": float(top.iloc[0]["BETA"]),
            "se": float(top.iloc[0]["SE"]), "p": float(top.iloc[0]["P"]),
            "allele0": str(top.iloc[0]["ALLELE0"]), "allele1": str(top.iloc[0]["ALLELE1"]),
            "a1freq": float(top.iloc[0]["A1FREQ"]), "n": int(top.iloc[0]["N"]),
        },
        "lead_snp_odds_ratio_per_allele0": float(np.exp(-float(top.iloc[0]["BETA"]))),
        "n_genes_scored": int(len(gscore)),
        "immune_gene_set_test_all_autosomes": full,
        "immune_gene_set_test_chr19_excluded": no19,
        "n_suggestive_non_chr19_loci": int(len(sug)),
        "top15_genes_by_score": gscore.sort_values("score", ascending=False)
            .head(15)[["gene", "chrom", "n_snps", "min_p", "score"]]
            .to_dict(orient="records"),
        "note": ("Lead-SNP BETA is for ALLELE1 as coded in the released file; "
                 "OR for the ALLELE0 allele is exp(-BETA). rs429358 ALLELE0=C is the "
                 "APOE epsilon4-defining allele."),
    }
    with open(outdir / "layerF_real_gwas_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ("immune_gene_set_test_all_autosomes",
                                   "immune_gene_set_test_chr19_excluded",
                                   "top15_genes_by_score")}, indent=2, ensure_ascii=False))
    print("\n[set test / all autosomes]", json.dumps(
        {k: v for k, v in full.items() if k != "set_genes"}, ensure_ascii=False))
    print("[set test / chr19 excluded]", json.dumps(
        {k: v for k, v in no19.items() if k != "set_genes"}, ensure_ascii=False))
    print("\nTOP15 GENES:")
    print(gscore.sort_values("score", ascending=False)
          .head(15)[["gene", "chrom", "n_snps", "min_p", "score"]].to_string(index=False))


if __name__ == "__main__":
    main()
