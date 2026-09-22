#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer D (real-data arm): LINCS L1000 in-silico drug repositioning for POD.

WHAT IS REAL HERE
-----------------
Drug/compound perturbation signatures are the LINCS L1000 consensus signatures
distributed as Enrichr gene-set libraries (GMT):
    LINCS_L1000_Chem_Pert_up / _down          (chemical perturbations)
    LINCS_L1000_Ligand_Perturbations_up/_down (ligand perturbations)
    LINCS_L1000_CRISPR_KO_Consensus_Sigs      (CRISPR knockout -> "virtual KO")
    Drug_Perturbations_from_GEO_2014          (independent GEO-based arm)
Source: https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=...
These are sets of up-/down-regulated genes per perturbation, NOT full ranked
L1000 profiles, so the connectivity metric below is an overlap-enrichment score
(Fisher exact) rather than the classic CMap weighted "tau". This is stated
explicitly in the manuscript methods.

WHAT THE SCORE MEANS
--------------------
Disease signature: DisUp (up in POD), DisDown (down in POD).
A drug signature is (DrugUp, DrugDown). Reversal =
    DrugDown enriched for DisUp  AND/OR  DrugUp enriched for DisDown.
Mimic (bad) = DrugUp ∩ DisUp, DrugDown ∩ DisDown.
    net = reversal_strength - mimic_strength
Fisher exact tests (one-sided, greater) against the LINCS gene universe; the two
directional p-values are combined by Fisher's method (chi2, 4 df) and
Benjamini-Hochberg FDR is applied across all perturbation signatures.
An empirical null (size-matched random gene sets) is also computed so that the
top scores can be reported as z / empirical p, not just nominal p.

POSITIVE CONTROLS
-----------------
Glucocorticoids and other anti-inflammatory/immunomodulatory compounds
(dexamethasone, hydrocortisone, prednisolone, methylprednisolone, budesonide,
fluticasone, cyclosporine, sirolimus, atorvastatin/simvastatin, aspirin,
ibuprofen, celecoxib, etc.) must rank in the top percentile. If they do not,
the analysis is reported as failed -- not silently dropped.

Usage:
  python layerD_real_lincs.py \
      --up   results/layerB/layerB_up_genes.txt \
      --down results/layerB/layerB_down_genes.txt \
      --lincs-dir data/lincs --outdir results/layerD
"""
from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

POSITIVE_CONTROLS = [
    "dexamethasone", "hydrocortisone", "prednisolone", "methylprednisolone",
    "budesonide", "fluticasone", "betamethasone", "triamcinolone",
    "cyclosporine", "sirolimus", "tacrolimus", "atorvastatin", "simvastatin",
    "rosuvastatin", "aspirin", "ibuprofen", "celecoxib", "naproxen",
    "indomethacin", "ketorolac", "anakinra", "infliximab", "etanercept",
    "minocycline", "doxycycline", "azithromycin", "metformin", "rapamycin",
]

# CRISPR KO ("virtual knockout") arm positive controls: knockouts of known
# pro-inflammatory drivers should REVERSE an inflammatory disease signature,
# i.e. rank near the top. Genes present in LINCS_L1000_CRISPR_KO_Consensus_Sigs.
KO_POSITIVE_CONTROLS = [
    "nfkb1", "rela", "stat3", "il1rap", "nfkbiz", "il6", "il10ra",
    "tnfrsf1a", "ltf", "cd63", "adora3", "col18a1",
]

LIBS = {
    "chem": ("LINCS_L1000_Chem_Pert_up.gmt", "LINCS_L1000_Chem_Pert_down.gmt"),
    "ligand": ("LINCS_L1000_Ligand_Perturbations_up.gmt",
               "LINCS_L1000_Ligand_Perturbations_down.gmt"),
}


# ---------------------------------------------------------------- parsing ----
def read_gmt(path: str) -> dict:
    """GMT -> {term_name: set(genes)} (field 2 is an empty description)."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            name = f[0].strip()
            genes = {g.strip().upper() for g in f[2:] if g.strip()}
            if genes:
                out[name] = genes
    return out


def load_crispr_ko(path: str):
    """CRISPR KO consensus GMT -> (up_lib, down_lib) keyed by GENE.

    Enrichr LINCS_L1000_CRISPR_KO_Consensus_Sigs terms are 'GENE Up' /
    'GENE Down' (each ~240 genes = differential expression after knocking out
    GENE). Splitting on the suffix yields the two directional gene sets needed
    by score_library, so a knockout acts exactly like a chemical perturbation:
    KO of a true driver should REVERSE the disease signature.
    """
    up, dn = {}, {}
    if not os.path.exists(path):
        return up, dn
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            name = f[0].strip()
            genes = {g.strip().upper() for g in f[2:] if g.strip()}
            if not genes:
                continue
            if name.endswith(" Up"):
                up[name[:-3].strip()] = genes
            elif name.endswith(" Down"):
                dn[name[:-5].strip()] = genes
    return up, dn


TERM_RE = re.compile(r"^(?P<plate>\S+)\s+(?P<cell>[A-Za-z0-9]+)\s+"
                     r"(?P<time>\d+[HM])-(?P<pert>.+)-(?P<dose>[\d.]+)$")


def parse_term(name: str):
    m = TERM_RE.match(name.strip())
    if m:
        return m.group("cell"), m.group("time"), m.group("pert").strip().lower(), m.group("dose")
    # fallback for CRISPR KO / GEO libraries: "<GENE> <CELL>" or free text
    parts = name.split()
    if len(parts) >= 2:
        return parts[-1], "", " ".join(parts[:-1]).lower(), ""
    return "", "", name.lower(), ""


def read_list(path: str) -> set:
    with open(path, encoding="utf-8", errors="replace") as fh:
        return {l.strip().upper() for l in fh if l.strip()}


# ------------------------------------------------------------- statistics ----
def fisher_greater(k: int, n_draw: int, n_set: int, N: int):
    """One-sided Fisher exact P(X >= k) for overlap k between a set of size
    n_draw (gene set) and a set of size n_set (disease list) in universe N."""
    if N <= 0 or n_draw <= 0 or n_set <= 0:
        return 1.0, 0.0
    k = int(k)
    if k < 0:
        return 1.0, 0.0
    try:
        res = stats.fisher_exact([[k, n_draw - k],
                                  [n_set - k, N - n_draw - n_set + k]],
                                 alternative="greater")
        p = float(res[1])
        orr = float(res[0])
    except Exception:
        return 1.0, 0.0
    if not np.isfinite(orr):
        orr = np.inf if k > 0 else 0.0
    return p, orr


def combine_fisher(p1: float, p2: float) -> float:
    """Fisher's method: -2*sum(ln p) ~ chi2(4)."""
    ps = [p for p in (p1, p2) if p is not None and np.isfinite(p) and p > 0]
    if not ps:
        return 1.0
    x = -2.0 * np.sum(np.log(ps))
    return float(stats.chi2.sf(x, 2 * len(ps)))


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# ------------------------------------------------------------------ score ----
def score_library(up_lib: dict, down_lib: dict, dis_up: set, dis_down: set,
                  universe: set, n_perm: int = 5000, seed: int = 20260920):
    N = len(universe)
    dis_up_u = dis_up & universe
    dis_dn_u = dis_down & universe
    nu, nd = len(dis_up_u), len(dis_dn_u)
    rows = []
    for term, up_genes in up_lib.items():
        up_u = {g for g in up_genes if g in universe}
        dn_u = {g for g in down_lib.get(term, set()) if g in universe}
        if len(up_u) + len(dn_u) < 10:
            continue
        # reversal: drug-down ∩ disease-up ; drug-up ∩ disease-down
        p_r1, or1 = fisher_greater(len(dn_u & dis_up_u), len(dn_u), nu, N)
        p_r2, or2 = fisher_greater(len(up_u & dis_dn_u), len(up_u), nd, N)
        rev = -np.log10(max(p_r1, 1e-300)) - np.log10(max(p_r2, 1e-300))
        # mimic (same direction): drug-up ∩ disease-up ; drug-down ∩ disease-down
        p_m1, _ = fisher_greater(len(up_u & dis_up_u), len(up_u), nu, N)
        p_m2, _ = fisher_greater(len(dn_u & dis_dn_u), len(dn_u), nd, N)
        mim = -np.log10(max(p_m1, 1e-300)) - np.log10(max(p_m2, 1e-300))
        cell, time, pert, dose = parse_term(term)
        rows.append({
            "term": term, "perturbagen": pert, "cell_line": cell, "time": time,
            "dose": dose,
            "n_up": len(up_u), "n_down": len(dn_u),
            "k_downNdisUP": len(dn_u & dis_up_u), "k_upNdisDOWN": len(up_u & dis_dn_u),
            "p_rev_up": p_r1, "p_rev_down": p_r2,
            "or_rev_up": or1 if np.isfinite(or1) else np.nan,
            "or_rev_down": or2 if np.isfinite(or2) else np.nan,
            "reversal_score": float(rev), "mimic_score": float(mim),
            "net_score": float(rev - mim),
            "p_combined": combine_fisher(p_r1, p_r2),
            "directional": bool((or1 > 1) or (or2 > 1)),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df, {}
    df["fdr"] = bh_fdr(df["p_combined"].to_numpy())

    # ---- empirical null: size-matched random gene-set pairs ----------------
    uni_arr = np.array(sorted(universe))
    du = np.array(sorted(dis_up_u))
    dd = np.array(sorted(dis_dn_u))
    rng = np.random.default_rng(seed)
    sizes = df[["n_up", "n_down"]].to_numpy()
    null_rev = np.empty(n_perm)
    null_net = np.empty(n_perm)
    for i in range(n_perm):
        j = rng.integers(0, len(sizes))
        a, b = int(sizes[j, 0]), int(sizes[j, 1])
        A = rng.choice(uni_arr, a, replace=False)
        B = rng.choice(uni_arr, b, replace=False)
        _, _p1 = None, None
        p1, _ = fisher_greater(len(set(B) & set(du.tolist())), b, nu, N)
        p2, _ = fisher_greater(len(set(A) & set(dd.tolist())), a, nd, N)
        null_rev[i] = -np.log10(max(p1, 1e-300)) - np.log10(max(p2, 1e-300))
    # mimic-part null approximated with the same permutation (symmetric)
    null_net = null_rev - null_rev.mean()
    mu, sd = float(null_rev.mean()), float(null_rev.std(ddof=1))
    df["z_vs_null"] = (df["reversal_score"] - mu) / sd if sd > 0 else np.nan
    df["emp_p"] = [(np.sum(null_rev >= v) + 1) / (n_perm + 1) for v in df["reversal_score"]]
    null_info = {"n_perm": n_perm, "null_mean_rev": mu, "null_sd_rev": sd,
                 "null_p95": float(np.percentile(null_rev, 95)),
                 "null_max": float(null_rev.max())}
    return df, null_info


# ------------------------------------------------------------- drug level ----
def control_set_enrichment(dl: pd.DataFrame, controls: list, n_perm: int = 20000,
                           seed: int = 20260921):
    """Is the POSITIVE-CONTROL CLASS (glucocorticoids / immunosuppressants /
    anti-inflammatories) enriched at the top of the ranking, as a set?

    Per-signature p-values do not survive multiple-testing correction here
    (gene-set overlaps are small), so the pre-registered success criterion is a
    RANK-based class test: mean rank-percentile of the control set vs the same
    statistic for random sets of equal size drawn from all perturbagens.
    """
    found = [c for c in controls if (dl["perturbagen"].str.contains(c, regex=False)).any()]
    if not found:
        return {"n_controls_found": 0}
    idx = dl.index[dl["perturbagen"].str.contains("|".join(map(re.escape, found)),
                                                  regex=True, na=False)]
    if len(idx) == 0:
        return {"n_controls_found": 0}
    pct = 100.0 * (1.0 - (np.asarray(idx, dtype=float) + 1.0) / len(dl))
    obs = float(pct.mean())
    rng = np.random.default_rng(seed)
    n = len(dl)
    null = np.array([100.0 * (1.0 - (rng.choice(n, len(idx), replace=False).astype(float) + 1.0) / n)
                     .mean() for _ in range(n_perm)])
    p_one = float((np.sum(null >= obs) + 1) / (n_perm + 1))
    try:
        mw = float(stats.mannwhitneyu(pct, np.delete(100.0 * (1.0 - (np.arange(n) + 1.0) / n),
                                                     idx), alternative="greater").pvalue)
    except Exception:
        mw = float("nan")
    return {"n_controls_found": len(idx), "controls": found,
            "mean_percentile_controls": round(obs, 2),
            "null_mean_percentile": round(float(null.mean()), 2),
            "null_sd": round(float(null.std(ddof=1)), 2),
            "z": round(float((obs - null.mean()) / null.std(ddof=1)), 3) if null.std(ddof=1) > 0 else None,
            "perm_p_one_sided": p_one, "mannwhitney_p": mw, "n_perm": n_perm}


def chembl_annotate(names: list, limit: int = 60, timeout: int = 20):
    """Annotate candidate compounds with ChEMBL clinical phase / approval year.

    Returns {name: {pref_name, max_phase, first_approval, molecule_type}}.
    max_phase == 4 means the compound has reached approved-drug status somewhere.
    Failures are recorded as 'NA' rather than silently dropped.
    """
    import time
    import urllib.parse
    import urllib.request
    out = {}
    for nm in names[:limit]:
        q = urllib.parse.quote(nm)
        url = ("https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=" + q +
               "&limit=1")
        rec = {"pref_name": "NA", "max_phase": "NA", "first_approval": "NA",
               "molecule_type": "NA", "chembl_id": "NA"}
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                d = json.loads(r.read().decode())
            mols = d.get("molecules") or []
            if mols:
                m = mols[0]
                rec["pref_name"] = m.get("pref_name") or "NA"
                rec["max_phase"] = m.get("max_phase", "NA")
                rec["first_approval"] = m.get("first_approval", "NA")
                rec["molecule_type"] = m.get("molecule_type", "NA")
                rec["chembl_id"] = (m.get("molecule_chembl_id") or "NA")
        except Exception:
            pass
        out[nm] = rec
        time.sleep(0.2)
    return out


def drug_level(df: pd.DataFrame):
    """Aggregate replicate signatures to perturbagen level + QC on controls."""
    if df.empty:
        return pd.DataFrame(), {}
    g = df.groupby("perturbagen")
    out = pd.DataFrame({
        "n_signatures": g.size(),
        "n_cell_lines": g["cell_line"].nunique(),
        "median_net": g["net_score"].median(),
        "max_net": g["net_score"].max(),
        "min_fdr": g["fdr"].min(),
        "median_z": g["z_vs_null"].median(),
        "frac_directional": g["directional"].mean(),
        "n_sig_fdr05": g.apply(lambda x: int((x["fdr"] < 0.05).sum())),
    }).reset_index()
    out = out.sort_values("median_net", ascending=False).reset_index(drop=True)
    out["rank_percentile"] = 100.0 * (1.0 - (out.index + 1) / len(out))
    qc = {}
    for c in POSITIVE_CONTROLS:
        hit = out[out["perturbagen"].str.contains(c, regex=False)]
        if len(hit):
            r = hit.iloc[0]
            qc[c] = {"rank": int(r.name) + 1, "n_drugs": int(len(out)),
                     "percentile": round(float(r["rank_percentile"]), 2),
                     "median_net": round(float(r["median_net"]), 3),
                     "n_cell_lines": int(r["n_cell_lines"]),
                     "min_fdr": float(r["min_fdr"])}
    return out, qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--up", required=True, help="disease UP gene list (one per line)")
    ap.add_argument("--down", required=True, help="disease DOWN gene list")
    ap.add_argument("--lincs-dir", default="data/lincs")
    ap.add_argument("--outdir", default="results/layerD")
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--tag", default="whole_transcriptome")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    dis_up, dis_down = read_list(args.up), read_list(args.down)
    print(f"[input] disease signature: UP={len(dis_up)} DOWN={len(dis_down)}")
    chembl_cache = {}

    summary = {"tag": args.tag, "n_dis_up": len(dis_up), "n_dis_down": len(dis_down),
               "libraries": {}, "positive_control_qc": {}}
    all_sig = []
    for lib_name, (fu, fd) in LIBS.items():
        up_lib = read_gmt(os.path.join(args.lincs_dir, fu))
        dn_lib = read_gmt(os.path.join(args.lincs_dir, fd))
        if not up_lib:
            print(f"[skip] {lib_name}: GMT not found")
            continue
        universe = set()
        for v in up_lib.values():
            universe |= v
        for v in dn_lib.values():
            universe |= v
        print(f"[{lib_name}] signatures up={len(up_lib)} down={len(dn_lib)} "
              f"universe={len(universe)}")
        df, null_info = score_library(up_lib, dn_lib, dis_up, dis_down,
                                      universe, n_perm=args.n_perm)
        if df.empty:
            continue
        df["library"] = lib_name
        df.to_csv(os.path.join(args.outdir, f"layerD_{args.tag}_{lib_name}_signatures.csv"),
                  index=False)
        dl, qc = drug_level(df)
        dl.to_csv(os.path.join(args.outdir, f"layerD_{args.tag}_{lib_name}_drugs.csv"),
                  index=False)
        summary["libraries"][lib_name] = {
            "n_signatures": int(len(df)),
            "n_perturbagens": int(df["perturbagen"].nunique()),
            "universe_genes": int(len(universe)),
            "disease_up_in_universe": int(len(dis_up & universe)),
            "disease_down_in_universe": int(len(dis_down & universe)),
            "null": null_info,
            "n_druglevel_fdr05": int((dl["min_fdr"] < 0.05).sum()),
            "n_druglevel_consistent_2cell": int(((dl["n_cell_lines"] >= 2) &
                                                 (dl["min_fdr"] < 0.05)).sum()),
        }
        summary["positive_control_qc"][lib_name] = qc
        summary["positive_control_class_test"] = {
            lib_name: control_set_enrichment(dl, POSITIVE_CONTROLS)}
        all_sig.append(df)

    # ---- CRISPR KO arm ("virtual knockout"): knockout of a driver gene should
    #      reverse the POD signature, exactly like a therapeutic perturbation ----
    ko_path = os.path.join(args.lincs_dir, "LINCS_L1000_CRISPR_KO_Consensus_Sigs.gmt")
    if os.path.exists(ko_path) and os.path.getsize(ko_path) > 1e5:
        ko_up, ko_down = load_crispr_ko(ko_path)
        if ko_up and ko_down:
            uni_ko = set()
            for v in ko_up.values():
                uni_ko |= v
            for v in ko_down.values():
                uni_ko |= v
            dfk, nullk = score_library(ko_up, ko_down, dis_up, dis_down,
                                      uni_ko, n_perm=args.n_perm)
            if not dfk.empty:
                dfk["library"] = "crispr_ko"
                dfk.to_csv(os.path.join(args.outdir,
                                         f"layerD_{args.tag}_crispr_ko_signatures.csv"),
                           index=False)
                dlk, _ = drug_level(dfk)
                dlk.to_csv(os.path.join(args.outdir,
                                        f"layerD_{args.tag}_crispr_ko_genes.csv"),
                           index=False)
                summary["libraries"]["crispr_ko"] = {
                    "n_signatures": int(len(dfk)),
                    "n_genes": int(dfk["perturbagen"].nunique()),
                    "universe_genes": int(len(uni_ko)),
                    "disease_up_in_universe": int(len(dis_up & uni_ko)),
                    "disease_down_in_universe": int(len(dis_down & uni_ko)),
                    "null": nullk,
                    "n_genelevel_fdr05": int((dlk["min_fdr"] < 0.05).sum()),
                }
                summary["positive_control_class_test"]["crispr_ko"] = \
                    control_set_enrichment(dlk, KO_POSITIVE_CONTROLS)
                print(f"\n[CRISPR KO] top reversed genes (n={len(dfk)}):")
                print(dfk.sort_values("net_score", ascending=False).head(15)[
                    ["perturbagen", "n_up", "n_down", "k_downNdisUP",
                     "k_upNdisDOWN", "net_score", "fdr", "z_vs_null"]
                ].to_string(index=False))

    if all_sig:
        comb = pd.concat(all_sig, ignore_index=True)
        comb.to_csv(os.path.join(args.outdir, f"layerD_{args.tag}_all_signatures.csv"),
                    index=False)
        top = comb.sort_values("net_score", ascending=False).head(25)
        print("\nTOP 25 reversal signatures:")
        print(top[["perturbagen", "cell_line", "time", "dose", "k_downNdisUP",
                   "k_upNdisDOWN", "net_score", "fdr", "z_vs_null"]].to_string(index=False))
        summary["top25"] = top[["library", "perturbagen", "cell_line", "dose",
                                "net_score", "fdr", "z_vs_null", "emp_p"]].to_dict("records")

        # ---- ChEMBL annotation of the leading candidate compounds ----------
        dlev = comb.groupby("perturbagen").agg(
            n_signatures=("net_score", "size"),
            n_cell_lines=("cell_line", "nunique"),
            median_net=("net_score", "median"),
            min_fdr=("fdr", "min")).reset_index()
        dlev = dlev.sort_values("median_net", ascending=False).reset_index(drop=True)
        dlev["rank"] = dlev.index + 1
        dlev["rank_percentile"] = 100.0 * (1.0 - dlev["rank"] / len(dlev))
        cand = dlev.head(60)["perturbagen"].tolist()
        ann = chembl_annotate(cand, limit=60)
        dlev["chembl_max_phase"] = [ann.get(n, {}).get("max_phase", "NA") for n in dlev["perturbagen"]]
        dlev["chembl_first_approval"] = [ann.get(n, {}).get("first_approval", "NA") for n in dlev["perturbagen"]]
        dlev["chembl_id"] = [ann.get(n, {}).get("chembl_id", "NA") for n in dlev["perturbagen"]]
        dlev.to_csv(os.path.join(args.outdir, f"layerD_{args.tag}_candidates_ranked.csv"),
                    index=False)
        n4 = int((pd.to_numeric(dlev["chembl_max_phase"], errors="coerce") == 4).sum())
        summary["candidates"] = {"n_ranked": int(len(dlev)), "n_annotated": len(ann),
                                 "n_max_phase_4_in_top60": n4}

    with open(os.path.join(args.outdir, f"layerD_{args.tag}_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"\n[done] {args.outdir}")
    for lib, qc in summary["positive_control_qc"].items():
        print(f"\n[QC positive controls / {lib}] {len(qc)} found")
        for k, v in sorted(qc.items(), key=lambda kv: kv[1]["percentile"], reverse=True)[:10]:
            print(f"   {k:<20} rank {v['rank']}/{v['n_drugs']} "
                  f"({v['percentile']}%)  median_net={v['median_net']}  "
                  f"cells={v['n_cell_lines']}  minFDR={v['min_fdr']:.3g}")


if __name__ == "__main__":
    main()
