#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Layer E - SUPPLEMENTARY molecular docking / ADMET triage.
==========================================================
NOT on the critical path. Per prior experience (CPSP virtual-screening project)
docking frequently fails (ion channels, no holo structure, no pocket, scoring
correlation ceiling rho~0.78). Failure here must NOT block the manuscript.

This script therefore does two things:
  (A) DOCKABILITY AUDIT  -- for each hub/anchor target, query UniProt + RCSB
      and report whether an experimental structure with a ligandable pocket
      exists. This is a REAL, reportable supplementary table (honest negative
      included: "target X has no experimental structure -> not dockable").
  (B) LIGAND TRIAGE      -- for each candidate drug from Layer D, pull PubChem
      identity/properties (MW, XLogP, TPSA, HBD/HBA, rotatable bonds),
      apply Lipinski/Veber/CNS-MPO-style rules and flag BBB permeability risk.
      Emits AutoDock Vina config + a runner script for the user's machine
      (Vina/obabel are not available in every environment).

Usage
-----
  # (A) audit only - works anywhere with internet
  python layerE_docking.py --audit --targets ADORA3,TMIGD3,COL13A1,SPATA13,COL18A1,CD63,LTF \
      --outdir results/layerE

  # (B) ligand triage from Layer D output
  python layerE_docking.py --candidates results/layerD/layerD_perturbagen_consensus.csv \
      --top 30 --outdir results/layerE

  # both
  python layerE_docking.py --audit --candidates ... --outdir results/layerE
"""
from __future__ import annotations

import argparse
import json
import os
import time
import urllib.parse
import urllib.request

import pandas as pd

UA = {"User-Agent": "pod-epigenetic-repositioning/1.0"}
TIMEOUT = 60


def _get(url: str, headers=None, data=None):
    req = urllib.request.Request(url, headers={**UA, **(headers or {})}, data=data)
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read().decode("utf-8", "replace")


# ---------------------------------------------------------------- (A) targets
def uniprot_accession(gene: str):
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote(f"gene:{gene} AND organism_id:9606 AND reviewed:true")
           + "&format=json&size=1&fields=accession,protein_name,id")
    for attempt in range(3):
        try:
            res = json.loads(_get(url))
            break
        except Exception as exc:
            if attempt == 2:
                print("[warn] uniprot", gene, exc)
                return None, None
            time.sleep(1.0 + attempt)
    try:
        if not res.get("results"):
            return None, None
        r = res["results"][0]
        name = r.get("proteinDescription", {}).get("recommendedName", {}) \
                .get("fullName", {}).get("value", "")
        return r.get("primaryAccession"), name
    except Exception as exc:
        print("[warn] uniprot", gene, exc)
        return None, None


def rcsb_structures(uniprot_acc: str, retries: int = 3):
    """Return list of PDB entries for a UniProt accession (experimental).

    NOTE: RCSB text search is only enabled on
    '...reference_sequence_identifiers.database_accession', NOT on
    '...uniprot_ids' (that one returns HTTP 400). Verified 2026-09-20.
    """
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": uniprot_acc}},
          "return_type": "entry",
          "request_options": {"paginate": {"start": 0, "rows": 25}}}
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(q))
    for attempt in range(retries):
        try:
            res = json.loads(_get(url))
            return [x.get("identifier") for x in res.get("result_set", [])]
        except Exception as exc:
            if attempt == retries - 1:
                print("[warn] rcsb", uniprot_acc, exc)
            time.sleep(1.0 + attempt)
    return []


def audit_targets(genes: list[str]) -> pd.DataFrame:
    rows = []
    for g in genes:
        acc, name = uniprot_accession(g)
        if not acc:
            rows.append({"gene": g, "uniprot": None, "protein": name,
                         "n_pdb_entries": 0, "dockable_experimental": False,
                         "note": "no reviewed human UniProt entry"})
            continue
        pdb = rcsb_structures(acc)
        rows.append({
            "gene": g, "uniprot": acc, "protein": name,
            "n_pdb_entries": len(pdb),
            "pdb_ids": ";".join(pdb[:10]),
            "dockable_experimental": len(pdb) > 0,
            "note": "" if pdb else "no experimental structure -> use AlphaFold model "
                                   "or skip docking (report as limitation)",
        })
        time.sleep(0.3)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- (B) ligands
PUBCHEM_PROPS = ["MolecularFormula", "MolecularWeight", "XLogP", "HBondDonorCount",
                 "HBondAcceptorCount", "TPSA", "RotatableBondCount", "CanonicalSMILES",
                 "IsomericSMILES"]


def pubchem_compound(name: str):
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    url = (base + urllib.parse.quote(name)
           + "/property/" + ",".join(PUBCHEM_PROPS) + "/JSON")
    try:
        res = json.loads(_get(url))
        return res["PropertyTable"]["Properties"][0]
    except Exception:
        return None


def triage(cands: pd.DataFrame, top: int, outdir: str):
    names = cands.sort_values("mean_score").head(top)["perturbagen"].astype(str).tolist()
    rows = []
    for n in names:
        info = pubchem_compound(n)
        if not info:
            rows.append({"candidate": n, "pubchem_resolved": False})
            continue
        mw = float(info.get("MolecularWeight", 0) or 0)
        logp = float(info.get("XLogP", 0) or 0)
        tpsa = float(info.get("TPSA", 0) or 0)
        hbd = int(info.get("HBondDonorCount", 0) or 0)
        hba = int(info.get("HBondAcceptorCount", 0) or 0)
        rot = int(info.get("RotatableBondCount", 0) or 0)
        lipinski_hits = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
        veber_ok = (tpsa <= 140) and (rot <= 10)
        # crude CNS/BBB heuristic (NOT a validated model - flag as such)
        bbb_likely = (mw <= 450) and (tpsa <= 90) and (logp >= 1 and logp <= 4) and hbd <= 3
        rows.append({"candidate": n, "pubchem_resolved": True,
                     "CID": info.get("CID"), "MW": mw, "XLogP": logp, "TPSA": tpsa,
                     "HBD": hbd, "HBA": hba, "RotB": rot,
                     "lipinski_violations": lipinski_hits, "veber_ok": veber_ok,
                     "bbb_heuristic_likely": bbb_likely,
                     "SMILES": info.get("IsomericSMILES") or info.get("CanonicalSMILES")})
        time.sleep(0.2)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "layerE_ligand_triage.csv"), index=False)

    # emit Vina runner (executed on a machine with vina + obabel)
    run = ["#!/usr/bin/env bash",
           "# AutoDock Vina runner generated by layerE_docking.py",
           "# Prereqs: vina in PATH; openbabel (obabel) in PATH; receptor PDBQT prepared.",
           "set -euo pipefail", ""]
    ok = df[df.get("pubchem_resolved", pd.Series(False)).astype(bool)] if len(df) else df
    for _, r in ok.iterrows():
        safe = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in str(r["candidate"]))
        run += [f"# --- {r['candidate']} ---",
                f"obabel -:'{r['SMILES']}' -O lig_{safe}.pdbqt --gen3d -p 7.4 2>/dev/null \\",
                f"  && vina --receptor receptor.pdbqt --ligand lig_{safe}.pdbqt \\",
                f"     --center_x CX --center_y CY --center_z CZ \\",
                f"     --size_x 22 --size_y 22 --size_z 22 --exhaustiveness 16 \\",
                f"     --seed 12345 --num_modes 9 --energy_range 3 \\",
                f"     --out out_{safe}.pdbqt --log log_{safe}.txt || echo '[skip] {safe}'",
                ""]
    with open(os.path.join(outdir, "run_vina.sh"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(run))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", action="store_true")
    ap.add_argument("--targets", default="ADORA3,TMIGD3,COL13A1,SPATA13,COL18A1,CD63,LTF")
    ap.add_argument("--candidates", help="layerD consensus CSV")
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--outdir", default="results/layerE")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = {}

    if args.audit:
        genes = [g.strip() for g in args.targets.split(",") if g.strip()]
        aud = audit_targets(genes)
        aud.to_csv(os.path.join(args.outdir, "layerE_dockability_audit.csv"), index=False)
        out["dockability_audit"] = aud.to_dict(orient="records")
        print(aud.to_string(index=False))

    if args.candidates:
        cands = pd.read_csv(args.candidates)
        tri = triage(cands, args.top, args.outdir)
        out["n_candidates_triaged"] = int(len(tri))
        out["n_resolved_in_pubchem"] = int(tri.get("pubchem_resolved", pd.Series(dtype=bool)).sum())
        print(tri.head(20).to_string(index=False))

    with open(os.path.join(args.outdir, "layerE_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
