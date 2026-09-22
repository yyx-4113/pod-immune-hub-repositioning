#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pre-specified positive-control SUBCLASS rank tests for Layer D.

Why this exists: the pooled control list (29 matched compounds) mixes strong
prior classes (glucocorticoids) with weak ones (NSAIDs, statins), so a single
class test is diluted. The subclasses below were fixed a priori (before looking
at the ranking) and are reported separately, honestly.
"""
import argparse
import re

import numpy as np
import pandas as pd

CLASSES = {
    "glucocorticoids": ["dexamethasone", "hydrocortisone", "prednisolone",
                        "methylprednisolone", "budesonide", "betamethasone",
                        "triamcinolone", "fluticasone", "cortisone", "prednisone",
                        "beclomethasone", "mometasone", "ciclesonide", "halcinonide",
                        "fluocinolone", "desonide", "alclometasone", "diflorasone",
                        "flurandrenolide", "flunisolide", "prednicarbate",
                        "clobetasol", "amcinonide"],
    "immunosuppressants": ["cyclosporine", "sirolimus", "tacrolimus", "everolimus",
                           "mycophenolate", "azathioprine", "methotrexate",
                           "rapamycin"],
    "nsaids": ["aspirin", "ibuprofen", "naproxen", "celecoxib", "indomethacin",
               "ketorolac", "diclofenac", "piroxicam", "sulindac", "flurbiprofen",
               "mefenamic"],
    "statins": ["atorvastatin", "simvastatin", "rosuvastatin", "lovastatin",
                "pravastatin", "fluvastatin", "mevastatin"],
}


def class_test(dl, names, n_perm=20000, seed=7):
    pat = "|".join(re.escape(x) for x in names)
    idx = dl.index[dl["perturbagen"].str.contains(pat, regex=True, na=False)]
    n = len(dl)
    if len(idx) == 0:
        return {"n": 0}
    pct = 100.0 * (1.0 - (np.asarray(idx, dtype=float) + 1.0) / n)
    obs = float(pct.mean())
    rng = np.random.default_rng(seed)
    null = np.array([100.0 * (1.0 - (rng.choice(n, len(idx), replace=False).astype(float)
                                     + 1.0) / n).mean() for _ in range(n_perm)])
    p = float((np.sum(null >= obs) + 1) / (n_perm + 1))
    z = float((obs - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) > 0 else None
    members = [(dl.loc[i, "perturbagen"], round(100 * (1 - (i + 1) / n), 1),
                round(float(dl.loc[i, "median_net"]), 3), int(dl.loc[i, "n_cell_lines"]))
               for i in idx]
    return {"n": len(idx), "mean_percentile": round(obs, 2),
            "null_mean": round(float(null.mean()), 2),
            "null_sd": round(float(null.std(ddof=1)), 2),
            "z": round(z, 3) if z is not None else None,
            "perm_p_one_sided": p, "n_perm": n_perm, "members": members}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs-csv", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    dl = pd.read_csv(args.drugs_csv)
    print(f"[file] {args.drugs_csv}  n_perturbagens={len(dl)}")
    rows = []
    for cls, names in CLASSES.items():
        r = class_test(dl, names)
        if not r.get("n"):
            print(f"  {cls:<20} : none found")
            continue
        print(f"  {cls:<20} n={r['n']:>2}  mean_pct={r['mean_percentile']:5.2f}  "
              f"null={r['null_mean']}±{r['null_sd']}  z={r['z']:+.3f}  "
              f"perm_p={r['perm_p_one_sided']:.4f}")
        for m in r["members"]:
            print(f"        - {m[0]:<38} {m[1]:>5}%  median_net={m[2]:>6}  cells={m[3]}")
        rows.append({"class": cls, **{k: v for k, v in r.items() if k != "members"},
                     "members": "; ".join(f"{a}({b}%)" for a, b, _, _ in r["members"])})
    df = pd.DataFrame(rows)
    out = args.out or args.drugs_csv.replace("_chem_drugs.csv", "_control_class_tests.csv")
    df.to_csv(out, index=False)
    print("[wrote]", out)


if __name__ == "__main__":
    main()
