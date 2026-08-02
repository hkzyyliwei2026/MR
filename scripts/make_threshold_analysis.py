# -*- coding: utf-8 -*-
"""
Instrument-threshold instability analysis.

Outputs:
  Table 1  cross-threshold comparison: where the traits reaching FDR significance under the genome-wide threshold (P<5e-8) sit in the primary analysis (P<1e-5)
  Table 2  instrument-count distribution under the genome-wide threshold
  threshold_instability.csv  trait-level comparison under both thresholds (used for plotting)

All inputs are read-only.
"""
import os
from pathlib import Path
import sys

import pandas as pd

# Input tables: local copy by default; an external copy may be supplied via MR_UPSTREAM_TABLES
HERE = Path(__file__).resolve().parent
LOCAL = HERE.parent / "results" / "tables"
# Optional upstream source; set MR_UPSTREAM_TABLES to re-run against an external copy
UPSTREAM = Path(os.environ.get("MR_UPSTREAM_TABLES", "results/tables"))
SRC = LOCAL if LOCAL.exists() else UPSTREAM

OUT = HERE.parent / "derived"
OUT.mkdir(parents=True, exist_ok=True)

OUTCOMES = {"CONDUCTIO": "Cardiac conduction disorders", "AVBLOCK": "Atrioventricular block"}
FDR_CUT = 0.05


def load(outcome):
    main = pd.read_csv(SRC / f"MR_immune_{outcome}.csv")
    gws = pd.read_csv(SRC / f"MR_immune_GWS_{outcome}.csv")
    main = main.sort_values("p").reset_index(drop=True)
    main["rank"] = main.index + 1
    gws = gws.sort_values("p").reset_index(drop=True)
    gws["grank"] = gws.index + 1
    return main, gws


def sparsity(gws):
    """Instrument-count distribution under the genome-wide threshold."""
    n = len(gws)
    rows = []
    for label, mask in [
        ("1 (Wald ratio only)", gws["nIV"] == 1),
        ("2", gws["nIV"] == 2),
        ("3", gws["nIV"] == 3),
        (">=4", gws["nIV"] >= 4),
    ]:
        rows.append({"nIV": label, "n_phenotypes": int(mask.sum()),
                     "pct": round(100 * mask.mean(), 1)})
    return pd.DataFrame(rows), n


def crosswalk(main, gws):
    """Rank, P and FDR in the primary analysis for traits reaching FDR significance under the genome-wide threshold."""
    sig = gws[gws["FDR"] < FDR_CUT]
    rows = []
    for _, r in sig.iterrows():
        q = main[main["id"] == r["id"]]
        rec = {
            "id": r["id"],
            "trait": r["trait"],
            "gws_nIV": int(r["nIV"]),
            "gws_method": r["method"],
            "gws_OR": round(float(r["OR"]), 3),
            # confidence bounds are required by build_medicine_docs.py::add_table1
            "gws_OR_L": round(float(r["OR_L"]), 3),
            "gws_OR_U": round(float(r["OR_U"]), 3),
            "gws_p": float(r["p"]),
            "gws_FDR": round(float(r["FDR"]), 4),
        }
        if len(q):
            q = q.iloc[0]
            rec.update({
                "main_nIV": int(q["nIV"]),
                "main_rank": int(q["rank"]),
                "main_p": float(q["p"]),
                "main_FDR": round(float(q["FDR"]), 4),
            })
        else:
            rec.update({"main_nIV": None, "main_rank": None,
                        "main_p": None, "main_FDR": None})
        rows.append(rec)
    return pd.DataFrame(rows)


def overlap(main, gws):
    """Overlap of the two thresholds' hit sets, and the trait-wide rank correlation."""
    main_nom = set(main.loc[main["p"] < 0.05, "id"])
    main_fdr = set(main.loc[main["FDR"] < FDR_CUT, "id"])
    gws_fdr = set(gws.loc[gws["FDR"] < FDR_CUT, "id"])

    merged = main.merge(gws, on="id", suffixes=("_main", "_gws"))
    # Spearman = Pearson correlation of ranks; computed directly to avoid a scipy dependency
    if len(merged) > 2:
        rho = merged["p_main"].rank().corr(merged["p_gws"].rank(), method="pearson")
    else:
        rho = float("nan")

    return {
        "main_nominal_n": len(main_nom),
        "main_FDR_sig_n": len(main_fdr),
        "gws_FDR_sig_n": len(gws_fdr),
        "gws_FDR_sig_also_main_FDR_sig": len(gws_fdr & main_fdr),
        "gws_FDR_sig_also_main_nominal": len(gws_fdr & main_nom),
        "phenotypes_in_both_analyses": len(merged),
        "spearman_rho_of_p": round(rho, 3) if rho == rho else None,
    }, merged


def main_():
    print(f"input tables: {SRC}\n")
    all_cross, all_merged = [], []

    for oc, label in OUTCOMES.items():
        m, g = load(oc)
        print("=" * 78)
        print(f"{label}  ({oc})   primary n={len(m)}   genome-wide n={len(g)}")
        print("=" * 78)

        sp, n = sparsity(g)
        print("\n[Table 2] instrument counts under the genome-wide threshold (P<5e-8)")
        print(sp.to_string(index=False))

        st, merged = overlap(m, g)
        print("\n[overlap statistics]")
        for k, v in st.items():
            print(f"  {k:38s} {v}")

        cw = crosswalk(m, g)
        print(f"\n[Table 1] the {len(cw)} traits with FDR<0.05 under the genome-wide threshold, and their position in the primary analysis")
        if len(cw):
            print(cw.to_string(index=False))
        else:
            print("  (none)")

        cw.insert(0, "outcome", label)
        merged.insert(0, "outcome", label)
        all_cross.append(cw)
        all_merged.append(merged)
        print()

    pd.concat(all_cross).to_csv(OUT / "table1_threshold_crosswalk.csv", index=False)
    md = pd.concat(all_merged)
    md["trait"] = md["trait_main"]
    cols = ["outcome", "id", "trait", "nIV_main", "p_main", "FDR_main", "nIV_gws", "p_gws", "FDR_gws"]
    md[[c for c in cols if c in md.columns]].to_csv(OUT / "threshold_instability.csv", index=False)
    print(f"written: {OUT}/table1_threshold_crosswalk.csv")
    print(f"written: {OUT}/threshold_instability.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main_())
