#!/usr/bin/env python3
"""
23: Post-hoc rare-variant diagnostics for the genome-wide-threshold findings
    (Supplementary Tables S14 and S15).

These analyses were run after the two threshold screens (scripts 12b and 20),
in response to the discordance between them. They are post-hoc; they were not
prespecified, and they do not redefine the primary analysis.

What it does
    S14  For each phenotype reaching FDR < 0.05 under the genome-wide instrument
         threshold, tabulate the contributing variants (effect-allele frequency,
         F-statistic), the share of inverse-variance weight each carries under
         both thresholds, leave-one-out estimates where at least three
         instruments are available, and the estimate after excluding
         rs17583875.
    S15  Repeat both threshold screens after excluding instruments with a minor
         allele frequency below 1%, and (S16) across a family of frequency and
         minor-allele-count cut-offs, since any single cut-off is arbitrary. Records with no reported effect-allele
         frequency are also excluded, because their MAF cannot be confirmed;
         there are 107 such records at P < 1e-5 and 27 at P < 5e-8.

Estimator
    Reproduces the convention used by scripts 12b and 20: a single-variant Wald
    ratio when one instrument is available, otherwise inverse-variance
    weighting, fixed-effect at three or fewer variants and random-effects above
    that (the default of MendelianRandomization::mr_ivw). Verified against the
    shipped result tables: 597 of 607 phenotype-level estimates for the
    genome-wide screen agree to 1e-6 on both the estimate and its standard
    error, and all aggregate counts reproduce exactly.

Inputs
    results/tables/immune_ivs.csv        harmonised exposure-side instruments
    <finngen_dir>/finngen_R11_I9_CONDUCTIO.gz
    <finngen_dir>/finngen_R11_I9_AVBLOCK.gz
        FinnGen release 11 summary statistics; not redistributed here, see
        https://www.finngen.fi/en/access_results

Usage
    python3 scripts/23_rare_variant_diagnostics.py <finngen_dir> [outdir]
"""
from __future__ import annotations

import collections
import csv
import gzip
import math
import os
import sys

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
FOCAL = "rs17583875"
Z = 1.959963985
MAF_MIN = 0.01
# Any single frequency cut-off is arbitrary, so the exclusion is repeated across a
# family of thresholds: minor allele frequency, and minor-allele count, which scales
# with the analysis sample size of each trait.
MAF_GRID = (0.005, 0.01, 0.02, 0.05)
MAC_GRID = (20, 50, 100)
GWS = 5e-8


def load_instruments(path):
    with open(path, newline="", encoding="utf-8") as f:
        return [dict(r) for r in csv.DictReader(f)]


def load_outcome(path, wanted):
    """Return {rsid: (ref, alt, beta, se)} for the requested rsIDs."""
    out = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            for rs in p[col["rsids"]].split(","):
                if rs not in wanted:
                    continue
                try:
                    beta = float(p[col["beta"]])
                    se = float(p[col["sebeta"]])
                except ValueError:
                    break
                if se > 0:
                    out[rs] = (p[col["ref"]].upper(), p[col["alt"]].upper(), beta, se)
                break
    return out


def harmonise(instruments, outcome):
    """Allele harmonisation, palindromic filter and F > 10, as in script 20.

    An absent effect-allele frequency leaves the palindromic condition
    undecidable. data.table does not select rows whose `i` expression is NA, so
    the R code keeps those records; this reproduces that behaviour.
    """
    kept = []
    for r in instruments:
        rs = r["rsid"]
        if rs not in outcome:
            continue
        ref, alt, beta_o, se_o = outcome[rs]
        ea, nea = r["ea"].upper(), r["nea"].upper()
        if ea == alt and nea == ref:
            aligned = beta_o
        elif ea == ref and nea == alt:
            aligned = -beta_o
        else:
            continue
        eaf = float(r["eaf"]) if r["eaf"] else float("nan")
        if COMP.get(ea) == nea and eaf == eaf and min(eaf, 1 - eaf) > 0.42:
            continue
        beta_e, se_e = float(r["beta"]), float(r["se"])
        f_stat = (beta_e / se_e) ** 2
        if f_stat <= 10:
            continue
        kept.append(
            dict(id=r["id"], trait=r["trait"], rsid=rs, eaf=eaf, beta_e=beta_e,
                 se_e=se_e, beta_o=aligned, se_o=se_o, F=f_stat)
        )
    return kept


def weights(group):
    return [(g["beta_e"] ** 2) / (g["se_o"] ** 2) for g in group]


def estimate(group):
    """Wald ratio for one instrument, otherwise IVW (fixed at n <= 3)."""
    if len(group) == 1:
        g = group[0]
        return g["beta_o"] / g["beta_e"], abs(g["se_o"] / g["beta_e"]), "Wald"
    w = weights(group)
    beta = sum(g["beta_e"] * g["beta_o"] / g["se_o"] ** 2 for g in group) / sum(w)
    se = math.sqrt(1 / sum(w))
    if len(group) > 3:
        q = sum((g["beta_o"] - beta * g["beta_e"]) ** 2 / g["se_o"] ** 2 for g in group)
        se *= max(1.0, math.sqrt(q / (len(group) - 1)))
    return beta, se, "IVW"


def pvalue(beta, se):
    return 2 * 0.5 * math.erfc(abs(beta / se) / math.sqrt(2))


def screen(instruments, outcome):
    """One phenotype-wide screen with Benjamini-Hochberg FDR within the outcome."""
    by_id = collections.defaultdict(list)
    for x in harmonise(instruments, outcome):
        by_id[x["id"]].append(x)
    rows = []
    for pid, group in by_id.items():
        beta, se, method = estimate(group)
        rows.append(dict(id=pid, trait=group[0]["trait"], nIV=len(group), method=method,
                         beta=beta, se=se, p=pvalue(beta, se)))
    rows.sort(key=lambda r: r["p"])
    running = 1.0
    for i in range(len(rows) - 1, -1, -1):
        running = min(running, rows[i]["p"] * len(rows) / (i + 1))
        rows[i]["FDR"] = running
    return rows


def maf_of(record):
    """min(EAF, 1 - EAF), or None when no effect-allele frequency is reported."""
    if not record["eaf"]:
        return None
    eaf = float(record["eaf"])
    return min(eaf, 1 - eaf)


def maf_pass(record, cutoff=MAF_MIN):
    maf = maf_of(record)
    return maf is not None and maf >= cutoff


def mac_pass(record, cutoff):
    """Minor-allele count, 2 * n * MAF, using the analysis n reported per trait."""
    maf = maf_of(record)
    return maf is not None and 2 * int(record["n"]) * maf >= cutoff


def weight_share(group, rsid):
    w = weights(group)
    focal = [wi for wi, g in zip(w, group) if g["rsid"] == rsid]
    return focal[0] / sum(w) * 100 if focal else 0.0


def fmt_or(beta, se):
    return (f"{math.exp(beta):.3f} "
            f"({math.exp(beta - Z * se):.3f}-{math.exp(beta + Z * se):.3f})")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    finngen = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else "results/tables"
    os.makedirs(outdir, exist_ok=True)

    ivs = load_instruments(os.path.join("results", "tables", "immune_ivs.csv"))
    # Published estimates come from the shipped result tables, so that the columns
    # this table shares with Table 1 cannot drift from it. Only the new diagnostic
    # columns (weight share, leave-one-out, estimate without rs17583875) are
    # recomputed here.
    published_gws = {r["id"]: r for r in load_instruments(
        os.path.join("results", "tables", "MR_immune_GWS_CONDUCTIO.csv"))}
    published_primary = {r["id"]: r for r in load_instruments(
        os.path.join("results", "tables", "MR_immune_CONDUCTIO.csv"))}
    primary_rank = {r["id"]: i + 1 for i, r in enumerate(
        sorted(published_primary.values(), key=lambda r: float(r["p"])))}
    gws = [r for r in ivs if float(r["p"]) < GWS]
    wanted = {r["rsid"] for r in ivs}

    outcomes = {}
    for tag in ("CONDUCTIO", "AVBLOCK"):
        path = os.path.join(finngen, f"finngen_R11_I9_{tag}.gz")
        print(f"reading {path}", flush=True)
        outcomes[tag] = load_outcome(path, wanted)

    # ---- Supplementary Table S14 -------------------------------------------
    conduct = outcomes["CONDUCTIO"]
    gws_rows = screen(gws, conduct)

    gws_groups = collections.defaultdict(list)
    for x in harmonise(gws, conduct):
        gws_groups[x["id"]].append(x)
    primary_groups = collections.defaultdict(list)
    for x in harmonise(ivs, conduct):
        primary_groups[x["id"]].append(x)

    s14 = []
    for row in sorted([r for r in gws_rows if r["FDR"] < 0.05], key=lambda r: r["FDR"]):
        group = sorted(gws_groups[row["id"]], key=lambda g: g["rsid"])
        remainder = [g for g in group if g["rsid"] != FOCAL]
        if remainder:
            b, se, _ = estimate(remainder)
            without = f"{fmt_or(b, se)}, P = {pvalue(b, se):.3g}"
        else:
            without = "no instrument remains; not estimable"
        if len(group) >= 3:
            loo = "; ".join(
                f"omit {g['rsid']}: OR {math.exp(estimate([x for x in group if x is not g])[0]):.3f}, "
                f"P {pvalue(*estimate([x for x in group if x is not g])[:2]):.3g}"
                for g in group)
        else:
            loo = "not informative (<3 instruments)"
        pub_g = published_gws[row["id"]]
        pub_p = published_primary[row["id"]]
        s14.append(dict(
            Immunophenotype=row["trait"],
            GWAS_Catalog_accession=row["id"].replace("ebi-a-", ""),
            Instruments_genome_wide=int(pub_g["nIV"]),
            Instrument_rsIDs="; ".join(
                f"{g['rsid']} (EAF {g['eaf']:.4f}, F {g['F']:.0f})" for g in group),
            OR_95CI=f"{float(pub_g['OR']):.3f} "
                    f"({float(pub_g['OR_L']):.3f}-{float(pub_g['OR_U']):.3f})",
            P=f"{float(pub_g['p']):.3g}",
            FDR=f"{float(pub_g['FDR']):.3f}".lstrip("0"),
            Estimator=pub_g["method"],
            rs17583875_weight_genome_wide_pct=f"{weight_share(group, FOCAL):.1f}",
            Instruments_primary=len(primary_groups[row["id"]]),
            rs17583875_weight_primary_pct=f"{weight_share(primary_groups[row['id']], FOCAL):.1f}",
            Primary_rank=f"{primary_rank[row['id']]}/731",
            Primary_P=f"{float(pub_p['p']):.2f}".lstrip("0"),
            Primary_FDR=f"{float(pub_p['FDR']):.2f}".lstrip("0"),
            Leave_one_out=loo,
            Estimate_excluding_rs17583875=without,
        ))

    # ---- Supplementary Table S15 -------------------------------------------
    s15 = []
    for threshold, subset, label in (
        ("P < 1 x 10-5", ivs, "Primary (instrument availability)"),
        ("P < 5 x 10-8", gws, "Genome-wide"),
    ):
        kept = [r for r in subset if maf_pass(r)]
        missing = sum(1 for r in subset if not r["eaf"])
        for tag in ("CONDUCTIO", "AVBLOCK"):
            base = screen(subset, outcomes[tag])
            filtered = screen(kept, outcomes[tag])
            s15.append(dict(
                Instrument_threshold=threshold, Analysis=label, Outcome=tag,
                IV_records_all=len(subset), IV_records_MAF_ge_1pct=len(kept),
                IV_records_MAF_missing_excluded=missing,
                Phenotypes_all=len(base), Phenotypes_MAF_ge_1pct=len(filtered),
                Nominal_all=sum(1 for r in base if r["p"] < 0.05),
                Nominal_MAF_ge_1pct=sum(1 for r in filtered if r["p"] < 0.05),
                Expected_nominal_MAF_ge_1pct=math.floor(len(filtered) * 0.05 + 0.5),
                FDR_significant_all=sum(1 for r in base if r["FDR"] < 0.05),
                FDR_significant_MAF_ge_1pct=sum(1 for r in filtered if r["FDR"] < 0.05),
                Smallest_FDR_all=f"{min(r['FDR'] for r in base):.3f}".lstrip("0"),
                Smallest_FDR_MAF_ge_1pct=f"{min(r['FDR'] for r in filtered):.3f}".lstrip("0"),
            ))

    # ---- Supplementary Table S16: cut-point family -------------------------
    s16 = []
    for label, subset in (("P < 1 x 10-5", ivs), ("P < 5 x 10-8", gws)):
        for kind, grid in (("MAF", MAF_GRID), ("MAC", MAC_GRID)):
            for cutoff in grid:
                keep = ([r for r in subset if maf_pass(r, cutoff)] if kind == "MAF"
                        else [r for r in subset if mac_pass(r, cutoff)])
                for tag in ("CONDUCTIO", "AVBLOCK"):
                    rows_ = screen(keep, outcomes[tag])
                    s16.append(dict(
                        Instrument_threshold=label, Outcome=tag,
                        Exclusion=(f"MAF >= {cutoff * 100:g}%" if kind == "MAF" else f"MAC >= {cutoff:g}"),
                        IV_records_retained=len(keep),
                        Phenotypes=len(rows_),
                        Nominal=sum(1 for r in rows_ if r["p"] < 0.05),
                        Expected_nominal=math.floor(len(rows_) * 0.05 + 0.5),
                        FDR_significant=sum(1 for r in rows_ if r["FDR"] < 0.05),
                        Smallest_FDR=f"{min(r['FDR'] for r in rows_):.3f}".lstrip("0"),
                    ))

    for name, rows in (("S14_rare_variant_diagnostics", s14),
                       ("S15_maf_sensitivity", s15),
                       ("S16_cutpoint_family", s16)):
        path = os.path.join(outdir, f"{name}.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
