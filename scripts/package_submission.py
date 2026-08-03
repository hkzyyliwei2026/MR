#!/usr/bin/env python3
"""Rebuild every submission deliverable from its generator, in dependency order.

The workbook and the code archive used to be assembled by hand, so a change to
make_supplementary.py could reach the manuscript without reaching the .xlsx that
reviewers open, and a change to any script could ship inside a stale zip. Running
this file is the only supported way to rebuild the package.

  1. make_supplementary.py         -> upload/supplementary/Supplementary_Tables.xlsx
  2. copy                          -> supplemental/Supplementary_Tables.xlsx
  3. zip the upload tree           -> supplemental/Supplementary_Code.zip
  4. build_medicine_docs.py        -> manuscript / title page / STROBE checklist

Usage: python3 source/package_submission.py [--upload DIR]
"""
import argparse
import os
import shutil
import subprocess
import sys
import zipfile

from PIL import Image

SUB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # overridden by --submission
DEFAULT_UPLOAD = "/Users/primihub/Desktop/MR_github_upload"

# Everything under the upload tree ships except build artefacts and the raw GWAS
# downloads, which are gigabytes and are already public at their original sources.
# The code archive freezes the analysis code at submission time, so it carries only what is
# needed to read and re-run it: the scripts, the session information, the README, and the two
# small derived inputs that the figure and document builders read. Result tables and figures are
# deliberately excluded - the results are uploaded separately as Supplementary Tables S1-S23, and
# the raw GWAS summary statistics are public at their original sources and far too large to ship.
# This is an allowlist rather than a skip list, so a large file added to the repository later
# cannot silently inflate the archive past the journal's 10 MB per-file limit.
ARCHIVE_DIRS = ("scripts", "derived")
ARCHIVE_FILES = ("README.md", "sessionInfo.txt")
SKIP_DIRS = {"__pycache__", ".git", ".ipynb_checkpoints"}
SKIP_FILES = {".DS_Store"}
# 16_flowchart.R's tiff device ignores its own compression setting and writes an 80 MB
# uncompressed file; the submission TIFFs are built from the PNGs by build_figures().
SKIP_SUFFIXES = (".tif", ".tiff")


def archive_provenance(upload):
    """Record which commit the archived code came from, so the frozen copy can be matched
    against the repository even after the repository moves on."""
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=upload,
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=upload,
                               capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "Source repository: https://github.com/hkzyyliwei2026/MR\nCommit: not available\n"
    state = "with uncommitted local changes" if dirty else "clean working tree"
    return (
        "This archive is the analysis code as submitted, frozen at the commit below.\n"
        "The repository may have moved on since; use this copy to reproduce the article.\n\n"
        "Source repository: https://github.com/hkzyyliwei2026/MR\n"
        f"Commit: {head}\n"
        f"Working tree at archive time: {state}\n\n"
        "Contents: analysis and plotting scripts (scripts/), R session information\n"
        "(sessionInfo.txt), the repository README, and the two derived inputs read by the\n"
        "figure and document builders (derived/).\n\n"
        "Not included: result tables and figures, which are uploaded separately as\n"
        "Supplementary Tables S1-S23 and Supplementary Figures S1-S2; and the GWAS summary\n"
        "statistics, which are public at the GWAS Catalog and FinnGen and are too large to ship.\n"
    )


def run(script, cwd):
    print(f"--- {os.path.basename(script)}")
    r = subprocess.run([sys.executable, script], cwd=cwd, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write(r.stdout + r.stderr)
        raise SystemExit(f"{script} failed with exit code {r.returncode}")
    tail = [l for l in r.stdout.strip().splitlines() if l.strip()][-1:]
    for l in tail:
        print("   ", l)


def build_zip(upload, out):
    # Deterministic member order and a fixed timestamp, so an unchanged tree produces an
    # archive that differs only where the content differs.
    members = []
    for name in ARCHIVE_FILES:
        full = os.path.join(upload, name)
        if os.path.exists(full):
            members.append((full, name))
    for d in ARCHIVE_DIRS:
        base = os.path.join(upload, d)
        for root, dirs, files in os.walk(base):
            dirs[:] = sorted(x for x in dirs if x not in SKIP_DIRS)
            for f in sorted(files):
                if f in SKIP_FILES or f.lower().endswith(SKIP_SUFFIXES):
                    continue
                full = os.path.join(root, f)
                members.append((full, os.path.relpath(full, upload).replace(os.sep, "/")))
    members.sort(key=lambda m: m[1])

    skipped = sorted(
        e.name for e in os.scandir(upload)
        if not e.name.startswith(".")
        and e.name not in ARCHIVE_DIRS
        and e.name not in ARCHIVE_FILES
    )
    tmp = out + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for full, rel in members:
            info = zipfile.ZipInfo("Supplementary_Code/" + rel, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            with open(full, "rb") as fh:
                z.writestr(info, fh.read())
        prov = zipfile.ZipInfo("Supplementary_Code/ARCHIVE_INFO.txt", date_time=(2026, 1, 1, 0, 0, 0))
        prov.compress_type = zipfile.ZIP_DEFLATED
        prov.external_attr = 0o644 << 16
        z.writestr(prov, archive_provenance(upload))
    os.replace(tmp, out)
    if skipped:
        print(f"    not archived (uploaded separately or public at source): {', '.join(skipped)}")
    return len(members) + 1


# Figures flow from the plotting scripts' output directory, not from hand-placed copies in the
# submission folder, so a re-run of a plotting script always reaches the upload. Medicine accepts
# TIF, EPS, PPT or DOC and requires grayscale or RGB; the scripts emit RGBA PNG, so the upload
# copies are converted here. LZW keeps 600-dpi files near 1.5 MB, inside the 10 MB per-file
# limit, so nothing is downsampled. Supplemental figures may be any format and stay PNG.
FIGURES = [
    ("Fig1_flowchart.png", "primary", "Figure_1_Flowchart.png", "Figure 1.tif"),
    ("Figure_2_Volcano.png", "primary", "Figure_2_Volcano.png", "Figure 2.tif"),
    ("Figure_3_ThresholdInstability.png", "primary", "Figure_3_ThresholdInstability.png", "Figure 3.tif"),
    ("Fig2_forest.png", "supplemental", "Supplementary_Figure_S1_Forest.png", None),
    ("Fig3_leaveoneout.png", "supplemental", "Supplementary_Figure_S2_Leaveoneout.png", None),
]


def build_figures(upload):
    out = []
    for src, folder, png_name, tif_name in FIGURES:
        src_path = os.path.join(upload, "results", "figures", src)
        if not os.path.exists(src_path):
            raise SystemExit(f"missing plotting-script output: {src_path}")
        shutil.copy2(src_path, os.path.join(SUB, folder, png_name))
        if tif_name is None:
            continue
        im = Image.open(src_path)
        dpi = im.info.get("dpi", (600, 600))
        if im.mode == "RGBA":
            alpha = im.split()[3]
            if alpha.getextrema() != (255, 255):
                raise SystemExit(f"{src} has real transparency; flattening would change it")
            flat = Image.new("RGB", im.size, (255, 255, 255))
            flat.paste(im, mask=alpha)
        else:
            flat = im.convert("RGB")
        target = os.path.join(SUB, folder, tif_name)
        flat.save(target, format="TIFF", compression="tiff_lzw", dpi=dpi)
        out.append((tif_name, os.path.getsize(target)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--upload", default=DEFAULT_UPLOAD)
    ap.add_argument("--submission", default=None,
                    help="submission tree containing primary/ and supplemental/. Defaults to the "
                         "parent of this file's directory, which is correct when run from "
                         "<submission>/source/ but not from the repository's scripts/ copy.")
    a = ap.parse_args()
    global SUB
    if a.submission:
        SUB = os.path.abspath(a.submission)
    if not os.path.isdir(os.path.join(SUB, "primary")):
        raise SystemExit(
            f"not a submission tree: {SUB}\nPass --submission <dir> when running the archived "
            "copy from the repository.")
    upload = os.path.abspath(a.upload)
    if not os.path.isdir(os.path.join(upload, "scripts")):
        raise SystemExit(f"not an upload tree: {upload}")

    run(os.path.join(upload, "scripts", "make_supplementary.py"), upload)

    src = os.path.join(upload, "supplementary", "Supplementary_Tables.xlsx")
    dst = os.path.join(SUB, "supplemental", "Supplementary_Tables.xlsx")
    shutil.copy2(src, dst)
    print(f"--- workbook copied ({os.path.getsize(dst):,} bytes)")

    for name, size in build_figures(upload):
        print(f"--- {name} ({size/1e6:.2f} MB, TIFF/RGB)")

    run(os.path.join(SUB, "source", "build_medicine_docs.py"), SUB)
    print("package rebuilt")


if __name__ == "__main__":
    main()
