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
SKIP_DIRS = {"__pycache__", ".git", ".ipynb_checkpoints", "data"}
SKIP_FILES = {".DS_Store"}
# TIFFs are upload artefacts, not analysis code. 16_flowchart.R's tiff() device also ignores its
# own compression setting and writes an 80 MB uncompressed RGBA file, which alone would push the
# archive against the journal's 10 MB per-file limit. The submission TIFFs are rebuilt from the
# PNGs by build_figures(), so nothing is lost by excluding them here.
SKIP_SUFFIXES = (".tif", ".tiff")
# Paths excluded because the archive would otherwise carry the same bytes twice, or carry a
# stale rename of a file the plotting scripts regenerate under their own name. The workbook is
# uploaded separately as Supplemental Digital Content and is rebuilt from results/tables by
# make_supplementary.py, so dropping it costs no reproducibility.
SKIP_RELPATHS = {
    "supplementary/Supplementary_Tables.xlsx",
}


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
    # archive that differs only where the content differs. The one exception is the workbook:
    # openpyxl does not write reproducible bytes, so it is compared cell-wise, not by hash.
    members, dropped = [], []
    for root, dirs, files in os.walk(upload):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
        for f in sorted(files):
            if f in SKIP_FILES or f.lower().endswith(SKIP_SUFFIXES):
                continue
            full = os.path.join(root, f)
            rel = os.path.relpath(full, upload).replace(os.sep, "/")
            if rel in SKIP_RELPATHS:
                dropped.append((rel, os.path.getsize(full)))
                continue
            members.append((full, rel))
    tmp = out + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for full, rel in members:
            info = zipfile.ZipInfo("Supplementary_Code/" + rel, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            with open(full, "rb") as fh:
                z.writestr(info, fh.read())
    os.replace(tmp, out)
    for rel, size in dropped:
        print(f"    excluded from archive: {rel} ({size/1e6:.2f} MB)")
    return len(members)


# Medicine accepts TIF, EPS, PPT or DOC for figures, and requires grayscale or RGB. The
# plotting scripts emit RGBA PNG, so the upload copies are converted here rather than by hand.
# LZW keeps 600-dpi files near 1.5 MB, well inside the 10 MB per-file limit, so the figures
# are not downsampled.
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

    n = build_zip(upload, os.path.join(SUB, "supplemental", "Supplementary_Code.zip"))
    print(f"--- code archive rebuilt ({n} files)")

    run(os.path.join(SUB, "source", "build_medicine_docs.py"), SUB)
    print("package rebuilt")


if __name__ == "__main__":
    main()
