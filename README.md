# Supplementary Code

Manuscript: Circulating immune-cell traits and cardiac conduction disease: a phenotype-wide Mendelian randomization study

This folder contains the scripts used to reproduce the Mendelian randomization screens, sensitivity analyses, reverse-direction analyses, power calculations, figures, and supplementary tables reported in the manuscript.

The analysis code is also available at:

https://github.com/hkzyyliwei2026/MR

## Software

The analyses were run in R 4.6.1. Main R packages used by the MR and figure workflow include:

- MendelianRandomization
- ieugwasr
- data.table
- ggplot2
- dplyr
- readr

MR-PRESSO is computed with the official `MRPRESSO` package
(`remotes::install_github("rondolab/MR-PRESSO")`), which `scripts/17_presso_steiger.R` requires. The
Steiger directionality test in the same script is a self-contained implementation of Hemani et al and
needs no additional package.

MR-PRESSO is used **only** as a global heterogeneity diagnostic. The global test was non-significant
for all five analyses reported in Supplementary Table S8, so the package did not proceed to the
variant-level outlier test and no outlier-corrected estimates exist; no reported result depends on
them.

Supplementary table and figure assembly also uses Python 3 with:

- pandas
- matplotlib
- openpyxl

`sessionInfo.txt` records the R session available in the final manuscript-generation and checking environment. Users reproducing the full MR workflow should install the packages listed above and authenticate OpenGWAS access where required.

## Public Input Data

Raw GWAS summary statistics are not redistributed in this code supplement because they are available from their original public resources.

Immune-cell GWAS summary statistics:

- GWAS Catalog accessions GCST90001391-GCST90002121
- Study: Orrù et al., Nature Genetics, 2020

Lymphocyte-count GWAS summary statistics (targeted re-analysis only):

- GWAS Catalog accession GCST90002316
- Study: Chen MH et al., Cell, 2020 (Blood Cell Consortium), European ancestry, 524,923 participants
- Instruments are read from the GWAS Catalog REST endpoint by `scripts/22_chen_reanalysis.R`

Outcome GWAS summary statistics:

- FinnGen release 11
- Cardiac conduction disorders: I9_CONDUCTIO
- Atrioventricular block: I9_AVBLOCK
- Access: https://www.finngen.fi/en/access_results

## Directory layout

```
scripts/            analysis scripts (run from the package root)
results/tables/     result tables produced by the scripts and used as inputs downstream
results/figures/    figure files (PNG, PDF, SVG)
derived/            cross-threshold comparison tables
supplementary/      Supplementary_Tables.xlsx
```

`data/outcome/` is referenced by several scripts but is not redistributed here; the FinnGen
release 11 files must be obtained separately (see Public Input Data above).

The sub-endpoint and specificity analyses need three further FinnGen release 11 files beyond the
two primary outcomes: `finngen_R11_I9_LBBB.gz`, `finngen_R11_I9_RBBB.gz` and
`finngen_R11_I9_AF.gz`, obtained from the same source.

The 600-dpi TIFF versions of Figure 1, the forest plot and the leave-one-out plot are not
included: they add roughly 216 MB uncompressed while carrying no information beyond the PNG,
PDF and SVG files shipped alongside them. `scripts/15_figures.R` and `scripts/16_flowchart.R`
write the TIFFs on each run, so they can be regenerated locally if a print-resolution raster
is needed.

## Run Order

Run scripts from the project root unless otherwise noted.

```bash
Rscript scripts/12a_immune_extract_ivs.R
Rscript scripts/12b_immune_mr.R
Rscript scripts/13_immune_sensitivity.R
Rscript scripts/14_reverse_mr.R
Rscript scripts/17_presso_steiger.R
Rscript scripts/20_gws_threshold_sensitivity.R
python3 scripts/make_threshold_analysis.py                     # writes derived/threshold_instability.csv
Rscript scripts/21_nominal_excess_and_concordance.R            # reads that file; must run after it
python3 scripts/23_rare_variant_diagnostics.py <finngen_dir>   # post-hoc; Tables S14-S16
Rscript scripts/25_phenotype_power.R <finngen_dir>             # Table S18; finngen_dir optional
Rscript scripts/26_diagnostic_coverage.R <finngen_dir>         # Tables S19-S20; finngen_dir optional
Rscript scripts/27_instrument_qc_sensitivity.R <finngen_dir>   # Table S21; finngen_dir optional
Rscript scripts/29_subendpoint_screen.R <finngen_dir>          # Table S23; needs I9_LBBB and I9_RBBB
Rscript scripts/30_endpoint_case_counts.R <finngen_dir>        # recovers component-endpoint case counts
Rscript scripts/18_volcano_power.R
Rscript scripts/15_figures.R
Rscript scripts/16_flowchart.R
Rscript scripts/regen_figure2.R
python3 scripts/make_figure3.py
python3 scripts/make_supplementary.py
```

The order matters in one place: `21_nominal_excess_and_concordance.R` reads
`derived/threshold_instability.csv`, which `make_threshold_analysis.py` writes, so the Python
step has to come first. Earlier versions of this file listed them the other way round.

`scripts/21_nominal_excess_and_concordance.R` re-derives the nominal-finding counts, exploratory chance-calculation outputs, the cross-endpoint concordance summaries and the effect of the FDR correction denominator, all from the result tables above. The manuscript reports the expected nominal counts descriptively and does not use a formal binomial or hypergeometric tail-probability test, because the immunophenotypes and endpoints are correlated.

`scripts/22_chen_reanalysis.R` is run separately because it needs external inputs: it retrieves the lymphocyte-count associations for GWAS Catalog accession GCST90002316 over the network and requires a PLINK binary, a PLINK-format 1000 Genomes European reference panel and the FinnGen release 11 outcome files. It writes the targeted re-analysis results reported as Supplementary Table S13. Its output is also included in this package as `results/tables/S13_chen_reanalysis.csv`, together with the run log, so that Table S13 can be checked and rebuilt without re-running the network and PLINK steps; `scripts/make_supplementary.py` reads that file.

```bash
Rscript scripts/22_chen_reanalysis.R <outdir> <plink> <ref_prefix> <finngen_dir>
Rscript scripts/28_af_specificity.R <outdir> <finngen_dir>     # Table S22; reuses that instrument set
cp <outdir>/chen_reanalysis_summary.csv results/tables/S13_chen_reanalysis.csv
cp <outdir>/chen_reanalysis_log.txt     results/tables/S13_chen_reanalysis_log.txt
```

The copy step is required: `22_chen_reanalysis.R` writes to the output directory given on the
command line, while `make_supplementary.py` reads the two files under `results/tables/`.

By default `make_supplementary.py` writes the workbook to `supplementary/` inside this package. To write it elsewhere, set `SUBMISSION_DIR`:

```bash
SUBMISSION_DIR=/path/to/output python3 scripts/make_supplementary.py
```

`12a_immune_extract_ivs.R` queries OpenGWAS through `ieugwasr` and requires a valid `OPENGWAS_JWT` environment variable if OpenGWAS authentication is enabled for the user account. The token itself is not included.

## Outputs

Main tabular output files are stored or generated under `results/tables/` and `derived/`. Key files include:

- `immune_ivs.csv`
- `MR_immune_CONDUCTIO.csv`
- `MR_immune_AVBLOCK.csv`
- `immune_sensitivity.csv`
- `immune_leaveoneout.csv`
- `MR_reverse.csv`
- `immune_presso_steiger.csv`
- `power_analysis.csv`
- `power_hits.csv`
- `MR_immune_GWS_CONDUCTIO.csv`
- `MR_immune_GWS_AVBLOCK.csv`
- `threshold_sensitivity_summary.csv`
- `derived/threshold_instability.csv`
- `derived/table1_threshold_crosswalk.csv`
- `S13_chen_reanalysis.csv` and `S13_chen_reanalysis_log.txt`
- `S14_rare_variant_diagnostics.csv`, `S15_maf_sensitivity.csv`, `S16_cutpoint_family.csv`
- `S18_phenotype_power.csv`
- `S19_coverage_summary.csv` (workbook Table S19, diagnostic computability summary) and `S19_diagnostic_coverage.csv` (workbook Table S20, per-phenotype Cochran's Q output)
- `S21_qc_summary.csv` and `S21_instrument_qc_sensitivity.csv`
- `S22_af_specificity.csv`
- `S23_subendpoint_summary.csv` and `S23_subendpoint_screen.csv`
- `S23_endpoint_case_counts.csv`

Main figure outputs include:

- `Fig1_flowchart.png`
- `Figure_2_Volcano.png`
- `Figure_3_ThresholdInstability.png`
- `Figure_3_ThresholdInstability.pdf`
- `Supplementary_Figure_S1_Forest.png`
- `Supplementary_Figure_S2_Leaveoneout.png`

`results/figures/` mixes two naming conventions, because the figure order changed after the plotting scripts were written. `Fig1_flowchart.*` is the manuscript's Figure 1, `Fig2_forest.*` is Supplementary Figure S1 and `Fig3_leaveoneout.*` is Supplementary Figure S2, while `Figure_2_Volcano.png` and `Figure_3_ThresholdInstability.*` already carry the manuscript numbering. None of these files is uploaded directly. `scripts/package_submission.py` maps each one to its manuscript number and writes the submission copies, converting the three main figures to `Figure 1.tif`, `Figure 2.tif` and `Figure 3.tif` (RGB, LZW, 600 dpi) because the journal does not accept PNG for main figures.

The supplementary workbook generated by `make_supplementary.py` writes `supplementary/Supplementary_Tables.xlsx` containing Supplementary Tables S1-S23. Every sheet is built from a file in `results/tables/`, so the workbook can be rebuilt from a clean checkout: S1-S12 from the MR result tables, S13 from `S13_chen_reanalysis.csv`, S14-S16 from the post-hoc rare-variant outputs, S17 from the endpoint definitions, S18 from `S18_phenotype_power.csv`, S19 from `S19_coverage_summary.csv`, S20 from `S19_diagnostic_coverage.csv`, S21 from `S21_qc_summary.csv` (27_instrument_qc_sensitivity.R), S22 from `S22_af_specificity.csv` (28_af_specificity.R), and S23 from `S23_subendpoint_summary.csv` plus `S23_endpoint_case_counts.csv` (29_subendpoint_screen.R and 30_endpoint_case_counts.R). Earlier versions carried S13-S16 forward from an existing copy of the workbook, which meant those sheets could not be regenerated from scratch. S18 is read from `S18_phenotype_power.csv` rather than recomputed, because instrument counts and R2 must come from the harmonised instrument set used by `12b_immune_mr.R` in order to match Supplementary Tables S1 and S2.

`scripts/build_medicine_docs.py` builds the manuscript, cover letter and STROBE-MR checklist
DOCX files from `manuscript_medicine.md`, the cover-letter Markdown and
`table1_threshold_crosswalk.csv`. It also holds the text of Table 1, the three figure legends
and every STROBE-MR checklist entry, so it is included here for completeness; it expects the
submission-package directory layout rather than this package's layout.

## Reproducibility Notes

The code is provided to document and reproduce the analytic workflow. Minor differences may occur if upstream GWAS files, OpenGWAS endpoints, or package versions are updated after publication. No individual-level participant data are required or included.
