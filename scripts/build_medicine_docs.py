from __future__ import annotations

import csv
import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source"
PRIMARY = ROOT / "primary"
SUPP = ROOT / "supplemental"
COVER = ROOT / "cover"

MANUSCRIPT_MD = SOURCE / "manuscript_medicine.md"
COVER_MD = COVER / "Cover_Letter_Medicine.md"
TABLE1_CSV = SOURCE / "table1_threshold_crosswalk.csv"

TITLE = (
    "Circulating immune-cell traits and cardiac conduction disease: "
    "a phenotype-wide Mendelian randomization study"
)

# AMA-style in-text reference numbers, e.g. [1], [3,4], [12-16], rendered superscript
CITATION_PATTERN = r"\[\d+(?:[,\-]\d+)*\]"


def setup_document(line_numbers: bool = False) -> Document:
    doc = Document()
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(11)
    for section in doc.sections:
        section.top_margin = Pt(72)
        section.bottom_margin = Pt(72)
        section.left_margin = Pt(72)
        section.right_margin = Pt(72)
        if not line_numbers:
            continue
        sect_pr = section._sectPr
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1")
        ln.set(qn("w:restart"), "continuous")
        ln.set(qn("w:distance"), "360")
        # CT_SectPr is a sequence: lnNumType must precede pgNumType/cols/docGrid.
        # Appending puts it last, which is schema-invalid and can be dropped by a
        # strict consumer, so insert it ahead of the first element that must follow it.
        successors = ("w:pgNumType", "w:cols", "w:formProt", "w:vAlign", "w:docGrid")
        anchor = next((c for c in sect_pr if c.tag in {qn(t) for t in successors}), None)
        if anchor is None:
            sect_pr.append(ln)
        else:
            anchor.addprevious(ln)
    return doc


def add_runs(paragraph, text: str) -> None:
    text = text.replace(" x 10", " × 10")
    text = text.replace("<=", "≤").replace(">=", "≥")
    parts = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text)
    for part in parts:
        if not part:
            continue
        bold = part.startswith("**") and part.endswith("**")
        italic = part.startswith("*") and part.endswith("*") and not bold
        inner = part[2:-2] if bold else part[1:-1] if italic else part
        for token in re.split(r"(\^-?\d+|" + CITATION_PATTERN + r")", inner):
            if not token:
                continue
            superscript = token.startswith("^") or re.fullmatch(CITATION_PATTERN, token)
            run = paragraph.add_run(token[1:] if token.startswith("^") else token)
            run.bold = bold
            run.italic = italic
            if superscript:
                run.font.superscript = True


def add_paragraph(doc: Document, text: str = "", style: str | None = None):
    p = doc.add_paragraph(style=style) if style else doc.add_paragraph()
    add_runs(p, text)
    return p


def split_refs(md: str) -> tuple[str, str]:
    body, refs = md.split("\n## References\n", 1)
    return body.rstrip(), refs.strip()


def section_text(md: str, header: str) -> str:
    match = re.search(rf"## {re.escape(header)}\n\n(.*?)(?=\n## |\Z)", md, re.S)
    return match.group(1).strip() if match else ""


def add_heading(doc: Document, title: str, level: int = 1) -> None:
    # Medicine numbers body sections as "1.", "2.1."; the in-text cross-references
    # ("Section 3.1") only resolve if the numbers survive into the heading.
    clean = re.sub(r"^(\d+(?:\.\d+)*)\s+", r"\1. ", title.strip())
    doc.add_heading(clean, level=level)


def add_table1(doc: Document) -> None:
    rows = []
    with TABLE1_CSV.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["outcome"] == "Cardiac conduction disorders":
                rows.append(row)

    add_paragraph(
        doc,
        "Table 1. Genome-wide-threshold signals in relation to the primary screen. "
        "The five threshold-specific immunophenotypes reaching FDR significance for "
        "cardiac conduction disorders under the stricter genome-wide instrument threshold "
        "(P < 5 x 10^-8) are shown against their position in the primary screen "
        "(P < 1 x 10^-5). None reached nominal significance in the primary analysis. "
        "Variant-level composition of these five instrument sets is given in "
        "Supplementary Table S14. GWS, genome-wide significant (P < 5 x 10^-8); "
        "IV, instrumental variable; IVW, inverse-variance weighted; Wald, single-variant "
        "Wald ratio; OR, odds ratio; CI, confidence interval; FDR, false discovery rate.",
    )
    columns = [
        "Immunophenotype",
        "GWS N IVs",
        "Method",
        "GWS OR (95% CI)",
        "GWS P",
        "GWS FDR",
        "Primary rank",
        "Primary P",
        "Primary FDR",
    ]
    table = doc.add_table(rows=1, cols=len(columns))
    table.style = "Table Grid"
    for i, col in enumerate(columns):
        run = table.rows[0].cells[i].paragraphs[0].add_run(col)
        run.bold = True
    for row in rows:
        cells = table.add_row().cells
        values = [
            row["trait"],
            str(int(float(row["gws_nIV"]))),
            row["gws_method"],
            f"{float(row['gws_OR']):.3f} ({float(row['gws_OR_L']):.3f}-{float(row['gws_OR_U']):.3f})",
            f"{float(row['gws_p']):.2e}",
            f"{float(row['gws_FDR']):.3f}".lstrip("0"),
            f"{int(float(row['main_rank']))}/731",
            f"{float(row['main_p']):.2f}".lstrip("0"),
            f"{float(row['main_FDR']):.2f}".lstrip("0"),
        ]
        for i, value in enumerate(values):
            cells[i].text = value
    for table_row in table.rows:
        for cell in table_row.cells:
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(7)


def add_figure_legend(doc: Document, text: str) -> None:
    p = add_paragraph(doc, text)
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT


def add_frontmatter(doc: Document, md: str) -> None:
    title = md.splitlines()[0].lstrip("# ").strip()
    p = add_paragraph(doc, title)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in p.runs:
        run.bold = True
        run.font.size = Pt(15)

    p = add_paragraph(doc, "Wei Li, Yanfu Wang*, Fang Liu and Yuanzheng Li")
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_paragraph(doc, "Department of Cardiovascular Medicine, Aviation General Hospital, Beijing 100012, China")
    add_paragraph(doc, "*Correspondence: Yanfu Wang, Department of Cardiovascular Medicine, Aviation General Hospital, Beijing 100012, China; Email: wangfy328@email.hkzyy.com.cn")

    abstract = re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", md, re.S).group(1).strip()
    add_heading(doc, "Abstract")
    for paragraph in re.split(r"\n\s*\n", abstract):
        add_paragraph(doc, paragraph.strip())
    keywords = re.search(r"\*\*Keywords:\*\*\s*(.*)", md).group(1).strip()
    add_paragraph(doc, "Keywords: " + keywords)


# Display items are inserted by matching a trigger phrase in the manuscript text. If the text
# is reworded the match silently stops firing and the item vanishes from the DOCX, so every
# trigger is registered here and verified after the body is written.
DISPLAY_TRIGGERS = {
    "Figure 1": "The analytical workflow is summarized in Figure 1",
    "Figure 2": "did not identify any FDR-significant association for either conduction outcome (Figure 2)",
    "Table 1": "Table 1 summarizes the contrast between the two thresholds",
    "Figure 3": "(Figure 3A)",
}


def add_body(doc: Document, md: str) -> None:
    emitted: set[str] = set()
    before_refs, _ = split_refs(md)
    body = before_refs.split("## 1 Introduction", 1)[1].split("## Author contributions", 1)[0]
    for raw_line in ("## 1 Introduction\n" + body).splitlines():
        line = raw_line.rstrip()
        if not line or line == "---":
            continue
        heading = re.match(r"^(#{2,4})\s+(.*)", line)
        if heading:
            level = max(1, min(3, len(heading.group(1)) - 1))
            add_heading(doc, heading.group(2), level)
            continue
        add_paragraph(doc, line)
        if DISPLAY_TRIGGERS["Figure 1"] in line:
            emitted.add("Figure 1")
            add_figure_legend(
                doc,
                "Figure 1. Study flowchart of the bidirectional Mendelian randomization analysis, "
                "including the genome-wide instrument-threshold sensitivity analysis.",
            )
        if DISPLAY_TRIGGERS["Figure 2"] in line:
            emitted.add("Figure 2")
            add_figure_legend(
                doc,
                "Figure 2. Volcano plots from the primary P < 1 x 10^-5 instrument-threshold analysis of 731 immune-cell phenotypes and cardiac conduction disorders or atrioventricular block. No phenotype survived false-discovery-rate correction in the primary screen. Green labeled points indicate the primary-screen diagnostic subset, the five directionally concordant phenotypes that passed the sensitivity-stability criteria of Section 2.5; none was FDR-significant. Gold diamonds mark the genome-wide-threshold FDR signals, the five phenotype-level signals identified only under the stricter P < 5 x 10^-8 threshold. The two sets are non-overlapping and share no phenotype. The genome-wide-threshold FDR signals are plotted in both panels to show where they fall in the primary screen; they reached FDR significance for cardiac conduction disorders only.",
            )
        if DISPLAY_TRIGGERS["Table 1"] in line:
            emitted.add("Table 1")
            add_table1(doc)
        if DISPLAY_TRIGGERS["Figure 3"] in line:
            emitted.add("Figure 3")
            add_figure_legend(
                doc,
                "Figure 3. Threshold dependence of the phenotype-wide results for cardiac conduction disorders. (A) Rank-rank plot of the phenotype-wide P-value ordering under the primary (P < 1 x 10^-5) and genome-wide (P < 5 x 10^-8) instrument thresholds, across the 607 phenotype records analyzable under both (Spearman rho = .11). Gold points mark the genome-wide-threshold FDR signals, none of which was even nominally significant in the primary screen. The LOESS curve is shown for visual orientation only. (B) Distribution of the number of instrumental variants per phenotype under the two thresholds, shown as a percentage of analyzable phenotypes because the denominators differ (607 at the genome-wide threshold, 731 at the primary threshold). At the genome-wide threshold, 26.5% of phenotypes were instrumented by a single variant and 68.0% by three or fewer; the genome-wide-threshold FDR signals rested on one to four variants. At the primary threshold, no phenotype was instrumented by a single variant, only 3 of 731 had three or fewer, and the median instrument count was 22. Post-hoc diagnostics for those signals, including the share of inverse-variance weight carried by a single rare variant and the effect of excluding instruments with a minor allele frequency below 1%, are given in Supplementary Tables S14 and S15.",
            )
    missing = sorted(set(DISPLAY_TRIGGERS) - emitted)
    if missing:
        raise RuntimeError(
            "display items were not inserted because their trigger phrase no longer appears "
            f"in the manuscript: {', '.join(missing)}. Update DISPLAY_TRIGGERS to match the text."
        )


def add_backmatter(doc: Document, md: str) -> None:
    headers = [
        "Author contributions",
        "Acknowledgments",
        "Funding",
        "Ethics approval",
        "Conflicts of interest",
        "Data availability",
        "Supplemental Digital Content",
        "Abbreviations",
    ]
    for header in headers:
        content = section_text(md, header)
        if content:
            add_heading(doc, header)
            # Author contributions lists one CRediT role per line; keep the breaks.
            for line in content.splitlines():
                if line.strip():
                    add_paragraph(doc, line.strip())

    _, refs = split_refs(md)
    add_heading(doc, "References")
    for line in refs.splitlines():
        if line.strip():
            add_paragraph(doc, line.strip())


def build_manuscript() -> None:
    md = MANUSCRIPT_MD.read_text(encoding="utf-8")
    doc = setup_document(line_numbers=True)
    add_frontmatter(doc, md)
    add_body(doc, md)
    add_backmatter(doc, md)
    PRIMARY.mkdir(parents=True, exist_ok=True)
    doc.save(PRIMARY / "Manuscript_Medicine.docx")


def build_cover_letter() -> None:
    md = COVER_MD.read_text(encoding="utf-8")
    doc = setup_document()
    # Keep the letter on one page. The default template adds 10pt after every
    # paragraph at 1.15 line spacing, and blank source lines were emitted as
    # empty paragraphs carrying the same spacing, which together pushed the
    # letter onto a second page. Paragraph separation now comes from space_after.
    normal = doc.styles["Normal"].paragraph_format
    normal.space_after = Pt(5)
    normal.line_spacing = 1.0
    for line in md.splitlines():
        if not line.strip():
            continue
        if line.startswith("# "):
            p = add_paragraph(doc, line[2:].strip())
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in p.runs:
                run.bold = True
            continue
        add_paragraph(doc, line.rstrip())
    doc.save(COVER / "Cover_Letter_Medicine.docx")


def build_strobe() -> None:
    items = [
        ("Title and abstract", "1",
         "Indicate Mendelian randomization (MR) as the study's design in the title and/or the abstract if that is a main purpose of the study.",
         "Title; Abstract (title and Methods paragraph)"),
        ("Introduction", "2",
         "Explain the scientific background and rationale for the reported study. What is the exposure? Is a causal effect plausible? Why use MR?",
         "Introduction, paragraphs 1-3"),
        ("Introduction", "3",
         "State specific objectives clearly, including prespecified causal hypotheses (if any).",
         "Introduction, paragraphs 4-6"),
        ("Methods", "4",
         "Study design and data sources. (4a) Study design and underlying population. (4b) Eligibility criteria and the sources and methods of selection of participants. (4c) Measurement, quality control and selection of genetic variants. (4d) Methods of assessment and diagnostic criteria for the exposures, outcomes and other relevant variables. (4e) Ethics committee approval and participant informed consent.",
         "4a: Methods 2.1. 4b and 4d: Methods 2.2, 2.3 (outcome case and control definitions, including the contributing diagnostic codes, are referenced to the FinnGen endpoint browser and summarized in Supplementary Table S17). 4c: Methods 2.4. 4e: Methods 2.1 and the Ethics approval statement."),
        ("Methods", "5",
         "Explicitly state the three core instrumental-variable assumptions (relevance, independence/exchangeability, exclusion restriction) and the conditions under which MR estimates a causal effect.",
         "Methods 2.1 (relevance, independence and exclusion-restriction assumptions, followed by the further identifying conditions - linearity and effect homogeneity, or a local interpretation in their absence - under which the ratio estimate recovers a causal effect)"),
        ("Methods", "6",
         "Statistical methods, main analysis. (6a) Handling of quantitative variables (scale, units, transformation). (6b) Handling of genetic variants and, if applicable, selection of their weights. (6c) MR estimator(s) and their assumptions. (6d) Handling of missing data. (6e) Multiple testing / multiplicity.",
         "6a: Methods 2.2 (traits inverse-normal transformed in the source GWAS; estimates per 1 standard deviation). 6b: Methods 2.4, and Methods 2.8 for the post-hoc minor-allele-frequency filter. 6c: Methods 2.5. 6d: Methods 2.4 (variants absent from the outcome dataset were dropped; no linkage-disequilibrium proxies were substituted). 6e: Methods 2.6."),
        ("Methods", "7",
         "Describe any methods or prior knowledge used to assess the instrumental-variable assumptions or to justify their validity.",
         "Methods 2.4 (F-statistic filtering, winner's curse), 2.5 (MR-Egger intercept, Cochran's Q, MR-PRESSO, leave-one-out; MR-PRESSO computed with the MRPRESSO package and used as a global heterogeneity diagnostic only), 2.6 (Steiger directionality), 2.8 (post-hoc variant-level diagnostics)"),
        ("Methods", "8",
         "Describe any sensitivity analyses or additional analyses performed (for example, comparison of different MR estimators, MR-Egger, weighted median, MR-PRESSO, assessment of sample overlap, negative controls).",
         "Methods 2.5, 2.6 (reverse-direction MR and multiple testing), 2.7 (targeted replication), 2.8 (post-hoc rare-variant diagnostics and minor-allele-frequency sensitivity analysis), and the genome-wide instrument-threshold sensitivity analysis in Methods 2.6"),
        ("Results", "9",
         "Descriptive data. (9a) Numbers at each stage of the study and reasons for exclusion; summary statistics for the genetic variants, the exposure and the outcome. (9b) If the data sources include meta-analyses of previous studies, assessments of heterogeneity across those studies. (9c) Data sources, participants, and whether exposure and outcome data come from the same or different (non-overlapping) samples.",
         "9a: Figure 1; Methods 2.3, 2.4; Results 3.4; Supplementary Tables S6 and S11. 9b: Methods 2.1 (neither source dataset is a meta-analysis of separately published studies, so no across-study heterogeneity assessment applies; instrument-level heterogeneity is reported in Methods 2.5 and Results 3.2, 3.6). 9c: Methods 2.1, 2.2, 2.3."),
        ("Results", "10",
         "Main results. (10a) Associations between genetic variant and exposure, and between genetic variant and outcome. (10b) MR estimates of the exposure-outcome association with measures of uncertainty, on an interpretable scale. (10c) Where relevant, translation of relative into absolute risk. (10d) Plots to visualise results.",
         "10a: Supplementary Table S6. 10b: Results 3.1-3.6; Table 1; Supplementary Tables S1, S2, S11, S13. 10d: Figures 2 and 3; Supplementary Figures S1 and S2. 10c: not applicable, because the primary screen returned no false-discovery-rate-significant association and therefore no relative risk to translate into absolute risk."),
        ("Results", "11",
         "Assessment of assumptions. (11a) Report the assessment of the validity of the assumptions. (11b) Report any additional statistics, such as assessments of heterogeneity across instruments (I-squared or Cochran's Q).",
         "11a: Results 3.2, 3.3, 3.4, 3.6; Supplementary Tables S4, S8 and S14. 11b: Results 3.2, 3.4, 3.6 (Cochran's Q, MR-Egger intercept); Methods 2.5 (I-squared GX); Supplementary Tables S4, S8 and S14."),
        ("Results", "12",
         "Sensitivity analyses and additional analyses. (12a) Robustness of the main results to violations of the assumptions. (12b) Assessment of the direction of the causal effect. (12c) Additional analyses. (12d) Where relevant, comparison with estimates from non-MR analyses. (12e) Indications of any other potential sources of bias, such as selection bias, sample overlap or winner's curse.",
         "12a: Results 3.2, 3.4, 3.5; Supplementary Tables S4, S11, S12, S14, S15. 12b: Results 3.3 (reverse-direction MR and Steiger filtering); the scope of that analysis is stated at the end of Results 3.4. 12c: Results 3.4-3.6. 12e: Methods 2.1 (no participant overlap between samples) and 2.4 (winner's curse); Results 3.4 and 3.5 and Supplementary Tables S14-S15 (weight carried by a single rare-variant instrument); Discussion, limitations. 12d: observational associations are summarized in the Introduction, but no formal quantitative comparison with non-MR estimates was performed."),
        ("Discussion", "13",
         "Summarise key results with reference to the study objectives.",
         "Discussion, paragraphs 1-2; Conclusion"),
        ("Discussion", "14",
         "Discuss limitations, taking into account the validity of the instrumental-variable assumptions, other sources of potential bias, and imprecision. Discuss both the direction and the magnitude of any potential bias.",
         "Discussion, paragraphs 12-18 (limitations First to Ninth); Results 3.1 and Discussion paragraph 3 (power and detectable effect sizes); Supplementary Table S18 (phenotype-level instrument strength and power)"),
        ("Discussion", "15",
         "Interpretation. (15a) A cautious overall interpretation in the context of the limitations and in comparison with other studies. (15b) Underlying biological mechanisms that could drive a potential causal effect, and whether the gene-environment equivalence assumption is reasonable. (15c) Whether the results have clinical or public-policy relevance, and to what extent they inform effect sizes of possible interventions.",
         "15a: Discussion, paragraphs 2-10. 15b: Discussion, paragraphs 6 and 9 (tissue-local and stage-specific mechanisms; gene-environment equivalence). 15c: Discussion, paragraph 10; Conclusion."),
        ("Discussion", "16",
         "Discuss the generalizability (external validity) of the results, including the relevance of the population(s) and ancestry to which the findings apply.",
         "Methods 2.1; Discussion, paragraph 13 (third limitation, founder-population transportability) and the closing sentence of paragraph 18"),
        ("Other information", "17",
         "Describe sources of funding and the role of funders in the present study and, if applicable, for the databases and original studies on which the present study is based.",
         "Funding"),
        ("Other information", "18",
         "Provide the data used to perform all analyses, or report where and how the data can be accessed, and reference these sources. Provide the statistical code needed to reproduce the results, or report where and how it can be accessed.",
         "Data availability; Supplemental Digital Content"),
        ("Other information", "19",
         "All authors should declare all potential conflicts of interest.",
         "Conflicts of interest"),
        ("Other information", "20",
         "Where applicable, report other study registration, protocol availability, or supplementary reporting such as a completed STROBE-MR checklist.",
         "Methods 2.1 states that the analysis was not prospectively registered and that no protocol was deposited in a public registry. Methods 2.8 and the seventh limitation identify the rare-variant diagnostics as post-hoc. This completed STROBE-MR checklist is provided as Supplemental Digital Content."),
    ]
    doc = setup_document()
    add_heading(doc, "STROBE-MR Checklist")
    add_paragraph(doc, "Manuscript: " + TITLE)
    add_paragraph(doc, "Reporting follows the STROBE-MR guideline, which extends STROBE to Mendelian randomization studies (Skrivankova VW, Richmond RC, Woolf BAR, et al. JAMA 2021;326:1614-1621). The 20 items below follow the section grouping of the STROBE-MR statement. Page numbers are to be completed after typesetting.")
    table = doc.add_table(rows=1, cols=4)
    table.style = "Table Grid"
    for i, heading in enumerate(["Section", "Item", "Checklist item", "Reported in (section)"]):
        run = table.rows[0].cells[i].paragraphs[0].add_run(heading)
        run.bold = True
    last_section = None
    for section, num, item, where in items:
        cells = table.add_row().cells
        cells[0].text = section if section != last_section else ""
        last_section = section
        cells[1].text = num
        cells[2].text = item
        cells[3].text = where
    doc.save(SUPP / "STROBE-MR_checklist.docx")


def main() -> None:
    build_manuscript()
    build_cover_letter()
    build_strobe()
    print("Generated Medicine submission documents in", ROOT)


if __name__ == "__main__":
    main()
