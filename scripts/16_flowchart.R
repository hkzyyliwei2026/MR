# =============================================================
# 16: Figure 1 study flowchart
# Outputs: results/figures/Fig1_flowchart.pdf/png/svg/tiff
# =============================================================
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
FIG  <- file.path(PROJ, "results/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

open_svg <- function(path, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(path, width = width, height = height)
  } else {
    svg(path, width = width, height = height, family = "Arial")
  }
}

draw <- function() {
  BLUE <- "#1F5A85"; BLUE_FILL <- "#EEF5FA"
  ORANGE <- "#B96A2C"; ORANGE_FILL <- "#FCF1E8"
  GREY <- "#4F4F4F"; GRID <- "#E8ECEF"
  par(mar = c(0.8, 1.0, 1.2, 1.0), family = "Arial", xpd = NA)
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "")
  # The figure title and the rule beneath it are deliberately absent: Medicine requires that
  # digital art carry no embedded figure title or legend text. The title now lives only in the
  # Figure 1 legend in the manuscript. Panel labels and box text stay, since they are what
  # makes the flowchart self-explanatory.

  box_ <- function(x, y, w, h, txt, fill = BLUE_FILL, border = BLUE, cex = 0.78, font = 1) {
    rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = border, lwd = 1.4)
    text(x, y, txt, cex = cex, font = font, col = "#222222")
  }
  arr <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, length = 0.085, lwd = 1.4, col = GREY)
  section <- function(x, y, label) text(x, y, label, cex = 0.68, font = 2, col = GREY)

  section(10, 88.5, "Input")
  box_(31, 85, 54, 10, "Exposure\n731 immune cell phenotypes\nImmune-cell GWAS, n = 3,757")
  box_(79, 85, 34, 10, "Outcomes: FinnGen R11\nConduction disorders: 12,371/342,690\nAV block: 6,935/342,690", cex = 0.68)

  section(10, 69.5, "MR design")
  box_(31, 69, 54, 12, "Instrument selection\nP < 1 x 10^-5; r^2 < 0.001; 10 Mb; F > 10\nAllele harmonization and palindromic SNP filtering")
  box_(31, 50, 54, 10, "Two-sample MR\nIVW primary; Wald ratio for single-IV traits")
  box_(79, 50, 34, 12, "Instrument-threshold sensitivity\nP < 5 x 10^-8; 607/731 analyzable\n5 FDR signals (conduction); 0 (AV block)", fill = ORANGE_FILL, border = ORANGE, cex = 0.60)

  section(10, 34.5, "Findings")
  box_(31, 34, 58, 12, "Primary result\n0/731 FDR-significant phenotypes\n33 nominal per outcome; 15 directionally concordant")
  box_(31, 16, 54, 10, "Sensitivity analyses\nWeighted median, MR-Egger, Cochran Q, LOO\n5 sensitivity-stable nominal phenotypes", fill = ORANGE_FILL, border = ORANGE, cex = 0.68)
  box_(79, 16, 34, 10, "Reverse MR\nConduction to immune\nNo support for reverse causation", fill = ORANGE_FILL, border = ORANGE, cex = 0.68)

  arr(58, 85, 62, 85)
  arr(31, 80, 31, 75)
  arr(31, 63, 31, 55)
  arr(31, 45, 31, 40)
  arr(31, 28, 31, 21)
  arr(58, 50, 62, 50)
  arr(58, 16, 62, 16)
}

cairo_pdf(file.path(FIG, "Fig1_flowchart.pdf"), width = 8.0, height = 7.0, family = "Arial")
draw()
dev.off()
open_svg(file.path(FIG, "Fig1_flowchart.svg"), width = 8.0, height = 7.0)
draw()
dev.off()
png(file.path(FIG, "Fig1_flowchart.png"), width = 4800, height = 4200, res = 600, pointsize = 13)
draw()
dev.off()
tiff(file.path(FIG, "Fig1_flowchart.tiff"), width = 8.0, height = 7.0, units = "in", res = 600, pointsize = 13, compression = "lzw")
draw()
dev.off()
cat("Fig1 output complete\n")
