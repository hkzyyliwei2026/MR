#!/usr/bin/env Rscript
# =============================================================================
# 21: Nominal-finding excess and cross-endpoint concordance across the two
#     instrument thresholds (reported in Results, "Quantifying instrument-
#     threshold instability", and in the Methods note on the FDR denominator).
#
# Inputs : results/tables/MR_immune_CONDUCTIO.csv   (731 traits, P < 1e-5)
#          results/tables/MR_immune_AVBLOCK.csv     (731 traits, P < 1e-5)
#          derived/threshold_instability.csv     (607 traits x 2 outcomes,
#                                                 both thresholds side by side)
# Output : results/tables/nominal_excess_concordance.txt
#
# No MR is re-run here; this script only re-derives counting statistics from
# the result tables so that every number quoted in the text is reproducible.
# =============================================================================
suppressMessages(library(data.table))

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results", "tables")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

ALPHA <- 0.05
out   <- file.path(RES, "nominal_excess_concordance.txt")
con   <- file(out, open = "wt")
say   <- function(...) { msg <- paste0(...); cat(msg, "\n"); writeLines(msg, con) }

# ---- 1. Nominal findings vs chance, primary threshold (all 731 traits) ------
main <- list(
  `Cardiac conduction disorders` = fread(file.path(PROJ, "results/tables/MR_immune_CONDUCTIO.csv")),
  `Atrioventricular block`       = fread(file.path(PROJ, "results/tables/MR_immune_AVBLOCK.csv")))

say("== Nominal findings vs chance ==")
say(sprintf("%-30s %-14s %5s %5s %8s %10s", "outcome", "threshold", "k", "n", "expected", "binom P"))

for (nm in names(main)) {
  d <- main[[nm]]
  k <- sum(d$p < ALPHA, na.rm = TRUE); n <- nrow(d)
  bt <- binom.test(k, n, ALPHA, alternative = "greater")
  say(sprintf("%-30s %-14s %5d %5d %8.1f %10.3g", nm, "P < 1e-5", k, n, ALPHA * n, bt$p.value))
}

ti <- fread(file.path(PROJ, "derived/threshold_instability.csv"))
# Both thresholds restricted to the same 607 traits, so only the threshold varies.
for (thr in c("p_main", "p_gws")) {
  lab <- if (thr == "p_main") "P < 1e-5 (607)" else "P < 5e-8 (607)"
  for (nm in unique(ti$outcome)) {
    d <- ti[outcome == nm & !is.na(get(thr))]
    k <- sum(d[[thr]] < ALPHA); n <- nrow(d)
    bt <- binom.test(k, n, ALPHA, alternative = "greater")
    say(sprintf("%-30s %-14s %5d %5d %8.1f %10.3g", nm, lab, k, n, ALPHA * n, bt$p.value))
  }
}

# ---- 2. Cross-endpoint concordance, restricted to the shared 607 traits -----
# Atrioventricular block cases are nested within the conduction-disorder
# definition and both endpoints share controls, so overlap is expected at
# either threshold; the comparison of interest is primary vs genome-wide.
say("")
say("== Cross-endpoint concordance (same trait set under both thresholds) ==")
say(sprintf("%-14s %5s %5s %8s %10s %10s %10s", "threshold", "nCD", "nAV", "overlap", "expected", "fold", "phyper P"))

wide <- dcast(ti, id ~ outcome, value.var = c("p_main", "p_gws"))
cd <- grep("Cardiac", names(wide), value = TRUE)
av <- grep("Atrioventricular", names(wide), value = TRUE)

concord <- function(pcd, pav, label) {
  keep <- !is.na(pcd) & !is.na(pav)
  pcd <- pcd[keep]; pav <- pav[keep]; N <- length(pcd)
  a <- pcd < ALPHA; b <- pav < ALPHA
  ov <- sum(a & b); expd <- sum(a) * sum(b) / N
  # P(overlap >= observed) given the two marginal counts
  ph <- phyper(ov - 1, sum(b), N - sum(b), sum(a), lower.tail = FALSE)
  say(sprintf("%-14s %5d %5d %8d %10.1f %10.2f %10.3g",
              label, sum(a), sum(b), ov, expd, ov / expd, ph))
}
concord(wide[[grep("p_main", cd, value = TRUE)]], wide[[grep("p_main", av, value = TRUE)]], "P < 1e-5")
concord(wide[[grep("p_gws",  cd, value = TRUE)]], wide[[grep("p_gws",  av, value = TRUE)]], "P < 5e-8")

# ---- 3. FDR denominator for the sensitivity analysis -----------------------
# The genome-wide-threshold screen was corrected across the 607 traits that
# remained analyzable; correcting across 731 is shown for comparison.
say("")
say("== Sensitivity-analysis FDR: effect of the correction denominator ==")
say(sprintf("%-30s %12s %10s %8s", "outcome", "denominator", "min FDR", "n < .05"))
for (nm in unique(ti$outcome)) {
  p <- sort(ti[outcome == nm & !is.na(p_gws), p_gws]); m <- length(p)
  for (den in c(m, 731L)) {
    q <- cummin(rev(pmin(1, p * den / seq_along(p))))
    q <- rev(q)
    say(sprintf("%-30s %12d %10.4f %8d", nm, den, min(q), sum(q < ALPHA)))
  }
}

close(con)
cat("\nwrote", out, "\n")
