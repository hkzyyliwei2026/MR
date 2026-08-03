# =============================================================
# 30: Recover component-endpoint case counts from the summary-statistics allele frequencies
#
# The downloaded FinnGen release 11 summary files carry no case-count field, and the release 11
# Risteys pages are the authoritative source. For the two component endpoints used in Section 3.4
# the counts can nevertheless be recovered from the files themselves, because each variant is
# reported with an overall alternate-allele frequency and with the case-only and control-only
# frequencies. The overall frequency is the case/control weighted mean, so
#
#     af_alt = (n_case * af_case + n_control * af_control) / (n_case + n_control)
#  => n_case / n_control = (af_alt - af_control) / (af_case - af_alt)
#
# The ratio is estimated as the median across common variants, which is insensitive to the
# rounding of individual frequencies. Multiplying by the control count gives the case count.
#
# The method is checked against the two endpoints whose case counts are known and reported in
# Section 2.3: it recovers 12,371 and 6,935 to within 0.1 of a case. It is applied only to the
# endpoints that share the same control group; the atrial-fibrillation endpoint of Section 3.4
# does not, so only its case-to-control ratio is recoverable and it is reported as such.
#
# Input : <finngen_dir>/finngen_R11_I9_{CONDUCTIO,AVBLOCK,LBBB,RBBB,AF}.gz
# Output: results/tables/S23_endpoint_case_counts.csv
# Usage : Rscript scripts/30_endpoint_case_counts.R [finngen_dir]
# =============================================================
suppressMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")
args <- commandArgs(TRUE)
FG <- if (length(args)) args[1] else file.path(PROJ, "data/outcome")

CONTROLS <- 342690L   # shared by the conduction-disorder endpoint family, Section 2.3
KNOWN <- c(I9_CONDUCTIO = 12371L, I9_AVBLOCK = 6935L)
SHARED_CONTROLS <- c("I9_CONDUCTIO", "I9_AVBLOCK", "I9_LBBB", "I9_RBBB")

ratio_for <- function(tag){
  f <- file.path(FG, sprintf("finngen_R11_%s.gz", tag))
  if (!file.exists(f)) return(NA_real_)
  d <- fread(f, select = c("af_alt", "af_alt_cases", "af_alt_controls"), nrows = 400000L)
  setnames(d, c("a", "c", "o"))
  d <- d[!is.na(a) & !is.na(c) & !is.na(o) & a > 0.05 & a < 0.95 & abs(c - a) > 1e-9]
  r <- (d$a - d$o) / (d$c - d$a)
  median(r[r > 0 & r < 1])
}

rows <- rbindlist(lapply(c(SHARED_CONTROLS, "I9_AF"), function(tag){
  r <- ratio_for(tag)
  shared <- tag %in% SHARED_CONTROLS
  data.table(
    Endpoint = tag,
    `Case-to-control ratio` = round(r, 8),
    `Controls` = if (shared) CONTROLS else NA_integer_,
    `Cases (recovered)` = if (shared) round(r * CONTROLS) else NA_real_,
    `Cases (reported in Section 2.3)` = if (tag %in% names(KNOWN)) KNOWN[[tag]] else NA_integer_)
}))
rows[, `Recovery error (cases)` := `Cases (recovered)` - `Cases (reported in Section 2.3)`]
fwrite(rows, file.path(RES, "S23_endpoint_case_counts.csv"))
print(rows)
cat("\nwrote", file.path(RES, "S23_endpoint_case_counts.csv"), "\n")
