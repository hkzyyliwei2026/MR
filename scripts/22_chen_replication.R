#!/usr/bin/env Rscript
# =============================================================================
# 22: Targeted replication of Chen Y et al (Front Immunol 2023;14:1041591)
#
#     Chen Y et al reported that genetically predicted lymphocyte count raises
#     the risk of atrioventricular block (OR 1.46, 95% CI 1.11-1.93, P = .0065),
#     using Blood Cell Consortium exposures and FinnGen *release 2* outcomes.
#     This script tests the same exposure against FinnGen *release 11*, which
#     has substantially more atrioventricular-block cases, and against the
#     broader conduction-disorder endpoint used in the present study.
#
# Exposure instruments
#     Source     : GWAS Catalog accession GCST90002316 (lymphocyte count,
#                  European ancestry, 524,923 participants; Chen MH et al,
#                  Cell 2020, the Blood Cell Consortium report cited by
#                  Chen Y et al as their exposure source).
#     Retrieval  : instruments are the study-level associations curated in the
#                  GWAS Catalog for this accession, i.e. the associations
#                  reported by the original investigators, read from the REST
#                  endpoint below. Repeated queries of the endpoint return an
#                  identical set, so the instrument set is reproducible.
#     Filtering  : P < 5e-8, then LD clumping at r2 < 0.001 within 10 Mb against
#                  a PLINK-format 1000 Genomes European reference panel, then
#                  F > 10 - the same rules used for the primary screen.
#     Alleles    : the curated records give the effect allele only; the other
#                  allele is taken from the reference-panel .bim file, and any
#                  variant whose effect allele matches neither reference allele
#                  is dropped.
#
# Outcomes : FinnGen R11 I9_AVBLOCK and I9_CONDUCTIO (not redistributed here;
#            available from https://www.finngen.fi/en/access_results)
#
# Usage: Rscript 22_chen_replication.R <outdir> <plink_binary> <ref_prefix> <finngen_dir>
# =============================================================================
suppressMessages({ library(data.table); library(jsonlite); library(MendelianRandomization) })

args   <- commandArgs(trailingOnly = TRUE)
OUTDIR <- args[1]; PLINK <- args[2]; REF <- args[3]; FG <- args[4]
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

ACC   <- "GCST90002316"
URL   <- paste0("https://www.ebi.ac.uk/gwas/rest/api/studies/", ACC,
                "/associations?size=2000&projection=associationByStudy")
ALPHA <- 5e-8; R2 <- 0.001; KB <- 10000; F_MIN <- 10

log <- c(); say <- function(...) { m <- paste0(...); cat(m, "\n"); log <<- c(log, m) }

# ---- 1. Curated exposure associations --------------------------------------
say("== Exposure: lymphocyte count, GWAS Catalog ", ACC, " ==")
js <- fromJSON(URL, simplifyVector = FALSE)
as <- js$`_embedded`$associations
say(sprintf("   curated associations retrieved: %d", length(as)))

rows <- rbindlist(lapply(as, function(x) {
  ra <- unlist(lapply(x$loci, function(l) l$strongestRiskAlleles), recursive = FALSE)
  if (length(ra) != 1 || is.null(x$betaNum) || is.null(x$standardError)) return(NULL)
  nm <- ra[[1]]$riskAlleleName
  m  <- regmatches(nm, regexec("^(rs[0-9]+)-([ACGT])$", nm))[[1]]
  if (length(m) != 3) return(NULL)
  b <- as.numeric(x$betaNum)
  if (!is.null(x$betaDirection) && tolower(x$betaDirection) == "decrease") b <- -b
  data.table(SNP = m[2], EA = m[3], BETA = b,
             SE = as.numeric(x$standardError), P = as.numeric(x$pvalue),
             EAF = suppressWarnings(as.numeric(ra[[1]]$riskFrequency)))
}), fill = TRUE)
rows <- unique(rows, by = "SNP")[P < ALPHA]
say(sprintf("   usable variants with P < %s: %d", ALPHA, nrow(rows)))

# ---- 2. LD clumping ---------------------------------------------------------
fwrite(rows[, .(SNP, P)], file.path(OUTDIR, "clump_in.txt"), sep = " ")
system(paste(shQuote(PLINK), "--bfile", shQuote(REF),
             "--clump", shQuote(file.path(OUTDIR, "clump_in.txt")),
             "--clump-p1", ALPHA, "--clump-r2", R2, "--clump-kb", KB,
             "--out", shQuote(file.path(OUTDIR, "clumped")), "--allow-no-sex"),
       ignore.stdout = TRUE, ignore.stderr = TRUE)
cl <- fread(file.path(OUTDIR, "clumped.clumped"))
ex <- rows[SNP %in% cl$SNP]
say(sprintf("   after LD clumping (r2 < %s, %d kb): %d", R2, KB, nrow(ex)))

# ---- 3. Other allele from the reference panel -------------------------------
bim <- fread(paste0(REF, ".bim"), col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
ex  <- merge(ex, bim[, .(SNP, CHR, BP, A1 = toupper(A1), A2 = toupper(A2))], by = "SNP")
ex[, OA := fifelse(EA == A1, A2, fifelse(EA == A2, A1, NA_character_))]
ex <- ex[!is.na(OA)]
say(sprintf("   with reference-panel alleles resolved: %d", nrow(ex)))

# ---- 4. Outcome, harmonisation, MR ------------------------------------------
comp <- c(A = "T", T = "A", C = "G", G = "C")
run <- function(tag) {
  say(""); say(sprintf("== Outcome: FinnGen R11 %s ==", tag))
  oc <- fread(cmd = paste0("gzcat ", shQuote(file.path(FG, paste0("finngen_R11_I9_", tag, ".gz")))),
              select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta"),
              col.names = c("CHR_o", "BP_o", "REF", "ALT", "rsids", "P_o", "BETA_o", "SE_o"))
  oc <- oc[, .(SNP = unlist(strsplit(rsids, ","))), by = .(REF, ALT, BETA_o, SE_o, P_o)][SNP %in% ex$SNP]
  m  <- merge(ex, unique(oc, by = "SNP"), by = "SNP")
  say(sprintf("   instruments present in outcome: %d", nrow(m)))

  m[, flip := toupper(EA) == toupper(REF) & toupper(OA) == toupper(ALT)]
  m <- m[(toupper(EA) == toupper(ALT) & toupper(OA) == toupper(REF)) | flip]
  m[flip == TRUE, BETA_o := -BETA_o]
  pal <- !is.na(comp[m$EA]) & comp[m$EA] == m$OA
  m <- m[!(pal & pmin(EAF, 1 - EAF) > 0.42)]
  m[, F := (BETA / SE)^2]; m <- m[F > F_MIN]
  say(sprintf("   after harmonisation, palindrome and F > %d filters: %d (min F %.0f, median F %.0f)",
              F_MIN, nrow(m), min(m$F), median(m$F)))

  o  <- mr_input(bx = m$BETA, bxse = m$SE, by = m$BETA_o, byse = m$SE_o, snps = m$SNP)
  iv <- mr_ivw(o, model = "random"); wm <- mr_median(o, weighting = "weighted"); eg <- mr_egger(o)
  f  <- function(b, l, u, p) sprintf("OR %.3f (95%% CI %.3f-%.3f), P = %.3g", exp(b), exp(l), exp(u), p)
  say(sprintf("   IVW (random)      %s", f(iv@Estimate, iv@CILower, iv@CIUpper, iv@Pvalue)))
  say(sprintf("   Weighted median   %s", f(wm@Estimate, wm@CILower, wm@CIUpper, wm@Pvalue)))
  say(sprintf("   MR-Egger          %s", f(eg@Estimate, eg@CILower.Est, eg@CIUpper.Est, eg@Pvalue.Est)))
  say(sprintf("   MR-Egger intercept %.5f, P = %.3g", eg@Intercept, eg@Pvalue.Int))
  say(sprintf("   Cochran's Q %.1f (df %d), P = %.3g", iv@Heter.Stat[1], nrow(m) - 1, iv@Heter.Stat[2]))
  fwrite(m, file.path(OUTDIR, paste0("harmonised_", tag, ".csv")))
  data.table(outcome = tag, nIV = nrow(m), minF = min(m$F), medF = median(m$F),
             OR = exp(iv@Estimate), L = exp(iv@CILower), U = exp(iv@CIUpper), P = iv@Pvalue,
             WM_OR = exp(wm@Estimate), WM_P = wm@Pvalue,
             EG_OR = exp(eg@Estimate), EG_P = eg@Pvalue.Est,
             EG_int = eg@Intercept, EG_int_P = eg@Pvalue.Int,
             Q = iv@Heter.Stat[1], Q_df = nrow(m) - 1, Q_P = iv@Heter.Stat[2])
}
res <- rbindlist(list(run("AVBLOCK"), run("CONDUCTIO")))

say(""); say("== Comparison with the original report ==")
say("   Chen Y et al 2023 (FinnGen release 2): OR 1.460 (1.110-1.930), P = .0065")
a <- res[outcome == "AVBLOCK"]
say(sprintf("   This analysis (FinnGen release 11):     OR %.3f (%.3f-%.3f), P = %.3g",
            a$OR, a$L, a$U, a$P))
say("   The comparison is descriptive only. FinnGen release 11 contains the release 2")
say("   participants, so the two estimates are positively correlated and a two-sample test")
say("   of their difference, which assumes independence, is not valid. An earlier version of")
say("   this script reported such a test; it has been removed.")

fwrite(res, file.path(OUTDIR, "chen_replication_summary.csv"))
writeLines(log, file.path(OUTDIR, "chen_replication_log.txt"))
