# =============================================================
# 28: Specificity control for the lymphocyte-count re-analysis (Supplementary Table S22)
#
# The conduction endpoints are defined on registry codes, and atrial fibrillation is treated with
# rate-control drugs, antiarrhythmic therapy, catheter ablation, atrioventricular-node ablation and
# pacemaker implantation, each of which can generate an I44 or I45 code. An immune-cell trait that
# acts on atrial fibrillation could therefore appear to act on conduction disease. This script tests
# the same lymphocyte-count instrument set against atrial fibrillation itself: if the association
# with atrial fibrillation is of similar size or larger, the conduction finding is not specific.
#
# The exposure instruments are those written by 22_chen_reanalysis.R after clumping (the file
# exposure_clumped_alleles.csv in its output directory); harmonisation follows the same rules.
#
# Input : <chen_outdir>/exposure_clumped_alleles.csv ; <finngen_dir>/finngen_R11_I9_AF.gz
#         and the two conduction endpoints for comparison
# Output: results/tables/S22_af_specificity.csv
# Usage : Rscript scripts/28_af_specificity.R <chen_outdir> [finngen_dir]
# =============================================================
suppressMessages({library(data.table); library(MendelianRandomization)})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")
args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript 28_af_specificity.R <chen_outdir> [finngen_dir]")
CHEN <- args[1]
FG   <- if (length(args) > 1) args[2] else file.path(PROJ, "data/outcome")

EX <- fread(file.path(CHEN, "exposure_clumped_alleles.csv"))
EX[, `:=`(EA = toupper(EA), OA = toupper(OA))]
EX[, F := (BETA/SE)^2]
EX <- EX[F > 10]

is_palindromic <- function(a1, a2){ p <- c("A"="T","T"="A","C"="G","G"="C"); !is.na(p[a1]) & p[a1]==a2 }

run <- function(tag){
  f <- file.path(FG, sprintf("finngen_R11_%s.gz", tag))
  if (!file.exists(f)) { cat("missing:", f, "\n"); return(NULL) }
  OUT <- fread(f, select = c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
  setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
                c("rsids","REF","ALT","BETA_o","SE_o","P_o","EAF_o"))
  OUT <- OUT[!is.na(BETA_o) & SE_o > 0]
  # FinnGen packs multiple rsIDs into one field; 22_chen_reanalysis.R expands them and then
  # deduplicates by SNP. Matching on the raw field instead loses multi-rsID rows and can also
  # duplicate a variant that appears twice, so the same expansion is applied here.
  OUT <- OUT[, .(SNP = unlist(strsplit(rsids, ","))),
             by = .(REF, ALT, BETA_o, SE_o, P_o, EAF_o)][SNP %in% EX$SNP]
  OUT <- unique(OUT, by = "SNP")
  OUT[, `:=`(REF = toupper(REF), ALT = toupper(ALT))]
  m <- merge(EX, OUT, by = "SNP")
  m[, keep := TRUE]; m[, BETA_al := BETA_o]
  m[EA == REF & OA == ALT, BETA_al := -BETA_o]
  m[!((EA == ALT & OA == REF) | (EA == REF & OA == ALT)), keep := FALSE]
  m[is_palindromic(EA, OA) & pmin(EAF, 1-EAF) > 0.42, keep := FALSE]
  m <- m[keep == TRUE]
  o  <- mr_input(bx = m$BETA, bxse = m$SE, by = m$BETA_al, byse = m$SE_o, snps = m$SNP)
  iv <- mr_ivw(o, model = "random"); wm <- mr_median(o, weighting = "weighted"); eg <- mr_egger(o)
  data.table(Outcome = tag, `N instruments` = nrow(m),
             `IVW OR` = round(exp(iv@Estimate), 3),
             `IVW 95%CI lower` = round(exp(iv@CILower), 3),
             `IVW 95%CI upper` = round(exp(iv@CIUpper), 3),
             `IVW P` = signif(iv@Pvalue, 3),
             `Weighted-median OR` = round(exp(wm@Estimate), 3),
             `Weighted-median P` = signif(wm@Pvalue, 3),
             `MR-Egger OR` = round(exp(eg@Estimate), 3),
             `MR-Egger intercept P` = signif(eg@Pvalue.Int, 3),
             `Cochran Q` = round(iv@Heter.Stat[1], 1),
             `Q P` = signif(iv@Heter.Stat[2], 3))
}

res <- rbindlist(lapply(c("I9_AF", "I9_AVBLOCK", "I9_CONDUCTIO"), run))
res[, Outcome := c("I9_AF" = "Atrial fibrillation",
                   "I9_AVBLOCK" = "Atrioventricular block",
                   "I9_CONDUCTIO" = "Cardiac conduction disorders")[Outcome]]
fwrite(res, file.path(RES, "S22_af_specificity.csv"))
print(res)
cat("\nwrote", file.path(RES, "S22_af_specificity.csv"), "\n")
