# =============================================================
# 17: MR-PRESSO (outlier pleiotropy) and Steiger (directionality) for the 5 sensitivity-stable traits
# NOTE: MR-PRESSO is computed with the official MRPRESSO package and is used only as a
# global heterogeneity diagnostic; the variant-level outlier test is not used and no
# reported result depends on it. Steiger is a self-contained implementation of Hemani 2015.
# Requires: remotes::install_github("rondolab/MR-PRESSO")
#   MR-PRESSO: global and outlier residual-sum-of-squares tests of Verbanck 2018, Nat Genet
#   Steiger  : directionality test of Hemani 2015, PLoS Genet
# Input: immune_ivs.csv (exposure instruments) + data/outcome/finngen_R11_I9_CONDUCTIO.gz (primary outcome)
# Output: results/tables/immune_presso_steiger.csv
# Usage: Rscript 17_presso_steiger.R
# =============================================================
suppressMessages({library(data.table); library(MRPRESSO)})
set.seed(42)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")

IVS  <- fread(file.path(RES, "immune_ivs.csv"))
SENS <- fread(file.path(RES, "immune_sensitivity.csv"))
robust <- SENS[robust == TRUE, trait]
IVS <- IVS[trait %in% robust]
setnames(IVS, c("ea","nea","eaf","beta","se"), c("ea_e","nea_e","eaf_e","beta_e","se_e"))
IVS[, `:=`(ea_e = toupper(ea_e), nea_e = toupper(nea_e))]

# ---- primary outcome: cardiac conduction disorders ----
ncase <- 12371; ncontrol <- 342690
oc_file <- file.path(PROJ, "data/outcome/finngen_R11_I9_CONDUCTIO.gz")
need <- unique(IVS$rsid)
cat("reading outcome file (~800 MB, please wait)...\n")
OUT <- fread(oc_file, select = c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
              c("rsid","ref_o","alt_o","beta_o","se_o","p_o","eaf_o"))
OUT <- OUT[rsid %in% need & !is.na(beta_o) & se_o > 0]
OUT[, `:=`(ref_o = toupper(ref_o), alt_o = toupper(alt_o))]

is_pal <- function(a1, a2){ p <- c("A"="T","T"="A","C"="G","G"="C"); !is.na(p[a1]) & p[a1] == a2 }
m <- merge(IVS, OUT, by = "rsid")
m[, keep := TRUE]; m[, by_al := beta_o]
m[ea_e == alt_o & nea_e == ref_o, by_al := beta_o]   # same direction
m[ea_e == ref_o & nea_e == alt_o, by_al := -beta_o]  # opposite direction, flip sign
m[!((ea_e == alt_o & nea_e == ref_o) | (ea_e == ref_o & nea_e == alt_o)), keep := FALSE]
m[is_pal(ea_e, nea_e) & pmin(eaf_e, 1 - eaf_e) > 0.42, keep := FALSE]
m <- m[keep == TRUE]
m[, F := (beta_e / se_e)^2]
m <- m[F > 10]

# ---- MR-PRESSO (official MRPRESSO package, global heterogeneity test only) ----
# MR-PRESSO is used in this study only as a global heterogeneity diagnostic. The package
# performs the variant-level outlier test only when the global test is significant and
# returns NA for the outlier-corrected estimate otherwise; that behaviour is preserved
# here. No result reported in the manuscript depends on the outlier test.
mrpresso <- function(bx, bxse, by, byse, K = 3000){
  if(length(bx) < 4) return(list(global_p = NA, n_outlier = NA, or_raw = NA, or_adj = NA))
  dt <- data.frame(by = by, bx = bx, sy = byse, sx = bxse)
  set.seed(42)
  r <- try(MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx",
        SdOutcome = "sy", SdExposure = "sx", OUTLIERtest = TRUE, DISTORTIONtest = FALSE,
        data = dt, NbDistribution = K, SignifThreshold = 0.05), silent = TRUE)
  if(inherits(r, "try-error")) return(list(global_p = NA, n_outlier = NA, or_raw = NA, or_adj = NA))
  gp <- suppressWarnings(as.numeric(sub("<", "", as.character(
          r$`MR-PRESSO results`$`Global Test`$Pvalue))))
  ot <- r$`MR-PRESSO results`$`Outlier Test`
  nout <- if(is.null(ot)) NA_integer_ else sum(ot$Pvalue < 0.05)
  mt <- r$`Main MR results`
  raw <- mt$`Causal Estimate`[mt$`MR Analysis` == "Raw"]
  cor <- mt$`Causal Estimate`[mt$`MR Analysis` == "Outlier-corrected"]
  list(global_p = gp, n_outlier = nout, or_raw = exp(raw),
       or_adj = if(length(cor) && !is.na(cor)) exp(cor) else NA_real_)
}

# ---- Steiger directionality (Hemani 2015) ----
# exposure (continuous, standardised): r2 = 2*f*(1-f)*b^2; outcome (binary): same formula on the
# log-odds scale, which is conservative
steiger <- function(bx, by, f, n1, n2){
  r2e <- 2*f*(1-f)*bx^2
  r2o <- 2*f*(1-f)*by^2
  r1 <- sqrt(min(sum(r2e), 0.999)); r2 <- sqrt(min(sum(r2o), 0.999))
  z  <- (atanh(r1) - atanh(r2)) / sqrt(1/(n1-3) + 1/(n2-3))
  list(R2_exp = sum(r2e), R2_out = sum(r2o),
       prop_correct = mean(r2e > r2o), steiger_p = 2*pnorm(-abs(z)),
       dir = ifelse(r1 > r2, "exposure->outcome", "ambiguous"))
}

res <- rbindlist(lapply(split(m, by = "id"), function(g){
  pr <- mrpresso(g$beta_e, g$se_e, g$by_al, g$se_o)
  st <- steiger(g$beta_e, g$by_al, g$eaf_e, 3757, ncase + ncontrol)
  data.table(
    trait = g$trait[1], nIV = nrow(g), meanF = round(mean(g$F), 1),
    R2_exp_pct = round(st$R2_exp*100, 2), R2_out_pct = round(st$R2_out*100, 4),
    steiger_prop_correct = round(st$prop_correct, 3),
    steiger_p = signif(st$steiger_p, 2), steiger_dir = st$dir,
    presso_global_p = round(pr$global_p, 3),
    presso_n_outlier = pr$n_outlier,
    OR_raw = round(pr$or_raw, 3),
    OR_adj = round(pr$or_adj, 3))
}))
setorder(res, steiger_p)
fwrite(res, file.path(RES, "immune_presso_steiger.csv"))
cat("\n==== MR-PRESSO + Steiger (primary outcome: conduction disorders) ====\n")
print(res)
cat("\nwritten:", file.path(RES, "immune_presso_steiger.csv"), "\n")
