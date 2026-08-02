# =============================================================
# 20: Genome-wide-significant instrument threshold sensitivity
# Re-runs the immune-cell MR screen after retaining only exposure IVs with
# P < 5e-8. This evaluates sensitivity to the primary P < 1e-5 threshold.
# Outputs:
#   results/tables/MR_immune_GWS_CONDUCTIO.csv
#   results/tables/MR_immune_GWS_AVBLOCK.csv
#   results/tables/threshold_sensitivity_summary.csv
# =============================================================
suppressMessages({library(data.table); library(MendelianRandomization)})

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")

IVS0 <- fread(file.path(RES, "immune_ivs.csv"))
IVS0 <- IVS0[p < 5e-8]
setnames(IVS0, c("ea", "nea", "eaf", "beta", "se"),
               c("ea_e", "nea_e", "eaf_e", "beta_e", "se_e"))
IVS0[, `:=`(ea_e = toupper(ea_e), nea_e = toupper(nea_e))]

OUTCOMES <- list(
  CONDUCTIO = file.path(PROJ, "data/outcome/finngen_R11_I9_CONDUCTIO.gz"),
  AVBLOCK   = file.path(PROJ, "data/outcome/finngen_R11_I9_AVBLOCK.gz"))

is_palindromic <- function(a1, a2){
  p <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  !is.na(p[a1]) & p[a1] == a2
}

run_outcome <- function(oc_name, oc_file){
  cat(sprintf("\n==== Strict threshold outcome: %s ====\n", oc_name))
  IVS <- copy(IVS0)
  need <- unique(IVS$rsid)
  OUT <- fread(oc_file, select = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"))
  setnames(OUT, c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"),
                c("rsid", "ref_o", "alt_o", "beta_o", "se_o", "p_o", "eaf_o"))
  OUT <- OUT[rsid %in% need & !is.na(beta_o) & se_o > 0]
  OUT[, `:=`(ref_o = toupper(ref_o), alt_o = toupper(alt_o))]

  m <- merge(IVS, OUT, by = "rsid")
  m[, keep := TRUE]
  m[, beta_o_al := beta_o]
  m[ea_e == ref_o & nea_e == alt_o, beta_o_al := -beta_o]
  m[!((ea_e == alt_o & nea_e == ref_o) | (ea_e == ref_o & nea_e == alt_o)), keep := FALSE]
  m[is_palindromic(ea_e, nea_e) & pmin(eaf_e, 1 - eaf_e) > 0.42, keep := FALSE]
  m <- m[keep == TRUE]
  m[, F := (beta_e / se_e)^2]
  m <- m[F > 10]

  res <- rbindlist(lapply(split(m, by = "id"), function(g){
    nIV <- nrow(g)
    if(nIV < 1) return(NULL)
    if(nIV == 1){
      b <- g$beta_o_al / g$beta_e
      se <- abs(g$se_o / g$beta_e)
      method <- "Wald"
      p <- 2 * pnorm(-abs(b / se))
    } else {
      fit <- tryCatch(
        mr_ivw(mr_input(bx = g$beta_e, bxse = g$se_e, by = g$beta_o_al, byse = g$se_o)),
        error = function(e) NULL)
      if(is.null(fit)) return(NULL)
      b <- fit$Estimate
      se <- fit$StdError
      method <- "IVW"
      p <- fit$Pvalue
    }
    data.table(
      id = g$id[1], trait = g$trait[1], nIV = nIV, method = method,
      mean_F = mean(g$F), min_F = min(g$F),
      b = b, se = se, OR = exp(b), OR_L = exp(b - 1.96 * se),
      OR_U = exp(b + 1.96 * se), p = p)
  }))

  if(nrow(res) > 0){
    res[, FDR := p.adjust(p, "BH")]
    res <- res[order(p)]
  }
  out <- file.path(RES, paste0("MR_immune_GWS_", oc_name, ".csv"))
  fwrite(res, out)
  cat(sprintf("traits=%d | FDR<0.05=%d | P<0.05=%d | output=%s\n",
              nrow(res), sum(res$FDR < 0.05, na.rm = TRUE), sum(res$p < 0.05, na.rm = TRUE), out))
  res
}

res_list <- lapply(names(OUTCOMES), function(nm) run_outcome(nm, OUTCOMES[[nm]]))
names(res_list) <- names(OUTCOMES)

main_C <- fread(file.path(RES, "MR_immune_CONDUCTIO.csv"))
main_A <- fread(file.path(RES, "MR_immune_AVBLOCK.csv"))

summ <- data.table(
  threshold = c("P < 1e-5", "P < 5e-8"),
  exposure_iv_records = c(nrow(fread(file.path(RES, "immune_ivs.csv"))), nrow(IVS0)),
  traits_with_CONDUCTIO_results = c(nrow(main_C), nrow(res_list$CONDUCTIO)),
  CONDUCTIO_nominal = c(sum(main_C$p < 0.05, na.rm = TRUE), sum(res_list$CONDUCTIO$p < 0.05, na.rm = TRUE)),
  CONDUCTIO_FDR = c(sum(main_C$FDR < 0.05, na.rm = TRUE), sum(res_list$CONDUCTIO$FDR < 0.05, na.rm = TRUE)),
  CONDUCTIO_min_FDR = c(min(main_C$FDR, na.rm = TRUE), min(res_list$CONDUCTIO$FDR, na.rm = TRUE)),
  traits_with_AVBLOCK_results = c(nrow(main_A), nrow(res_list$AVBLOCK)),
  AVBLOCK_nominal = c(sum(main_A$p < 0.05, na.rm = TRUE), sum(res_list$AVBLOCK$p < 0.05, na.rm = TRUE)),
  AVBLOCK_FDR = c(sum(main_A$FDR < 0.05, na.rm = TRUE), sum(res_list$AVBLOCK$FDR < 0.05, na.rm = TRUE)),
  AVBLOCK_min_FDR = c(min(main_A$FDR, na.rm = TRUE), min(res_list$AVBLOCK$FDR, na.rm = TRUE))
)

fwrite(summ, file.path(RES, "threshold_sensitivity_summary.csv"))
cat("\n==== Threshold sensitivity summary ====\n")
print(summ)
