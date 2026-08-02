# =============================================================
# 27: Per-instrument credibility filter applied to the primary screen (Supplementary Table S21)
#
# The aggregate quality-control rule of 25_phenotype_power.R (R^2 > 1 or sum(F)/n > 1) acts on
# phenotype-level quantities, and aggregate quantities cannot see a failure that is confined to
# individual variants: a very-low-frequency variant contributes little to aggregate R^2 however
# implausible its own effect estimate. This script therefore repeats the primary screen after
# removing individual instruments on two variant-level criteria:
#
#   |beta| > 1 SD   for an inverse-normal transformed exposure, an effect exceeding one standard
#                   deviation per allele is not credible
#   MAC < 20        expected minor-allele count 2 * n * MAF; below this the variant-exposure
#                   estimate rests on a handful of carriers regardless of its F statistic
#
# The harmonisation and estimation blocks are copied from 12b_immune_mr.R so that the only
# difference from the primary analysis is the added filter.
#
# Input : results/tables/immune_ivs.csv ; <finngen_dir>/finngen_R11_I9_{CONDUCTIO,AVBLOCK}.gz
# Output: results/tables/S21_instrument_qc_sensitivity.csv (per-phenotype)
#         results/tables/S21_qc_summary.csv                (screen-level comparison)
# Usage : Rscript scripts/27_instrument_qc_sensitivity.R [finngen_dir]
# =============================================================
suppressMessages({library(data.table); library(MendelianRandomization)})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")
args <- commandArgs(TRUE)
FG <- if (length(args)) args[1] else file.path(PROJ, "data/outcome")

BETA_MAX <- 1      # SD per allele
MAC_MIN  <- 20     # expected minor-allele count

IVS <- fread(file.path(RES, "immune_ivs.csv"))
setnames(IVS, c("ea","nea","eaf","beta","se"), c("ea_e","nea_e","eaf_e","beta_e","se_e"))
IVS[, `:=`(ea_e=toupper(ea_e), nea_e=toupper(nea_e))]
IVS[, maf := pmin(eaf_e, 1 - eaf_e)]
IVS[, mac := 2 * n * maf]
IVS[, qc_pass := abs(beta_e) <= BETA_MAX & mac >= MAC_MIN]
is_palindromic <- function(a1,a2){ p<-c("A"="T","T"="A","C"="G","G"="C"); !is.na(p[a1]) & p[a1]==a2 }

run <- function(nm, ivs){
  OUT <- fread(file.path(FG, sprintf("finngen_R11_I9_%s.gz", nm)),
               select=c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
  setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
                c("rsid","ref_o","alt_o","beta_o","se_o","p_o","eaf_o"))
  OUT <- OUT[rsid %in% unique(ivs$rsid) & !is.na(beta_o) & se_o>0]
  OUT[, `:=`(ref_o=toupper(ref_o), alt_o=toupper(alt_o))]
  m <- merge(ivs, OUT, by="rsid")
  m[, keep := TRUE]; m[, beta_o_al := beta_o]
  m[ea_e==ref_o & nea_e==alt_o, beta_o_al := -beta_o]
  m[!( (ea_e==alt_o & nea_e==ref_o) | (ea_e==ref_o & nea_e==alt_o) ), keep := FALSE]
  m[is_palindromic(ea_e,nea_e) & pmin(eaf_e,1-eaf_e)>0.42, keep := FALSE]
  m <- m[keep==TRUE]
  m[, F := (beta_e/se_e)^2]; m <- m[F>10]
  rbindlist(lapply(split(m, by="id"), function(g){
    nIV <- nrow(g); if (nIV < 1) return(NULL)
    if (nIV == 1){
      b <- g$beta_o_al/g$beta_e; se <- abs(g$se_o/g$beta_e)
      method <- "Wald"; p <- 2*pnorm(-abs(b/se))
    } else {
      fit <- tryCatch(mr_ivw(mr_input(bx=g$beta_e,bxse=g$se_e,by=g$beta_o_al,byse=g$se_o)),
                      error=function(e) NULL)
      if (is.null(fit)) return(NULL)
      b <- fit@Estimate; se <- fit@StdError; method <- "IVW"; p <- fit@Pvalue
    }
    data.table(id=g$id[1], trait=g$trait[1], nIV=nIV, method=method,
               b=b, se=se, OR=exp(b), OR_L=exp(b-1.96*se), OR_U=exp(b+1.96*se), p=p)
  }))
}

acc <- list(); summ <- list()
for (nm in c("CONDUCTIO","AVBLOCK")) {
  full <- run(nm, IVS)
  filt <- run(nm, IVS[qc_pass == TRUE])
  for (tag in c("all instruments","QC-filtered")) {
    r <- if (tag == "all instruments") copy(full) else copy(filt)
    r[, FDR := p.adjust(p, method="BH")]
    r[, `:=`(outcome=nm, instrument_set=tag)]
    acc[[paste(nm,tag)]] <- r
    summ[[paste(nm,tag)]] <- data.table(
      Outcome=nm, `Instrument set`=tag, `Phenotypes analyzable`=nrow(r),
      `Median instruments`=as.integer(median(r$nIV)),
      `Nominally significant`=sum(r$p < 0.05),
      `Expected nominal`=round(nrow(r)*0.05, 1),
      `FDR significant`=sum(r$FDR < 0.05),
      `Smallest FDR`=round(min(r$FDR), 4),
      `Median OR upper 95% CI`=round(median(r$OR_U), 4))
  }
}
A <- rbindlist(acc)
A[, `GWAS Catalog accession` := sub("^ebi-a-","",id)]
setnames(A, c("trait","nIV","method","b","se","p"),
            c("Immunophenotype","Instruments","Method","beta","SE","P"))
keep <- c("outcome","instrument_set","Immunophenotype","GWAS Catalog accession","Instruments",
          "Method","beta","SE","OR","OR_L","OR_U","P","FDR")
setnames(A, c("outcome","instrument_set","OR_L","OR_U"),
            c("Outcome","Instrument set","OR 95%CI lower","OR 95%CI upper"))
keep <- c("Outcome","Instrument set","Immunophenotype","GWAS Catalog accession","Instruments",
          "Method","beta","SE","OR","OR 95%CI lower","OR 95%CI upper","P","FDR")
fwrite(A[, ..keep], file.path(RES,"S21_instrument_qc_sensitivity.csv"))
S <- rbindlist(summ)
fwrite(S, file.path(RES,"S21_qc_summary.csv"))
cat(sprintf("instruments retained: %d of %d (%.1f%%)\n",
            sum(IVS$qc_pass), nrow(IVS), 100*mean(IVS$qc_pass)))
print(S)
