# =============================================================
# 25: Phenotype-level instrument strength and statistical power (Supplementary Table S18)
#
# The harmonisation block below is copied verbatim from scripts/12b_immune_mr.R so that
# the instrument set used here is the same set the primary MR analysis used. Computing
# R^2 from the unharmonised exposure-side table would overstate both the instrument
# count and the variance explained.
#
# Input : results/tables/immune_ivs.csv ; <finngen_dir>/finngen_R11_I9_{CONDUCTIO,AVBLOCK}.gz
# Output: results/tables/S18_phenotype_power.csv
# Usage : Rscript scripts/25_phenotype_power.R [finngen_dir]
#         finngen_dir defaults to data/outcome inside the project, as in 12b_immune_mr.R
# =============================================================
suppressMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")
args <- commandArgs(TRUE)
FG <- if (length(args)) args[1] else file.path(PROJ, "data/outcome")

IVS <- fread(file.path(RES, "immune_ivs.csv"))
setnames(IVS, c("ea","nea","eaf","beta","se"), c("ea_e","nea_e","eaf_e","beta_e","se_e"))
IVS[, `:=`(ea_e=toupper(ea_e), nea_e=toupper(nea_e))]
is_palindromic <- function(a1,a2){ p<-c("A"="T","T"="A","C"="G","G"="C"); !is.na(p[a1]) & p[a1]==a2 }

OC <- list(CONDUCTIO=c(N=12371+342690, K=12371/(12371+342690)),
           AVBLOCK  =c(N= 6935+342690, K= 6935/( 6935+342690)))
za <- qnorm(1-(0.05/731)/2); zb <- qnorm(0.80)

pw <- function(r2,N,K,orv){ if(is.na(r2)||r2<=0||r2>1) return(NA_real_)
  se<-1/sqrt(N*r2*K*(1-K)); b<-abs(log(orv)); pnorm(b/se-za)+pnorm(-b/se-za) }
mdo <- function(r2,N,K){ if(is.na(r2)||r2<=0||r2>1) return(NA_real_)
  exp((za+zb)/sqrt(N*r2*K*(1-K))) }

# Effect-size bounds are reported both phenotype-specifically and simultaneously across the
# screen, because inference and power in this study are multiplicity-adjusted and an
# unadjusted interval cannot support a study-wide exclusion claim.
MR <- list(CONDUCTIO=fread(file.path(RES,'MR_immune_CONDUCTIO.csv')),
           AVBLOCK  =fread(file.path(RES,'MR_immune_AVBLOCK.csv')))
z95 <- qnorm(0.975); zsim <- qnorm(1-(0.05/731)/2)

acc <- list()
for (nm in names(OC)) {
  need <- unique(IVS$rsid)
  OUT <- fread(file.path(FG, sprintf("finngen_R11_I9_%s.gz", nm)),
               select=c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
  setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
                c("rsid","ref_o","alt_o","beta_o","se_o","p_o","eaf_o"))
  OUT <- OUT[rsid %in% need & !is.na(beta_o) & se_o>0]
  OUT[, `:=`(ref_o=toupper(ref_o), alt_o=toupper(alt_o))]
  m <- merge(IVS, OUT, by="rsid")
  m[, keep := TRUE]; m[, beta_o_al := beta_o]
  m[ea_e==alt_o & nea_e==ref_o, beta_o_al := beta_o]
  m[ea_e==ref_o & nea_e==alt_o, beta_o_al := -beta_o]
  m[!( (ea_e==alt_o & nea_e==ref_o) | (ea_e==ref_o & nea_e==alt_o) ), keep := FALSE]
  m[is_palindromic(ea_e,nea_e) & pmin(eaf_e,1-eaf_e)>0.42, keep := FALSE]
  m <- m[keep==TRUE]
  m[, F := (beta_e/se_e)^2]; m <- m[F>10]

  N <- OC[[nm]]["N"]; K <- OC[[nm]]["K"]
  # sumF_n = sum(F)/n is the scale-invariant counterpart of the aggregate R^2: rescaling the
  # reported effect sizes changes R^2 but leaves F = (beta/se)^2 untouched. Phenotypes whose
  # sumF_n exceeds 1 are implausible on a metric that no reporting-scale error can produce.
  g <- m[, .(nIV=.N, R2=sum(2*eaf_e*(1-eaf_e)*beta_e^2, na.rm=TRUE),
             sumF_n=sum(F, na.rm=TRUE)/median(n, na.rm=TRUE)), by=.(id,trait)]
  g[, `:=`(outcome=nm,
    p11=sapply(R2,pw,N,K,1.1), p12=sapply(R2,pw,N,K,1.2), p13=sapply(R2,pw,N,K,1.3),
    mdoh=sapply(R2,mdo,N,K))]
  acc[[nm]] <- g
}
res <- rbindlist(acc)
res <- merge(res, rbindlist(lapply(names(MR), function(nm){
  x <- MR[[nm]][, .(id, b, se)]; x[, outcome := nm]; x}), fill=TRUE),
  by=c('id','outcome'), all.x=TRUE)
res[, `:=`(ci95_lo=exp(b-z95*se), ci95_hi=exp(b+z95*se),
           sim_lo =exp(b-zsim*se), sim_hi =exp(b+zsim*se))]
# Only aggregate R2 > 1 is excluded from the power summaries: for a variance-standardised
# exposure it is arithmetically impossible. sum(F)/n > 1 is an empirical flag for unusually
# large cumulative instrument strength, not an estimate of joint variance explained and not a
# formal validity criterion, so those phenotypes are flagged but retained.
res[, bad := R2 > 1]
res[, `:=`(
  `Outcome`=outcome,
  `Immunophenotype`=trait,
  `GWAS Catalog accession`=sub("^ebi-a-","",id),
  `Instruments`=nIV,
  # the R^2 value is kept for flagged phenotypes as well: Section 3.1 quotes their range, and a
  # blank cell would leave that claim uncheckable. The QC flag carries the interpretation caveat.
  `R2 (%)`=sprintf("%.2f",R2*100),
  `sum(F)/n`=sprintf("%.3f",sumF_n),
  `Power at OR 1.1`=ifelse(bad|is.na(p11),"not interpretable",sprintf("%.3f",p11)),
  `Power at OR 1.2`=ifelse(bad|is.na(p12),"not interpretable",sprintf("%.3f",p12)),
  `Power at OR 1.3`=ifelse(bad|is.na(p13),"not interpretable",sprintf("%.3f",p13)),
  `Minimum detectable protective OR at 80% power`=ifelse(bad|is.na(mdoh),"not interpretable",sprintf("%.3f",1/mdoh)),
  `Minimum detectable harmful OR at 80% power`=ifelse(bad|is.na(mdoh),"not interpretable",sprintf("%.3f",mdoh)),
  `OR 95%CI lower`=sprintf('%.4f',ci95_lo),
  `OR 95%CI upper`=sprintf('%.4f',ci95_hi),
  `OR simultaneous CI lower (alpha=.05/731)`=sprintf('%.4f',sim_lo),
  `OR simultaneous CI upper (alpha=.05/731)`=sprintf('%.4f',sim_hi),
  `Entire 95% CI within 0.83-1.20`=fifelse(ci95_lo>1/1.20 & ci95_hi<1.20,'yes','no'),
  `Entire simultaneous CI within 0.83-1.20`=fifelse(sim_lo>1/1.20 & sim_hi<1.20,'yes','no'),
  `Entire 95% CI within 0.91-1.10`=fifelse(ci95_lo>1/1.10 & ci95_hi<1.10,'yes','no'),
  `Entire simultaneous CI within 0.91-1.10`=fifelse(sim_lo>1/1.10 & sim_hi<1.10,'yes','no'),
  `QC flag`=fifelse(R2>1 & sumF_n>1, "aggregate R2 > 1 (excluded from power summaries); sum(F)/n > 1",
             fifelse(R2>1, "aggregate R2 > 1 (excluded from power summaries)",
             fifelse(sumF_n>1, "sum(F)/n > 1 (empirical flag; retained in summaries)", ""))))]
keep <- c("Outcome","Immunophenotype","GWAS Catalog accession","Instruments","R2 (%)","sum(F)/n",
          "Power at OR 1.1","Power at OR 1.2","Power at OR 1.3",
          "Minimum detectable protective OR at 80% power",
          "Minimum detectable harmful OR at 80% power",
          "OR 95%CI lower","OR 95%CI upper",
          "OR simultaneous CI lower (alpha=.05/731)","OR simultaneous CI upper (alpha=.05/731)",
          "Entire 95% CI within 0.83-1.20","Entire simultaneous CI within 0.83-1.20",
          "Entire 95% CI within 0.91-1.10","Entire simultaneous CI within 0.91-1.10","QC flag")
fwrite(res[, ..keep], file.path(RES,"S18_phenotype_power.csv"))
cat("wrote", file.path(RES,"S18_phenotype_power.csv"), nrow(res), "rows\n")
