# =============================================================
# 13: sensitivity analyses for the 15 cross-outcome concordant traits (primary outcome CONDUCTIO)
# IVW / weighted median / MR-Egger (intercept pleiotropy, I2GX-NOME) / Cochran's Q / leave-one-out
# Input: results/tables/immune_ivs.csv + MR_immune_{CONDUCTIO,AVBLOCK}.csv + full outcome file
# Output: results/tables/immune_sensitivity.csv, immune_leaveoneout.csv
# =============================================================
suppressMessages({library(data.table); library(MendelianRandomization)})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES <- file.path(PROJ,"results/tables")
OC   <- file.path(PROJ,"data/outcome/finngen_R11_I9_CONDUCTIO.gz")

# --- identify the 15 cross-outcome concordant ids ---
C <- fread(file.path(RES,"MR_immune_CONDUCTIO.csv"))
A <- fread(file.path(RES,"MR_immune_AVBLOCK.csv"))
m <- merge(C[,.(id,trait,OR_C=OR,p_C=p)], A[,.(id,OR_A=OR,p_A=p)], by="id")
hits <- m[p_C<0.05 & p_A<0.05 & sign(OR_C-1)==sign(OR_A-1)][order(p_C)]
cat("concordant traits:", nrow(hits), "\n")

# --- harmonisation (as in 12b) ---
IVS <- fread(file.path(RES,"immune_ivs.csv"))
setnames(IVS,c("ea","nea","eaf","beta","se"),c("ea_e","nea_e","eaf_e","beta_e","se_e"))
IVS[,`:=`(ea_e=toupper(ea_e),nea_e=toupper(nea_e))]
IVS <- IVS[id %in% hits$id]
OUT <- fread(OC, select=c("rsids","ref","alt","beta","sebeta"))
setnames(OUT,c("rsids","ref","alt","beta","sebeta"),c("rsid","ref_o","alt_o","beta_o","se_o"))
OUT <- OUT[rsid %in% unique(IVS$rsid) & !is.na(beta_o) & se_o>0]
OUT[,`:=`(ref_o=toupper(ref_o),alt_o=toupper(alt_o))]
d <- merge(IVS, OUT, by="rsid")
pal <- function(a1,a2){p<-c(A="T",T="A",C="G",G="C"); !is.na(p[a1]) & p[a1]==a2}
d[, beta_o_al := beta_o]
d[ea_e==ref_o & nea_e==alt_o, beta_o_al := -beta_o]
d <- d[(ea_e==alt_o & nea_e==ref_o)|(ea_e==ref_o & nea_e==alt_o)]
d <- d[!(pal(ea_e,nea_e) & pmin(eaf_e,1-eaf_e)>0.42)]
d <- d[(beta_e/se_e)^2 > 10]

# --- per-trait sensitivity analyses ---
sens <- list(); loo <- list()
for(hid in hits$id){
  g <- d[id==hid]; if(nrow(g)<3){ next }
  inp <- mr_input(bx=g$beta_e,bxse=g$se_e,by=g$beta_o_al,byse=g$se_o)
  ivw <- tryCatch(mr_ivw(inp), error=function(e)NULL)
  wm  <- tryCatch(mr_median(inp, weighting="weighted"), error=function(e)NULL)
  eg  <- tryCatch(mr_egger(inp), error=function(e)NULL)
  row <- data.table(
    trait=hits[id==hid,trait], nIV=nrow(g),
    IVW_OR = if(!is.null(ivw)) exp(ivw$Estimate) else NA,
    IVW_p  = if(!is.null(ivw)) ivw$Pvalue else NA,
    WM_OR  = if(!is.null(wm))  exp(wm$Estimate) else NA,
    WM_p   = if(!is.null(wm))  wm$Pvalue else NA,
    Egger_OR = if(!is.null(eg)) exp(eg$Estimate) else NA,
    Egger_p  = if(!is.null(eg)) eg$Pvalue.Est else NA,
    Egger_intercept   = if(!is.null(eg)) eg$Intercept else NA,
    Egger_intercept_p = if(!is.null(eg)) eg$Pvalue.Int else NA,   # >0.05 = no directional pleiotropy
    Q      = if(!is.null(ivw)) ivw$Heter.Stat[1] else NA,
    Q_p    = if(!is.null(ivw)) ivw$Heter.Stat[2] else NA,          # >0.05 = no heterogeneity
    I2_GX  = if(!is.null(eg)) eg$I.sq else NA,                     # MR-Egger NOME I2GX, not the heterogeneity I2
    mean_F = mean((g$beta_e/g$se_e)^2),
    min_F  = min((g$beta_e/g$se_e)^2),
    # R2 is calculated on the harmonized SNP list. A small number of exposure
    # records have missing EAF; those SNPs are retained for MR but excluded from
    # the R2 summation because 2*f*(1-f)*beta^2 is not defined without f.
    R2_pct = sum(2*g$eaf_e*(1-g$eaf_e)*g$beta_e^2, na.rm = TRUE) * 100)
  # concordant direction across the three estimators + significant weighted median
    # + no directional pleiotropy + no Cochran's Q heterogeneity.
    # Leave-one-out is computed below but does NOT enter this classification.
  row[, robust := (!is.na(WM_OR) && !is.na(Egger_OR) &&
                   sign(IVW_OR-1)==sign(WM_OR-1) &&
                   sign(IVW_OR-1)==sign(Egger_OR-1) &&
                   WM_p<0.05 &&
                   (is.na(Egger_intercept_p) || Egger_intercept_p>0.05) &&
                   (is.na(Q_p) || Q_p>0.05))]
  sens[[hid]] <- row
  # leave-one-out (reported separately, not part of the robust flag)
  for(i in seq_len(nrow(g))){
    sub <- g[-i]; if(nrow(sub)<2) next
    f <- tryCatch(mr_ivw(mr_input(bx=sub$beta_e,bxse=sub$se_e,by=sub$beta_o_al,byse=sub$se_o)),error=function(e)NULL)
    if(!is.null(f)) loo[[length(loo)+1]] <- data.table(trait=hits[id==hid,trait], drop_rsid=g$rsid[i], OR=exp(f$Estimate), p=f$Pvalue)
  }
}
S <- rbindlist(sens)
fwrite(S, file.path(RES,"immune_sensitivity.csv"))
fwrite(rbindlist(loo), file.path(RES,"immune_leaveoneout.csv"))
cat("\n=== sensitivity summary (robust = concordant direction + significant WM + no pleiotropy/heterogeneity) ===\n")
print(S[,.(trait,nIV,IVW_OR=round(IVW_OR,3),IVW_p=signif(IVW_p,2),WM_p=signif(WM_p,2),
           Egg_int_p=signif(Egger_intercept_p,2),Q_p=signif(Q_p,2),robust)])
cat("\nrobust traits:", sum(S$robust,na.rm=TRUE), "\n")
