# =============================================================
# 12b: two-sample MR of 731 immune-cell traits on conduction outcomes (local, no network required)
# Input: results/tables/immune_ivs.csv (instruments from step 12a)
#       data/outcome/finngen_R11_I9_CONDUCTIO.gz (primary), _AVBLOCK.gz (related secondary endpoint)
# Method: allele harmonisation -> F>10 -> IVW (>=2 IV) / Wald ratio (1 IV) -> FDR across traits
# Output: results/tables/MR_immune_<outcome>.csv
# Usage: Rscript 12b_immune_mr.R          (runs both outcomes)
# =============================================================
suppressMessages({library(data.table); library(MendelianRandomization)})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ,"results/tables")
IVS  <- fread(file.path(RES,"immune_ivs.csv"))
setnames(IVS, c("ea","nea","eaf","beta","se"), c("ea_e","nea_e","eaf_e","beta_e","se_e"))
IVS[, `:=`(ea_e=toupper(ea_e), nea_e=toupper(nea_e))]

OUTCOMES <- list(
  CONDUCTIO = file.path(PROJ,"data/outcome/finngen_R11_I9_CONDUCTIO.gz"),
  AVBLOCK   = file.path(PROJ,"data/outcome/finngen_R11_I9_AVBLOCK.gz"))

is_palindromic <- function(a1,a2){ p<-c("A"="T","T"="A","C"="G","G"="C"); !is.na(p[a1]) & p[a1]==a2 }

run_outcome <- function(oc_name, oc_file){
  cat(sprintf("\n==== outcome: %s ====\n", oc_name))
  need <- unique(IVS$rsid)
  OUT <- fread(oc_file, select=c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
  setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
                c("rsid","ref_o","alt_o","beta_o","se_o","p_o","eaf_o"))
  OUT <- OUT[rsid %in% need & !is.na(beta_o) & se_o>0]
  OUT[, `:=`(ref_o=toupper(ref_o), alt_o=toupper(alt_o))]        # FinnGen: effect allele = alt

  m <- merge(IVS, OUT, by="rsid")
  # harmonisation: align effects to the exposure effect allele ea_e
  m[, keep := TRUE]; m[, beta_o_al := beta_o]
  m[ea_e==alt_o & nea_e==ref_o, beta_o_al := beta_o]            # same direction
  m[ea_e==ref_o & nea_e==alt_o, beta_o_al := -beta_o]           # opposite direction, flip sign
  m[!( (ea_e==alt_o & nea_e==ref_o) | (ea_e==ref_o & nea_e==alt_o) ), keep := FALSE] # allele mismatch
  m[is_palindromic(ea_e,nea_e) & pmin(eaf_e,1-eaf_e)>0.42, keep := FALSE]            # ambiguous palindromic
  m <- m[keep==TRUE]
  m[, F := (beta_e/se_e)^2]
  m <- m[F>10]                                                   # weak-instrument filter

  res <- rbindlist(lapply(split(m, by="id"), function(g){
    nIV <- nrow(g)
    if(nIV<1) return(NULL)
    if(nIV==1){
      b <- g$beta_o_al/g$beta_e; se <- abs(g$se_o/g$beta_e)
      method<-"Wald"; p<-2*pnorm(-abs(b/se))
    } else {
      fit <- tryCatch(mr_ivw(mr_input(bx=g$beta_e,bxse=g$se_e,by=g$beta_o_al,byse=g$se_o)),
                      error=function(e) NULL)
      if(is.null(fit)) return(NULL)
      b<-fit$Estimate; se<-fit$StdError; method<-"IVW"; p<-fit$Pvalue
    }
    data.table(id=g$id[1], trait=g$trait[1], nIV=nIV, method=method,
               b=b, se=se, OR=exp(b), OR_L=exp(b-1.96*se), OR_U=exp(b+1.96*se), p=p)
  }))
  res[, FDR := p.adjust(p, "BH")]
  res <- res[order(p)]
  fwrite(res, file.path(RES, paste0("MR_immune_",oc_name,".csv")))
  cat(sprintf("traits=%d | FDR<0.05: %d | p<0.05: %d\n",
              nrow(res), sum(res$FDR<0.05,na.rm=TRUE), sum(res$p<0.05,na.rm=TRUE)))
  cat("Top 10 (by p):\n"); print(head(res[,.(trait,nIV,method,OR=round(OR,3),p=signif(p,2),FDR=signif(FDR,2))],10))
  res
}

for(nm in names(OUTCOMES)) run_outcome(nm, OUTCOMES[[nm]])
cat("\n==== done; see results/tables/MR_immune_*.csv ====\n")
