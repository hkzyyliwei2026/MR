# =============================================================
# 29: Phenotype-wide screen repeated on the component sub-endpoints (Supplementary Table S23)
#
# I9_CONDUCTIO is a composite: it applies the whole of ICD-10 I44-I45 and ICD-9 426 with no
# sub-code exclusions, so its case group mixes degenerative atrioventricular and fascicular block
# with pre-excitation, the long QT syndrome and sinoatrial block (Section 2.3). Repeating the
# primary screen on the narrower component endpoints tests whether the null result is a property
# of the immune-cell traits or an artefact of endpoint heterogeneity: I9_AVBLOCK, I9_LBBB and
# I9_RBBB are progressively more homogeneous phenotypes drawn from the same cohort.
#
# The harmonisation and estimation blocks are copied from 12b_immune_mr.R so that the only
# difference from the primary analysis is the outcome file.
#
# Input : results/tables/immune_ivs.csv ; <finngen_dir>/finngen_R11_I9_{LBBB,RBBB}.gz
# Output: results/tables/S23_subendpoint_screen.csv (per-phenotype)
#         results/tables/S23_subendpoint_summary.csv (screen-level)
# Usage : Rscript scripts/29_subendpoint_screen.R [finngen_dir]
# =============================================================
suppressMessages({library(data.table); library(MendelianRandomization)})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")
args <- commandArgs(TRUE)
FG <- if (length(args)) args[1] else file.path(PROJ, "data/outcome")

IVS <- fread(file.path(RES, "immune_ivs.csv"))
setnames(IVS, c("ea","nea","eaf","beta","se"), c("ea_e","nea_e","eaf_e","beta_e","se_e"))
IVS[, `:=`(ea_e=toupper(ea_e), nea_e=toupper(nea_e))]
is_palindromic <- function(a1,a2){ p<-c("A"="T","T"="A","C"="G","G"="C"); !is.na(p[a1]) & p[a1]==a2 }

run <- function(tag){
  f <- file.path(FG, sprintf("finngen_R11_%s.gz", tag))
  if (!file.exists(f)) { cat("missing:", f, "\n"); return(NULL) }
  cat("reading", tag, "...\n")
  OUT <- fread(f, select=c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
  setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
                c("rsid","ref_o","alt_o","beta_o","se_o","p_o","eaf_o"))
  OUT <- OUT[rsid %in% unique(IVS$rsid) & !is.na(beta_o) & se_o>0]
  OUT[, `:=`(ref_o=toupper(ref_o), alt_o=toupper(alt_o))]
  m <- merge(IVS, OUT, by="rsid")
  m[, keep := TRUE]; m[, beta_o_al := beta_o]
  m[ea_e==ref_o & nea_e==alt_o, beta_o_al := -beta_o]
  m[!( (ea_e==alt_o & nea_e==ref_o) | (ea_e==ref_o & nea_e==alt_o) ), keep := FALSE]
  m[is_palindromic(ea_e,nea_e) & pmin(eaf_e,1-eaf_e)>0.42, keep := FALSE]
  m <- m[keep==TRUE]
  m[, F := (beta_e/se_e)^2]; m <- m[F>10]
  r <- rbindlist(lapply(split(m, by="id"), function(g){
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
               OR=exp(b), OR_L=exp(b-1.96*se), OR_U=exp(b+1.96*se), p=p)
  }))
  r[, FDR := p.adjust(p, method="BH")]
  r[, outcome := tag]
  r
}

acc <- lapply(c("I9_LBBB","I9_RBBB"), run)
A <- rbindlist(acc[!vapply(acc, is.null, logical(1))])
A[, `GWAS Catalog accession` := sub("^ebi-a-","",id)]

S <- A[, .(`Phenotypes analyzable`=.N,
           `Median instruments`=as.integer(median(nIV)),
           `Nominally significant`=sum(p<0.05),
           `Expected nominal`=round(.N*0.05,1),
           `FDR significant`=sum(FDR<0.05),
           `Smallest FDR`=round(min(FDR),4),
           `Median OR upper 95% CI`=round(median(OR_U),4)), by=outcome]
setnames(S, "outcome", "Outcome")

setnames(A, c("outcome","trait","nIV","method","OR_L","OR_U","p"),
            c("Outcome","Immunophenotype","Instruments","Method",
              "OR 95%CI lower","OR 95%CI upper","P"))
keep <- c("Outcome","Immunophenotype","GWAS Catalog accession","Instruments","Method",
          "OR","OR 95%CI lower","OR 95%CI upper","P","FDR")
fwrite(A[, ..keep], file.path(RES,"S23_subendpoint_screen.csv"))
fwrite(S, file.path(RES,"S23_subendpoint_summary.csv"))
print(S)
cat("\nwrote", file.path(RES,"S23_subendpoint_screen.csv"), nrow(A), "rows\n")
