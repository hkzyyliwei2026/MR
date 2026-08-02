# =============================================================
# 14: reverse-direction MR: cardiac conduction disorders (CONDUCTIO) -> the 5 sensitivity-stable immune traits
# Exposure: FinnGen CONDUCTIO p<5e-8, LD-clumped (r2<0.001); Outcome: immune-cell traits (OpenGWAS)
# Reverse-direction analyses are used to assess whether prioritized forward
# candidates show evidence compatible with reverse causation.
# Output: results/tables/MR_reverse.csv
# =============================================================
suppressMessages({library(data.table); library(ieugwasr); library(MendelianRandomization)})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT"))>0)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES <- file.path(PROJ,"results/tables")

# --- ids of the 5 sensitivity-stable traits ---
S <- fread(file.path(RES,"immune_sensitivity.csv"))[robust==TRUE]
C <- fread(file.path(RES,"MR_immune_CONDUCTIO.csv"))
hits <- merge(S[,.(trait)], unique(C[,.(id,trait)]), by="trait")
cat("reverse-MR targets (sensitivity-stable traits):\n"); print(hits)

# --- exposure: CONDUCTIO instruments ---
d <- fread(file.path(PROJ,"data/outcome/finngen_R11_I9_CONDUCTIO.gz"),
           select=c("rsids","ref","alt","beta","sebeta","pval","af_alt","#chrom","pos"))
setnames(d, c("rsids","ref","alt","beta","sebeta","pval","af_alt","#chrom","pos"),
            c("rsid","nea_e","ea_e","beta_e","se_e","p_e","eaf_e","chrom","pos"))  # FinnGen effect allele = alt
d <- d[p_e<5e-8 & rsid!="" & !is.na(rsid) & se_e>0]
d[, `:=`(ea_e=toupper(ea_e), nea_e=toupper(nea_e))]
# LD clumping (OpenGWAS reference, r2<0.001, 10 Mb)
cl <- tryCatch(ld_clump(data.frame(rsid=d$rsid, pval=d$p_e, id="CONDUCTIO"),
                        clump_r2=0.001, clump_kb=10000, pop="EUR"),
               error=function(e){ cat("ld_clump failed:",conditionMessage(e),"; falling back to distance-based pruning (+/-10 Mb)\n"); NULL })
if(!is.null(cl)){ d <- d[rsid %in% cl$rsid] } else {
  d <- d[order(p_e)]; kept <- d[0]
  while(nrow(d)>0){ top<-d[1]; kept<-rbind(kept,top)
    d <- d[!(chrom==top$chrom & abs(pos-top$pos)<1e7)] }
  d <- kept
}
d[, F := (beta_e/se_e)^2]; d <- d[F>10]
cat(sprintf("CONDUCTIO independent instruments (after clumping, F>10) = %d\n", nrow(d)))

# --- per immune trait: query outcome effects from OpenGWAS -> reverse MR ---
out <- rbindlist(lapply(seq_len(nrow(hits)), function(k){
  id<-hits$id[k]; tr<-hits$trait[k]
  a <- tryCatch(associations(variants=d$rsid, id=id), error=function(e) NULL)
  if(is.null(a)||nrow(a)==0) return(data.table(trait=tr,nIV=0,method=NA,OR=NA,p=NA,note="no outcome data"))
  a <- as.data.table(a)[, .(rsid, ea_o=toupper(ea), nea_o=toupper(nea), beta_o=beta, se_o=se)]
  g <- merge(d, a, by="rsid")
  g[, beta_o_al := beta_o]
  g[ea_e==nea_o & nea_e==ea_o, beta_o_al := -beta_o]           # opposite direction, flip sign
  g <- g[(ea_e==ea_o & nea_e==nea_o)|(ea_e==nea_o & nea_e==ea_o)]
  if(nrow(g)<1) return(data.table(trait=tr,nIV=0,method=NA,OR=NA,p=NA,note="no matching SNP"))
  if(nrow(g)==1){ b<-g$beta_o_al/g$beta_e; se<-abs(g$se_o/g$beta_e); mth<-"Wald"; p<-2*pnorm(-abs(b/se)) }
  else { fit<-mr_ivw(mr_input(bx=g$beta_e,bxse=g$se_e,by=g$beta_o_al,byse=g$se_o)); b<-fit$Estimate;se<-fit$StdError;mth<-"IVW";p<-fit$Pvalue }
  data.table(trait=tr, nIV=nrow(g), method=mth, OR=round(exp(b),3),
             OR_L=round(exp(b-1.96*se),3), OR_U=round(exp(b+1.96*se),3), p=signif(p,3), note="ok")
}), fill=TRUE)
fwrite(out, file.path(RES,"MR_reverse.csv"))
cat("\n=== reverse MR: conduction disorders -> immune traits ===\n"); print(out)
cat("\nInterpretation: p>0.05 does not support reverse-direction evidence for the prioritized candidates.\n")
