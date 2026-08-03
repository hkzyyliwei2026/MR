# =============================================================
# 26: Which sensitivity diagnostics are computable at each instrument threshold, and
#     Cochran's Q for every phenotype for which it is defined (Supplementary Table S19)
#
# Each diagnostic has a minimum instrument count: Cochran's Q and MR-Egger need at least
# three variants, MR-PRESSO at least four, and leave-one-out is uninformative below three.
# Reporting the proportion of phenotypes that clear each requirement separates a threshold
# at which the screen can be checked from one at which it cannot.
#
# The harmonisation block is copied from 12b_immune_mr.R, as in 25_phenotype_power.R.
#
# Input : results/tables/immune_ivs.csv ; <finngen_dir>/finngen_R11_I9_{CONDUCTIO,AVBLOCK}.gz
# Output: results/tables/S19_diagnostic_coverage.csv  (per-phenotype Q)
#         results/tables/S19_coverage_summary.csv     (runnability by threshold)
# Usage : Rscript scripts/26_diagnostic_coverage.R [finngen_dir]
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

THRESH <- list(`P<1e-5` = 1e-5, `P<5e-8` = 5e-8)
OUTCOMES <- c("CONDUCTIO", "AVBLOCK")

# Cochran's Q for the IVW ratio estimate: Q = sum w_j (beta_j - beta_IVW)^2 on df = nIV - 1.
#
# Two weightings are reported. The first-order weight w = (bx/byse)^2 assumes the variant-exposure
# coefficients are measured without error (the NOME assumption). That assumption is not satisfied
# here: instruments are selected at P < 1e-5, so their F statistics are bounded below by 19.51 by
# construction and are not large. Under violated NOME the first-order statistic overstates
# heterogeneity, so the modified second-order weight of Bowden et al is reported alongside it; it
# propagates the variance of bx and is evaluated at the corresponding IVW estimate.
cochran_q <- function(bx, bxse, by, byse){
  r  <- by / bx
  w1 <- (bx / byse)^2
  b1 <- sum(w1 * r) / sum(w1)
  q1 <- sum(w1 * (r - b1)^2)
  # Modified second-order weight of Bowden et al (Int J Epidemiol 2019): the variance of the
  # ratio is evaluated at the IVW estimate, not at each variant's own ratio. Substituting the
  # variant's own by_j downweights precisely the variants that drive heterogeneity, which
  # suppresses Q by construction; an earlier version of this script made that error.
  bh <- b1
  for (i in 1:100) {
    w  <- 1 / (byse^2 / bx^2 + bh^2 * bxse^2 / bx^2)
    nb <- sum(w * r) / sum(w)
    if (abs(nb - bh) < 1e-12) { bh <- nb; break }
    bh <- nb
  }
  w2 <- 1 / (byse^2 / bx^2 + bh^2 * bxse^2 / bx^2)
  q2 <- sum(w2 * (r - bh)^2)
  df <- length(bx) - 1L
  list(Q = q1, df = df, P = pchisq(q1, df, lower.tail = FALSE),
       Q2 = q2, P2 = pchisq(q2, df, lower.tail = FALSE))
}

per_pheno <- list(); summary_rows <- list(); harm <- list()
for (nm in OUTCOMES) {
  OUT <- fread(file.path(FG, sprintf("finngen_R11_I9_%s.gz", nm)),
               select=c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
  setnames(OUT, c("rsids","ref","alt","beta","sebeta","pval","af_alt"),
                c("rsid","ref_o","alt_o","beta_o","se_o","p_o","eaf_o"))
  OUT <- OUT[rsid %in% unique(IVS$rsid) & !is.na(beta_o) & se_o>0]
  OUT[, `:=`(ref_o=toupper(ref_o), alt_o=toupper(alt_o))]

  for (tn in names(THRESH)) {
    ivs <- IVS[p < THRESH[[tn]]]
    m <- merge(ivs, OUT, by="rsid")
    m[, keep := TRUE]; m[, beta_o_al := beta_o]
    m[ea_e==ref_o & nea_e==alt_o, beta_o_al := -beta_o]
    m[!( (ea_e==alt_o & nea_e==ref_o) | (ea_e==ref_o & nea_e==alt_o) ), keep := FALSE]
    m[is_palindromic(ea_e,nea_e) & pmin(eaf_e,1-eaf_e)>0.42, keep := FALSE]
    m <- m[keep==TRUE]
    m[, F := (beta_e/se_e)^2]; m <- m[F>10]

    g <- split(m, by="id")
    rows <- rbindlist(lapply(names(g), function(k){
      d <- g[[k]]
      if (nrow(d) < 3) return(data.table(id=k, trait=d$trait[1], nIV=nrow(d),
                                         Q=NA_real_, Q_df=NA_integer_, Q_P=NA_real_,
                                         Q2=NA_real_, Q2_P=NA_real_))
      q <- cochran_q(d$beta_e, d$se_e, d$beta_o_al, d$se_o)
      data.table(id=k, trait=d$trait[1], nIV=nrow(d), Q=q$Q, Q_df=q$df, Q_P=q$P,
                 Q2=q$Q2, Q2_P=q$P2)
    }))
    rows[, `:=`(outcome=nm, threshold=tn)]
    per_pheno[[paste(nm,tn)]] <- rows

    # STROBE-MR items 9a and 10a require variant-outcome associations alongside the
    # variant-exposure ones. Dump the harmonised outcome side of the primary instrument
    # set so that Supplementary Table S6 can report both.
    if (tn == "P<1e-5") {
      h <- m[, .(rsid, id, beta_outcome=beta_o_al, se_outcome=se_o,
                 p_outcome=p_o, eaf_outcome=eaf_o)]
      setnames(h, "beta_outcome", paste0("beta_", nm))
      setnames(h, "se_outcome",   paste0("se_", nm))
      setnames(h, "p_outcome",    paste0("p_", nm))
      setnames(h, "eaf_outcome",  paste0("eaf_", nm))
      harm[[nm]] <- h
    }

    n <- nrow(rows)
    summary_rows[[paste(nm,tn)]] <- data.table(
      Outcome=nm, `Instrument threshold`=tn, `Analyzable phenotypes`=n,
      `Median instruments`=as.integer(median(rows$nIV)),
      `Single-variant (Wald only), n`=sum(rows$nIV==1),
      `Single-variant (Wald only), %`=round(100*mean(rows$nIV==1),1),
      `Cochran Q / MR-Egger computable (>=3 IVs), n`=sum(rows$nIV>=3),
      `Cochran Q / MR-Egger computable (>=3 IVs), %`=round(100*mean(rows$nIV>=3),1),
      `MR-PRESSO computable (>=4 IVs), n`=sum(rows$nIV>=4),
      `MR-PRESSO computable (>=4 IVs), %`=round(100*mean(rows$nIV>=4),1),
      `Q computed, n`=sum(!is.na(rows$Q)),
      `Q significant at P<.05, n`=sum(rows$Q_P < 0.05, na.rm=TRUE),
      `Q significant at P<.05, %`=round(100*mean(rows$Q_P < 0.05, na.rm=TRUE),1),
      `Q significant, modified weights, n`=sum(rows$Q2_P < 0.05, na.rm=TRUE),
      `Q significant, modified weights, %`=round(100*mean(rows$Q2_P < 0.05, na.rm=TRUE),1),
      `Median Q/df`=round(median(rows$Q/rows$Q_df, na.rm=TRUE),3),
      `Median Q/df, modified weights`=round(median(rows$Q2/rows$Q_df, na.rm=TRUE),3))
  }
}

PP <- rbindlist(per_pheno)
PP[, `GWAS Catalog accession` := sub("^ebi-a-","",id)]
setnames(PP, c("trait","nIV","Q","Q_df","Q_P","Q2","Q2_P","outcome","threshold"),
             c("Immunophenotype","Instruments","Cochran Q","Q df","Q P",
               "Cochran Q (modified weights)","Q P (modified weights)","Outcome","Instrument threshold"))
keep <- c("Outcome","Instrument threshold","Immunophenotype","GWAS Catalog accession",
          "Instruments","Cochran Q","Q df","Q P",
          "Cochran Q (modified weights)","Q P (modified weights)")
for (cc in c("Cochran Q (modified weights)","Q P (modified weights)")) PP[[cc]] <- round(PP[[cc]], 4)
fwrite(PP[, ..keep], file.path(RES,"S19_diagnostic_coverage.csv"))
fwrite(rbindlist(summary_rows), file.path(RES,"S19_coverage_summary.csv"))
H <- Reduce(function(x,y) merge(x,y,by=c("rsid","id"),all=TRUE), harm)
fwrite(H, file.path(RES,"harmonised_outcome_associations.csv"))
cat("wrote harmonised_outcome_associations.csv", nrow(H), "rows\n")
cat("wrote S19_diagnostic_coverage.csv", nrow(PP), "rows\n")
print(rbindlist(summary_rows))
