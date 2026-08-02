# =============================================================
# 12a: extract instruments for the 731 immune-cell phenotypes (OpenGWAS tophits, p<1e-5, server-side LD clumping)
# Resumable: completed IDs are recorded in immune_done.txt and skipped on re-run
# Output: results/tables/immune_ivs.csv (all instruments), immune_done.txt (completed IDs)
# =============================================================
suppressMessages({library(ieugwasr); library(data.table)})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT"))>0)   # requires an OpenGWAS token in .Renviron
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ,"results/tables"); dir.create(RES,showWarnings=FALSE,recursive=TRUE)
IV_CSV   <- file.path(RES,"immune_ivs.csv")
DONE_TXT <- file.path(RES,"immune_done.txt")

ids <- sprintf("ebi-a-GCST%d", 90001391:90002121)   # 731 consecutive accessions
done <- if(file.exists(DONE_TXT)) readLines(DONE_TXT) else character(0)
done <- unique(done[nzchar(done)])
todo <- setdiff(ids, done)
cat(sprintf("%d traits total, %d done, %d remaining\n", length(ids), length(done), length(todo)))

get_iv <- function(id){
  for(try in 1:4){
    r <- tryCatch(tophits(id, pval=1e-5, clump=TRUE),
                  error=function(e){ msg<-conditionMessage(e)
                    if(grepl("429|rate|Too Many", msg, ignore.case=TRUE)) Sys.sleep(20) else Sys.sleep(3)
                    NULL })
    if(!is.null(r)) return(r)
  }
  NULL
}

n_ok<-0; n_empty<-0; n_fail<-0
for(i in seq_along(todo)){
  id <- todo[i]
  r <- get_iv(id)
  if(is.null(r)){ n_fail<-n_fail+1; next }          # failure: not recorded as done, retried next run
  if(nrow(r)>0){
    dt <- as.data.table(r)[, .(id,trait,chr,position,rsid,ea,nea,eaf,beta,se,p,n)]
    fwrite(dt, IV_CSV, append=file.exists(IV_CSV))
    n_ok<-n_ok+1
  } else n_empty<-n_empty+1
  write(id, DONE_TXT, append=TRUE)
  if(i %% 25 == 0) cat(sprintf("  [%d/%d] ok=%d empty=%d fail=%d\n", i,length(todo),n_ok,n_empty,n_fail))
  Sys.sleep(0.15)
}
cat(sprintf("\nthis run: ok=%d empty(0 IV)=%d fail(retry)=%d\n", n_ok,n_empty,n_fail))
cat("if fail>0, re-run this script to resume the remainder\n")
