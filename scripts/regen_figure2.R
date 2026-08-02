# Regenerate the Figure 2 volcano plot with the manuscript's two group names
# (consistent with the manuscript and legend); reads results/tables, writes results/figures/Figure_2_Volcano.png
suppressMessages(library(data.table))
par_family <- "Arial"
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables")
OUT  <- file.path(PROJ, "results/figures", "Figure_2_Volcano.png")
C <- fread(file.path(RES, "MR_immune_CONDUCTIO.csv"))
A <- fread(file.path(RES, "MR_immune_AVBLOCK.csv"))
robust <- fread(file.path(RES, "immune_sensitivity.csv"))[robust == TRUE, trait]
gws_hits <- fread(file.path(RES, "MR_immune_GWS_CONDUCTIO.csv"))[order(FDR)][FDR < 0.05, trait]
short_lab <- c(
  "HLA DR+ Natural Killer %Natural Killer" = "HLA-DR+ NK",
  "CD45 on granulocyte" = "CD45 granulocyte",
  "CD45RA on CD39+ resting CD4 regulatory T cell" = "CD45RA+ Treg",
  "CD25++ CD45RA+ CD4 not regulatory T cell %T cell" = "CD25++ non-Treg",
  "CD62L- CD86+ myeloid Dendritic Cell Absolute Count" = "CD62L-CD86+ mDC")
label_meta <- data.table(
  trait = names(short_lab), lab = unname(short_lab),
  c_dx = c(-0.008,-0.008,-0.008, 0.008, 0.008), c_dy = c(0.16,-0.18,0.14,0.17,0.16), c_adj = c(1,1,1,0,0),
  a_dx = c(-0.008,-0.028,-0.008, 0.010, 0.008), a_dy = c(0.16,-0.20,0.16,0.23,0.20), a_adj = c(1,1,1,0,0))

draw_volc <- function(D, ttl){
  D <- copy(D); D[, y := -log10(p)]; D[, x := log(OR)]
  D[, grp := "ns"]; D[p < 0.05, grp := "nom"]; D[trait %in% robust, grp := "rob"]
  cols <- c(ns="#B8BDC2", nom="#0072B2", rob="#009E73")
  fills <- c(ns="#DDE1E5", nom="#D6ECF7", rob="#DDF0E8")
  strict_col <- "#E69F00"
  strict_fill <- "#F6D77A"
  xr <- max(abs(D$x), na.rm=TRUE)*1.08; yr <- max(D$y, na.rm=TRUE)*1.16
  plot(NA, xlim=c(-xr,xr), ylim=c(0,yr),
       xlab="ln(OR) per SD increase in immune trait", ylab=expression(-log[10](italic(P))),
       main=ttl, cex.main=0.90, las=1, axes=FALSE)
  axis(1, las=1, col="#333333", col.axis="#333333"); axis(2, las=1, col="#333333", col.axis="#333333")
  box(col="#333333", lwd=0.8); grid(col="#E8ECEF", lty=1, lwd=0.8)
  abline(v=0, col="#5D646B", lwd=1.0); abline(h=-log10(0.05), lty=2, col="#5D646B", lwd=1.0)
  ord <- order(factor(D$grp, levels=c("ns","nom","rob")))
  points(D$x[ord], D$y[ord], pch=21, bg=fills[D$grp[ord]], col=cols[D$grp[ord]],
         lwd=ifelse(D$grp[ord]=="rob",1.0,0.5), cex=ifelse(D$grp[ord]=="rob",1.15,0.60))
  G <- D[trait %in% gws_hits]
  if(nrow(G)>0) points(G$x, G$y, pch=23, bg=strict_fill, col=strict_col, lwd=0.85, cex=0.95)
  R <- D[grp=="rob"]
  if(nrow(R)>0){
    panel <- if(grepl("Atrioventricular", ttl)) "a" else "c"
    R <- merge(R, label_meta, by="trait", all.x=TRUE, sort=FALSE); R[is.na(lab), lab := trait]
    if(panel=="a") R[, `:=`(lx=x+a_dx, ly=y+a_dy, ladj=a_adj)] else R[, `:=`(lx=x+c_dx, ly=y+c_dy, ladj=c_adj)]
    segments(R$x, R$y, R$lx, R$ly, col="#9AA1A8", lwd=0.55)
    text(R$lx, R$ly, labels=R$lab, adj=c(R$ladj,0.5), cex=0.48, col="#333333")
  }
  text(-xr*0.98, -log10(0.05), "P = 0.05", pos=4, cex=0.65, col="#4F4F4F")
  legend("topright",
         c("Not significant", "Nominal P < 0.05", "Primary-screen diagnostic subset",
           "Genome-wide-threshold FDR signal"),
         pch=c(21,21,21,23), pt.bg=c(fills,strict_fill), col=c(cols,strict_col),
         pt.cex=c(0.75,0.75,1.05,0.95), cex=0.62, bty="n")
}
png(OUT, width=6600, height=3180, res=600, pointsize=13)
par(mfrow=c(1,2), mar=c(4.5,4.8,3,1), oma=c(0,0,0,0), family=par_family)
draw_volc(C, "Cardiac conduction disorders (731 phenotypes)")
draw_volc(A, "Atrioventricular block (731 phenotypes)")
dev.off()
cat("wrote", OUT, "\n")
