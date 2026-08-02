# =============================================================
# 15: Supplementary figures
# Fig2 = forest plot of the 5 sensitivity-stable traits (forward vs reverse); Fig3 = leave-one-out for the same 5
# Output: results/figures/Fig2_forest.pdf/png, Fig3_leaveoneout.pdf/png
# =============================================================
suppressMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES <- file.path(PROJ,"results/tables")
FIG  <- file.path(PROJ,"results/figures"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

open_svg <- function(path, width, height){
  if(requireNamespace("svglite", quietly=TRUE)){
    svglite::svglite(path, width=width, height=height)
  } else {
    svg(path, width=width, height=height, family="Arial")
  }
}

S  <- fread(file.path(RES,"immune_sensitivity.csv"))[robust==TRUE]
Cf <- fread(file.path(RES,"MR_immune_CONDUCTIO.csv"))
Rv <- fread(file.path(RES,"MR_reverse.csv"))
LOO<- fread(file.path(RES,"immune_leaveoneout.csv"))

# short display labels
lab <- c("HLA DR+ Natural Killer %Natural Killer"              = "HLA-DR+ NK (%NK)",
         "CD45 on granulocyte"                                 = "CD45 on granulocyte",
         "CD45RA on CD39+ resting CD4 regulatory T cell"       = "CD45RA on CD39+ resting Treg",
         "CD25++ CD45RA+ CD4 not regulatory T cell %T cell"    = "CD25++CD45RA+ CD4 non-Treg (%T)",
         "CD62L- CD86+ myeloid Dendritic Cell Absolute Count"  = "CD62L- CD86+ myeloid DC (AC)")
fwd <- merge(S[,.(trait)], Cf[,.(trait,OR,OR_L,OR_U,p)], by="trait")
fwd <- fwd[order(match(trait, S$trait))]
rev <- merge(S[,.(trait)], Rv[,.(trait,OR,OR_L,OR_U,p)], by="trait")
rev <- rev[order(match(trait, S$trait))]
fwd[, lab := lab[trait]]; rev[, lab := lab[trait]]

# ---------------- Fig 2: forest plot (forward vs reverse) ----------------
draw_forest <- function(){
  n <- nrow(fwd); y <- n:1
  BLUE <- "#1F5A85"; GREY <- "#5D646B"; GRID <- "#E8ECEF"
  par(mar=c(4.7,14.5,2.8,1.2), family="Arial", xpd=FALSE)
  xr <- range(c(fwd$OR_L,fwd$OR_U,rev$OR_L,rev$OR_U,0.88,1.14), na.rm=TRUE)
  plot(NA, xlim=xr, ylim=c(0.5,n+0.7), log="x", axes=FALSE,
       xlab="Odds ratio (95% CI)", ylab="")
  rect(par("usr")[1], y - 0.5, par("usr")[2], y + 0.5,
       col=rep(c("#FFFFFF","#F7F9FA"), length.out=n), border=NA)
  abline(h=y, col=GRID, lwd=0.8)
  abline(v=1, lty=2, col="#4F4F4F", lwd=1.0)
  axis(1, las=1, col="#333333", col.axis="#333333")
  axis(2, at=y, labels=fwd$lab, las=1, cex.axis=0.78, tick=FALSE, line=-0.25)
  box(col="#333333", lwd=0.8)
  # forward (blue, offset up); reverse (grey, offset down)
  off<-0.15
  segments(fwd$OR_L, y+off, fwd$OR_U, y+off, col=BLUE, lwd=2.1)
  points(fwd$OR, y+off, pch=22, bg="#D8E8F2", col=BLUE, lwd=1.0, cex=1.25)
  segments(rev$OR_L, y-off, rev$OR_U, y-off, col=GREY, lwd=1.9)
  points(rev$OR, y-off, pch=24, bg="#ECEFF1", col=GREY, lwd=1.0, cex=1.05)
  legend("topright", c("Forward: immune to conduction","Reverse: conduction to immune"),
         col=c(BLUE,GREY), pt.bg=c("#D8E8F2","#ECEFF1"), pch=c(22,24), bty="n", cex=0.75)
  title("Primary-screen diagnostic subset", cex.main=0.95)
}
cairo_pdf(file.path(FIG,"Fig2_forest.pdf"), width=9, height=4.6, family="Arial"); draw_forest(); dev.off()
open_svg(file.path(FIG,"Fig2_forest.svg"), width=9, height=4.6); draw_forest(); dev.off()
png(file.path(FIG,"Fig2_forest.png"), width=5400, height=2760, res=600, pointsize=13); draw_forest(); dev.off()
tiff(file.path(FIG,"Fig2_forest.tiff"), width=9, height=4.6, units="in", res=600, pointsize=13, compression="lzw"); draw_forest(); dev.off()

# ---------------- Fig 3: leave-one-out (5 traits, small multiples) ----------------
draw_loo <- function(){
  BLUE <- "#1F5A85"; GRID <- "#E8ECEF"
  par(mfrow=c(2,3), mar=c(3.2,3.2,2.4,1.0), oma=c(0.2,0.2,0.2,0.2), family="Arial")
  for(tr in S$trait){
    sub <- LOO[trait==tr]; if(nrow(sub)==0) next
    xrg <- range(c(sub$OR,1), na.rm=TRUE)
    pad <- diff(xrg) * 0.08; if(!is.finite(pad) || pad == 0) pad <- 0.02
    plot(NA, xlim=xrg + c(-pad, pad), ylim=c(0.5, nrow(sub)+0.5), axes=FALSE,
         ylab="", xlab="IVW OR (leave-one-out)", main=lab[tr], cex.main=0.76)
    axis(1, las=1, cex.axis=0.78, col="#333333", col.axis="#333333")
    axis(2, at=pretty(seq_len(nrow(sub)), n=3), labels=FALSE, col="#333333")
    box(col="#333333", lwd=0.8)
    grid(col=GRID, lty=1, lwd=0.8)
    abline(v=1, lty=2, col="#4F4F4F", lwd=1.0)
    points(sub$OR, seq_len(nrow(sub)), pch=21, cex=0.72, bg="#D8E8F2", col=BLUE, lwd=0.8)
  }
  plot.new()
  text(0.5, 0.60, "All SNPs removed one at a time", cex=0.75, font=2, col="#333333")
  text(0.5, 0.47, "Dashed line marks OR = 1", cex=0.68, col="#4F4F4F")
}
cairo_pdf(file.path(FIG,"Fig3_leaveoneout.pdf"), width=10, height=6, family="Arial"); draw_loo(); dev.off()
open_svg(file.path(FIG,"Fig3_leaveoneout.svg"), width=10, height=6); draw_loo(); dev.off()
png(file.path(FIG,"Fig3_leaveoneout.png"), width=6000, height=3600, res=600, pointsize=13); draw_loo(); dev.off()
tiff(file.path(FIG,"Fig3_leaveoneout.tiff"), width=10, height=6, units="in", res=600, pointsize=13, compression="lzw"); draw_loo(); dev.off()

cat("figures written:\n"); print(list.files(FIG, pattern="Fig[23]"))
