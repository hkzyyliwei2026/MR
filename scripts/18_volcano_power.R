# =============================================================
# 18: 731-phenotype volcano plot + power/minimum detectable effect analysis
# Volcano plots: both outcomes, x = ln(OR), y = -log10(P), highlighting the 5 candidates and the P=0.05 line
# Power: binary-outcome MR power (Brion 2013 approximation), var(bhat) ~ 1/(N*R2*K*(1-K))
# Output: results/figures/Fig2_volcano.pdf/png ; results/tables/power_analysis.csv
# =============================================================
suppressMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
PROJ <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
RES  <- file.path(PROJ, "results/tables"); FIG <- file.path(PROJ, "results/figures")
C <- fread(file.path(RES, "MR_immune_CONDUCTIO.csv"))
A <- fread(file.path(RES, "MR_immune_AVBLOCK.csv"))
robust <- fread(file.path(RES, "immune_sensitivity.csv"))[robust == TRUE, trait]
gws_hits <- fread(file.path(RES, "MR_immune_GWS_CONDUCTIO.csv"))[order(FDR)][FDR < 0.05, trait]
PS <- fread(file.path(RES, "immune_presso_steiger.csv"))   # contains R2_exp_pct for each candidate
short_lab <- c(
  "HLA DR+ Natural Killer %Natural Killer" = "HLA-DR+ NK",
  "CD45 on granulocyte" = "CD45 granulocyte",
  "CD45RA on CD39+ resting CD4 regulatory T cell" = "CD45RA+ Treg",
  "CD25++ CD45RA+ CD4 not regulatory T cell %T cell" = "CD25++ non-Treg",
  "CD62L- CD86+ myeloid Dendritic Cell Absolute Count" = "CD62L-CD86+ mDC")
label_meta <- data.table(
  trait = names(short_lab),
  lab = unname(short_lab),
  c_dx = c(-0.008, -0.008, -0.008,  0.008,  0.008),
  c_dy = c( 0.16, -0.18,  0.14,  0.17,  0.16),
  c_adj = c(1, 1, 1, 0, 0),
  a_dx = c(-0.008, -0.028, -0.008,  0.010,  0.008),
  a_dy = c( 0.16, -0.20,  0.16,  0.23,  0.20),
  a_adj = c(1, 1, 1, 0, 0)
)

open_svg <- function(path, width, height){
  if(requireNamespace("svglite", quietly=TRUE)){
    svglite::svglite(path, width=width, height=height)
  } else {
    svg(path, width=width, height=height, family="Arial")
  }
}

# ---------------- Fig 2: volcano plots (two panels) ----------------
draw_volc <- function(D, ttl){
  D <- copy(D); D[, y := -log10(p)]; D[, x := log(OR)]
  D[, grp := "ns"]; D[p < 0.05, grp := "nom"]; D[trait %in% robust, grp := "rob"]
  cols <- c(ns = "#C8CDD2", nom = "#1F5A85", rob = "#B94A3A")
  fills <- c(ns = "#DDE1E5", nom = "#D8E8F2", rob = "#F5D9D3")
  xr <- max(abs(D$x), na.rm = TRUE) * 1.08
  yr <- max(D$y, na.rm = TRUE) * 1.16
  plot(NA, xlim = c(-xr, xr), ylim = c(0, yr),
       xlab = "ln(OR) per SD increase in immune trait", ylab = expression(-log[10](italic(P))),
       main = ttl, cex.main = 0.90, las = 1, axes = FALSE)
  axis(1, las = 1, col = "#333333", col.axis = "#333333")
  axis(2, las = 1, col = "#333333", col.axis = "#333333")
  box(col = "#333333", lwd = 0.8)
  grid(col = "#E8ECEF", lty = 1, lwd = 0.8)
  abline(v = 0, col = "#5D646B", lwd = 1.0)
  abline(h = -log10(0.05), lty = 2, col = "#5D646B", lwd = 1.0)
  ord <- order(factor(D$grp, levels = c("ns", "nom", "rob")))
  points(D$x[ord], D$y[ord],
         pch = 21, bg = fills[D$grp[ord]], col = cols[D$grp[ord]],
         lwd = ifelse(D$grp[ord] == "rob", 1.0, 0.5),
         cex = ifelse(D$grp[ord] == "rob", 1.15, 0.60))
  G <- D[trait %in% gws_hits]
  if(nrow(G) > 0){
    points(G$x, G$y, pch = 23, bg = "#F2C94C", col = "#7A5A00",
           lwd = 0.85, cex = 0.95)
  }
  R <- D[grp == "rob"]
  if(nrow(R) > 0){
    panel <- if(grepl("Atrioventricular", ttl)) "a" else "c"
    R <- merge(R, label_meta, by = "trait", all.x = TRUE, sort = FALSE)
    R[is.na(lab), lab := trait]
    if(panel == "a"){
      R[, `:=`(lx = x + a_dx, ly = y + a_dy, ladj = a_adj)]
    } else {
      R[, `:=`(lx = x + c_dx, ly = y + c_dy, ladj = c_adj)]
    }
    segments(R$x, R$y, R$lx, R$ly, col = "#9AA1A8", lwd = 0.55)
    text(R$lx, R$ly, labels = R$lab, adj = c(R$ladj, 0.5),
         cex = 0.48, col = "#333333")
  }
  text(-xr * 0.98, -log10(0.05), "P = 0.05", pos = 4, cex = 0.65, col = "#4F4F4F")
  legend("topright",
         c("Not significant", "Nominal P < 0.05", "Primary-screen diagnostic subset",
           "Genome-wide-threshold FDR signal"),
         pch = c(21, 21, 21, 23),
         pt.bg = c(fills, "#F2C94C"),
         col = c(cols, "#7A5A00"),
         pt.cex = c(0.75, 0.75, 1.05, 0.95),
         cex = 0.62, bty = "n")
}
for(dev_fun in list(function() cairo_pdf(file.path(FIG, "Fig2_volcano.pdf"), width = 11, height = 5.3, family = "Arial"),
                    function() open_svg(file.path(FIG, "Fig2_volcano.svg"), width = 11, height = 5.3),
                    function() png(file.path(FIG, "Fig2_volcano.png"), width = 6600, height = 3180, res = 600, pointsize = 13),
                    function() tiff(file.path(FIG, "Fig2_volcano.tiff"), width = 11, height = 5.3, units = "in", res = 600, pointsize = 13, compression = "lzw"))){
  dev_fun(); par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3, 1), oma = c(0, 0, 0, 0), family = "Arial")
  draw_volc(C, "Cardiac conduction disorders (731 phenotypes)")
  draw_volc(A, "Atrioventricular block (731 phenotypes)")
  dev.off()
}
cat("volcano plots written: Fig2_volcano.pdf/png\n")

# ---------------- Power analysis ----------------
za  <- qnorm(1 - 0.05/2)              # nominal α=0.05
zsw <- qnorm(1 - (0.05/731)/2)        # study-wide α=0.05/731
zb  <- qnorm(0.80)                    # 80% power
oc  <- list(CONDUCTIO = c(N = 12371 + 342690, K = 12371/(12371 + 342690)),
            AVBLOCK   = c(N = 6935 + 342690,  K = 6935/(6935 + 342690)))

# (a) minimum detectable OR at 80% power across values of R2
R2v <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20)
gridpow <- rbindlist(lapply(names(oc), function(nm){
  N <- oc[[nm]]["N"]; K <- oc[[nm]]["K"]
  rbindlist(lapply(R2v, function(r2){
    se <- 1/sqrt(N * r2 * K * (1 - K))
    data.table(outcome = nm, R2_pct = r2*100,
               minDetectableOR_nominal   = round(exp((za  + zb) * se), 3),
               minDetectableOR_studywide = round(exp((zsw + zb) * se), 3))
  }))
}))
fwrite(gridpow, file.path(RES, "power_analysis.csv"))
cat("\n== minimum detectable OR at 80% power (>1 direction; <1 is the reciprocal) ==\n"); print(gridpow)

# (b) the 5 candidates: realised power given each trait's R2 and observed OR
K <- oc$CONDUCTIO["K"]; N <- oc$CONDUCTIO["N"]
hit <- merge(PS[, .(trait, R2_exp_pct)],
             C[trait %in% robust, .(trait, OR)], by = "trait")
hit[, se := 1/sqrt(N * (R2_exp_pct/100) * K * (1 - K))]
hit[, b := abs(log(OR))]
hit[, power_nominal   := round(pnorm(b/se - za)  + pnorm(-b/se - za), 3)]
hit[, power_studywide := round(pnorm(b/se - zsw) + pnorm(-b/se - zsw), 3)]
hit[, `:=`(OR = round(OR, 3), R2_exp_pct = round(R2_exp_pct, 1), se = NULL, b = NULL)]
setorder(hit, -power_studywide)
fwrite(hit, file.path(RES, "power_hits.csv"))
cat("\n== power of the 5 candidates at their observed effects ==\n"); print(hit)
cat("\ndone.\n")
