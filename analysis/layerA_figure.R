# Layer A publication figure (base graphics, no extra deps).
# Recomputes the limma fit on the lean rds (fast) and draws a 2-panel figure:
#   (A) volcano of within-POD pre->post t vs -log10 P, immune-gene probes in red,
#       the 6 DMPs labelled;
#   (B) density of t-statistics, immune vs background, with the directional MWU p.
.libPaths(c("C:/Users/Administrator/R/win-library/4.4", .libPaths()))
suppressMessages({ library(data.table); library(limma) })
ROOT <- "D:/podprj"; GSE <- file.path(ROOT, "analysis/results/GSE330869")
SEED <- 12345; set.seed(SEED)

beta <- readRDS(file.path(GSE, "GSE330869_beta_matrix_lean.rds"))
raw_cols <- colnames(beta)
beta_cols <- raw_cols[grep("_beta$", raw_cols)]
sheet <- data.frame(raw = beta_cols,
                    clean = sub("_beta$", "", beta_cols),
                    patient = sub("_(Pre|Post)$", "", sub("_beta$", "", beta_cols)),
                    time = ifelse(grepl("_Pre_beta$", beta_cols), "Pre", "Post"),
                    stringsAsFactors = FALSE)
beta <- beta[, beta_cols, drop = FALSE]; colnames(beta) <- sheet$clean

anno <- read.csv(file.path(GSE, "EPICv2_anno_small.csv"), stringsAsFactors = FALSE)
anno_full <- anno; rownames(anno_full) <- anno_full$ID
anno_full <- anno_full[rownames(anno_full) %in% rownames(beta), ]
chr_ok <- anno_full$chr %in% c(as.character(1:22), paste0("chr", 1:22))
gene_ok <- !is.na(anno_full$UCSC_RefGene_Name) & anno_full$UCSC_RefGene_Name != ""
keep <- chr_ok & gene_ok
anno <- anno_full[keep, , drop = FALSE]; beta <- beta[rownames(anno), , drop = FALSE]

bb <- pmin(pmax(beta, 1e-4), 1 - 1e-4); M <- log2(bb / (1 - bb))
paired <- split(colnames(M), sheet$patient[match(colnames(M), sheet$clean)])
paired <- paired[sapply(paired, function(x) {
  tt <- sheet$time[match(x, sheet$clean)]; all(c("Pre", "Post") %in% tt) })]
use <- unlist(paired); sheet_u <- sheet[match(use, sheet$clean), ]
M_u <- M[, use, drop = FALSE]
patient_fac <- factor(sheet_u$patient); time_fac <- factor(sheet_u$time, levels = c("Pre", "Post"))
design <- model.matrix(~ time_fac + patient_fac)
fit <- eBayes(lmFit(M_u, design))
coef_post <- "time_facPost"
top <- topTable(fit, coef = coef_post, number = Inf, sort.by = "P")
post_cols <- sheet_u$clean[sheet_u$time == "Post"]
pre_cols  <- sheet_u$clean[sheet_u$time == "Pre"]
top$deltaBeta <- rowMeans(beta[rownames(top), post_cols, drop = FALSE], na.rm = TRUE) -
                 rowMeans(beta[rownames(top), pre_cols, drop = FALSE], na.rm = TRUE)
top$gene <- anno[rownames(top), "UCSC_RefGene_Name"]

IMMUNE_GENES <- c("TNF","TNFRSF1A","TNFRSF1B","IL1B","IL1A","IL6","IL6R","IL10",
  "NFKB1","NFKB2","RELA","RELB","IKBKB","IKBKE","NFKBIA","TLR2","TLR4","TLR7","TLR9",
  "MYD88","TRAF6","CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CCL2","CCL3","CCL5",
  "IFNG","STAT1","STAT3","IRF1","IRF3","IRF7","CD68","ITGAM","LYZ","FCGR1A",
  "HLA-DRA","HLA-DRB1","PTGS2","NOS2","IL6ST","MAP3K7","MYC","BCL3","TNFAIP3","RELN")
probe_genes <- strsplit(anno[rownames(top), "UCSC_RefGene_Name"], "[;/]+")
immune_mask <- sapply(probe_genes, function(g) any(na.omit(g) %in% IMMUNE_GENES))
mw <- wilcox.test(top$t[immune_mask], top$t[!immune_mask], alternative = "two.sided")
mw_p <- mw$p.value
message("[fig] immune probes=", sum(immune_mask), " MWU p=", format.pval(mw_p, digits = 3))

# save full top + mask for any later tweaks (no recompute needed)
saveRDS(list(top = top, immune_mask = immune_mask, mw_p = mw_p),
        file.path(GSE, "layerA_top_full.rds"))

dmp <- top[top$adj.P.Val < 0.05 & abs(top$deltaBeta) > 0.05, ]

png(file.path(GSE, "layerA_figure.png"), width = 7, height = 3.4, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 2.6, 1), mgp = c(2.2, 0.7, 0))
# (A) volcano
plot(top$t, -log10(top$P.Value), pch = 20, cex = 0.22,
     col = ifelse(immune_mask, "#C0392B", "#C7CBD1"),
     xlim = c(-12, 12), ylim = c(0, max(-log10(top$P.Value), na.rm = TRUE)),
     xlab = "limma t (Post - Pre)", ylab = expression(-log[10]~P),
     main = "(A) Within-POD methylation volcano", cex.main = 0.85, cex.lab = 0.8, cex.axis = 0.7)
points(top$t[immune_mask], -log10(top$P.Value)[immune_mask], pch = 20, cex = 0.5, col = "#C0392B")
abline(v = c(-2, 2), h = -log10(0.05), lty = 2, col = "steelblue")
if (nrow(dmp) > 0) {
  dmp <- dmp[order(dmp$t), ]
  lab_y <- -log10(dmp$P.Value)
  for (s in c(-1, 1)) {                       # stagger same-side labels >= 0.32 apart
    idx <- which(sign(dmp$t) == s)
    if (length(idx) > 1) {
      oo <- idx[order(lab_y[idx])]
      for (k in 2:length(oo)) if (lab_y[oo[k]] - lab_y[oo[k - 1]] < 0.32)
        lab_y[oo[k]] <- lab_y[oo[k - 1]] + 0.32
    }
  }
  x_end <- dmp$t + ifelse(dmp$t < 0, -0.35, 0.35)
  segments(dmp$t, -log10(dmp$P.Value), x_end, lab_y, col = "gray45", lwd = 0.5)
  text(x_end, lab_y, labels = sub(";.*", "", dmp$gene),
       pos = ifelse(dmp$t < 0, 2, 4), cex = 0.55, col = "black")
}
legend("topleft", legend = c("immune-gene probe", "background"), pch = 20,
       col = c("#C0392B", "#C7CBD1"), cex = 0.6, bty = "n")
# (B) t-statistic density
di <- density(top$t[immune_mask], na.rm = TRUE, from = -12, to = 12)
db <- density(top$t[!immune_mask], na.rm = TRUE, from = -12, to = 12)
plot(db, col = "gray50", lwd = 1.6, xlim = c(-12, 12),
     main = "(B) t-statistic distribution", xlab = "limma t (Post - Pre)",
     ylab = "density", cex.main = 0.85, cex.lab = 0.8, cex.axis = 0.7)
lines(di, col = "#C0392B", lwd = 1.6)
abline(v = 0, lty = 2, col = "gray40")
legend("topleft", legend = c("immune", "background"), col = c("#C0392B", "gray50"),
       lwd = 1.6, cex = 0.6, bty = "n")
mtext(paste0("immune shift (two-sided MWU) p = ", format.pval(mw_p, digits = 3)),
      side = 3, line = 0.2, adj = 1, cex = 0.7)
dev.off()
cat("FIGURE_WRITTEN\n")
