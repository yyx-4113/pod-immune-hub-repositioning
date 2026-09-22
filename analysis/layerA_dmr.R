.libPaths(c("C:/Users/Administrator/R/win-library/4.4", .libPaths()))
suppressMessages({ library(limma); library(DMRcate) })
ROOT <- "D:/podprj"; GSE <- file.path(ROOT, "analysis/results/GSE330869")
set.seed(12345)

cat("loading lean beta rds ...\n")
beta <- readRDS(file.path(GSE, "GSE330869_beta_matrix_lean.rds"))
raw_cols <- colnames(beta)
beta_cols <- raw_cols[grep("_beta$", raw_cols)]
sheet <- data.frame(raw = beta_cols,
                    clean = sub("_beta$", "", beta_cols),
                    patient = sub("_(Pre|Post)$", "", sub("_beta$", "", beta_cols)),
                    time = ifelse(grepl("_Pre_beta$", beta_cols), "Pre", "Post"),
                    stringsAsFactors = FALSE)
beta <- beta[, beta_cols, drop = FALSE]; colnames(beta) <- sheet$clean

bb <- pmin(pmax(beta, 1e-4), 1 - 1e-4); M <- log2(bb / (1 - bb)); rm(beta); gc()
paired <- split(colnames(M), sheet$patient[match(colnames(M), sheet$clean)])
paired <- paired[sapply(paired, function(x) {
  tt <- sheet$time[match(x, sheet$clean)]; all(c("Pre", "Post") %in% tt) })]
use <- unlist(paired); sheet_u <- sheet[match(use, sheet$clean), ]
M_u <- M[, use, drop = FALSE]; rm(M); gc()
patient_fac <- factor(sheet_u$patient); time_fac <- factor(sheet_u$time, levels = c("Pre", "Post"))
design <- model.matrix(~ time_fac + patient_fac)
cat("M dim:", nrow(M_u), "x", ncol(M_u), "\n")

# Small annotation (chr,pos,strand) from minfi EPICv2 annotation -> avoid loading huge S4
anno_pkg <- "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
cat("loading annotation package:", anno_pkg, "\n")
ann <- minfi::getAnnotation(anno_pkg)
keep <- rownames(ann) %in% rownames(M_u)
ann_s <- data.frame(chr = ann$chr[keep], pos = ann$pos[keep], strand = ann$strand[keep])
rownames(ann_s) <- rownames(ann)[keep]; rm(ann); gc()
ann_s <- ann_s[rownames(M_u), , drop = FALSE]
cat("all.annot rows:", nrow(ann_s), "\n")

cat("cpg.annotate ...\n")
myAnno <- tryCatch(
  cpg.annotate(object = M_u, datatype = "array", design = design,
               coef = "time_facPost", arraytype = "EPICv2", all.annot = ann_s),
  error = function(e) {
    cat("arraytype EPICv2 failed, retrying with all.annot only:", conditionMessage(e), "\n")
    cpg.annotate(object = M_u, datatype = "array", design = design,
                 coef = "time_facPost", all.annot = ann_s) })
cat("cpg.annotate OK, annotated probes:", nrow(myAnno), "\n")

cat("dmrcate ...\n")
dmrs <- dmrcate(myAnno, lambda = 1000, C = 2)
res <- extractRanges(dmrs, genome = "hg38")
res_df <- as.data.frame(res)
cat("n_DMR:", nrow(res_df), "\n")
write.csv(res_df, file.path(GSE, "layerA_DMR.csv"), row.names = FALSE)

# summary: DMRs with |mean meth diff| thresholds + immune-gene overlap
res_df$abs_mean_diff <- abs(res_df$meanmeth.diff)
n_dmr <- nrow(res_df)
n_dmr_01 <- sum(res_df$abs_mean_diff > 0.01, na.rm = TRUE)
n_dmr_05 <- sum(res_df$abs_mean_diff > 0.05, na.rm = TRUE)
# immune gene overlap using anno gene names
IG <- c("TNF","TNFRSF1A","TNFRSF1B","IL1B","IL1A","IL6","IL6R","IL10","NFKB1","NFKB2","RELA",
  "RELB","IKBKB","IKBKE","NFKBIA","TLR2","TLR4","TLR7","TLR9","MYD88","TRAF6","CXCL1","CXCL2",
  "CXCL8","CXCL9","CXCL10","CCL2","CCL3","CCL5","IFNG","STAT1","STAT3","IRF1","IRF3","IRF7",
  "CD68","ITGAM","LYZ","FCGR1A","HLA-DRA","HLA-DRB1","PTGS2","NOS2","IL6ST","MAP3K7","MYC",
  "BCL3","TNFAIP3","RELN")
anno2 <- read.csv(file.path(GSE, "EPICv2_anno_small.csv"), stringsAsFactors = FALSE)
anno2 <- anno2[anno2$ID %in% rownames(M_u), ]
ov_immune <- 0
if ("overlapping.genes" %in% colnames(res_df) || "genes" %in% colnames(res_df)) {
  gcol <- if ("overlapping.genes" %in% colnames(res_df) "overlapping.genes" else "genes"
  ov_immune <- sum(sapply(res_df[[gcol]], function(s) any(na.omit(strsplit(as.character(s), "[;/]+")[[1]]) %in% IG)))
}
cat("n_DMR>0.01:", n_dmr_01, " n_DMR>0.05:", n_dmr_05, " immune_overlap_DMR:", ov_immune, "\n")
saveRDS(list(myAnno = myAnno, dmrs = dmrs, res = res_df),
        file.path(GSE, "layerA_dmr_results.rds"))
cat("DMR_DONE\n")
