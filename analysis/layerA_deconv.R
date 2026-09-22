# Layer A — EPIC blood cell deconvolution (EpiDISH / Houseman-RPC)
# TARGET: an Rtools-equipped machine with a COMPLETE FlowSorted.Blood.EPIC binary
#         (the source-only install available on the sandbox lacks the sorted-sample
#          reference dataset, so this script cannot run there; see Rtools_handoff.md).
# Correct epidish signature: epidish(beta.m, ref.m, method = c("RPC","CBS","CP"))
#   ref.m must be a matrix: probes x cell-types (aggregated per-cell-type profile).
.libPaths(c("C:/Users/Administrator/R/win-library/4.4", .libPaths()))
suppressMessages({ library(limma); library(EpiDISH); library(FlowSorted.Blood.EPIC) })
ROOT <- "D:/podprj"; GSE <- file.path(ROOT, "analysis/results/GSE330869")

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
paired <- split(colnames(beta), sheet$patient[match(colnames(beta), sheet$clean)])
paired <- paired[sapply(paired, function(x) {
  tt <- sheet$time[match(x, sheet$clean)]; all(c("Pre", "Post") %in% tt) })]
use <- unlist(paired); sheet_u <- sheet[match(use, sheet$clean), ]
beta_u <- beta[, use, drop = FALSE]; rm(beta); gc()

# --- resolve the FlowSorted.Blood.EPIC reference object (present only in a complete binary) ---
ref.obj <- tryCatch(FlowSorted.Blood.EPIC, error = function(e) NULL)
if (is.null(ref.obj))
  ref.obj <- tryCatch(get("FlowSorted.Blood.EPIC", envir = asNamespace("FlowSorted.Blood.EPIC")),
                      error = function(e) NULL)
if (is.null(ref.obj)) stop("FlowSorted.Blood.EPIC reference object not available; run on Rtools machine")

ref.m <- NULL; ct <- NULL
if (inherits(ref.obj, "SummarizedExperiment")) {
  ref.m <- as.matrix(assay(ref.obj, assayNames(ref.obj)[1]))
  ct <- as.character(colData(ref.obj)$CellType)
} else if (is.list(ref.obj)) {
  ref.m <- ref.obj$beta; ct <- as.character(ref.obj$cellTypes)
}
if (is.null(ref.m) || is.null(ct)) stop("could not extract ref.m / cellTypes from reference")
cat("reference sorted-sample beta dim:", paste(dim(ref.m), collapse = "x"), "\n")

# aggregate to per-cell-type mean profile (probes x cell types) -> ref.m for RPC
cts <- unique(ct)
prof <- do.call(cbind, lapply(cts, function(c) rowMeans(ref.m[, ct == c, drop = FALSE], na.rm = TRUE)))
colnames(prof) <- cts; rownames(prof) <- rownames(ref.m)
cat("ref profile (probes x cell types):", paste(dim(prof), collapse = "x"), "\n")

out <- epidish(beta.m = beta_u, ref.m = prof, method = "RPC")
est <- out$est
# number of EPIC v1/v2 overlapping probes actually used by RPC
n_ref_probes <- if (!is.null(out$ref) && !is.null(out$ref$beta)) nrow(out$ref$beta) else NA
cat("cell types:", paste(colnames(est), collapse = ","), "\n")
cat("n samples:", nrow(est), " ref probes used (EPIC v1/v2 overlap):", n_ref_probes, "\n")

# paired Pre vs Post per cell type (Wilcoxon signed-rank within patient)
summ <- data.frame(cell_type = colnames(est),
                   mean_pre = NA, mean_post = NA, delta = NA, p_paired_wilcox = NA)
for (i in seq_along(colnames(est))) {
  ct.i <- colnames(est)[i]
  pre <- est[sheet_u$clean[sheet_u$time == "Pre"], ct.i]
  post <- est[sheet_u$clean[sheet_u$time == "Post"], ct.i]
  names(pre) <- sheet_u$patient[sheet_u$time == "Pre"]
  names(post) <- sheet_u$patient[sheet_u$time == "Post"]
  common <- intersect(names(pre), names(post))
  summ$mean_pre[i] <- mean(pre, na.rm = TRUE)
  summ$mean_post[i] <- mean(post, na.rm = TRUE)
  summ$delta[i] <- summ$mean_post[i] - summ$mean_pre[i]
  summ$p_paired_wilcox[i] <- if (length(common) >= 2) wilcox.test(pre[common], post[common], paired = TRUE)$p.value else NA
}
write.csv(est, file.path(GSE, "layerA_cell_proportions.csv"), row.names = TRUE)
write.csv(summ, file.path(GSE, "layerA_deconv_summary.csv"), row.names = FALSE)
cat("DECONV_DONE\n")
print(summ)
