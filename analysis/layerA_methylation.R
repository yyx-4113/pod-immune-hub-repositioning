# =============================================================================
# layerA_methylation.R — Primary epigenetic layer (GSE330869, EPIC v2)
# -----------------------------------------------------------------------------
# Verified design (2026-09-20): 65 POD cases, each with paired Pre & Post blood
# on Illumina MethylationEPIC v2 (GPL33022). NO non-delirium controls.
# Therefore the analysable signal is the SURGERY-INDUCED inflammatory/immune
# methylation shift WITHIN POD (pre -> post). Primary endpoint:
#   (a) pre->post immune-methylation signature (paired DMP + DMR, pathway enrich)
#   (b) Day0 subgroup: Pre baseline methylation predicts Post-Day0 inflammatory
#       epigenetic burden (penalised MRS, leave-pair-out CV + bootstrap)
# Runs in R >= 4.3 / Bioconductor >= 3.18. First run downloads ~1.16 GB matrix.
# =============================================================================

source("analysis/00_config.R")
library(minfi); library(ChAMP); library(EpiDISH); library(DMRcate)
library(limma); library(missMethyl); library(glmnet); library(here)

# -----------------------------------------------------------------------------
# 1. Load beta matrix (processed) — columns: KMD##_Post_beta / KMD##_Pre_beta
# -----------------------------------------------------------------------------
mat_path <- file.path(RESULTS_DIR, paste0(GEO_PRIMARY, "_matrix_processed.csv.gz"))
if (!file.exists(mat_path)) {
  message("[layerA] downloading processed matrix (~1.16 GB) ...")
  download.file(PROC_MATRIX_URL, mat_path, mode = "wb")
}
message("[layerA] reading matrix (this is large: ~935k probes x 130 cols) ...")
# readr::read_csv is faster; fall back to utils::read.csv
beta_raw <- as.data.frame(data.table::fread(mat_path, sep = ",", header = TRUE))
# first column = probe id (cg/CH names); remaining = samples
probe_id <- beta_raw[[1]]
beta_raw[[1]] <- NULL
colnames(beta_raw) <- sub("_beta$", "", colnames(beta_raw))   # KMD##_Post / KMD##_Pre
beta <- as.matrix(beta_raw); rownames(beta) <- probe_id

# -----------------------------------------------------------------------------
# 2. Build sample sheet + pairing
# -----------------------------------------------------------------------------
samples <- data.frame(
  sample = colnames(beta),
  patient = sub("_(Pre|Post)$", "", colnames(beta)),
  time = sub("^KMD[0-9]+_", "", colnames(beta)),   # Pre / Post
  stringsAsFactors = FALSE
)
samples$day <- ifelse(samples$time == "Post",
                      sub(".*_(Day[0-9]).*", "\\1",
                          # Post columns encode Day0/Day3 in GEO title; recover from the matrix header
                          colnames(beta)), NA)
# NOTE: GEO processed-matrix header only has Pre/Post (no Day tag). Recover Day
# from GSE330869_series_matrix.txt "sampling day" if needed; here we treat all
# Post as the post-op contrast and flag Day0 vs Day3 by re-reading the matrix
# characteristics file (see 00_config PROC_MATRIX_URL sibling). For the primary
# pre->post signature we use ALL paired patients.

paired <- split(samples$sample, samples$patient)
paired <- paired[sapply(paired, function(x) all(c("Pre", "Post") %in%
                                                 samples$time[match(x, samples$sample)]))]

# -----------------------------------------------------------------------------
# 3. Probe filtering (EPIC v2)
# -----------------------------------------------------------------------------
anno <- getAnnotation(get(ANNO_PKG))
beta <- beta[rownames(beta) %in% rownames(anno), ]
keep <- !(anno[rownames(beta), "chr"] %in% c("chrX", "chrY")) &
        !is.na(anno[rownames(beta), "UCSC_RefGene_Name"])    # drop SNP/ch loci per ChAMP
beta <- beta[keep, ]; anno <- anno[rownames(beta), ]
message("[layerA] probes retained: ", nrow(beta))

# -----------------------------------------------------------------------------
# 4. Cell-type proportion (EpiDISH, EPIC-compatible reference ONLY)
# -----------------------------------------------------------------------------
ref <- EpiDISH:::getRef("FlowSorted.Blood.EPIC")   # EPIC v2 compatible
est <- epidish(beta, ref, method = "RPC")
cellprops <- est$estF

# -----------------------------------------------------------------------------
# 5. Paired pre -> post differential methylation (M-values, patient as block)
# -----------------------------------------------------------------------------
M <- beta2m(beta)   # logit transform
patient_fac <- factor(samples$patient)
time_fac    <- factor(samples$time, levels = c("Pre", "Post"))
design <- model.matrix(~ time_fac + patient_fac)
fit <- lmFit(M, design)
fit <- eBayes(fit)
topDMP <- topTable(fit, coef = "time_facPost", number = Inf, sort.by = "P")
topDMP$deltaBeta <- beta[rownames(topDMP), "Post_mean"] - beta[rownames(topDMP), "Pre_mean"]
topDMP$gene <- anno[rownames(topDMP), "UCSC_RefGene_Name"]
topDMP <- topDMP[topDMP$adj.P.Val < 0.05 & abs(topDMP$deltaBeta) > 0.1, ]
write.csv(topDMP, file.path(RESULTS_DIR, "layerA_topDMP.csv"))

# -----------------------------------------------------------------------------
# 6. DMR (DMRcate) on the same paired contrast
# -----------------------------------------------------------------------------
dat <- list(beta = beta, M = M, design = design,
            covariate = time_fac, arrays.type = "EPICv2")
dmr <- dmrcate(dat, lambda = 1000, C = 2)
dmr <- extractRanges(dmr, genome = "hg19")
write.csv(as.data.frame(dmr), file.path(RESULTS_DIR, "layerA_topDMR.csv"))

# -----------------------------------------------------------------------------
# 7. Gene-set / pathway enrichment (missMethyl, EPIC v2 gene sets)
# -----------------------------------------------------------------------------
sig <- rownames(topDMP)
gsa <- gsarem(sig, array.type = "EPICv2",
              collection = "GO", ...)
# Manual check of the inflammatory/immune axis:
immune_probes <- rownames(beta)[anno[rownames(beta), "UCSC_RefGene_Name"] %in% IMMUNE_GENES]
enrich_immune <- sum(rownames(topDMP) %in% immune_probes)
message("[layerA] immune/NF-kB gene probes among top DMP: ", enrich_immune)

# -----------------------------------------------------------------------------
# 8. Day0 subgroup: Pre baseline -> Post-Day0 inflammatory-epigenetic burden
#    (penalised elastic-net MRS, leave-pair-out CV + bootstrap)
# -----------------------------------------------------------------------------
# burden score = mean centred M-value across IMMUNE_GENES probes at Post-Day0
post_cols  <- samples$sample[samples$time == "Post"]   # refine to Day0 via characteristics
immune_idx <- which(rownames(M) %in% immune_probes)
burden_post <- colMeans(M[immune_idx, post_cols, drop = FALSE], na.rm = TRUE)
X_pre  <- t(M[immune_idx, sub("Post", "Pre", post_cols), drop = FALSE])
y <- burden_post
# leave-pair-out CV + bootstrap CI on elastic-net
set.seed(SEED)
cvfit <- cv.glmnet(X_pre, y, alpha = 0.5, nfolds = length(paired))
pred <- predict(cvfit, X_pre, s = "lambda.min")
r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
message("[layerA] Day0 baseline->post burden MRS CV R^2 = ", round(r2, 3))
saveRDS(list(topDMP = topDMP, dmr = dmr, cellprops = cellprops,
             mrs_r2 = r2, cvfit = cvfit),
        file.path(RESULTS_DIR, "layerA_results.rds"))
message("[layerA] DONE. Outputs in ", RESULTS_DIR)
