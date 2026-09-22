# T1-3 supplemental sensitivity analysis (Round-8 panel request).
#
# WHY: the manuscript's "46-61% of DMRs are composition-attributable" band mixes two
# DIFFERENT covariate specifications (DMR level used ~aNeu alone -> -61%, and
# ~aNeu+aEos+aBaso -> -46%). The probe-level composition adjustment, by contrast, was
# run with the three-lineage set (Neu + Lymph + Mono). Here we re-run the DMR level
# with that SAME three-lineage covariate set, so probe level and DMR level become
# directly comparable under one specification.
#
# Self-check: the unadjusted model is re-run from scratch here; it must reproduce
# 13,357 DMRs (L>=3), 79 immune-overlapping DMRs, 48 distinct immune genes. If it
# does not, the new 3-lineage numbers must not be reported.
#
# NOTE: bumphunter runs with B = 0 (no permutation); all DMR counts are descriptive.
.libPaths(c("C:/Users/Administrator/R/win-library/4.4", .libPaths()))
ROOT <- "D:/podprj"; GSE <- file.path(ROOT, "analysis/results/GSE330869")
set.seed(12345)
logf <- file.path(GSE, "layerA_bumphunter_3lineage.log")
sink(logf, split = TRUE)

suppressMessages({ library(limma); library(bumphunter); library(GenomicRanges); library(minfi) })

cat("=== building per-subject paired M-value deltas ===\n")
beta <- readRDS(file.path(GSE, "GSE330869_beta_matrix_lean.rds"))
beta_cols <- grep("_beta$", colnames(beta), value = TRUE)
sheet <- data.frame(
  clean   = sub("_beta$", "", beta_cols),
  patient = sub("_(Pre|Post)$", "", sub("_beta$", "", beta_cols)),
  time    = ifelse(grepl("_Pre_beta$", beta_cols), "Pre", "Post"),
  stringsAsFactors = FALSE)

bb <- pmin(pmax(beta[, beta_cols, drop = FALSE], 1e-4), 1 - 1e-4)
M  <- log2(bb / (1 - bb)); rm(beta, bb); invisible(gc())
colnames(M) <- sheet$clean

fr <- read.csv("D:/podprj/analysis/deconv_ref_ready/subject_delta_fractions.csv",
               stringsAsFactors = FALSE, check.names = FALSE)
subjects <- fr[, 1]
available <- intersect(subjects, unique(sheet$patient))
for (s in available) {
  need <- c(paste0(s, "_Post"), paste0(s, "_Pre"))
  if (!all(need %in% colnames(M))) available <- setdiff(available, s)
}
cat("subjects usable:", length(available), "of", length(subjects), "\n")

Dmat <- matrix(NA_real_, nrow = nrow(M), ncol = length(available),
               dimnames = list(rownames(M), available))
for (s in available) Dmat[, s] <- M[, paste0(s, "_Post")] - M[, paste0(s, "_Pre")]
rm(M); invisible(gc())
cat("Delta matrix:", paste(dim(Dmat), collapse = " x "), "\n")

fr2 <- fr[fr[, 1] %in% available, ]
rownames(fr2) <- fr2[, 1]
fr2 <- fr2[available, , drop = FALSE]
for (cc in setdiff(colnames(fr2), fr2[, 1][1]))
  fr2[[cc]] <- as.numeric(fr2[[cc]])

# three-lineage covariate set, spelled exactly as in the probe-level MB model
COV3 <- c("aNeu", "lymphoid", "aMono")
stopifnot(all(COV3 %in% colnames(fr2)))
cat("3-lineage covariates present:", paste(COV3, collapse = " + "), "\n")

cat("\n=== annotation (chr/pos/gene) ===\n")
ann <- minfi::getAnnotation("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
kp  <- rownames(ann) %in% rownames(Dmat)
chr <- as.character(ann$chr[kp]); pos <- as.numeric(ann$pos[kp])
names(chr) <- rownames(ann)[kp]; names(pos) <- rownames(ann)[kp]
gname <- as.character(ann$UCSC_RefGene_Name[kp]); names(gname) <- rownames(ann)[kp]
rm(ann); invisible(gc())
chr <- chr[rownames(Dmat)]; pos <- pos[rownames(Dmat)]; gname <- gname[rownames(Dmat)]
keep2 <- !is.na(pos) & !is.na(chr)
Dmat <- Dmat[keep2, ]; chr <- chr[keep2]; pos <- pos[keep2]; gname <- gname[keep2]
cat("final matrix:", paste(dim(Dmat), collapse = " x "), "\n")

ig <- readLines("D:/podprj/analysis/signatures/immune_genes.txt")
ig <- trimws(ig); ig <- ig[nzchar(ig) & !grepl("^#", ig)]
cat("immune gene set:", length(ig), "\n")

pr_gr <- GRanges(seqnames = sub("^chr", "", chr), IRanges(pos, pos))

immune_by_dmr <- function(tb) {
  dgr <- GRanges(seqnames = sub("^chr", "", as.character(tb$chr)),
                 IRanges(tb$start, tb$end))
  ov  <- findOverlaps(dgr, pr_gr)
  qh  <- queryHits(ov); sh <- subjectHits(ov)
  sym <- list(); flag <- logical(nrow(tb))
  for (qi in unique(qh)) {
    ss <- unique(unlist(strsplit(as.character(gname[sh[qh == qi]]), ";")))
    ss <- ss[nzchar(ss) & !is.na(ss)]
    sym[[as.character(qi)]] <- ss
    flag[qi] <- any(ss %in% ig)
  }
  list(flag = flag, sym = sym)
}

runModel <- function(form, label) {
  cat("\n=================", label, "=================\n")
  design <- if (is.null(form)) matrix(1, nrow = length(available), ncol = 1,
                                      dimnames = list(NULL, "(Intercept)"))
            else model.matrix(as.formula(form), data = fr2)
  cat("design:", paste(colnames(design), collapse = " | "), "\n")
  bh  <- bumphunter(Dmat, design = design, chr = chr, pos = pos,
                    coef = 1, cutoff = 0.05, B = 0)
  tb  <- bh$table
  cat("raw bumps:", nrow(tb), "\n")
  tb3 <- tb[tb$L >= 3, ]
  cat("DMRs with L>=3:", nrow(tb3), "\n")
  res <- immune_by_dmr(tb3)
  cat("immune-overlapping DMRs:", sum(res$flag), "\n")
  hits <- table(unlist(lapply(res$sym[res$flag], function(s) intersect(s, ig))))
  cat("distinct immune genes hit:", length(hits), "\n")
  if (length(hits)) print(head(sort(hits, decreasing = TRUE), 12))
  write.csv(tb3, file.path(GSE, paste0("bumphunter_DMR_", label, ".csv")), row.names = FALSE)
  invisible(list(tb = tb3, flag = res$flag, hits = hits))
}

# 1) self-check: must reproduce the published unadjusted figures
r0 <- runModel(NULL, "unadjusted")
# 2) the requested unified three-lineage specification
rC <- runModel("~ aNeu + lymphoid + aMono", "adj_3lineage")

cat("\n\n############ SUMMARY ############\n")
summ <- data.frame(
  model        = c("unadjusted", "adj_3lineage"),
  covariates   = c("none (intercept only)", "aNeu + lymphoid + aMono"),
  n_DMR_L3     = c(nrow(r0$tb), nrow(rC$tb)),
  immune_DMR   = c(sum(r0$flag), sum(rC$flag)),
  immune_genes = c(length(r0$hits), length(rC$hits)),
  stringsAsFactors = FALSE)
summ$pct_immune    <- round(100 * summ$immune_DMR / summ$n_DMR_L3, 2)
summ$pct_DMR_left  <- round(100 * summ$n_DMR_L3 / nrow(r0$tb), 1)
summ$pct_DMR_removed <- round(100 - summ$pct_DMR_left, 1)
print(summ)
write.csv(summ, file.path(GSE, "bumphunter_3lineage_summary.csv"), row.names = FALSE)

cat("\n--- SELF-CHECK (unadjusted must equal 13357 / 79 / 48) ---\n")
cat("n_DMR_L3   :", nrow(r0$tb),   "(expect 13357)\n")
cat("immune_DMR :", sum(r0$flag),  "(expect 79)\n")
cat("imm_genes  :", length(r0$hits), "(expect 48)\n")
ok <- (nrow(r0$tb) == 13357) && (sum(r0$flag) == 79) && (length(r0$hits) == 48)
cat("SELF-CHECK PASSED:", ok, "\n")

shared <- intersect(names(r0$hits), names(rC$hits))
cat("\nimmune genes still hit after 3-lineage adjustment:",
    length(shared), "of", length(r0$hits), "\n")
cat("shared:", paste(sort(shared), collapse = ", "), "\n")
lost <- setdiff(names(r0$hits), names(rC$hits))
cat("lost  :", paste(sort(lost), collapse = ", "), "\n")
write.csv(data.frame(gene = names(r0$hits),
                     n_DMR_unadj   = as.integer(r0$hits[names(r0$hits)]),
                     n_DMR_adj3lin = as.integer(rC$hits[names(r0$hits)])),
          file.path(GSE, "bumphunter_immune_gene_before_after_3lineage.csv"),
          row.names = FALSE)

cat("\n=== DONE ===\n")
sink()
