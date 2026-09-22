# Provenance: per-DMR immune-gene overlap tables for the composition-adjusted DMR sets.
# Reproduces EXACTLY the mapping used by layerA_bumphunter_adjusted.R
#   (annotation: minfi::getAnnotation("IlluminaHumanMethylationEPICv2anno.20a1.hg38"),
#    immune set: analysis/signatures/immune_genes.txt, seed 12345)
# so that the manuscript's cited counts (79/24/36 DMRs; 48/21/29 genes) are traceable
# to a committed artifact, not just to bumphunter_adjusted_summary.csv.
.libPaths(c("C:/Users/Administrator/R/win-library/4.4", .libPaths()))
suppressMessages({ library(minfi); library(GenomicRanges) })

ROOT <- "D:/podprj"
GSE  <- file.path(ROOT, "analysis/results/GSE330869")
set.seed(12345)

# --- analyzed probe set: exactly the probes the bumphunter run saw (lean beta matrix) ---
cat("loading analyzed probe set from lean rds ...\n")
rds <- readRDS(file.path(GSE, "GSE330869_beta_matrix_lean.rds"))
probes <- rownames(rds); rm(rds); invisible(gc())
cat("analyzed probes:", length(probes), "\n")

# --- annotation identical to original pipeline ---
ann <- minfi::getAnnotation("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
ann <- ann[rownames(ann) %in% probes, ]
kp  <- !is.na(ann$pos) & !is.na(ann$chr)
ann <- ann[kp, ]
cat("annotated analyzed probes with chr/pos:", nrow(ann), "\n")

pr_gr <- GRanges(seqnames = sub("^chr", "", as.character(ann$chr)),
                 IRanges(as.numeric(ann$pos), as.numeric(ann$pos)))
names(pr_gr) <- rownames(ann)
gname <- as.character(ann$UCSC_RefGene_Name); names(gname) <- rownames(ann)

ig <- readLines(file.path(ROOT, "analysis/signatures/immune_genes.txt"))
ig <- trimws(ig); ig <- ig[nzchar(ig) & !grepl("^#", ig)]
cat("immune gene set (Set-2, 83-gene):", length(ig), "\n")

immune_by_dmr <- function(tb) {
  dgr <- GRanges(seqnames = sub("^chr", "", as.character(tb$chr)),
                 IRanges(tb$start, tb$end))
  ov  <- findOverlaps(dgr, pr_gr)
  qh  <- queryHits(ov); sh <- subjectHits(ov)
  n   <- nrow(tb)
  sym <- vector("list", n); flag <- logical(n)
  for (qi in unique(qh)) {
    ss <- unique(unlist(strsplit(as.character(gname[sh[qh == qi]]), ";")))
    ss <- ss[nzchar(ss) & !is.na(ss)]
    sym[[qi]] <- ss
    flag[qi]  <- any(ss %in% ig)
  }
  list(flag = flag, sym = sym)
}

sets <- c("unadjusted"            = "bumphunter_DMR_unadjusted.csv",
          "adj_dNeu"              = "bumphunter_DMR_adj_dNeu.csv",
          "adj_dNeu_Eos_Baso"     = "bumphunter_DMR_adj_dNeu_Eos_Baso.csv",
          "adj_3lineage"          = "bumphunter_DMR_adj_3lineage.csv")

for (nm in names(sets)) {
  tb  <- read.csv(file.path(GSE, sets[[nm]]), stringsAsFactors = FALSE)
  res <- immune_by_dmr(tb)
  tb$direction <- ifelse(tb$value > 0, "hyper", "hypo")
  imm_genes <- sapply(res$sym, function(s) paste(intersect(s, ig), collapse = ";"))
  n_imm     <- sapply(res$sym, function(s) length(intersect(s, ig)))
  out <- data.frame(chr = tb$chr, start = tb$start, end = tb$end,
                    value = round(tb$value, 4), L = tb$L, direction = tb$direction,
                    n_immune_genes = n_imm, immune_genes = imm_genes,
                    stringsAsFactors = FALSE)
  out <- out[res$flag, ]                 # keep only immune-overlapping DMRs
  fn  <- paste0("bumphunter_DMR_", nm, "_immune_overlap.csv")
  write.csv(out, file.path(GSE, fn), row.names = FALSE)
  genes_hit <- sort(unique(unlist(strsplit(imm_genes[imm_genes != ""], ";"))))
  cat(sprintf("%-18s -> %4d immune-overlapping DMRs, %3d distinct immune genes: %s\n",
              nm, nrow(out), length(genes_hit), paste(genes_hit, collapse = ",")))
}
cat("=== DONE ===\n")
