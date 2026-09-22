# =============================================================================
# 00_config.R — project configuration for POD epigenetic immune repositioning
# Author: Yongxin Yang  |  Verified 2026-09-20
# All paths are relative to the repository root.
# =============================================================================

ROOT          <- here::here()                       # repo root
ANALYSIS_DIR  <- file.path(ROOT, "analysis")
RESULTS_DIR   <- file.path(ANALYSIS_DIR, "results")
FIG_DIR       <- file.path(RESULTS_DIR, "figures")
dir.create(c(RESULTS_DIR, FIG_DIR), recursive = TRUE, showWarnings = FALSE)

# ---- Cohorts (verified via E-utilities / GEO FTP, 2026-09-20) --------------
GEO_PRIMARY   <- "GSE330869"   # EPIC v2, 65 POD cases, Pre/Post paired, NO controls
GEO_EXPR      <- "GSE163943"   # RNA-seq n=8 (hub anchor)
GEO_SCRNA     <- "GSE252572"   # scRNA PBMC n=8 (hub localisation)
GEO_CSF       <- "GSE242736"   # CSF proteomics n=48
POD_GWAS_DOI  <- "10.5523/bris.1m83zai2e26yq2lro3tixz9kqq"  # Armstrong 2026, APOE-driven

# GSE330869 processed matrix (beta values); ~1.16 GB on first download
PROC_MATRIX_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/",
  GEO_PRIMARY, "/suppl/", GEO_PRIMARY, "_matrix_processed.csv.gz")

# ---- Platform / annotation (EPIC v2 — do NOT use 450K refs) ----------------
PLATFORM      <- "EPICv2"      # GPL33022
ANNO_PKG      <- "IlluminaHumanMethylationEPICv2anno.ilmn12.hg19"  # or ...hg38
# EpiDISH reference MUST be EPIC-compatible:
EPIDISH_REF   <- "FlowSorted.Blood.EPIC"   # loads with EpiDISH; EPIC v2 compatible

# ---- Curated inflammatory/immune gene set for the "burden score" -----------
# SEED list (NF-kB / leukocyte-immunity axis). REFINE from the actual top DMPs
# of this cohort + Seki 2026 (PMID 42143058) + MSigDB before final run.
IMMUNE_GENES <- c(
  "TNF", "TNFRSF1A", "TNFRSF1B", "IL1B", "IL1A", "IL6", "IL6R", "IL10",
  "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "IKBKE", "NFKBIA",
  "TLR2", "TLR4", "TLR7", "TLR9", "MYD88", "TRAF6",
  "CXCL1", "CXCL2", "CXCL8", "CXCL9", "CXCL10", "CCL2", "CCL3", "CCL5",
  "IFNG", "STAT1", "STAT3", "IRF1", "IRF3", "IRF7",
  "CD68", "ITGAM", "LYZ", "FCGR1A", "HLA-DRA", "HLA-DRB1",
  "PTGS2", "NOS2", "IL6ST", "MAP3K7", "MYC", "BCL3", "TNFAIP3", "RELN"
)

# ---- Published POD-vs-control anchor genes (Seki 2026, Time x Group model) --
# These come from the ORIGINAL PAPER (which had non-POD controls, n~33 preop),
# NOT from GEO (GEO GSE330869 deposited POD cases only). FDR-significant locus:
#   cg01534316 -> TMIGD3 / ADORA3 locus (p=2.41E-08, FDR=0.021)
#   cg05325643 -> SPATA13 (upstream regulatory region, FDR=0.075)
#   cg20382094 -> COL13A1 (gene body, FDR=0.075)
# ADORA3 (adenosine A3 receptor) is a plausible anti-inflammatory/neuroprotective
# hub; use as published anchors for the consensus signature (not a training set).
POD_CASE_CONTROL_ANCHORS <- c("TMIGD3", "ADORA3", "SPATA13", "COL13A1")

# ---- Pre->post immune shift pathways (Seki 2026 main finding, within POD) --
# Strongest at postoperative Day0, fades by Day3 -> transient, diagnostic potential.
IMMUNE_SHIFT_PATHWAYS <- c(
  "leukocyte mediated immunity", "NF-kappa B signaling pathway",
  "inflammatory response", "cytokine-mediated signaling"
)

# ---- Reproducible random seed ----------------------------------------------
SEED <- 12345

# ---- Software versions (pin at submission) ---------------------------------
VERSIONS <- list(
  R = "4.3", Bioconductor = "3.18",
  minfi = "1.48", ChAMP = "2.32", EpiDISH = "1.99",
  DMRcate = "2.16", limma = "3.58", missMethyl = "1.36",
  glmnet = "4.1", Python = "3.11"
)

message("[config] primary cohort = ", GEO_PRIMARY, " (", PLATFORM, ", paired Pre/Post, no controls)")
