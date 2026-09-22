# pod-immune-hub-repositioning

Reproducible analysis package for:

> Yang Y. A two-axis model of postoperative delirium: APOE ε4 constitutive susceptibility and a peripheral immune state axis supported by multi-omics integration and in-silico drug repositioning. *Journal of Neuroinflammation* (submitted).

This repository contains the analysis code, derived result tables, the LINCS
gene-set libraries used, and a `MANIFEST.sha256` checksum file for the
multi-omics + in-silico drug-repositioning study of postoperative delirium
(POD). It accompanies the manuscript `manuscript/POD_immune_hub_manuscript_v1.md`
(and its rendered `manuscript/POD_immune_hub_manuscript_v1.docx`).

## Data sources (all public)

| Layer | Source | Accession / DOI |
|-------|--------|-----------------|
| Blood transcriptome | GSE163943 (GPL26963, 4 POD vs 4 non-POD) | NCBI GEO |
| PBMC scRNA-seq | GSE252572 (4 patients, pre/post 24 h paired) | NCBI GEO |
| Blood methylation | GSE330869 (EPIC v2, 65 paired POD cases, recomputed in-house) | NCBI GEO |
| POD GWAS | Armstrong et al., *PLOS Medicine* 2026 | DOI 10.1371/journal.pmed.1004963; summary-statistics data DOI 10.5523/bris.1m83zai2e26yq2lro3tixz9kqq |
| LINCS L1000 | Enrichr / maayanlab.cloud consensus gene sets | maayanlab.cloud |

## File map

- `analysis/` — layer scripts (Layer A methylation/deconvolution, B
  transcriptome, C scRNA-seq, D LINCS repositioning, E docking audit, F
  GWAS/MR) and `results/` derived tables + figures.
- `analysis/adjust/` — composition-adjustment sensitivity analyses.
- `analysis/signatures/` — published POD anchor genes.
- `manuscript/` — manuscript source (`.md`) and rendered (`.docx`).
- `MANIFEST.sha256` — checksums of all tracked files.

## Reproduce

- **R 4.4.3**: minfi, limma, glmnet, EpiDISH, bumphunter, AnnotationHub,
  ExperimentHub.
- **Python 3.13**: numpy, scipy, pandas.
- Run scripts per layer; see each script header for inputs/outputs. Layer A
  recomputes on the GEO-deposited GSE330869 EPIC v2 beta matrix (312,514
  probes after author filtering).

## Caveats (read before reuse)

- DMR calling used `bumphunter` with `B = 0` (no permutation FDR); DMR counts
  are **descriptive**, not FDR-controlled.
- LINCS enrichment uses **consensus overlap** (Enrichr LINCS GMT), not
  full-ranked CMap τ/ξ connectivity scores.
- The POD GWAS signal is **monolithic on chromosome 19** (APOE ε4);
  off-chromosome-19 immune-gene-set enrichment is null.
- Layer B (n = 8), Layer C (n = 4 patients), and Layer A (n = 65 paired) are
  underpowered for gene-level FDR; results are framed as hypothesis-generating.

## License

MIT — see `LICENSE`.

## How to cite

See `CITATION.cff`. A Zenodo archive with a versioned DOI will be linked here
upon acceptance.
