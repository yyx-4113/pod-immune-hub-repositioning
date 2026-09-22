# Round-8 Independent Multi-Expert Review — Panel Brief

**Manuscript:** `manuscript/POD_immune_hub_manuscript_v1.md` (version 1.1, Round-7 revision, 2026-09-22)
**Target venue:** *Journal of Neuroinflammation* (JNI)
**Article type:** original research / hypothesis-generating multi-omics integration + in-silico drug-repositioning
**Author:** single-author (Yongxin Yang, B.M.)
**Subject:** postoperative delirium (POD); two-axis hypothesis (constitutive APOE ε4 susceptibility + a candidate perioperative peripheral immune-activation state axis); LINCS L1000 repositioning screen.

This is a **fresh, first-submission-style** review. The manuscript has been through several prior rounds; you must NOT read any of them.

---

## Independence discipline (MANDATORY)

**Forbidden to read (any of these disqualifies your review):**
- `REVIEW_*.md`, `RESPONSE_*.md`, `REVISION_*.md`, `review_round*/`, `review_panel/`, any `*_round*.md`
- `PROGRESS.md`, `G1_gate_result.md`, `G2_gate_result.md`, `G3_gate_result.md`, `G4_gate_result.md`
- `方案一_*.md`, `方案重构_*.md`, `方案*.md`
- the `.workbuddy/` directory (memory / logs)
- `submission_notes.md`, `author_verification_statement.md`
- `DEPOSIT_FILELIST.txt`, the `zenodo_deposit_v1.0.0/` folder
- **the other four experts' output files inside `review_round8/`** — do not read them; form your own judgement.

Do **not** assume the manuscript is mature or has passed prior review. Treat it as a first submission.
Every judgement must come from text or source data you read yourself.
Any claim in the manuscript that you CAN verify, you MUST verify.

---

## Output contract (for EVERY item)

Four mandatory parts:
- **【Problem】** one sentence
- **【Evidence】** pinned to file:line, or section + exact numbers; numbers you cite must be ones you recomputed yourself
- **【Why it matters】** concrete effect on conclusions / credibility / acceptance
- **【Specific fix】** a paste-ready English replacement sentence, or an explicit spec for a new analysis (variables, strata, output columns)

"Consider strengthening the discussion" is banned.

---

## Required end-of-report sections

- **§ Stands up (≥3, with evidence)** — things you suspected but found correct. Deliverable, not filler.
- **§ Questions for the authors** — what you need to know; do not guess answers.
- **§ What I actually checked** — files read, commands run, values recomputed vs the manuscript, with the discrepancy stated.

---

## Allowed source files (read what your role needs)

Repository root: `pod-epigenetic-immune-repositioning/`
- Manuscript: `manuscript/POD_immune_hub_manuscript_v1.md`
- `README.md`, `CITATION.cff`, `cover_letter.md`, `preregistration_plan_draft.md`, `zenodo_metadata.json`
- `analysis/results/GSE330869/`: `bumphunter_adjusted_summary.csv`, `bumphunter_immune_gene_before_after.csv`, `bumphunter_DMR_{unadjusted,adj_dNeu,adj_dNeu_Eos_Baso}.csv`, `bumphunter_DMR_{unadjusted,adj_dNeu,adj_dNeu_Eos_Baso}_immune_overlap.csv`, `layerA_topDMP.csv`, `layerA_core_snapshot.txt`, `layerA_summary.txt`, `layerA_deconv_summary.csv`, `layerA_cell_proportions_python.csv`, `layerA_DMR_immune_gene_counts.csv`
- `analysis/adjust/`: `adjust_immune_probes.csv`, `adjust_summary.csv`, `set1_signed_t_adjusted.csv`
- `analysis/results/layerB/`, `layerC/`, `layerD/`, `layerE/`, `layerF/` (all CSV/JSON)
- `analysis/signatures/`: `immune_genes.txt`, `published_pod_anchors.csv`
- `analysis/pdb_audit.json`
- Scripts: `analysis/layerA_bumphunter_adjusted.R`, `analysis/layerA_dmr_immune_overlap_provenance.R`, `analysis/layerA_set1_signed_t_adjusted.py`, `analysis/layerD_real_lincs.py`, `analysis/layerF_real_gwas.py`, `analysis/layerA_deconv_python.py`

**Recompute with Python** (managed: `C:/Users/Administrator/.workbuddy/binaries/python/versions/3.13.12/python.exe`) reading the committed CSVs — R is not assumed available to you. For web checks use WebSearch/WebFetch (EuropePMC, PubMed, Springer, Clarivate JCR, RCSB, PLOS, GEO are reachable; api.clue.io needs a key and is NOT needed since Enrichr consensus CSVs are committed locally).

---

## Environment traps (do not re-discover these the hard way)

- The deposited GSE330869 beta matrix has only **312,514 of ~935,000** EPIC v2 probes (author-filtered); any reference intersect is ~30–36% by design — not an address bug.
- `layerA_cell_proportions.csv` is a **degenerate** centBloodSub NNLS output (Neu/CD4T collapsed to 0) and is explicitly NOT the source of the reported neutrophil Δ; the reported Δ comes from `layerA_deconv_summary.csv`.
- bumphunter was run with **B = 0** (no permutation) → DMR counts are descriptive, no empirical FDR.
- Layer C scRNA: pooled clustering (not per-sample); 4 patients (2 POD / 2 non) × pre/post, 7 of 8 libraries usable.
- LINCS: "chem.n_druglevel_fdr05 = 0" means the *chemical* sub-library has 0 drug-level FDR hits; "24" refers to *ligand* perturbations, a different sub-library — do not conflate.
- KO arm: each gene = single consensus signature, `cell_line` empty for 5208/5210 → `n_cell_lines` is structurally 1, not measured replication.
- GWAS: Armstrong file has **37 SNPs at P<5e-8, all chr19** (not "7 lead SNPs" — 7 are Armstrong's conditionally-independent leads); N effective = 134,310; λ = 1.0133.

---

## Verdict scale

Give an overall verdict: **Accept / Minor revision / Major revision / Desk-reject**, plus a tier grade per finding (T0 conclusion-invalidating, T1 analysis to add, T2 wording, T3 format). Be concrete.
