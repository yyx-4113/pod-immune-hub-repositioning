# OSF Preregistration — Postoperative Delirium Immune Hub & Virtual-Knockout Drug Repositioning

> **Status:** Preregistration plan (recommended, non-blocking). This document is the a-priori analytic blueprint for the manuscript *"Postoperative delirium two-axis hypothesis: APOE ε4 constitutive susceptibility and a peripheral immune state axis with in-silico drug repositioning"*.
> **Version:** 1.0.0 (2026-09-22)
> **Corresponding author:** Yongxin Yang (ORCID 0009-0004-9698-6552), B.M. — no graduate/professional title (e.g. MD/PhD) is asserted. Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China.
> **Code & data:** GitHub (public) `https://github.com/yyx-4113/pod-immune-hub-repositioning`; `MANIFEST.sha256` accompanies every tagged release; Zenodo DOI to be minted upon acceptance.
> **This is a computational reanalysis + in-silico repositioning study. No new human-subjects data are generated.**

---

## 1. Study identification

- **Title:** Postoperative delirium two-axis hypothesis: APOE ε4 constitutive susceptibility and a peripheral immune state axis with in-silico drug repositioning
- **Study type:** Observational reanalysis of public omics + in-silico perturbation (LINCS) + genetic epidemiology (MR/PRS). No randomization, no prospective enrollment in the analytic core.
- **Registration intention:** Preregister the a-priori hypotheses, analytic plan, decision gates, positive controls, and expected-null specifications **before** journal submission, to prevent post-hoc rationalization (p-hacking / HARKing).
- **Update rule:** Any deviation from this plan after registration will be logged in `PROGRESS.md` and in the GitHub release notes under an explicit "Deviations from preregistration" heading.

---

## 2. Background & a-priori hypotheses

Postoperative delirium (POD) has a strong but partly unexplained heritable component. A 2026 genome-wide association study (Armstrong et al., *PLOS Medicine*) localized essentially all genome-wide significant POD signal to the chromosome-19 *APOE* region, implying a **constitutive (trait) susceptibility axis**. Independently, perioperative systemic inflammation is a long-standing mechanistic candidate, implying a **state (candidate-modifiable) axis**. We pre-specify a **two-axis working hypothesis**:

- **H1 (constitutive / hypothesis-generating, to be confirmed via public GWAS):** The *APOE* ε4 allele is the dominant germline POD susceptibility locus; non-APOE immune loci do **not** survive genome-wide correction once the chr19 block is accounted for.
- **H2 (state / primary empirical hypothesis):** Major surgery induces a peripheral immune epigenomic shift in POD cases — specifically a neutrophil-predominant / lymphocyte-diminished leukocyte composition change accompanied by NF-κB / leukocyte-immunity CpG methylation movement — that is detectable as a **pre→post paired methylation signature** in whole blood.
- **H3 (repositioning / translational hypothesis):** In-silico knockdown (LINCS L1000) of genes in the pre→post immune signature can nominate **already-approved drugs** whose documented transcriptomic effect reverses the POD immune signature direction (i.e. virtual "knockout/normalization" of the activated state).

These hypotheses were fixed **before** finalizing effect sizes and were not derived from peeking at the primary endpoint results.

---

## 3. Data sources (pre-specified)

| Dataset | Modality | Role | Pre-specified use |
|---|---|---|---|
| **GSE330869** (EPIC v2) | DNA methylation, whole blood | **Primary discovery** | 65 POD cases, each with paired pre- and post-surgery samples. **No non-delirium control was deposited** (verified via GEO/E-utilities + GEO FTP, 2026-09-20). Analysis is strictly **within-POD pre→post**. |
| **GSE163943** (4 v 4) | Bulk transcriptomics, whole blood | Disease signature + anchoring | Immune gene-set competitive test; directional anchoring of methylation inference. |
| **GSE252572** (4 patients × pre/post PBMC) | scRNA-seq | Hub cell-type localization (descriptive) | Pooled clustering across samples (per-sample clustering makes cluster IDs non-comparable; 1 library corrupted/unavailable). |
| **Armstrong 2026 POD GWAS** | GWAS summary stats | Genetics core | UKB delirium, 1,016 cases / 139,148 controls; data DOI 10.5523/bris.1m83zai2e26yq2lro3tixz9kqq; publication DOI 10.1371/journal.pmed.1004963, PMID 41770756. |

All datasets are public; no individual-level data are accessed beyond what the repositories provide.

---

## 4. Variables

- **Exposure / state axis (Layer A):** perioperative status (pre vs post surgery), operationalized as the paired Δ (post − pre) in (i) leukocyte composition (EpiDISH/Houseman deconvolution) and (ii) CpG methylation (M-value).
- **Constitutive axis (Layer F):** *APOE* ε4 dosage (0/1/2 copies) and POD polygenic score.
- **Outcome (repositioning):** LINCS in-silico perturbation direction concordance — whether a drug's documented gene-expression signature opposes the POD immune signature.
- **Covariates / confounders (pre-specified):** estimated leukocyte composition (neutrophil, eosinophil, basophil, monocyte, lymphocyte fractions) entered as covariates in every methylation model, because the immune methylation signal is expected to be **partly explained by post-surgery leukocyte-composition shift** (stated a-priori, not discovered post-hoc).

---

## 5. Analysis plan (layer-by-layer, pre-specified)

### Layer A — Methylation (primary)
- Paired ΔM per CpG; DMP via moderated paired statistics; DMR via bump-hunting on paired Δ.
- **Cell-composition sensitivity (mandatory):** refit every model with estimated neutrophil (and eosinophil/basophil) fractions as covariates; report immune-probe |t| ratio and DMR count under unadjusted vs adjusted models. Per a-priori expectation, a substantial fraction of the immune methylation signal will attenuate after adjustment; this is reported honestly, not framed as cell-intrinsic reprogramming.
- Building block: subject is the pairing key, so the paired design cancels subject effects; models need only intercept + composition covariates.

### Layer B — Transcriptomics (anchoring)
- Competitive gene-set test (Mann–Whitney on immune-set vs background t-statistics). Full-transcriptome DE reported as exploratory (FDR likely null due to n=8 power); immune signal is claimed only at the gene-set level.

### Layer E — Docking (secondary, optional)
- Only genes with an available experimental structure (e.g. ADORA3 has 5 PDBs) enter docking/MD. Genes without structures (TMIGD3, COL13A1, SPATA13) are explicitly excluded from docking and noted as a limitation.

### Layer F — Genetics (MR / PRS)
- **Reverse MR (immune trait → POD):** bidirectional two-sample MR using Armstrong 2026 summary stats; pre-specified as **expected positive/feasible**.
- **Hub-gene → POD MR:** pre-specified as **expected null** (POD GWAS is APOE/chr19-monolithic); reported as a falsified/null result, not a killer.
- Immune gene-set enrichment in the GWAS: pre-specified as **expected null after chr19 removal** (permutation p).
- All p-interval boundaries: `5e-8 <= p < 1e-5` (strict lower-bound inclusive; written correctly to avoid silent zero-site returns). UCSC refGene spans >3 Mb are filtered before gene-level scoring (pseudo-record QC).

---

## 6. Decision gates (go/no-go) — pre-registered

| Gate | Trigger | Decision |
|---|---|---|
| **G1 Data verification** | Confirmed: EPIC v2, 65 POD paired, no control | Primary endpoint changed from "predict POD yes/no" (no negative class) to "pre→post immune signature + baseline-predicts-postoperative-inflammatory-epigenetic-burden (continuous)"; discovery set retained |
| **G2 Primary endpoint** | MRS internal CV AUC < 0.65 | Fall back to simpler immune score (EpiDISH fractions + logistic regression); remains positive |
| **G3 LINCS** | Zero candidates at strict threshold | Enable graded Top-N + positive-control-validated method → report hypothesis-generating list |
| **G4 Docking** | All candidates fail docking | Drop docking section; LINCS + genetics stand alone |
| **G5 Clinical** | In-house enrollment <40 or hub null | Use public-data direction + literature anchors for argument; clinical downgraded to exploratory |

These gates were fixed before analysis completion and are not revised retroactively to fit results.

---

## 7. Positive controls (pre-specified, mandatory)

1. **LINCS:** dexamethasone and statins (anti-inflammatory / immunomodulatory approved drugs) must appear among reversed-direction candidates; if they do not, the LINCS pipeline is judged failed and not reported as evidence.
2. **Genetics:** the *APOE* ε4 association (lead rs429358) must replicate the published direction and magnitude (OR ≈ 1.86 per ε4 allele; 1-copy OR ≈ 1.75, 2-copy OR ≈ 4.19). Failure to replicate indicates a data/QA problem.
3. **Methylation sentinel:** known surgery/inflammation-associated CpGs should move in the expected direction as a sanity check on the paired Δ pipeline.

---

## 8. Expected-null specifications (pre-registered expectations)

Stated up front so that null results are recognized as **confirmations of pre-specified expectations**, not disappointments:

- **E1.** Hub-gene → POD MR is null (chr19 monolithic, APOE-driven). Pre-registered; reported as such.
- **E2.** Immune gene-set enrichment in POD GWAS is null after chr19 removal (permutation p non-significant).
- **E3.** A meaningful portion of the Layer-A immune methylation signal is explained by post-surgery leukocyte-composition shift; the residual (composition-adjusted) signal is reported as hypothesis-generating, **not** as cell-intrinsic epigenetic reprogramming.
- **E4.** Subtype-level (CD4 naive/memory, Treg, B naive/memory) deconvolution is not identifiable due to lymphocyte-centroid collinearity; only lineage-level (neutrophil, total lymphocyte, monocyte) estimates are interpreted.

---

## 9. Statistical inference criteria (pre-specified)

- Genome-wide significance: `p < 5e-8`. Suggestive: `5e-8 ≤ p < 1e-5`.
- Transcriptomics / methylation gene-set: two-sided Mann–Whitney; nominal `p < 0.05` reported with effect direction; FDR q-value reported where computable; n-driven null FDR stated honestly.
- Repositioning ranking: cross-cell-line consensus + permutation p (>0 candidates required); positive-control co-occurrence required for pipeline validity.
- No peeking-based threshold adjustment after seeing results.

---

## 10. Reproducibility & data availability

- Public GitHub repository with tagged `v1.0.0` release, accompanied by `MANIFEST.sha256` (sorted SHA-256 of every tracked file minus itself).
- Zenodo versioned archive (DOI to be assigned upon acceptance) mirroring the tagged release.
- All scripts (R 4.4.3 + Python 3.13) and intermediate CSVs are committed; large raw matrices are excluded from the public package but their provenance is documented.
- Per journal EM policy, the manuscript selects "Yes / 是的" to the research-data question (public research data were analyzed and reported), with the repository URL + version tag + MANIFEST checksum as the data-availability statement.

---

## 11. Limitations (pre-noted, honesty)

- GSE330869 contains **POD cases only** — no POD classifier can be trained; the study is inherently within-POD (pre→post), not case-control.
- Layer B (n=8) and Layer C (n=4 patients) are underpowered; their results are exploratory/descriptive.
- `bumphunter` was run with B=0 (no bootstrap) in the sandbox environment (Rtools absent); DMR counts are therefore descriptive, not bootstrap-informed — stated explicitly.
- LINCS consensus overlap is a rank-based surrogate, **not** a Connectivity-Map τ score; drug nominations are hypothesis-generating.
- Cell-composition confounding is a central, pre-acknowledged limitation of whole-blood methylation.

---

## 12. Deviations & versioning

Any post-registration change to hypotheses, endpoints, gates, or inference criteria will be recorded under **"Deviations from preregistration"** in `PROGRESS.md` and the GitHub release notes, with date and rationale. Silent post-hoc changes are prohibited by this plan.

---

*Prepared 2026-09-22 by Yongxin Yang (ORCID 0009-0004-9698-6552). This preregistration is a living document mirrored to OSF upon registration.*
