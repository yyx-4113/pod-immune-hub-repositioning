# C3 — Implementation / Provenance Audit (independent first-submission review)

**Manuscript:** `manuscript/POD_immune_hub_manuscript_v1.md` (v1.1, Round-7 revision, 2026-09-22)
**Auditor role:** Verify every quantitative claim against committed artefacts; recompute; catch double-rounding, cross-section inconsistencies, broken tables, unreproducible claims.
**Verdict (headline):** The large majority of numbers are faithfully reproduced from the committed artefacts. There are **two substantive provenance defects** — (i) the immune-gene counts 48/21/29 are not contained in the overlap CSVs the manuscript cites as their provenance, and three committed files disagree (48 vs 75 vs 79); (ii) the lymphoid −1.67 pp leucocyte-shift figure is not reproducible from the same reference/method and is contradicted by the author's own cover letter (−1.6 pp). Neither invalidates the central argument, but both must be corrected before acceptance. Full detail below.

---

## 1. DMR counts (13,357 / 5,153 / 7,263) and immune-overlapping DMRs (79 / 24 / 36)

### 1a. DMR totals and immune-overlapping DMR counts — VERIFIED
【Problem】 None for the DMR totals and immune-overlapping DMR counts themselves.
【Evidence】 `analysis/results/GSE330869/bumphunter_adjusted_summary.csv` rows 2–4 report `unadjusted 13357 / 79 / 48 / 0.59`, `adj_dNeu 5153 / 24 / 21 / 0.47`, `adj_dNeu_Eos_Baso 7263 / 36 / 29 / 0.5`. Row counts of the three overlap CSVs (`bumphunter_DMR_{unadjusted,adj_dNeu,adj_dNeu_Eos_Baso}_immune_overlap.csv`) are exactly 79 / 24 / 36; full DMR CSVs contain 13,357 / 5,153 / 7,263 rows. Manuscript lines 102 and 106 state 13,357 / 5,153 / 7,263 DMRs and 79 / 24 / 36 immune-overlapping DMRs, immune proportion 0.59% → 0.47–0.50%. All match.
【Why it matters】 These are the spine of Layer A's composition-adjustment story and they are correct.
【Specific fix】 None.

### 1b. Immune-gene counts (48 / 21 / 29) — DISCREPANCY (does NOT match the cited overlap CSVs; three artefacts disagree)
【Problem】 The manuscript asserts (line 50) that the per-DMR overlap CSVs "fix the manuscript's cited counts of 79 / 24 / 36 immune-overlapping DMRs **(and 48 / 21 / 29 immune genes)**". The overlap CSVs do **not** contain 48/21/29 genes. Each of the 79/24/36 overlap rows carries exactly one immune gene (`n_immune_genes` max = 1 across all three files), so counting distinct genes in those very files yields **79 / 24 / 36**, not 48 / 21 / 29.
【Evidence】
- `bumphunter_DMR_unadjusted_immune_overlap.csv` etc.: 79 / 24 / 36 rows, every row `n_immune_genes = 1` → distinct genes = 79 / 24 / 36.
- `bumphunter_adjusted_summary.csv` (rows 2–4) reports `immune_genes = 48 / 21 / 29`.
- `bumphunter_immune_gene_before_after.csv`: 48 genes with non-NA `n_DMR_unadj` (and only 14 with non-NA `n_DMR_adj`), i.e. the "48" source.
- `layerA_DMR_immune_gene_counts.csv`: 75 rows (75 distinct genes overlapping unadjusted DMRs).
So for the *unadjusted* set alone, three committed artefacts disagree: **48** (summary + before_after), **75** (DMR_immune_gene_counts), **79** (overlap CSV). The same pattern holds for the adjusted sets (21/24, 29/36).
【Why it matters】 A reader following the manuscript's provenance pointer ("the overlap CSVs fix 48/21/29") opens the overlap CSV and finds 79/24/36 — a direct contradiction. The gene-count provenance is genuinely muddled, and the manuscript's parenthetical claim is simply false. This is the single most important provenance defect in the paper.
【Specific fix】
1. Remove the parenthetical "(and 48 / 21 / 29 immune genes)" from line 50, or replace it with the accurate statement: "the overlap CSVs contain 79/24/36 immune-overlapping DMRs, each annotated with its overlapping immune gene; the count of *distinct* immune genes is 79/24/36 in those files, while `bumphunter_adjusted_summary.csv` reports 48/21/29 under a stricter Set-2 probe→gene criterion."
2. Explain the 48 vs 75 vs 79 divergence explicitly: state which gene-mapping/overlap rule produces each number (the overlap CSVs appear to use a broader probe→gene annotation than the strict 83-gene Set-2 filter that yields 48). Pick ONE definition for the paper body and footnote the others, or reconcile them in code. Do not present 48 as if it came from the overlap CSVs.

---

## 2. Deconvolution (neutrophil +3.38 pp p=7.6e-4; lymphoid −1.67 pp p=0.003; mono −1.11 pp p=1.2e-4)

### 2a. Neutrophil and monocyte — VERIFIED (exactly)
【Problem】 None.
【Evidence】 `analysis/layerA_deconv_summary.csv` (root), config `centCAB100i_a12`, method `CP_nnls_norm`: row `aNeu` → `median_delta = 0.0338375` (=3.38 pp), `wilcoxon_p = 0.0007551` (=7.6×10⁻⁴); row `aMono` → `median_delta = −0.0110814` (=−1.11 pp), `wilcoxon_p = 0.0001202` (=1.2×10⁻⁴). Independent recomputation from `analysis/layerA_cell_proportions_python.csv` (centCAB100i_a12, CP_nnls_norm, 65 paired subjects) gives neutrophil median delta = **+3.384 pp** and monocyte median delta = **−1.108 pp**, matching the manuscript exactly. The CP_penalty_K1000 cross-check (neutrophil +3.39 pp, p=8.3e-4) also matches manuscript line 106 ("+3.39 pp, p = 8.3×10⁻⁴").
【Why it matters】 These are the paper's positive-control result and they are solid.
【Specific fix】 None. (Minor note: the *other* `layerA_deconv_summary.csv` under `analysis/results/GSE330869/` has a different schema — columns `cell_type, paired_wilcox_p` with `Neu = NA` — and is NOT the source cited in the Methods. The Methods citation "rows where celltype = aNeu/Neutro and config = centCAB100i_a12" describes the **root** file, which is correct. Clarify in the text that the cited file is the root `analysis/layerA_deconv_summary.csv`, not the GSE330869 copy, to avoid confusion.)

### 2b. Lymphoid −1.67 pp (p = 0.003) — NOT REPRODUCIBLE; contradicts the cover letter
【Problem】 The manuscript's "total lymphoid falls by −1.67 pp (p = 0.003)" (line 106) is not recoverable from the same reference/method used for neutrophil and monocyte, and it disagrees with the author's own cover letter.
【Evidence】
- Recomputation from `analysis/layerA_cell_proportions_python.csv` (centCAB100i_a12, CP_nnls_norm — the same config/method that yields the cited neutrophil +3.38 and mono −1.11) gives lymphoid median delta = **−1.599 pp** (n=65), i.e. −1.60, not −1.67.
- The cover letter (`manuscript/cover_letter.md` line 11) states "total lymphoid −1.6 pp" — consistent with the recompute, **not** with the manuscript's −1.67.
- Across *all* methods/configs in `centCAB100i_a12`, the lymphoid median ranges −1.30 (RPC_irls) to −1.73 (CBS_svr_nu0.5); only CBS_svr_nu0.5 approaches −1.67 and it is a different solver than the one used for the other two lineages. No single consistent reference reproduces −1.67.
- The p = 0.003 for lymphoid is also not traceable to a committed artefact row: `layerA_deconv_summary.csv` has no "lymphoid" (or "Lymph") cell-type row, so the −1.67 / 0.003 pair is the least well-provenanced of the three leucocyte-shift numbers.
【Why it matters】 A ~0.07 pp discrepancy is small in absolute terms, but it is a visible internal inconsistency (manuscript vs cover letter) and the number is not reproducible from the stated source. For a result that is presented as a positive control, the provenance should be airtight.
【Specific fix】 Either (a) change the manuscript lymphoid figure to −1.6 pp to match CP_nnls_norm and the cover letter, or (b) explicitly define "total lymphoid" (which cell types are summed) and state the config/method that yields −1.67, then confirm that config is applied consistently. Also add the lymphoid Wilcoxon p (0.003) to `layerA_deconv_summary.csv` as a derived composite row so it is traceable.

### 2c. Abstract vs Figure 2 rounding — CONSISTENT
【Problem】 None.
【Evidence】 Abstract line 17: "neutrophils +3.4 pp, p ≈ 10⁻³". Figure 2A line 106: "+3.38 pp … p = 7.6 × 10⁻⁴". 3.38 → 3.4 and 7.6×10⁻⁴ ≈ 10⁻³. Consistent rounding.
【Specific fix】 None.

---

## 3. Figure 2C |t|-enrichment numbers (0.968 vs 0.903, p=2.5e-4; 0.650 vs 0.656, p=0.72; 0.677 vs 0.657, p=0.55)

【Problem】 None — verified, with one source caveat.
【Evidence】 `analysis/results/GSE330869/adjust_summary.csv` rows 2–7 give median |t|: unadjusted immune 0.9679 / background 0.9027; adj_dNeu 0.6500 / 0.6564; adj_3lineage 0.6769 / 0.6572 — matching the manuscript's 0.968/0.903, 0.650/0.656, 0.677/0.657. The MWU p-values are not in `adjust_summary.csv` (only medians); they live in `analysis/adjust/adjust_report.txt` §2: "M0_unadjusted MWU p=0.0002463", "MA_adj_dNeu MWU p=0.7244", "MB_adj_3lineage MWU p=0.55" — i.e. 2.5×10⁻⁴ / 0.72 / 0.55, exactly as reported.
【Why it matters】 The composition-attenuation of the |t| enrichment is a key claim and it reproduces.
【Specific fix】 Minor: add the MWU p-values to `adjust_summary.csv` (or cite `adjust_report.txt` in the manuscript) so the Figure 2C p-values are traceable from the primary table, not only from a log file.

---

## 4. Set-1 signed-t MWU p = 0.041 vs |t| enrichment p = 2.5e-4 — two distinct tests, correctly kept separate

【Problem】 None — the manuscript keeps these as two distinct tests.
【Evidence】 `analysis/results/GSE330869/layerA_summary.txt` line 48–49: `directional_mwu_p = 0.04127092` (Set-1, 737 immune probes, mean t immune +0.0644 vs background −0.0348 — matching manuscript line 102). The |t| enrichment p = 0.0002463 is the Figure 2C test (Section 3 above). Manuscript lines 102, 106, and 141 each refer to "Set 1, p = 0.041" and "the |t|-based enrichment test … distinct from the signed-t directional test." No conflation; line 141 explicitly says the directional signal is "not-yet-composition-validated." Consistent.
【Why it matters】 Confusing these two tests would be a serious error; it is avoided.
【Specific fix】 None.

---

## 5. GWAS (37 SNPs chr19; N_eff 134,310; λ=1.0133; BETA sign logic; OR=1.86 per ε4; 7 conditionally-independent = Armstrong's)

【Problem】 None — fully verified, no internal contradiction.
【Evidence】
- `analysis/results/layerF/layerF_gwas_top_loci.csv`: 37 rows, every `CHROM = 19` → 37 genome-wide-significant SNPs, all chromosome 19. Matches manuscript lines 17, 116.
- `layerF_real_gwas_summary.json`: `n_variants = 12,353,257`, `genomic_inflation_lambda = 1.01325566` → manuscript λ = 1.0133 (line 72, 116) ✓; `lead_snp.beta = −0.621355`, `lead_snp_odds_ratio_per_allele0 = 1.8614486` → manuscript OR = exp(0.621) = 1.86 ✓; `lead_snp.n = 134310` → manuscript N_eff 134,310 ✓.
- BETA sign logic is internally consistent: the JSON note states "Lead-SNP BETA is for ALLELE1 … OR for the ALLELE0 allele is exp(−BETA). rs429358 ALLELE0=C is the APOE ε4-defining allele." Manuscript line 116: "BETA = −0.621 … is the effect of the T (ε3) allele; the C (ε4) allele therefore has log-OR = +0.621, OR = exp(0.621) = 1.86." Both say the ε4 (C) allele OR = exp(0.621) = 1.86. No contradiction between "BETA is effect of T allele" and "OR = 1.86 per ε4 allele."
- "7 conditionally-independent lead SNPs" (line 116) is explicitly attributed to Armstrong ("Armstrong reported 7 conditionally-independent lead SNPs … We did not re-derive conditional independence") and is correctly treated as distinct from the 37 genome-wide-significant SNPs. No conflation.
- Immune-set test off chr19: `immune_gene_set_test_chr19_excluded.perm_p_one_sided = 0.2769446` → manuscript p = 0.277 (lines 17, 116) ✓. All-autosomes p = 0.0812 → manuscript "p = 0.081" (line 116) ✓.
【Why it matters】 The GWAS layer is the genetic cornerstone; its numbers are clean.
【Specific fix】 None.

---

## 6. LINCS (glucocorticoid 0.041/0.080/0.056; statin 0.85/0.80; NSAID 0.78/0.93; 6 IFNγ FDR<0.05; KO 1,109 vs 0; chem FDR0 vs ligand 24 not conflated)

### 6a. Glucocorticoid / statin / NSAID class tests — VERIFIED (one minor rounding note)
【Problem】 Minor: the glucocorticoid two-sided p is reported as 0.080 but recomputes to 0.0814.
【Evidence】 `layerD_immune_control_class_tests.csv` row 2 (glucocorticoids): `mean_percentile 58.46`, `null_mean 49.92`, `z 1.751`, `perm_p_one_sided 0.0406979651` → manuscript one-sided 0.041 ✓; two-sided = 2 × 0.040698 = **0.0814**, manuscript reports 0.080 (lines 17, 110). `layerD_whole_tx_control_class_tests.csv` row 2: glucocorticoids `perm_p_one_sided 0.0555972` → manuscript whole-tx 0.056 ✓. Statins: immune 0.8532 → 0.85 ✓, whole-tx 0.7982 → 0.80 ✓. NSAIDs: immune 0.7839 → 0.78 ✓, whole-tx 0.9311 → 0.93 ✓. The high per-drug percentiles (dexamethasone 95.1%, hydrocortisone valerate 92.0%, betamethasone acetate 90.9%, prednisolone hemisuccinate 87.2% on whole-tx; dexamethasone acetate 76.6%, hydrocortisone 79.4%, betamethasone 27.4%, prednisolone 76.3% on immune-restricted) all match the members lists exactly (lines 110).
【Why it matters】 The positive-control class logic is the credibility anchor of Layer D; it reproduces.
【Specific fix】 Change "two-sided p = 0.080" to "two-sided p = 0.081" (or state the two-sided value was computed by a method giving 0.080). Trivial.

### 6b. Six IFNγ perturbations — VERIFIED
【Problem】 None.
【Evidence】 `layerD_immune_ligand_signatures.csv` (six `IFNG-*` rows): ifng-mcf7 FDR = 1.9248e-4 (strongest), ifng-mdamb231 FDR = 3.3581e-4 (second-strongest), ifng-bt20 6.2073e-4, ifng-skbr3 3.7948e-3, ifng-mcf10a 2.2459e-3, ifng-hs578t 2.4622e-2 (weakest). All six FDR < 0.05. Manuscript lines 17, 110: "six IFNγ perturbations … FDR < 0.05; cited FDR = 3.4×10⁻⁴ = ifng-mdamb231, second-strongest; strongest ifng-mcf7 FDR = 1.9×10⁻⁴; weakest ifng-hs578t FDR = 2.5×10⁻²." All match.
【Why it matters】 The IFNγ reversal is the headline pharmacologic finding; correct.
【Specific fix】 None.

### 6c. KO arm 1,109 vs 0 — VERIFIED
【Problem】 None.
【Evidence】 `layerD_immune_crispr_ko_genes.csv`: 1,109 genes with `min_fdr < 0.05` (matches `layerD_immune_summary.json` `n_genelevel_fdr05 = 1109`). `layerD_whole_tx_crispr_ko_genes.csv`: 0 genes with `min_fdr < 0.05` (matches whole-tx `n_genelevel_fdr05 = 0`). Top genes (Discussion line 112): nupl2 5.81e-8, hla-dpa1 2.21e-7, itgb4 4.24e-7, pde4a/emb/plk5/il4r 4.16e-6 — all match `min_fdr` exactly. KO positive-control z/p: immune z=+0.427 p=0.342 (`layerD_immune_summary.json`) → manuscript "+0.43, p = 0.34" ✓; whole-tx z=−2.472 p=0.9934 → manuscript "z = −2.47, p = 0.9934" ✓. The four anchor genes (adora3, cd63, ltf, col18a1) are present in the KO control list (`layerD_immune_summary.json` / `layerD_whole_tx_summary.json` `controls`) as the manuscript states (line 146).
【Why it matters】 The 1,109-vs-0 asymmetry and the KO positive-control failure are correctly reported.
【Specific fix】 None.

### 6d. chem `n_druglevel_fdr05 = 0` vs ligand 24 — NOT conflated
【Problem】 None in the manuscript.
【Evidence】 `layerD_immune_summary.json`: chem `n_druglevel_fdr05 = 0`, ligand `n_druglevel_fdr05 = 24`, crispr_ko `n_genelevel_fdr05 = 1109`. The manuscript never states "24 ligand hits" in the body and consistently keeps the chemical arm (drug-level FDR) separate from the ligand sub-library; Discussion line 135 "no compound cleared FDR at the drug level" refers to the chemical arm (0), which is correct and not conflated with the ligand 24. The two sub-libraries are distinct as the Methods (line 64) require.
【Why it matters】 Avoiding this conflation is exactly what the a-priori design demanded; it is respected.
【Specific fix】 None.

---

## 7. Docking audit (ADORA3 5; COL18A1 9; CD63 2; LTF 10; TMIGD3/COL13A1/SPATA13 none)

【Problem】 None for the manuscript numbers. Minor internal artefact note (UniProt IDs differ between the two docking files — not in the manuscript).
【Evidence】 `layerE_dockability_audit.csv`: ADORA3 5 (8X16,8X17,9EBH,9EBI,9EHS), COL18A1 9, CD63 2, LTF 10, TMIGD3/COL13A1/SPATA13 0 PDB entries. `analysis/pdb_audit.json`: identical PDB-ID lists and `dockable_experimental` flags. Manuscript line 122 matches exactly. (Cross-check: `pdb_audit.json` gives TMIGD3 UniProt Q9H6F9 and SPATA13 Q9P2Y7, whereas `layerE_dockability_audit.csv` lists TMIGD3 P0DMS9 and SPATA13 Q96N96 — a UniProt mismatch between two committed artefacts, but neither UniProt ID is cited in the manuscript, so it does not affect a reported number.)
【Why it matters】 Docking is an audit only; the structural feasibility claims are correct.
【Specific fix】 Optional: align the UniProt IDs between `pdb_audit.json` and `layerE_dockability_audit.csv` for data hygiene (low priority; not a manuscript defect).

---

## 8. Cross-section consistency (Abstract / Results / Discussion / Limitations)

【Problem】 The major cross-section numbers are consistent; the only cross-section mismatch is lymphoid (−1.67 manuscript vs −1.6 cover letter).
【Evidence】
- Neutrophil pp & p: Abstract "+3.4 pp, p≈10⁻³" vs Figure 2 "+3.38 pp, p=7.6×10⁻⁴" — consistent rounding (§2c).
- DMR 79/24/36: stated in Results (line 102) and Figure 2D (line 106); not restated inconsistently elsewhere. Consistent.
- Glucocorticoid p: Abstract "0.041 / 0.080 / 0.056" = Results line 110 "0.041 / 0.080 / 0.056" = Discussion "class level on the immune-restricted arm (whole-transcriptome arm marginal, p = 0.056)". Consistent (aside from the 0.080→0.081 rounding, §6a).
- IFNγ 6 perturbations: Abstract "six IFNγ … best FDR = 1.9×10⁻⁴" = Results line 110 = Discussion "six IFNγ perturbations … ifng-mdamb231 FDR = 3.4×10⁻⁴, second-strongest". Consistent.
- GWAS OR 1.86 / 37 SNPs: Abstract = Results line 116 = Discussion line 130 ("OR ≈ 1.86 per allele"; "37 genome-wide-significant SNPs, all in the APOE/TOMM40/APOC1/PVRL2 block"). Consistent.
- Lymphoid: manuscript Figure 2A "−1.67 pp" vs cover letter "−1.6 pp" — INCONSISTENT (see §2b). This is the only number that disagrees across the submitted package's own documents.
【Why it matters】 Cross-section drift is a classic sign of unreconciled revisions; here it is isolated to one leucocyte-shift figure.
【Specific fix】 Reconcile lymphoid to −1.6 pp (or redefine and document −1.67) across manuscript and cover letter.

---

## 9. Reproducibility — are cited CSVs present and matching? Any broken markdown?

【Problem】 All cited CSVs exist and match (except the gene-count attribution in §1b). No broken markdown tables detected.
【Evidence】
- Present and matching: `bumphunter_adjusted_summary.csv`, the three `bumphunter_DMR_*_immune_overlap.csv` (79/24/36 rows), `layerA_deconv_summary.csv` (root), `layerA_summary.txt`, `adjust_summary.csv`, `set1_signed_t_adjusted.csv` (0.296/0.098/0.141 — matches manuscript line 17/141), `layerF_gwas_top_loci.csv`, `layerF_real_gwas_summary.json`, `layerD_immune_control_class_tests.csv`, `layerD_whole_tx_control_class_tests.csv`, `layerD_immune_crispr_ko_genes.csv`, `layerD_immune_ligand_signatures.csv`, `layerE_dockability_audit.csv`, `pdb_audit.json`, `layerB_summary.json`, `CITATION.cff`, `zenodo_metadata.json`, `manuscript/cover_letter.md`.
- `set1_signed_t_adjusted.csv` (adjust): M0 0.2958 / MA 0.0979 / MB 0.1414 → manuscript "0.296 / 0.098 / 0.141" ✓ (line 17, 141). Confirms the Set-1 signed-t vs adjusted-analogue distinction.
- The datasets table (lines 39–44), figure captions, and the GWAS/Results prose render as valid markdown; no broken table delimiters or unclosed code fences were found in the manuscript body.
- The single provenance mismatch is §1b (immune-gene counts 48/21/29 not in the cited overlap CSVs).
【Why it matters】 Reproducibility is the paper's stated strength; the gap in §1b is the one place it breaks.
【Specific fix】 Address §1b; everything else is present and reconciles.

---

## 10. cover_letter.md + CITATION.cff vs manuscript title ("two-axis hypothesis") and Zenodo DOI 10.5281/zenodo.22896443

【Problem】 None — title framing and DOI are consistent; no "model" vs "hypothesis" overclaim.
【Evidence】
- `manuscript/cover_letter.md` line 7: title identical to manuscript ("A two-axis hypothesis for postoperative delirium: APOE ε4 constitutive susceptibility and a peripheral immune state axis …"). Line 15: Zenodo DOI "10.5281/zenodo.22896443" — matches manuscript Data availability (line 167).
- `CITATION.cff` line 3: title "Postoperative Delirium two-axis hypothesis: APOE epsilon4 constitutive susceptibility and peripheral immune state axis — multi-omics integration and in-silico, hypothesis-generating drug-repositioning reproducibility package" — consistent "two-axis hypothesis" framing, no "model" wording.
- `zenodo_metadata.json` line 3/4: same "two-axis hypothesis" title; description echoes the manuscript title verbatim. Neither `CITATION.cff` nor `zenodo_metadata.json` embeds the DOI string — expected, since Zenodo assigns the DOI on upload (the manuscript is the only place it is stated, and it is internally consistent).
- The manuscript's circumspect language ("working two-axis hypothesis", "candidate", "hypothesis-generating", "not directly testable") is consistent across Abstract/Discussion/Conclusions and does not overclaim a "model" or established mechanism. No overclaim mismatch vs the metadata.
【Why it matters】 Title/DOI consistency between submission and archive is a basic reproducibility check; it passes.
【Specific fix】 None. (Optional: add the Zenodo DOI to `CITATION.cff` `identifiers` and to `zenodo_metadata.json` `related_identifiers` for completeness, since the manuscript quotes it — low priority.)

---

# § Stands up (things the manuscript does well, provenance-wise)

1. **DMR totals and immune-overlapping DMR counts (13,357/5,153/7,263 and 79/24/36)** are exactly reproducible from the committed CSVs — the compositional attenuation story (46–61% of DMRs removed) is solid.
2. **Deconvolution neutrophil (+3.38 pp, p=7.6e-4) and monocyte (−1.11 pp, p=1.2e-4)** reproduce to the third decimal from the cited root summary and from raw per-sample fractions.
3. **GWAS layer is internally flawless**: 37 SNPs all chr19, λ=1.0133, N_eff=134,310, and the BETA/OR sign logic (OR = exp(0.621) = 1.86 per ε4 allele) is self-consistent with no contradiction between "BETA is effect of T allele" and "OR=1.86 per ε4 allele."
4. **LINCS class tests and IFNγ/KO findings** all match the committed JSON/CSV artefacts (glucocorticoid 0.041; 6 IFNγ FDRs; KO 1,109 vs 0; top-gene FDRs exact).
5. **Layer B transcriptomic numbers** (737/996, MWU 0.0079, mean t 0.253 vs −0.171, 44,824 probes, 83 immune genes, 48/35 input genes) all verify against `layerB_summary.json`.
6. **Set-1 signed-t (0.041) vs |t|-enrichment (2.5e-4)** are correctly presented as two distinct tests and never conflated.

# § Questions for the authors

1. For the immune-gene counts: which probe→gene mapping yields 48/21/29, and why do `bumphunter_DMR_unadjusted_immune_overlap.csv` (79 distinct genes), `layerA_DMR_immune_gene_counts.csv` (75), and `bumphunter_adjusted_summary.csv` (48) disagree? Please reconcile to a single definition and correct line 50.
2. For lymphoid: what exact cell-type aggregation and config/method produces −1.67 pp (p=0.003)? CP_nnls_norm/centCAB100i_a12 gives −1.60 pp, and the cover letter says −1.6 pp. Please align the manuscript and cover letter.
3. Glucocorticoid two-sided p: is 0.080 a rounded 0.0814, or computed differently? Please state the method.
4. The `pdb_audit.json` vs `layerE_dockability_audit.csv` UniProt IDs for TMIGD3 (Q9H6F9 vs P0DMS9) and SPATA13 (Q9P2Y7 vs Q96N96) differ — which is correct? (Not a manuscript defect, but worth fixing in the archive.)
5. The Methods (line 50) cites `analysis/layerA_deconv_summary.csv` for the neutrophil Δ, but a same-named file exists under `analysis/results/GSE330869/` with a different schema (Neu = NA). Please confirm the cited path is the root file to avoid reader confusion.

# § What I actually checked (files, commands, recomputed vs manuscript, discrepancies)

**Files read/used (all under the allowed repo root):**
- `analysis/results/GSE330869/bumphunter_adjusted_summary.csv` — DMR 13357/5153/7263, immune_DMR 79/24/36, immune_genes 48/21/29, pct 0.59/0.47/0.5.
- `analysis/results/GSE330869/bumphunter_DMR_{unadjusted,adj_dNeu,adj_dNeu_Eos_Baso}_immune_overlap.csv` — 79/24/36 rows; `n_immune_genes` max=1.
- `analysis/results/GSE330869/bumphunter_immune_gene_before_after.csv` — 48 genes w/ n_DMR_unadj.
- `analysis/results/GSE330869/layerA_DMR_immune_gene_counts.csv` — 75 rows.
- `analysis/results/GSE330869/layerA_summary.txt` — directional_mwu_p 0.04127; 6 DMPs; mrs_cv_r2 0.7939; immune overlap 0/49 Fisher 1.
- `analysis/layerA_deconv_summary.csv` (root) — aNeu CP_nnls_norm median_delta 0.0338 p 0.000755; aMono −0.01108 p 0.000120; CP_penalty_K1000 neutrophil +0.0339 p 0.000830.
- `analysis/layerA_cell_proportions_python.csv` — recomputed paired medians: neutrophil +3.384 pp, monocyte −1.108 pp, lymphoid −1.599 pp (centCAB100i_a12, CP_nnls_norm, n=65).
- `analysis/adjust/adjust_summary.csv` + `adjust_report.txt` — |t| medians 0.968/0.903, 0.650/0.656, 0.677/0.657; MWU p 0.0002463 / 0.7244 / 0.55.
- `analysis/adjust/set1_signed_t_adjusted.csv` — 0.2958 / 0.0979 / 0.1414.
- `analysis/results/layerF/layerF_gwas_top_loci.csv` — 37 rows, all chr19; rs429358 BETA −0.621355.
- `analysis/results/layerF/layerF_real_gwas_summary.json` — λ 1.01326, OR 1.86145, N 134310, immune-set p 0.081 (all aut) / 0.277 (chr19 excl).
- `analysis/results/layerD/layerD_immune_control_class_tests.csv`, `layerD_whole_tx_control_class_tests.csv` — gluc 0.0407/0.0556; statin 0.8532/0.7982; NSAID 0.7839/0.9311.
- `analysis/results/layerD/layerD_immune_ligand_signatures.csv` — 6 IFNG FDRs (mcf7 1.9e-4, mdamb231 3.4e-4, hs578t 2.5e-2, …).
- `analysis/results/layerD/layerD_immune_crispr_ko_genes.csv`, `layerD_whole_tx_crispr_ko_genes.csv` — 1109 / 0 min_fdr<0.05; top-gene FDRs exact.
- `analysis/results/layerD/layerD_immune_summary.json`, `layerD_whole_tx_summary.json` — n_druglevel_fdr05 chem 0 / ligand 24; KO z/p immune +0.43/0.34, whole-tx −2.47/0.9934.
- `analysis/results/layerE/layerE_dockability_audit.csv`, `analysis/pdb_audit.json` — ADORA3 5, COL18A1 9, CD63 2, LTF 10, TMIGD3/COL13A1/SPATA13 0.
- `analysis/results/layerB/layerB_summary.json` — 44824 probes, 737/996, MWU 0.00787, mean t 0.253/−0.171, 83 immune genes, 48/35 input.
- `manuscript/cover_letter.md`, `CITATION.cff`, `zenodo_metadata.json`.

**Recomputed vs manuscript — matched:** DMR 13357/5153/7263; immune-overlap DMR 79/24/36; neutrophil +3.38/7.6e-4; mono −1.11/1.2e-4; CP_penalty_K1000 +3.39/8.3e-4; |t| medians & p 0.968/0.903/2.5e-4, 0.650/0.656/0.72, 0.677/0.657/0.55; Set-1 0.041; adjusted 0.296/0.098/0.141; GWAS 37 chr19, λ 1.0133, OR 1.86, N 134310, immune-set 0.081/0.277; gluc 0.041/0.056, statin 0.85/0.80, NSAID 0.78/0.93; 6 IFNγ FDRs; KO 1109/0 + top genes; docking PDB counts; Layer B 737/996/0.0079/0.253/−0.171/83/48-35.

**Recomputed vs manuscript — DID NOT MATCH:**
- Immune-gene counts 48/21/29: the cited overlap CSVs contain 79/24/36 distinct genes; a third file gives 75. Manuscript line 50 falsely attributes 48/21/29 provenance to the overlap CSVs. (§1b)
- Lymphoid −1.67 pp (p=0.003): recompute = −1.60 pp (CP_nnls_norm); cover letter says −1.6 pp. No consistent config yields −1.67. (§2b)
- Glucocorticoid two-sided p: manuscript 0.080 vs recomputed 0.0814. (§6a, trivial)

**Commands:** Python (managed 3.13) used to (a) count overlap-CSV rows and distinct genes, (b) recompute paired medians of neutrophil/mono/lymphoid delta from `layerA_cell_proportions_python.csv`, (c) count KO `min_fdr<0.05` and extract top-gene FDRs, (d) tabulate GWAS chromosomes and recompute OR=exp(0.621). Row counts and medians were cross-checked against the committed summary CSVs/JSONs by direct file reads.

---

# Overall verdict

The manuscript is, by single-author standards, unusually disciplined about provenance: almost every headline number (DMR totals, neutrophil/monocyte shift, GWAS OR and SNP count, LINCS class tests, IFNγ FDRs, KO 1,109, docking PDB counts, Layer B transcriptomics) is faithfully and exactly reproduced from the committed artefacts. The two defects are localized and fixable:

- **T2 (moderate, must fix before acceptance):** the immune-gene counts 48/21/29 are not in the overlap CSVs the manuscript cites as their provenance, and three committed files disagree (48 vs 75 vs 79). The manuscript's parenthetical claim at line 50 is false as written.
- **T2 (moderate):** the lymphoid −1.67 pp leucocyte-shift figure is not reproducible from the cited reference/method and contradicts the cover letter's −1.6 pp.

**Tier allocation**
- **T0 (show-stopper / would block):** none found.
- **T1 (major, should block unless convincingly answered):** none found that undermine a conclusion; §1b is the closest but concerns a secondary descriptive count, not a primary result. (If the editors treat provenance attribution strictly, §1b could be elevated to T1 because the manuscript explicitly directs verification to a file that contradicts it.)
- **T2 (moderate, fix required):** §1b (immune-gene count provenance, 48/75/79 disagreement); §2b (lymphoid −1.67 vs −1.60 / cover-letter −1.6).
- **T3 (minor / cosmetic):** §6a glucocorticoid two-sided 0.080 vs 0.0814; §3 add MWU p to `adjust_summary.csv`; §7 UniProt mismatch between the two docking artefacts; §10 optionally add DOI to CITATION.cff/zenodo_metadata.json; §2a clarify which `layerA_deconv_summary.csv` is cited.

**Recommendation:** Minor revision. Correct the two T2 items (and the trivial T3 rounding) and the paper is ready. The scientific narrative — a working two-axis hypothesis with honest null/negative reporting and a composition-attributable methylation signal — is well-supported by the artefacts I could examine.
