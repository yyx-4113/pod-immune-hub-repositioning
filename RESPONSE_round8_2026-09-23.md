# Round-8 revision — point-by-point response

**Manuscript:** *A two-axis hypothesis for postoperative delirium: APOE ε4 constitutive susceptibility and a peripheral immune state axis supported by multi-omics integration, with an in-silico, hypothesis-generating drug-repositioning screen*
**Version:** v1.1 → **v1.2** (2026-09-23)
**Source of the requests:** independent five-expert panel (`REVIEW_round8_2026-09-22.md`; expert reports in `review_round8/C1–C5`).
**Verdict received:** Major revision (C1, C2), Revise (C4), Minor revision (C3, C5); no desk-reject.
**Scope of this revision:** text/claims/citations and one artefact amendment. **No new experiment, no new analysis layer.** One artefact (`analysis/adjust/adjust_summary.csv`) was extended with two columns so that the reported Mann–Whitney p-values become directly traceable.

---

## Tier 1 — hard gates (all completed)

| ID | Request | Action taken | Where |
|---|---|---|---|
| **T1-1** | Reference [41] Evered 2018 carried the wrong PMID (30179813 = an unrelated chlorinated-paraffin paper) and a fabricated author list ("Whittaker P, Pu Y"), while omitting DeKosky / Rasmussen / Oh / Crosby / Berger / Eckenhoff. | Entire entry replaced with the verified record: *Anesthesiology* 2018;129(5):872–879, doi:10.1097/ALN.0000000000002334, **PMID 30325806**, with the correct Nomenclature Consensus Working Group author list. | References [41] |
| **T1-2** | BMC/Springer Nature requires that large-language-model use be documented in the Methods (or an equivalent location); absent from manuscript and cover letter. | New subsection **"Use of artificial intelligence"** added under Declarations, stating the assistant was used for language polishing / literature organisation / drafting, that all analyses and interpretations were designed, executed and verified by the author (Y.Y.), that no AI tool generated primary data or ran autonomous inferential statistics, and that the author takes full responsibility. A matching paragraph was added to the cover letter. | Declarations; `manuscript/cover_letter.md` |
| **T1-3** | "46–61% of DMRs composition-attributable" was presented as a sensitivity envelope, but its two ends come from **two different covariate specifications** (neutrophil-only vs granulocyte-subtype); and the immune-overlapping proportion was actually **flat** (0.59% → 0.47–0.50%), so the removed burden is a *global* post-surgical shift, not a disproportionately immune one. | **(a) Text.** Re-framed at every occurrence (Abstract, Results, Figure 2 legend, Limitations, cover letter). The 46–61% band is now explicitly described as the **bracket produced by switching specifications, not a sensitivity interval for one model**; the flat immune proportion is stated wherever the percentage appears; the immune-specific conclusion is assigned to the probe-level \|t\| test (2.46 × 10⁻⁴ → 0.724 / 0.55), not to DMR counts; every DMR count carries "B = 0, no permutation, hence descriptive". **(b) Analysis.** The panel's suggested supplemental run was **also performed** — see the dedicated section below. | Abstract; Results (Layer A); Fig. 2 legend; Limitations; cover letter; new artefacts |
| **T1-4** | Set-1 (49 genes / 737 probes, limma moderated t, unadjusted, p = 0.041) and Set-2 (83 genes / 1,195 probes, OLS-on-delta t, 0.296 / 0.098 / 0.141) are **not comparable** (different gene set *and* different t); the hedge existed but was easy to miss. | Inline specification added at the first occurrence of p = 0.041 (Set-1, unadjusted, limma moderated t, no composition-adjusted analogue on that exact definition); **new Table 1** "Layer-A methylation signal — all reported tests (not interchangeable)" added, with columns *Test / Gene set / t-statistic definition / Composition-adjusted? / p*; Conclusions and Limitations repeat that Set-2 can neither confirm nor refute Set-1. | Results (Layer A); new **Table 1**; Limitations |

---

## Tier 2 — moderate revisions (all completed)

| ID | Request | Action taken |
|---|---|---|
| **T2-1** | Minimum one-sided Mann–Whitney p for n₁ = n₂ = 2 was written as 1/3; exact enumeration gives **1/6**. | Corrected to "1/6 ≈ 0.167" at patient level, with the library-level (4 vs 3) floor 1/35 ≈ 0.029 and the observed most-extreme p = 0.23. The "read as a floor, not biological evidence" framing is retained. |
| **T2-2** | Lymphoid change −1.67 pp could not be reproduced (recomputation −1.60 pp) and contradicted the cover letter's −1.6 pp. | Corrected to **−1.60 pp** in the Results and Figure 2 legend; cover letter already read −1.6 pp, so the two now agree. |
| **T2-3** | The 48 / 21 / 29 immune-gene counts were mis-attributed to the per-DMR overlap CSVs (which hold 79 / 24 / 36 rows); the 48 / 75 / 79 spread across three artefacts was unexplained. | Attribution corrected: the overlap CSVs contain 79 / 24 / 36 immune-overlapping DMRs (one immune gene per row); the distinct-gene counts 48 / 21 / 29 come from `bumphunter_adjusted_summary.csv` under the stricter Set-2 probe→gene rule; a sentence now states that the 48 / 75 / 79 spread reflects **different mapping rules, not a numerical error**, and that the Results use only 48 / 21 / 29. |
| **T2-4** | "No compound cleared FDR at the drug level" contradicted the ligand arm's own 24 drug-level hits. | Split explicitly: **0 of 1,069 small-molecule chemical perturbagens** cleared drug-level FDR in either signature arm, whereas **24 of 96 cytokine/ligand perturbations** did so in the immune-restricted arm (six IFNγ) and 0 in the whole-transcriptome arm. Ligand hits are endogenous signalling molecules, not repositionable drugs, and are read at pathway level. |
| **T2-5** | The IFNγ "inhibition is actionable / administration would exacerbate" inference was too strong and conflated type I and type II interferon pathways. | Downgraded to "consistent with exacerbation, pending experimental validation"; the tautology with the input ISG signature is stated; type II (IFNγ / IFNGR) and type I (IFNAR) arms are now distinguished, with the note that JAK–STAT blockade covers both whereas IFNAR antagonism is type-I-specific. |
| **T2-6** | Results subheading "chemical arm valid" over-promised relative to the body text. | Subheading changed to **"chemical arm partially supportive, virtual-knockout arm hypothesis-generating"**. |
| **T2-7** | Canonical neuroinflammation literature was missing from the mechanistic background. | Four references added and cited at the point where neuroinflammation and the periphery-to-brain route are introduced: **Cerejeira 2010** (*Acta Neuropathol* 119:737, PMID 20490615), **Cunningham & Maclullich 2013** (*Brain Behav Immun* 28:1, PMID 23088936), **Inouye 2014** (*Lancet* 383:911, PMID 24560607), **Terrando 2016** (*Front Immunol* 7:441, PMID 27822212, systemic HMGB1 neutralisation). |
| **T2-8** | "Peripheral immune activation state axis" was used as if it were a single measured endpoint. | Explicitly defined as a **working, integrative construct** spanning directional blood-transcriptomic immune upregulation, a largely leucocyte-composition-associated perioperative methylation shift, and in-silico immune-signature reversibility — not a single endpoint; each component flagged hypothesis-generating. |
| **T2-9** | "Genuine, modifiable target" overstated the evidence. | Softened: the state axis is a **candidate** therapeutic target whose modifiability **remains unproven and requires randomized interventional testing**; repositioning supplies hypotheses, not efficacy evidence. |
| **T2-10** | Reference [35] cited as "therapeutically actionable". | Softened to "a core, potentially tractable POD/PND mechanism", consistent with the cited review's own call for further clinical evidence. |
| **T2-11** | Ethics statement named no committee or waiver basis. | Expanded: GSE330869 under the submitting institutions' IRBs per its published ethics statement; Armstrong 2026 under UK Biobank's generic Research Tissue Bank approval, with summary statistics accessed via the University of Bristol Research Data Repository; secondary-analysis waiver restated. (An earlier draft asserted a "University of Bristol ethics committee", which we could not verify; it was replaced with the verifiable UK Biobank / repository formulation.) |
| **T2-12** | The post-hoc lower-tail p from the whole-transcriptome KO class test needed an explicit exploratory label, and the anchor-gene circularity needed restating next to z = −2.47. | Both done: the post-hoc lower-tail p is labelled "exploratory, not pre-specified" and the anchor-gene circularity is restated in the same paragraph; the "failed positive control / mimicry" framing is retained. |
| **T2-13** | bumphunter ran with B = 0, so every DMR count and every derived percentage is descriptive. | "(descriptive, B = 0, no permutation)" now accompanies each occurrence of 13,357 / 5,153 / 7,263 and each derived percentage, in the Abstract, Results, Figure 2 legend, Methods and Limitations. |
| **T2-14** | "one-sample t p = 0.18 is a lower bound on the directional signal" inverted the meaning of a p-value. | Reworded: the immune-gene t distribution is shifted upward relative to background (Mann–Whitney p = 0.0079) while the mean per-gene shift is not significant at n = 8 (p = 0.18), i.e. the signal is broad and directional rather than driven by a few large effects; the "lower bound" phrasing was deleted. |

---

## Tier 3 — minor / housekeeping (all completed)

| ID | Action taken |
|---|---|
| **T3-1** | [3] Hsiao 2023 given **PMID 36827351**; on verification the journal locator was also wrong (18(8) → **18(2)**:e0282214) and was corrected. |
| **T3-2** | [15] Li 2019 given **PMID 31253079** and PMCID PMC6599229. |
| **T3-3** | Glucocorticoid two-sided class p 0.080 → **0.081** (2 × 0.0407). |
| **T3-4** | `analysis/adjust/adjust_summary.csv` extended with `mwu_p_immune_vs_bg` (0.0002463 / 0.7244 / 0.5500) and `median_ratio_immune_over_bg` (1.072 / 0.990 / 1.030), so the reported Mann–Whitney p-values are directly traceable instead of living only in `adjust_report.txt`. Consequently the manuscript now quotes **2.46 × 10⁻⁴ → 0.724 / 0.55** instead of the rounded 2.5 × 10⁻⁴ → 0.72 / 0.55, in the Abstract, Results, Figure 2 legend, Limitations and Table 1. |
| **T3-6** | Manuscript footer bumped to **v1.2 (Round-8, 2026-09-23)** with an explicit version note: the Zenodo deposit is frozen at v1.0.0 and does not yet contain the Round-7/Round-8 provenance artefacts; the GitHub repository `yyx-4113/pod-immune-hub-repositioning` is authoritative and a v1.1.0 deposit will be cut before submission. (The repository URL in the footer was also corrected.) |
| **T3-7** | Abstract abbreviations minimised: MWU → Mann–Whitney (or removed), "pp" → "percentage points"; "scRNA" and "PBMC" expanded at first abstract use. **The abstract was then compressed from 453 to 342 words** to sit inside BMC's structured-abstract limit (the Correspondence page states 250 words; the Research-article cap could not be read directly because the Springer guideline page was served behind a verification challenge, so 350 was used as the working target). Journal-name citations were removed from the abstract, per the "do not cite references in the abstract" rule. |
| **T3-8** | `README.md` and `preregistration_plan_draft.md` aligned to the manuscript's softened term: "two-axis **model**" → "two-axis **hypothesis**"; "state (modifiable) axis" → "state (candidate-modifiable) axis". |
| **T3-10** | The referenced deconvolution artefact is now disambiguated: paired neutrophil estimates come from the repository-root `analysis/layerA_deconv_summary.csv` (celltype = aNeu/Neutro, config = centCAB100i_a12, CP_nnls_norm), with the degenerate `analysis/results/GSE330869/layerA_cell_proportions.csv` explicitly excluded as a source. |
| **T3-11** | Figure 2C now states that the |t|-based (magnitude) enrichment test is distinct from the signed-t directional test reported in the Results (Set-1, p = 0.041); Table 1 makes the distinction systematic. |
| **T3-12** | The all-autosome immune-set stratified permutation p = 0.081 (z = 1.26) is now reported in the Results, framed as "weak, chr19-driven", with the caveat that the stratification did not separately control SNP density. |
| **T3-13** | "nominal FDR" for the KO top genes corrected to **BH-adjusted FDR**. |
| **T3-14** | Footnote added explaining that `n_controls_found` = 18 counts every matched gene including aliases/multi-mapping whereas the committed `controls` array lists the 12 unique named controls actually scored. |
| **T3-15** | The "5208/5210 empty cell_line" figure could not be reproduced against the artefact, so it was replaced with the directly quotable values from `layerD_immune_summary.json → libraries.crispr_ko`: **n_signatures = 5,212, n_genes = 5,210** (i.e. one consensus signature per gene, no per-cell-line replication) rather than an unverifiable pair of integers. |
| **T3-16** | The drug-level FDR is now defined: BH-FDR over an empirical permutation p computed **one-tailed in the reversal direction** on each perturbagen's most-negative mean connectivity score — a *reversal-only* FDR, distinct from the signed `net_score` / rank-percentile used for the class-level tests, which carries no FDR. |
| **T3-17** | The one-sided p = 0.041 is explicitly flagged as borderline and not meeting a two-sided 0.05 threshold. |
| **T3-18** | Reference [22] (Mardani 2013) is already qualified in the text as a single-centre RCT with non-CAM delirium ascertainment, set against the null dexamethasone meta-analysis [15]; retained as-is. |
| **T3-19** | [23] Rebola 2011 (*J Neurochem* 117:100) and [24] Martí Navia 2020 (*Cells* 9:1739) strings checked and confirmed. |
| **T3-20** | [41] incidence range 10–50% and the POD / POCD / PND distinction confirmed as correctly cited; retained. |

---

## Two findings we did **not** act on, and why

1. **C4 reported that reference [6] (Seki 2026, PMID 42143058) could not be found.** We re-checked: the paper exists (*Translational Psychiatry* 2026;16(1):349, doi:10.1038/s41398-026-04067-6), consistent with the G1 gate evidence already in the project record. This was an index-lag false alarm; [6] is unchanged.
2. **C3 reported a contradiction among three DMR immune-gene artefacts (48 vs 75 vs 79).** On inspection all three numbers are correct and refer to different objects under different probe→gene mapping rules. The real defect was the mis-attribution on line 50, which is fixed under T2-3; no number was changed.

---

## Artefact changes in this round

| File | Change |
|---|---|
| `manuscript/POD_immune_hub_manuscript_v1.md` | v1.1 → v1.2; all T1–T3 text changes; new Table 1; abstract compressed to 342 words. |
| `manuscript/POD_immune_hub_manuscript_v1.docx` | Regenerated from the v1.2 markdown. |
| `manuscript/cover_letter.md` | Generative-AI paragraph added; DMR-attenuation paragraph re-framed (46–61% band, flat immune proportion, precise p-values). |
| `manuscript/convert_md_to_docx.py` | Fixed: the table splitter now splits on **unescaped** pipes only, so markdown-escaped `\|` inside a cell (e.g. `\|t\| enrichment`) no longer creates phantom columns. Without this fix the regenerated .docx crashed with `IndexError` on the new Table 1. |
| `analysis/adjust/adjust_summary.csv` | Two columns added (`mwu_p_immune_vs_bg`, `median_ratio_immune_over_bg`); existing values unchanged. |
| `README.md`, `preregistration_plan_draft.md` | Terminology aligned to "two-axis hypothesis" / "candidate-modifiable". |

---

## Supplemental analysis performed for T1-3: unified three-lineage DMR adjustment

The panel flagged that the DMR-level adjustment had been run with **granulocyte-subtype** covariates (ΔNeu alone, then ΔNeu+ΔEos+ΔBaso) while the probe-level adjustment used the **three-lineage** set (ΔNeu+Lymph+Mono), i.e. the two levels were never compared under one specification. We therefore re-ran bumphunter on the same 65 paired ΔM profiles with `~ aNeu + lymphoid + aMono`, in a dedicated script (`analysis/layerA_bumphunter_3lineage.R`, `set.seed(12345)`, B = 0) that **re-runs the unadjusted model first as a self-check**.

**Self-check:** unadjusted reproduced **13,357 DMRs / 79 immune-overlapping / 48 immune genes** — identical to the published figures, so the new run is not a silent pipeline change.

**Result under the single three-lineage specification:**

| Model | Covariates | DMRs (L ≥ 3) | Immune-overlapping DMRs | Distinct immune genes | % immune |
|---|---|---|---|---|---|
| unadjusted | none | 13,357 | 79 | 48 | 0.59 |
| adj_dNeu | ΔNeu | 5,153 | 24 | 21 | 0.47 |
| adj_dNeu_Eos_Baso | ΔNeu+ΔEos+ΔBaso | 7,263 | 36 | 29 | 0.50 |
| **adj_3lineage (new)** | **ΔNeu+lymphoid+ΔMono** | **6,628** | **34** | **27** | **0.51** |

**What this settles:**
1. Under one specification the composition-attributable share is **50.4%** — squarely inside the previously quoted 46–61% bracket, so the bracket did not misstate the magnitude; what was wrong was presenting it as a sensitivity envelope. The manuscript now quotes **~50% as the single-specification estimate** and 46–61% only as the spread across specifications.
2. The **immune proportion is flat under all three specifications** (0.59% → 0.47% / 0.50% / 0.51%). This is now a like-for-like result rather than an inference across mismatched covariate sets, and it is the load-bearing evidence for "adjustment removes a global post-surgical methylation shift, not a disproportionate immune component".
3. **21 of the 48** unadjusted immune genes survive the three-lineage adjustment (shared: BCL3, CCL5, CD14, CXCL2, CXCL8, FOS, IKBKE, IL6R, IRAK4, IRF8, MAP3K7, MAPK14, NFKBIA, NFKBIZ, S100A9, SOCS1, STAT1, STAT3, TNFAIP3, TNFRSF1A, TNIP1). Lost: IL6, IL10, IL1B, IL1RN, TNF, IFNG, IRF1, IRF7, TLR4, NLRP3 and 16 others — i.e. the canonical acute-cytokine genes are precisely the composition-driven ones.

**New / amended artefacts:** `analysis/layerA_bumphunter_3lineage.R`, `bumphunter_DMR_adj_3lineage.csv`, `bumphunter_3lineage_summary.csv`, `bumphunter_immune_gene_before_after_3lineage.csv`, `bumphunter_DMR_adj_3lineage_immune_overlap.csv`, `layerA_bumphunter_3lineage.log`; `bumphunter_adjusted_summary.csv` extended with the fourth row; Figure 2 regenerated so panel D shows all four fits and panel C carries the exact p-values.

---

## What was deliberately **not** done

- **No new data layer and no new experiment.** The only analysis added is the unified-covariate DMR re-run above, which reuses the existing inputs.
- **Zenodo not re-published.** DOI 10.5281/zenodo.22896443 is already published and cannot be appended to; a new version would mint a new DOI and needs the author's token. The FILES list is updated and the deposit collector validates at zero missing; the v1.1.0 deposit is scheduled for just before submission.
