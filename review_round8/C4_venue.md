# C4 — Venue & Reporting-Standard Audit (acting as *Journal of Neuroinflammation* editor)

**Manuscript reviewed:** `manuscript/POD_immune_hub_manuscript_v1.md` (single author, Y. Yang; target *Journal of Neuroinflammation* / JNI; ~226 lines; structured abstract; 41 Vancouver references; full Declarations block).
**Audit date:** 2026-09-22.
**Role scope:** journal-metric truth, JNI reporting/compliance, Declarations completeness/honesty, reference reality, cover-letter↔manuscript consistency, format hard-fails. I performed all journal-metric, policy and reference checks **independently via web** (Springer Nature submission guidelines, Clarivate JCR 2025 listing, EuropePMC/PubMed, the live Zenodo API) rather than trusting the manuscript's self-claim. This is a first-submission, fresh-review posture.
**Bottom line up front:** The headline journal-metric claim is *true and independently verified*. The submission is broadly well-aligned to JNI's structured abstract, Declarations and scope. However, two reference defects (one touching the foundational methylation layer) and one missing mandatory BMC policy statement (generative-AI use) must be resolved before acceptance. No fatal/desk-reject item was found, but two T1 items warrant a Revise decision.

---

## 1. Journal metrics claim — VERIFIED TRUE

**【Problem】** Manuscript line 226 asserts: "*Journal of Neuroinflammation* (2025 JCR JIF 11.5, Q1, SCIE, open access; scope includes peripheral immune–nervous system interactions). Journal metrics verified 2026-09-20 via Springer Nature and Clarivate JCR listings."

**【Evidence (web source + what I found)】**
- Clarivate JCR 2025 listing (retrieved via JCR data table, 2026-09-22): *Journal of Neuroinflammation* (J NEUROINFLAMM, BMC, eISSN 1742-2094), **2025 JIF = 11.5**, **Q1** in both categories — NEUROSCIENCES rank 13/330 (JIF percentile 96.2) and IMMUNOLOGY rank 12/183 (JIF percentile 93.7); **SCIE** edition; 5-yr JIF 12.9; % OA gold 100.00.
- Springer Nature "Journal of Neuroinflammation" home/scope (BMC): fully open access; scope text explicitly states it "focuses on interactions of the immune system … with the nervous system … roles of peripheral neuro-immune interactions, T cells, monocytes, complement proteins …" — i.e. **peripheral immune–nervous-system interactions are in scope** (confirmed).
- Cross-check: enterscholar profile also lists "JCR Q1, 影响因子 11.5 (2025)". (One secondary source, cqvip, listed a 2025 IF of 9.995; that figure is not the authoritative JCR value and contradicts the JCR 2025 table — the manuscript's 11.5 is the correct JCR number.)

**【Why it matters】** A false metric claim (e.g., inflated JIF or fabricated Q1/SCIE status) is a serious integrity/venue-mismatch problem. Here the claim is accurate, which supports the author's diligence and the appropriateness of the target.

**【Specific fix】** No change required to the metric claim. Optional polish: cite the verification date and source explicitly as already done; consider adding the JCR category (Neurosciences/Immunology) for precision. Keep the scope sentence — it is correct and helps the editor triage.

---

## 2. Abstract format & limit — COMPLIANT (one soft point)

**【Problem】** JNI requires a structured abstract ≤350 words. The manuscript uses Background / Methods / Results / Conclusions headings and is ~338 words.

**【Evidence (web + file)】**
- JNI (Springer Nature) submission guideline: "*The Abstract should not exceed 350 words and should be structured with a background, main body of the abstract and short conclusion. Please minimize the use of abbreviations and do not cite references in the abstract.*"
- Manuscript abstract (lines 11–21): four labelled subsections; I count ≈ 38 (Background) + 81 (Methods) + 165 (Results) + 54 (Conclusions) ≈ **338 words** — within the 350 limit. No references are cited in the abstract (compliant).
- The four-part (Background/Methods/Results/Conclusions) variant is an accepted BMC structured-abstract style; JNI's literal template is "background / main body / short conclusion," but BMC routinely accepts the finer-grained form.

**【Why it matters】** Abstract over-length or missing structure is a common desk-return reason at BMC. Here it passes.

**【Specific fix】** (T3) The abstract is abbreviation-heavy (POD, APOE ε4, LINCS L1000, GWAS, scRNA-seq, PBMC, EPIC v2, FDR, MWU, DMRs, IFNγ). JNI says "minimize the use of abbreviations." Suggest spelling out APOE ε4 → "APOE ε4 allele" is fine, but consider expanding LINCS L1000, scRNA-seq, EPIC v2, DMRs on first use within the abstract, or accept as-is since the journal guidance is soft.

---

## 3. Declarations — mostly complete; two gaps

### 3a. Ethics approval / IRB-exemption statement — defensible but under-specified
**【Problem】** Declarations (line 163) states the work is secondary analysis of public de-identified data and "exempt from separate Institutional Review Board review." It does not name the original approving committees or reference numbers.

**【Evidence (web + file)】**
- JNI/BMC policy: "*Manuscripts reporting studies involving human participants, human data or human tissue must: include a statement on ethics approval and consent (even where the need for approval was waived) include the name of the ethics committee that approved the study and the committee's reference number if appropriate.*"
- Manuscript line 163: names only that "the original GEO and GWAS studies obtained ethics approval and informed consent" but gives no committee names/IRB numbers for GSE330869 or the Armstrong GWAS. Consent for publication = "Not applicable" (line 165) — acceptable.

**【Why it matters】** Secondary analysis of public, de-identified omics is generally acceptable to BMC without a new IRB, but the guideline explicitly asks for the original committee name/reference "if appropriate." Reviewers/editors may ask for this; an unspecified exemption is a soft compliance gap.

**【Specific fix】** (T2) Either (i) name the source-study ethics committees/IRBs and reference numbers (GSE330869's approving board; the UKB/University of Bristol ethics approval underpinning Armstrong 2026), or (ii) state explicitly "Not applicable — this is a secondary analysis of publicly available, de-identified datasets; no new human-subjects data were generated, and the need for IRB approval was waived per [institution]/[journal] policy," with the legal basis. This removes ambiguity.

### 3b. Data availability — LIVE and matches (verified)
**【Problem】** None material. The Zenodo DOI must be live and title/author/license consistent.

**【Evidence (web: Zenodo API)】**
- `GET https://zenodo.org/api/records/22896443` → status "published", DOI 10.5281/zenodo.22896443, title "*Postoperative Delirium two-axis hypothesis: APOE epsilon4 constitutive susceptibility and peripheral immune state axis — multi-omics integration and in-silico, hypothesis-generating drug-repositioning reproducibility package*", creator "Yang, Yongxin" (ORCID 0009-0004-9698-6552, affiliation matches manuscript), license **MIT**, resource_type "software", 9.1 MB zip, publication_date 2026-09-22. This matches the manuscript Data-availability statement (line 167) and the GitHub repo (`yyx-4113/pod-immune-hub-repositioning`).
- Repo `CITATION.cff` also lists MIT license and the same title — consistent with the live record.

**【Why it matters】** A dead/mismatched data-DOI is a hard reproducibility fail at JNI, which "strongly encourages" public data availability.

**【Specific fix】** (T3) Two minor consistencies: (i) manuscript footer line 226 says "Manuscript version 1.1 (Round-7 revision, 2026-09-22)" while the Zenodo record version is **1.0.0** — align the Zenodo-version tag to the manuscript version, or vice-versa. (ii) The Zenodo record's `version`/`related_identifiers` fields contain the branch name in Chinese ("术后谵妄免疫枢纽"); cosmetic, but consider a clean semantic version. Neither blocks acceptance.

### 3c. Funding — "None" (line 169) — acceptable.
### 3d. Author contributions & degree claim — HONEST (strength)
**【Problem】** None. The single author claims all four ICMJE criteria and explicitly states "The author holds a Bachelor of Medicine (B.M.) degree; no postgraduate degree (MD/PhD/MS) is claimed" (line 171).

**【Evidence (file)】** Manuscript line 171. This is transparent and does **not** inflate credentials — the opposite of a problem.

**【Why it matters】** Credential misrepresentation is a common integrity flag; here the author proactively disclaims any graduate title, which is commendable and should be trusted at face value.

**【Specific fix】** None on honesty. (Soft note, not a defect: a single-author computational multi-omics + GWAS-permutation + LINCS study with no declared statistical/genetics co-author or supervisor is unusual; reviewers may probe whether methodological choices were independently sanity-checked. Not a venue-format issue — flag only as a substance point for the other reviewers.)

### 3e. Competing interests — "The author declares no competing interests" (line 173) — acceptable and complete.

---

## 4. Generative-AI disclosure — MISSING (mandatory BMC policy)

**【Problem】** The manuscript contains **no statement** about use (or non-use) of generative-AI / LLM tools, in Methods or anywhere. The cover letter is also silent.

**【Evidence (web)】**
- BMC / Springer Nature policy (stated verbatim in JNI submission guidelines and the SN AI policy): "*Large Language Models (LLMs), such as ChatGPT, do not currently satisfy our authorship criteria … Use of an LLM should be properly documented in the Methods section (and if a Methods section is not available, in a suitable alternative part) of the manuscript.*" (Note: pure 'AI-assisted copy editing' for readability/grammar need not be declared; substantive generation must be.)
- Manuscript Methods (lines 35–84) has no AI-use sentence. Given the highly polished, idiomatically native English prose of a sole B.M.-credentialed author across a 226-line multi-omics manuscript, the likelihood of at least assistive AI use is high, and the policy requires documentation either way.

**【Why it matters】** This is a **required declaration**, not optional. Omission is a policy breach and is routinely a "return-to-author" item at BMC; it can also undermine the accountability prong of the ICMJE statement the author asserts in §3d.

**【Specific fix】** (T1) Add a short statement to the Methods section (and ideally the cover letter), e.g.: "*No generative-AI tools were used in the conception, analysis, or writing of this manuscript; writing was performed by the author without AI assistance.*" — or, if AI was used only for language polishing, state that explicitly ("AI-assisted copy editing for grammar/style only; no content was AI-generated"). Choose the truthful version. This must be in the submitted manuscript, not just the letter.

---

## 5. References — 41 Vancouver; two defects found, several weak spots

I spot-checked 10+ references against EuropePMC/PubMed (PMID/DOI/title/author/year/journal).

**VERIFIED CORRECT (real, correctly attributed):**
- [1] Meza Monge 2025 *J Clin Med* 14(23):8418 — PMID 41375722 ✓
- [2] Mosharaf 2025 *Sci Rep* 15:7830 — PMID 40050293 ✓
- [4] Wu 2025 *Biomedicines* 13(12):2962 (GSE252572) — PMID 41462974 ✓
- [7] Armstrong 2026 *PLoS Med* 23(3):e1004963 — PMID 41770756, DOI 10.1371/journal.pmed.1004963 ✓ (exact match; foundational GWAS)
- [17] Nishizawa 2024 *Transl Psychiatry* 14:275 — PMID 38965205 ✓
- [26] MoDUS, Page 2017 *Lancet Respir Med* 5:727–737 — PMID 28734823 ✓
- [32] Huang 2025 *Hum Genomics* 19:110 — PMID 41029350 ✓
- [35] Cheng 2022 *J Neuroinflammation* 19:297 — PMID 36503642 ✓ (also confirms JNI scope fit)
- [3] Hsiao 2023 *PLoS One* 18(8):e0282214 (APOE ε4 POCD meta-analysis, OR 1.28) — PMID 36827351 ✓ (manuscript did not print the PMID, but the paper is real and the attribution is correct)
- [22] Mardani 2013 *J Res Med Sci* 18(2):137–143 (dexamethasone RCT) — PMID 23914217 ✓

### 5a. **[6] Seki 2026 — UNVERIFIABLE; foundational methylation layer depends on it (T1)**
**【Problem】** Reference [6] is cited as the source publication of GSE330869 and as the origin of the methylation anchors (cg01534316→TMIGD3/ADORA3, FDR 0.021) used throughout Layer A. It could not be located.

**【Evidence (web)】**
- EuropePMC query `EXT_ID:42143058` → **hitCount 0 (record does not exist)**.
- EuropePMC title query `TITLE:"Epigenetic biomarker of delirium risk"` → **hitCount 0**.
- Web search confirms the *Shinozaki/Seki/Nishizawa group* has genuine 2026 delirium-epigenetics output (e.g., Ishii, Shibata, Nishitani, **Seki**, Yamanishi, Nishiguchi, Shimura, Aoyama, Gorantla, Phuong, Santiago, Shinozaki, "*Neuroinflammation as molecular landscape of post-operative delirium revealed by live human brain multi-omics profiling*," *Molecular Psychiatry* / Nature 2026, DOI 10.1038/s41380-026-03804-z; and a 2026 *Translational Psychiatry* hip-fracture delirium EWAS described in a press release by "Seki, Nishitani, Nishizawa …"). But **none of these resolves to PMID 42143058 or the exact title "Epigenetic biomarker of delirium risk" in PubMed/EuropePMC.**

**【Why it matters】** Layer A (the recomputed-in-house blood-methylation layer) and the "Seki 2026 anchor" narrative (lines 43, 48, 100, 122) are load-bearing. If [6]'s PMID/title is a typo, the citation is merely wrong; if the record is not indexed, the authors must supply a verifiable pointer; if the paper is misattributed to GSE330869, the provenance of the anchor claims is compromised. A foundational reference that fails web verification is a serious reporting-standard defect.

**【Specific fix】** (T1) The authors MUST (i) confirm the exact source publication of GSE330869 (query GEO directly for the "Citation" field), (ii) correct the PMID if 42143058 is a typo, or (iii) if the 2026 *Transl Psychiatry* paper is not yet PubMed-indexed, cite it with its DOI and a preprint/available link, and clearly distinguish it from the group's other delirium-epigenetics papers (Nishizawa 2024 *Transl Psychiatry* PMID 38965205; Wahba 2022 *J Psychiatr Res* PMID 36270064; Nishizawa 2024 *J Psychiatr Res* PMID 39043004). Do not let the manuscript rest a key anchor attribution on an unresolvable reference.

### 5b. **[41] Evered 2018 — WRONG PMID (T2)**
**【Problem】** Reference [41] lists DOI 10.1097/ALN.0000000000002334 but **PMID 30179813**, which is incorrect.

**【Evidence (web)】**
- EuropePMC `EXT_ID:30179813` → returns an unrelated paper: "*Spatial distributions and transport implications of short- and medium-chain chlorinated paraffins in soils and sediments …*," *Sci Total Environ* 2019. That PMID does **not** belong to Evered.
- The correct record: title-search "Recommendations for the Nomenclature of Cognitive Change Associated with Anaesthesia and Surgery" returns the Anesthesiology 2018 paper at **PMID 30325806**, DOI 10.1097/aln.0000000000002334 (volume 129(5):872–879) — i.e., the manuscript's **DOI is right, the PMID is wrong**. Author list ("Evered L, Silbert B, Knopman DS, Scott DA, … Nomenclature Consensus Working Group") matches the manuscript. The content attribution (POD vs POCD/PND nomenclature) is correct.

**【Why it matters】** A misattributed PMID is a factual error that signals the reference list was not validated against PubMed; it also breaks DOI↔PMID linking in the typeset article.

**【Specific fix】** (T2) Change [41] PMID from 30179813 to **30325806**. (Note: the same DOI appears correctly, so only the numeric PMID needs correction.)

### 5c. Missing PMIDs / unverifiable 2026 refs (T2, completeness)
**【Problem】** Several references print no PMID, and a cluster of 2026 references could not be independently verified because no PMID/DOI is given: [15] (Li 2019, no PMID), [26] (MoDUS, no PMID printed — but verified above via DOI), [33] Lee 2026 *Anesth Pain Med*, [39] Paterno 2026 *J Clin Med*, [40] Qin 2026 *Front Neurol*.

**【Evidence (file)】** Reference list lines 196–222. The 2026 entries [33][39][40] carry only author/year/journal/DOI, no PMID; [15] carries only DOI.

**【Why it matters】** BMC encourages (not mandates) PMID/PMCID, but unverifiable forward-dated references are a known vector for citation error. Reviewers cannot check them.

**【Specific fix】** (T2) Add PMIDs where available (e.g., [15] Li 2019 BMC Anesthesiol has a PMID) and supply DOIs/PMID for [33][39][40]; if any 2026 paper is in-press/unindexed, mark it "in press" with a verifiable DOI rather than a bare citation.

### 5d. Format consistency (T3)
Vancouver style is used throughout (numbered, "Authors. Title. Journal. Year;vol(issue):pages. doi."). Mixed presence of PMID/PMCID is acceptable to BMC. One cosmetic point: [7] prints the GWAS summary-stats DOI inline with the article DOI — fine. No systematic formatting break found.

---

## 6. Cover letter vs manuscript — CONSISTENT (no overclaim)

**【Problem】** Verify the cover letter does not overclaim "two dissociable axes," is consistent on "model" vs "hypothesis," and represents the negative/failed-positive-control results faithfully.

**【Evidence (file: `manuscript/cover_letter.md`)】**
- Dissociation: cover letter line 9 explicitly states "*the independence (dissociation) of the two axes is a framework we could not test (no layer carried APOE ε4 genotypes) and is supported only indirectly by the transcriptomic and pharmacologic layers.*" This matches the manuscript's stance (lines 19, 118, 133, 138, 153: "working two-axis hypothesis," "dissociation is a framework, not a tested result," "genetic data are agnostic on dissociation"). **No overclaim of "two dissociable axes" as established** — good.
- "Model" vs "hypothesis": the cover letter and manuscript both use "two-axis hypothesis" / "working two-axis hypothesis." However, the repo `README.md` title and `preregistration_plan_draft.md` title use "two-axis **model**." This is a minor cross-document terminology drift (model vs hypothesis), not an internal manuscript/cover-letter contradiction.
- Negative/failed controls faithfully represented: cover letter line 13 lists "*a failed virtual-knockout class test on the whole-transcriptome signature (z = −2.47), no immune enrichment in the GWAS once chromosome 19 is excluded (p = 0.277), and no cell-type proportion shift in PBMC (minimum p = 0.229)*" — all match the manuscript (lines 112, 116, 96). The cover letter's self-falsification narrative (composition adjustment removing 61% of DMRs, p 2.5×10⁻⁴ → 0.72) matches manuscript lines 106, 141.
- Data availability / repository / no-competing-interests / no-funding in the letter (line 15) match the manuscript Declarations.

**【Why it matters】** A cover letter that overclaims (e.g., asserting dissociation as proven) would mislead the editor and breach the manuscript's own honesty posture. Here the letter is conservative and consistent.

**【Specific fix】** (T3) Align terminology across repo docs: use "hypothesis" (not "model") in README.md and preregistration title to match the manuscript, or explicitly note the two are used interchangeably. Add the generative-AI statement to the cover letter too (see §4).

---

## 7. Format hard-fails & figure standards

**【Problem】** Confirm figures meet JNI/BMC resolution/legend/colourblind standards; check for missing legends or embedded-text redundancy.

**【Evidence (web + file)】**
- JNI/BMC figure policy: figures should be supplied at **≥300 dpi** (JNI prefers high resolution; 600 dpi is compliant and above threshold). The manuscript cites `layerA_figure.png (600 dpi)` + `.pdf` and `layerA_figure2.png`/`.pdf` (lines 104, 106) — resolution standard met. Colour-encoding is claimed "colorblind-safe" (line 104) — acceptable but the actual files were not available in the review package for me to verify perceptual safety; this is a submission-package content check, not a manuscript-text fail.
- Figure legends: Figure 1 (line 104) and Figure 2 (lines 104–106) have full, informative legends with panel labels (A–D) and statistical annotations — compliant with BMC "captioned, cited in order" requirement.
- No redundant embedded-text duplication of figure content was found in the manuscript body beyond the legitimate legend text.

**【Why it matters】** Low-resolution or legend-less figures are a common BMC desk-return; here the text declarations are compliant, but the actual figure files must accompany the submission.

**【Specific fix】** (T2) Ensure the two figure files (PNG≥300 dpi + PDF) are actually uploaded in the Editorial Manager submission bundle (they are referenced but not present in the reviewed folder). Confirm colourblind-safe palettes with a checker (e.g., Coblis) and state the palette name in the legend. (T3) Submit as .docx/.tex per JNI's Editorial Manager route — the `.md` is a source artifact, not a submission format.

---

## 8. Other reporting-checklist observations

- **Reporting guidelines:** The manuscript correctly notes ARRIVE/STROBE/STROBE-MR do not apply and references MRE-MoSTre + Subramanian 2017 CMap convention (line 84). Reasonable for a public-data multi-omics study. No mandatory checklist (e.g., PRISMA) is triggered.
- **Preregistration honesty:** The manuscript (line 31) and preregistration both avoid the term "pre-registered" in its formal sense and describe a version-controlled a-priori plan — this is honest and aligns with ICMJE/journal expectations. Good.
- **Scope fit:** Confirmed in-scope (peripheral immune–nervous-system interactions; JNI publishes delirium EWAS/omics, e.g., Hirsch 2016 *J Neuroinflammation*; Cheng 2022 *J Neuroinflammation*). Edge case: the paper is purely computational with no primary neurobiology/wet-lab validation; JNI editors may weigh mechanistic depth, but it is within scope and not a format fail.

---

## § Stands up (strengths)
1. **Journal-metric claim is accurate and independently verified** (JCR 2025 JIF 11.5, Q1, SCIE, fully OA; scope statement correct) — line 226 is trustworthy.
2. **Data availability is real and matching**: Zenodo DOI 10.5281/zenodo.22896443 is live, title/author/affiliation/OCRID consistent, MIT-licensed, with a deposited code+results zip; GitHub mirror named. Reproducibility posture is strong.
3. **Honest, well-structured Declarations** including a transparent single-author B.M. credential statement (no false graduate-degree claim), explicit IRB-exemption stance, "None" funding, and no competing interests.
4. **Cover letter is conservative and consistent** with the manuscript — it does not overclaim dissociation, accurately reports failed positive controls (KO z = −2.47; chr19-excluded p = 0.277; PBMC min p = 0.229), and matches the data/repo statements.
5. **Abstract is compliant** (structured, ≈338 words < 350 limit, no in-abstract citations) and the negative-result/hypothesis-generating framing is commendably transparent throughout.

---

## § Questions for the authors
1. Please confirm the exact, PubMed-verifiable source publication of GSE330869. Is [6] "Seki 2026, *Transl Psychiatry*, PMID 42143058" correct, or is the PMID a typo / the paper not yet indexed? What is the GEO-listed "Citation" for GSE330869?
2. Did you use any generative-AI / LLM tool in writing, analysis, or figure preparation? Per BMC policy this must be documented in Methods (copy-editing-only use should also be stated). Please add the statement.
3. For the ethics statement: can you name the original IRBs/reference numbers for GSE330869 and the Armstrong GWAS (or state the explicit waiver basis)?
4. Please correct [41] PMID to 30325806 and add missing PMIDs/DOIs to [15], [33], [39], [40].
5. Will the two figure files (≥300 dpi PNG + PDF, colourblind-safe) be uploaded with the submission? Which palette was used?

---

## § What I actually checked
**Files read (allowed sources only):** `manuscript/POD_immune_hub_manuscript_v1.md`, `manuscript/cover_letter.md`, `README.md`, `CITATION.cff`, `preregistration_plan_draft.md`, `zenodo_metadata.json`. (Per the brief, I did **not** open REVIEW_*/RESPONSE_*/REVISION_*, review_round8 sibling files, PROGRESS.md, G1–G4 gates, submission_notes, author_verification_statement, DEPOSIT_FILELIST, or zenodo_deposit_v1.0.0.)
**Web checks performed:**
- JNI metrics/policy: Springer Nature JNI submission guidelines (abstract ≤350 words structured; LLM-use must be documented in Methods; ethics-committee name/reference required; Declarations subheadings; scope includes peripheral neuro-immune interactions); Clarivate JCR 2025 table (JIF 11.5, Q1, SCIE, 100% OA).
- Zenodo live API: `records/22896443` (published, MIT, Yang Y, title matches, 9.1 MB zip).
- EuropePMC/PubMed verification of refs [1],[2],[4],[6],[7],[17],[26],[32],[35],[41],[3] (and title-search for [41] and [6]): [1][2][4][7][17][26][32][35][3][22] verified correct; **[6] not found (hitCount 0 by PMID and by title)**; **[41] PMID 30179813 resolves to an unrelated paper (correct PMID = 30325806)**.
- Web search for GSE330869/Seki/Shinozaki 2026 delirium epigenetics to contextualise [6].

**Discrepancies found:**
- [6] Seki 2026 unresolvable in PubMed/EuropePMC (T1).
- [41] Evered 2018 wrong PMID (30179813 → 30325806) (T2).
- No generative-AI disclosure in manuscript or cover letter, contrary to mandatory BMC policy (T1).
- Ethics statement lacks original committee names/reference numbers (T2).
- Manuscript version "1.1" vs Zenodo version "1.0.0" (T3); repo `zenodo_metadata.json` declares cc-by-4.0 while live record + CITATION.cff are MIT (T3, repo-only).
- Several 2026 refs lack PMIDs (T2); abstract abbreviation density (T3); "model" vs "hypothesis" drift between repo docs (T3).

---

## Overall verdict
**Recommendation: Revise (moderate).** No desk-reject or fatal integrity defect was identified, and the central venue claim (JNI metrics/scope) is verified true. Acceptance is contingent on resolving two T1 items — (a) verify/correct the GSE330869/Seki 2026 citation, and (b) add the mandatory BMC generative-AI use statement — plus the straightforward T2 corrections (Evered PMID, ethics-committee naming, figure-file upload, missing PMIDs). The manuscript's transparency posture, reproducibility packaging, and conservative hypothesis framing are genuinely strong and align with JNI's editorial values.

**Tier classification**
- **T0 (fatal / desk-reject):** none.
- **T1 (major; must resolve before acceptance):** (1) [6] Seki 2026 / GSE330869 citation unverifiable — confirm exact source and fix; (2) missing mandatory generative-AI disclosure in Methods (+ cover letter).
- **T2 (moderate; correct in revision):** (3) [41] Evered 2018 PMID → 30325806; (4) ethics statement — name original IRBs/references or state explicit waiver basis; (5) supply/upload figures at ≥300 dpi with colourblind-safe palette named; (6) add PMIDs/DOIs to [15],[33],[39],[40].
- **T3 (minor/cosmetic):** (7) manuscript v1.1 vs Zenodo v1.0.0 alignment; (8) repo license metadata (zenodo_metadata.json cc-by-4.0 vs live MIT); (9) abstract abbreviation minimisation; (10) "model" vs "hypothesis" terminology consistency across README/preregistration; (11) submit as .docx/.tex per Editorial Manager.
