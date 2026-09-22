# C2 — Design & Statistics Review (independent, first-submission posture)

**Manuscript:** *A two-axis hypothesis for postoperative delirium: APOE ε4 constitutive susceptibility and a peripheral immune state axis supported by multi-omics integration, with an in-silico, hypothesis-generating drug-repositioning screen* (target *Journal of Neuroinflammation*).
**Reviewer code:** C2 (design-level bias, confounding, sparse-cell / permutation validity, gene-set enrichment validity, and hedge-vs-uncertainty matching).
**Declaration:** I read this manuscript fresh and recomputed every quantitative claim below from the committed CSV/JSON artefacts under `analysis/results/` and `analysis/adjust/` using managed Python 3.13. I did **not** consult any prior review, response, or revision material. Numbers I recomputed are marked **[recomputed]**.

---

## Overall verdict (short)

This is an unusually honest multi-omics integration: the null results, failed positive controls, and the "untested dissociation" are carried as load-bearing arguments, which I strongly welcome. The design is appropriate for a hypothesis-generating study. My concerns are mostly about **how precisely the uncertainty is quantified and hedged**, not about fabrication or fatal flaws. The single concrete arithmetic error I found is the Layer-C "minimum p = 1/3" (it is 1/6). The most substantive design issue is the "46–61% composition-attributable" band, which mixes two different covariate specifications and is not a sensitivity/confidence range in the usual sense. I recommend **major-revision (T1)** on the composition-attribution framing and the Set-1/Set-2 comparability hedge; the rest are **T2/T3**.

---

## CHECK 1 — EWAS design: paired-within-POD, no controls; "6 DMPs, 0/49 immune overlap, Fisher p=1"

**【Problem】** The manuscript reports only 6 DMPs (FDR<0.05, |Δβ|>0.05), 0 of 49 immune-gene overlaps (Fisher p=1), and interprets this as "individual immune loci are not genome-wide significant in this paired-within-POD design." This interpretation is **correct but under-powered to the point that the statement is almost tautological**, and the manuscript understates how little this design can ever detect.

**【Evidence】**
- `layerA_topDMP.csv`: exactly 6 rows; all satisfy `adj.P.Val<0.05` **and** `|deltaBeta|>0.05` (range −0.0648 to +0.0642) **[recomputed]**.
- `layerA_summary.txt`: `n_DMP_fdr05_dB05 = 6`, `n_DMP_fdr05_dB10 = 0`, `immune_gene_overlap_n = 0`, `immune_gene_overlap_total = 49`, `fisher_p = 1`.
- Fisher exact, gene-level 2×2 (6 DMPs vs 312,514 background probes; 49 immune genes; 0 overlap): two-sided p = **1.0**, one-sided (greater) p = **1.0** **[recomputed]**.
- Power: with n=65 paired, limma empirical-Bayes, the detectable |Δβ| at FDR 5% for a moderate effect is on the order of 0.05–0.08 depending on residual variance; the fact that only 11 probes reach FDR<0.05 array-wide (stated in Limitations) means the array is operating near its detection floor.

**【Why it matters】** "Individual immune loci are not genome-wide significant" is true, but a paired-within-POD design with no control group and only 65 subjects has effectively **zero power to detect immune-specific DMPs** against background. Fisher p=1 here is not evidence of *absence* of immune methylation signal — it is the expected outcome of a severely under-powered, control-free design where only 6 probes clear any bar at all. Presenting "0/49 overlap, Fisher p=1" as if it were a test that *could have* found an enrichment risks being read as negative evidence against the immune axis, when it is really a non-test.

**【Specific fix】**
1. Reframe the 0/49 Fisher p=1 sentence explicitly as "no power to detect immune-specific DMPs; the 6 DMPs are background-level noise and their non-immune identity is uninformative, not evidence against the immune hypothesis."
2. Add a one-line power note: "At n=65 paired, EPIC v2, this design detects only large |Δβ| (>0.05) effects; immune loci of the magnitude seen in Seki 2026's Time×Group model (which included controls) would not survive here." This is already partially in Limitations — promote it to the Results sentence.

---

## CHECK 2 — Composition adjustment: probe-level vs DMR-level use *different* covariate sets; the "46–61% composition-attributable" band is not a sensitivity range

**【Problem】** The headline claim "composition removed 46–61% of DMRs, so the methylation immune signal is substantially composition-associated" conflates two different analyses built on **different covariate definitions**, and presents a two-point spread as if it were a robust attribution band.

**【Evidence】**
- `bumphunter_adjusted_summary.csv`: unadjusted DMR = 13,357; adj_dNeu = 5,153; adj_dNeu_Eos_Baso = 7,263. Immune-overlapping DMRs: 79 → 24 → 36; immune proportion 0.59% → 0.47% → 0.50% **[recomputed]**.
- Reduction vs unadjusted **[recomputed]**: adj_dNeu = **61.4%**; adj_dNeu_Eos_Baso = **45.6%**. This is the source of "46–61%."
- Covariate sets are *not* matched across levels:
  - **Probe level** (`adjust_summary.csv` / `analysis/layerA_set1_signed_t_adjusted.py`): uses `MA_adj_dNeu` (neutrophil delta only) and `MB_adj_3lineage` (`dNeu + dLymph + dMono`).
  - **DMR level** (`analysis/layerA_bumphunter_adjusted.R`): uses `~ aNeu` (neutrophil only) and `~ aNeu + aEos + aBaso` (granulocyte-subtype). The 3-lineage (Neu+Lymph+Mono) set used at the probe level is **never** used at the DMR level.
- The "band" endpoints come from two *different* covariate choices at the DMR level: dNeu-only (61%) and dNeu+Eos+Baso (46%). Eosinophil/basophil deltas barely move (dEos mean = −0.025 pp, sd 1.55 pp; `adjust_report.txt`), so adding them as covariates injects near-zero-variance regressors that absorb degrees of freedom and *reduce* the apparent removal (61% → 46%). The spread is therefore an artefact of covariate selection, not a confidence interval on the composition-attributable fraction.

**【Why it matters】** Two issues:
1. **The band is not a sensitivity envelope.** A reader will read "46–61%" as "we are 95%-confident the composition share is in [46,61]%." It is not; it is two numbers from two covariate specifications. The true composition-attributable share could be higher (if e.g. a richer lineage model were used) or lower (regression dilution — see below).
2. **The immune-proportion is essentially flat (0.59% → 0.47–0.50%).** This means composition adjustment removes immune and non-immune DMRs *proportionally*. If composition were specifically driving the *immune* methylation signal, you would expect the immune proportion to **drop** preferentially. It does not. So the DMR-level evidence actually shows the *global* post-surgical methylation shift is composition-associated, with immune DMRs just a constant ~0.5% slice — the *immune-specific* composition claim rests almost entirely on the **probe-level |t| enrichment test** (unadjusted p=2.5e-4 → dNeu-adjusted p=0.72, `adjust_summary.csv`), which is the cleaner, more specific demonstration. The abstract's causal linkage ("removed 46–61% of DMRs, *so the methylation immune signal is substantially composition-associated*") slightly over-reaches from the DMR count alone.
3. **Regression dilution.** Composition fractions are estimated with measurement error (the manuscript acknowledges this in Limitations). Adjusting for a noisy covariate removes both composition signal *and* some genuine within-cell signal, so "composition-attributable" is an **upper bound**, not a point estimate. This compounds the band imprecision.

**【Specific fix】**
1. Re-state the DMR result as: "adjusting for the neutrophil shift alone removed 61% of DMRs; a broader granulocyte-subtype covariate removed 46% — i.e. at least ~46% and plausibly up to ~61% of the *overall* DMR burden is composition-explicable, and this is an upper bound given deconvolution measurement error." Drop the "band = precision" implication.
2. Explicitly note the flat immune proportion (0.59%→0.47–0.50%) and state that the *immune-specific* composition conclusion is carried by the probe-level |t| enrichment test (p=2.5e-4 → 0.72), not by the DMR proportion.
3. Use **one consistent covariate set** (recommend the 3-lineage Neu+Lymph+Mono, since it is the most interpretable lineage partition) at *both* probe and DMR levels, or present the two covariate choices as explicitly separate scenarios rather than as a single "46–61%" figure.

---

## CHECK 3 — Set-1 vs Set-2 are effectively incomparable tests; is the "not-yet-composition-validated" hedge sufficient?

**【Problem】** The directional methylation signal is carried by Set-1 (p=0.041, 49 genes / 737 probes, limma moderated-t signed-t, **unadjusted**). The only composition-adjusted analogue is Set-2 (p=0.296 / 0.098 / 0.141, 83 genes / 1,195 probes, OLS-on-delta t). These differ in **both** gene set **and** t-statistic definition, so Set-2 cannot validate or refute Set-1. The manuscript's hedge ("not-yet-composition-validated") is present but is easy to under-weight.

**【Evidence】**
- `layerA_summary.txt`: `directional_mwu_p = 0.04127092`, `n_immune_probes = 737` (Set-1, 49 genes, limma moderated t).
- `set1_signed_t_adjusted.csv`: Set-2 (1,195 probes, 83 genes), OLS-on-delta t: M0 p=0.296, dNeu p=0.098, 3-lineage p=0.141 **[recomputed against file]**.
- `analysis/layerA_set1_signed_t_adjusted.py` line 4–8: the script itself states it only computes the Set-2 OLS signed-t; it explicitly does **not** recompute the Set-1 limma moderated-t after adjustment.

**【Why it matters】** A reader who remembers "immune methylation directional signal p=0.041 (significant)" from the Abstract/Results has no adjusted counterpart for that exact test. The adjusted null (0.296/0.098/0.141) is a *different* test, so it neither confirms nor disproves p=0.041. The risk is that p=0.041 is carried forward as a "significant directional immune-methylation signal" while its only adjusted sibling is null — and the two are never placed side by side. This is the central statistical ambiguity of Layer A.

**【Specific fix】**
1. At the **first mention** of p=0.041 (Abstract and Results line 102), add in-line: "Set-1, unadjusted, limma moderated-t, 49-gene/737-probe set; no composition-adjusted analogue of this exact test exists."
2. Present a compact table of *all four* Layer-A methylation p-values with columns {test object, gene set, t-definition, adjusted?, p}: (a) Set-1 signed-t MWU 0.041 unadj; (b) Set-2 OLS signed-t 0.296/0.098/0.141 adj; (c) |t| enrichment MWU 2.5e-4→0.72/0.55 adj. This makes the incomparability visible rather than prose-hidden.
3. In the Abstract, soften "a directional probe-set test was hypothesis-generating (Set-1 … p=0.041, not composition-validated)" to "a directional probe-set test was hypothesis-generating and unvalidated (Set-1 … p=0.041; its only adjusted analogue, a different Set-2 OLS test, was null: 0.296/0.098/0.141)."

---

## CHECK 4 — Two distinct probe-level tests (signed-t MWU 0.041 vs |t| enrichment MWU 2.5e-4): are they conflated?

**【Problem】** The directional signed-t MWU (p=0.041, Set-1) and the |t| enrichment MWU (p=2.5e-4, Set-2) are different tests on different objects. I checked whether any reader could mistake 0.041 for the composition-adjusted number.

**【Evidence】**
- `manuscript` line 106 (Figure 2C): "the immune-versus-background |t| enrichment observed unadjusted (median 0.968 vs 0.903; two-sided MWU p = 2.5 × 10⁻⁴) is abolished after regressing … (0.650 vs 0.656; p = 0.72) … This |t|-based enrichment test is distinct from the signed-t directional test reported in the Results text (Set 1, p = 0.041)." — the manuscript **explicitly separates** them.
- `adjust_summary.csv`: |t| enrichment unadjusted MWU p = 0.0002463; dNeu p=0.7244; 3-lineage p=0.55 **[recomputed]**.
- The composition-adjusted numbers (0.296/0.098/0.141 for Set-2; 0.72/0.55 for |t|) are all clearly *different* from 0.041, and 0.041 is never presented as adjusted.

**【Why it matters】** On this specific point the manuscript is **well handled** — I found no conflation. The only residual risk is that a scanning reader sees "0.041" and "2.5e-4" both described as "immune methylation MWU p" near each other and assumes they are the same quantity. They are not (0.041 = directional/signed; 2.5e-4 = magnitude/|t|). The distinction is made in text but not in a single comparative table.

**【Specific fix】** T3 — add a one-line footnote or table caption at first appearance of each p-value labelling it "signed-t (directional)" vs "|t| (magnitude)" so the two are never ambiguous. No substantive change needed.

---

## CHECK 5 — bumphunter B=0 (no permutation): reporting "13,357 DMRs" as descriptive

**【Problem】** bumphunter was run with `B = 0` (no bootstrap permutation), so DMR counts carry no empirical FDR. Reporting "13,357 DMRs" prominently (Results line 102, Figure 2) is acceptable *as descriptive* but invites over-interpretation, especially because the "46–61% composition-attributable" reduction is computed *from these descriptive counts* (see Check 2).

**【Evidence】**
- `analysis/layerA_bumphunter_adjusted.R` line 92: `bumphunter(Dmat, design = design, chr = chr, pos = pos, coef = 1, cutoff = 0.05, B = 0)` — `B=0` confirmed **[recomputed from script]**.
- `bumphunter_adjusted_summary.csv`: 13,357 DMRs (L≥3, |ΔM|≥0.05) — committed as descriptive.
- The manuscript does label B=0 "descriptive, no empirical FDR" in Methods (line 50) and Limitations (line 140). Good.

**【Why it matters】** Two compounding risks: (i) a reader skimming the Abstract/Results may treat "13,357 DMRs" as a substantive finding rather than a region-call with no multiplicity control; (ii) because the composition-attribution percentage is derived from these uncontrolled counts, the "46–61%" magnitude inherits the B=0 lack of FDR — it is a *descriptive* reduction, not an inferential one. The manuscript hedges the DMR count but not the *derived* percentage.

**【Specific fix】** T2 — append "(descriptive, B=0)" to *every* occurrence of "13,357 DMRs" and to the composition-reduction percentage, not only at first mention. Consider a single sentence: "All DMR counts below are descriptive (bumphunter B=0, no permutation FDR); the 46–61% reduction is therefore a descriptive proportion, not an inferential attribution."

---

## CHECK 6 — Layer B gene-set MWU (p=0.0079) at n=8; the "one-sample t p=0.18 as a lower bound on directional signal"

**【Problem】** The immune gene-set MWU on moderated t (p=0.0079) is a reasonable competitive test, but the companion statement that the one-sample t (p=0.18) is a "lower bound on the directional signal" is statistically loose and potentially misleading.

**【Evidence】**
- `layerB_summary.json`: `immune_set_test`: `n_immune_genes_on_array=83`, `mean_t_immune=0.2533`, `mean_t_background=−0.1705`, `mannwhitney_U=2169503`, `mannwhitney_p=0.0078746`, `onesample_t_p=0.18107` **[recomputed against file]**.
- Manuscript line 92: "an unweighted one-sample test for mean non-zero shift was non-significant (p = 0.18, a conservative lower bound on the directional signal given heterogeneous per-gene variances at n = 4 v 4)."

**【Why it matters】**
1. **"Lower bound on the directional signal" is not a standard statistical statement.** A p-value is not a bound on an effect size. What the authors likely mean is: "the per-gene mean shift is itself non-significant at this n, so the detectable signal is broad-but-shallow, and the MWU (p=0.0079) is the stronger evidence for a *distributional* shift." But p=0.18 is *larger* (less significant) than p=0.0079, so calling it a "lower bound" inverts the intuition — a reader could wrongly infer "the true directional evidence is at least as strong as p=0.18," which is backwards.
2. **Validity of the MWU at n=8.** The moderated t already borrows strength across 44,824 probes via empirical-Bayes, so the 83 immune-gene t's are well-behaved and the competitive MWU vs 44,741 background genes is a legitimate gene-set enrichment test. The n=8 (4v4) affects per-gene *precision* but the moderated t stabilises it. This part is fine.
3. **The one-sample t is over 83 genes, not n=8.** Saying "heterogeneous per-gene variances at n=4 v 4" explains why individual genes are noisy, but the one-sample t aggregates 83 gene-level t's; its non-significance (p=0.18) reflects that the *mean* immune t (0.253) is small relative to its spread across 83 heterogeneous genes, not a pure n=8 limitation. The "lower bound" framing misattributes the non-significance.

**【Specific fix】** T2 — rephrase: "The immune-gene t-distribution is shifted toward upregulation versus background (MWU p=0.0079), but the *mean* per-gene shift is non-significant (one-sample t p=0.18) because individual immune genes show small, heterogeneous effects at n=8; the signal is therefore broad and directional rather than driven by a few large per-gene changes." Delete "lower bound on the directional signal."

---

## CHECK 7 — GWAS gene-set permutation: stratified on SNP count, density not controlled; off-chr19 null (p=0.277)

**【Problem】** The off-chromosome-19 immune-set null (p=0.277) is the genetic cornerstone of "the immune axis is genetically consistent with, but not established by, GWAS." The permutation is stratified on SNP count (gene-length control), but SNP *density* (SNP/kb, mappability/GC) is not separately controlled. I assessed whether residual density confounding could overturn the null.

**【Evidence】**
- `layerF_real_gwas_summary.json`: `n_variants = 12,353,257`, `genomic_inflation_lambda = 1.0133`, `n_genome_wide_significant = 37`, all on chr19 (`chromosome_distribution_of_hits: {"19":37}`); lead rs429358 BETA=−0.621, OR per ε4 = exp(0.621)=**1.861**; `n=134,310` per SNP **[recomputed against file]**.
- Immune-set test, all autosomes: `perm_p_one_sided = 0.0812`, `z_vs_null = 1.263` (manuscript "p=0.081, z=1.26" ✓).
- Immune-set test, chr19 excluded: `perm_p_one_sided = 0.2769` (manuscript "p=0.277" ✓); 76 of 80 immune genes remain (the 4 dropped are chr19: RELB, IRF3, BCL3, TICAM1).
- `analysis/layerF_real_gwas.py` `stratified_set_test`: strata by `pd.qcut(n_snps, 10)` — i.e. SNP **count** deciles, which controls gene length / marker-count confounding. Density (SNP/kb) is *not* a stratum.

**【Why it matters】**
1. **Robustness of the off-chr19 null.** The dominant GWAS signal is a single chr19 LD block (APOE). Excluding chr19 leaves 76 immune genes with no genome-wide-significant member; the permutation p=0.277 is far from significant. Residual SNP-density confounding could only matter if immune genes systematically sit in denser/sparser regions than matched background *within the same SNP-count decile* — a second-order effect unlikely to move p from 0.28 to <0.05. So the null is **robust**; the manuscript's hedging ("residual density confounding cannot be excluded") is appropriately conservative without undermining the conclusion.
2. **The all-autosome p=0.081 is the more interpretable number.** It is "borderline" and driven entirely by 4 chr19-adjacent immune genes (RELB ~0.37 Mb, IRF3 ~4.2 Mb, BCL3, TICAM1) that are inside or flanking the APOE LD block — i.e. "guilty by proximity," not independent immune signals. The manuscript states this correctly. I would foreground p=0.081 as "weak, chr19-driven" rather than letting "off-chr19 p=0.277" stand as the only GWAS gene-set number, because the former is the honest summary of the whole-autosome test.
3. **Minor**: the abstract reports "1,016 cases / 139,148 nominal controls" while the summary's effective per-SNP N is 134,310. This is expected post-QC and not an error, but a one-line reconciliation ("effective N per SNP = 134,310 after QC") would pre-empt a reviewer query.

**【Specific fix】** T3 — (a) keep the density caveat (it is correct); (b) add the all-autosome p=0.081 prominently alongside the off-chr19 p=0.277 in the Results, framed as "weak, chr19-driven"; (c) add the effective-N reconciliation sentence.

---

## CHECK 8 — Layer C: 2 POD vs 2 non patients; "minimum one-sided p = 1/3 (patient), ≈0.23 (library)"

**【Problem】** The stated patient-level minimum achievable one-sided MWU p ("exactly 1/3") is a **factor-of-two arithmetic error**; the correct exact minimum is 1/6. The qualitative conclusion (grossly underpowered, descriptive) is unaffected.

**【Evidence】**
- `layerC_summary.json`: proportion_test MWU p-values are 1.0, 0.857, 0.229 (Mono_classical), 0.229 (NK), 1.0, 1.0, 1.0 — i.e. the *observed* most-extreme library-level p is **0.2286** (manuscript "≈0.23" ✓, but this is the observed p, not the minimum achievable).
- Direct enumeration of the exact permutation distribution **[recomputed]**: for n1=n2=2, there are C(4,2)=6 equally-likely rank assignments; complete separation (all POD > all control) occurs in 1. The exact one-sided greater p = **1/6 = 0.1667**, *not* 1/3. For n1=4,n2=3, the minimum achievable one-sided greater p = 1/C(7,3) = **1/35 = 0.0286** (the observed 0.23 is far above this floor).
- The manuscript's "1/3" does not match either scipy's exact test (0.1667) or direct enumeration (0.1667). It is off by 2×.

**【Why it matters】** This is a small but concrete numerical error in a sentence whose *purpose* is to demonstrate that the Layer-C p-value is a floor, not evidence. If the floor itself is mis-stated, the pedagogical point is weakened (though the conclusion "underpowered, descriptive" is even *stronger* with the correct 1/6). It also slightly undercuts the manuscript's otherwise meticulous numeracy.

**【Specific fix】** T2 — correct "minimum achievable one-sided Mann–Whitney p is exactly 1/3" to "**1/6 ≈ 0.167**" (n1=n2=2), and note the library-level (4v3) floor is 1/35 ≈ 0.029 while the observed most-extreme p is 0.23. The "read as a floor, not biological evidence" framing is correct and should stay.

---

## CHECK 9 — Whole-transcriptome KO class test: z=−2.47; one-tailed upper p=0.9934; post-hoc lower p≈0.0066 (two-sided ≈0.013)

**【Problem】** The pre-specified one-tailed permutation logic is sound, but the post-hoc lower-tail p=0.0066 / two-sided ≈0.013 is presented as if it were a discovery, when it is an *exploratory* observation made after seeing the data.

**【Evidence】**
- `layerF`-adjacent KO class: manuscript line 112 states z=−2.47, pre-specified one-tailed upper p=0.9934, post-hoc lower-tail P(null≤obs)≈0.0066, two-sided≈0.013, "computed by hand … not by the script."
- Recomputation from z=−2.47 **[recomputed]**: one-sided upper-tail p = P(Z≥−2.47) = **0.9932** (manuscript 0.9934, rounding OK); one-sided lower-tail p = P(Z≤−2.47) = **0.0068** (manuscript ~0.0066, OK); two-sided = **0.0135** (manuscript ~0.013, OK). The arithmetic is internally consistent.
- Logic: the pre-specified hypothesis was "KO positive controls rank at the TOP" (upper tail). Observed z=−2.47 is at the *bottom* — i.e. KO of pro-inflammatory drivers *mimics* the disease signature (wrong direction for a "reverse" positive control). So p_upper=0.99 correctly shows the pre-specified directional test is non-significant; the lower-tail 0.0066 indicates significant mimicry.

**【Why it matters】**
1. The pre-specified upper-tail test and its non-significance (p=0.99) are correct and honestly reported as a failed positive control. Good.
2. The post-hoc lower-tail p=0.0066 is **selective** — it was chosen *because* the observed z was negative. Reporting it as "two-sided ≈0.013 … in the wrong direction (mimicry)" is a fair *descriptive* reading, but it must not be presented as a confirmatory result. The manuscript does call it a "failed positive control," which is honest, but the phrasing "the significant two-sided p ≈ 0.013" could be misread as a significant *finding*. The two-sided p only says "the KO control class is not at the top — it is significantly at the bottom," which is exactly the mimicry caveat, not evidence for the two-axis hypothesis.
3. **Circular-list caveat (also in Limitations line 146).** The KO positive-control class includes four anchor genes (adora3, cd63, ltf, col18a1) that are themselves objects of study; the class z is therefore partially circular, and the manuscript already says a re-run excluding them is needed. This should be restated next to the z=−2.47 result, not only in Limitations.

**【Specific fix】** T2 — (a) label p=0.0066 / 0.013 explicitly as "post-hoc, exploratory, not pre-specified"; (b) state next to z=−2.47 that the KO class also contains the four anchor genes and is therefore partly circular; (c) keep the "failed positive control / mimicry" framing — it is correct.

---

## Cross-cutting design comment — traceability gap on the neutrophil +3.4 pp claim

**【Problem】** The Abstract and Results foreground "neutrophils +3.4 pp, p≈10⁻³" as the positive, reproducible Layer-A result, but the committed `layerA_deconv_summary.csv` (cell-type level, `paired_wilcox_p`) has **Neu = NA** (only CD8T 0.0004, Bcell 0.0177, Mono 0.0057 are populated). The neutrophil delta is cited from "rows where celltype = aNeu/Neutro and config = centCAB100i_a12" in a *different* (lineage-level) representation that is not the committed `layerA_deconv_summary.csv`.

**【Evidence】** `layerA_deconv_summary.csv` shows `Neu,NA`. The manuscript (line 50, 106) cites +3.38 pp / p=7.6e-4 (CP_nnls_norm) and +3.39 pp / p=8.3e-4 (CP_penalty_K1000) — these numbers are not visible in the committed cell-type summary.

**【Why it matters】** The neutrophil shift is the *one* Layer-A result the manuscript treats as solid ("the leucocyte shift is itself a positive, reproducible result"). If its p-value is not traceable to a committed artefact in the form cited, that undercuts the strongest claim. (The shift is almost certainly real — multiple references agree directionally — but the specific p should be reproducible from a committed file.)

**【Specific fix】** T2 — commit the lineage-level paired-delta table that actually contains the neutrophil row (or correct the citation to the file that does), so "+3.4 pp, p≈10⁻³" is verifiable.

---

## § Stands up (what I think is genuinely strong)

1. **Honest null/negative reporting.** Failed positive controls (statins p=0.85, NSAIDs p=0.78; whole-transcriptome KO class z=−2.47 mimicry) and the untested two-axis dissociation are carried as load-bearing arguments. This is rare and commendable.
2. **Clear separation of the two probe-level methylation tests** (signed-t directional p=0.041 vs |t| magnitude p=2.5e-4) in text — I found no conflation (Check 4).
3. **Appropriate GWAS framing.** "Monolithic on chr19, no MR feasible, null off chr19" is correctly hedged, density caveat included, and the all-autosome p=0.081 is honestly "chr19-driven" (Check 7).
4. **Pre-specification posture.** A-priori positive controls and the "expected-null MR" are documented; the B=0 descriptive-DMR caveat is stated in Methods and Limitations.
5. **Reproducible artefacts.** Every other quantitative claim I recomputed (DMP counts, DMR reduction %, Fisher p=1, GWAS N/λ/OR, KO z-math) matched the committed files exactly.

---

## § Questions for the authors

1. For the composition-attribution band (Check 2): can you re-run the DMR-level adjustment with the *same* 3-lineage (Neu+Lymph+Mono) covariate used at the probe level, so probe and DMR levels share one covariate definition? How does the "composition-attributable %" change?
2. For Set-1 (Check 3): is it computationally feasible to recompute the *actual* Set-1 limma moderated-t signed-t MWU after composition adjustment (not the Set-2 OLS surrogate)? If not, can you state plainly that the headline p=0.041 has *no* adjusted analogue?
3. For Layer C (Check 8): what convention produced "1/3"? Our enumeration gives 1/6 for n1=n2=2 — please reconcile.
4. For the neutrophil shift (traceability gap): which committed file contains the neutrophil paired-delta p-value, and can it be made verifiable?
5. For the KO class (Check 9): will you provide the z=−2.47 result *excluding* the four anchor genes, as you flag is needed?

---

## § What I actually checked

**Files read (allowed set only):**
- `manuscript/POD_immune_hub_manuscript_v1.md` (full, 226 lines).
- `analysis/results/GSE330869/`: `layerA_topDMP.csv`, `layerA_summary.txt`, `bumphunter_adjusted_summary.csv`, `bumphunter_immune_gene_before_after.csv`, `layerA_DMR_immune_gene_counts.csv`, `layerA_deconv_summary.csv`.
- `analysis/adjust/`: `adjust_summary.csv`, `set1_signed_t_adjusted.csv`, `adjust_report.txt`.
- `analysis/results/layerB/layerB_summary.json`, `analysis/results/layerC/layerC_summary.json`, `analysis/results/layerF/layerF_real_gwas_summary.json`.
- Scripts: `analysis/layerA_bumphunter_adjusted.R` (confirmed `B=0`), `analysis/layerA_set1_signed_t_adjusted.py`, `analysis/layerF_real_gwas.py`.

**Commands run (managed Python 3.13.12):**
- Fisher exact on the 6-DMP / 49-immune-gene overlap → p=1.0 (two-sided and one-sided).
- DMR reduction: 13,357→5,153 (61.4%) and →7,263 (45.6%); immune proportion 0.59%→0.47%/0.50%.
- Layer-C exact MWU minimum p by rank-enumeration: n1=n2=2 → 1/6=0.1667; n1=4,n2=3 → 1/35=0.0286.
- KO class z=−2.47 → upper p=0.9932, lower p=0.0068, two-sided 0.0135.

**Recomputed-vs-manuscript agreement:** DMP count (6) ✓; Fisher p=1 ✓; DMR 13,357/5,153/7,263 ✓; immune 79/24/36 ✓; immune proportion 0.59/0.47/0.50% ✓; |t| enrichment 0.968 vs 0.903, p=2.5e-4 → 0.72/0.55 ✓; Set-2 OLS p=0.296/0.098/0.141 ✓; directional p=0.041 ✓; Layer B MWU p=0.0079, one-sample p=0.18 ✓; GWAS N=134,310, λ=1.0133, OR=1.86, 37 chr19 SNPs, off-chr19 p=0.277, all-autosome p=0.081 ✓; KO z-math ✓.
**Discrepancies found:** (a) Layer-C patient-level minimum p stated as 1/3 but recomputes to 1/6 (factor-of-2 error); (b) the neutrophil +3.4 pp / p≈10⁻³ is not traceable to the committed `layerA_deconv_summary.csv` (Neu row NA) — needs the lineage-level table; (c) the "46–61% composition-attributable" band mixes two covariate specifications and is descriptive (B=0) — see Check 2.

---

## Overall verdict & T0–T3 tiers

**Recommendation: Major revision (revise-and-resubmit).** No fatal/blocking (T0) flaw; the science is honest and the integration is coherent. The blocking-quality issues are about *precision of uncertainty quantification*, concentrated in Layer A's composition-attribution framing and the Set-1/Set-2 comparability hedge.

- **T0 (blocking):** none.
- **T1 (major, must address):**
  - Composition-attribution "46–61%" band (Check 2): not a sensitivity range; mixes covariate sets; DMR-level immune proportion is flat so the *immune-specific* claim rests on the probe-level |t| test. Re-frame as descriptive upper bound; use one consistent covariate set.
  - Set-1 vs Set-2 incomparability (Check 3): promote the "unadjusted, no adjusted analogue" hedge to the Abstract/Results first mention; add a compact all-p-values table.
- **T2 (minor, should address):**
  - Layer-C "1/3" → correct to 1/6 (Check 8).
  - Layer-B "one-sample t p=0.18 as lower bound" rephrase (Check 6).
  - KO post-hoc lower-tail p flagged explicitly as exploratory; restate anchor-gene circularity next to z=−2.47 (Check 9).
  - Neutrophil +3.4 pp traceability gap (cross-cutting).
  - bumphunter B=0 "(descriptive)" appended to every DMR-count and derived-% occurrence (Check 5).
- **T3 (trivial/editorial):**
  - Add comparative signed-t vs |t| caption (Check 4).
  - Foreground GWAS all-autosome p=0.081; add effective-N reconciliation sentence (Check 7).
