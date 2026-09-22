# C5 独立评审意见 — 药物重定位药理学 / LINCS & 分子对接方向
**评审专家代号：** C5（Repositioning pharmacology / LINCS & docking）
**角色定位：** LINCS L1000 重定位方法学效度、阳性对照诚实性、ISG 同义反复、KO 臂注意事项、对接方向性声明、以及任何可行动重定位候选是否被夸大。
**处理方式：** 视为首次投稿（fresh first-submission review），未读取任何前序评审材料；仅依据稿文、仓库 `analysis/results/layerD/`、`analysis/results/layerE/`、`analysis/pdb_audit.json`、`analysis/layerD_real_lincs.py` 及公开网页核验。
**评审日期：** 2026-09-22

---

## 总体判断（先给结论）

这是一份对方法学局限**异常诚实**的重定位筛选稿。Layer D 的核心结论——"无小分子化合物在药物水平通过 FDR、KO 臂为探索性、结论仅为假设生成"——与底层数据（immune summary `chem.n_druglevel_fdr05 = 0`、`crispr_ko.n_genelevel_fdr05 = 1109`）一致，且作者主动报告了失败阳性对照（他汀、NSAID）、KO 全转录本方向性失败（z = −2.47 模拟）、KO 阳性对照循环性、BH-FDR 过宽、consensus 集合而非 CMap τ 等缺陷。就"重定位可信度"而言，这是该领域少见的自我约束良好的稿件。

但存在 **两处准确性/诚实性缺陷**值得在修回前纠正（均属 T2），以及若干清晰度问题（T3）。**无 T0（致命）或 T1（重大必须重做）问题**；没有发现数据捏造或方法学根本性错误。建议**小修（minor revision）**后接受。

最关键的一条是：**稿件在讨论/结论中宣称"no compound cleared FDR at the drug level"（无化合物在药物水平通过 FDR），但该声明与自身 ligand 臂数据直接矛盾**——immune 臂 ligand 子库有 **24 个** perturbagen 在 FDR<0.05 通过（`layerD_immune_summary.json: ligand.n_druglevel_fdr05 = 24`），且 6 个 IFNγ 扰动正是这 24 个里的成员。该声明若按字面理解会误导读者以为整屏在药物水平全为阴性，而实际上 ligand（细胞因子）臂产出了显著的扰动信号。详见 Item 9 与 Item 2。

---

## MANDATORY CHECK 1 — LINCS 方法学效度（consensus 集合 + Fisher 合并 + BH + 规模匹配零分布）

**【Problem】** 方法学框架本身合理，但**显著性计算与排序指标之间存在不一致**，且稿件未向读者讲清，影响对"命中"定义的理解。

**【Evidence (file:line / CSV / web)】**
- `layerD_real_lincs.py:171-178`：`combine_fisher` 用 `stats.chi2.sf(x, 2*len(ps))`，对 2 个 p 值即 χ²(4 df)——与稿件 Methods（line 62）"combined by Fisher's method (χ², 4 df)"完全一致 ✓。
- `layerD_real_lincs.py:207-228`：`p_combined = combine_fisher(p_r1, p_r2)`，即 **FDR 仅基于"reversal"（逆转）两个方向 p 值合并**；`net_score = rev − mim` 是另外构造的。
- `layerD_real_lincs.py:232`：`df["fdr"] = bh_fdr(df["p_combined"])`——**BH-FDR 作用对象是 reversal-only 的 p_combined**。
- 但 `layerD_real_lincs.py:465-470` 与稿件 top25/排序：**排序（net_score）使用的是 rev−mim，而 FDR/emp_p/z 基于 reversal-only**。两者并非同一量。
- 稿件 Methods line 62 仅写"combined by Fisher's method … with BH-FDR across all signatures"，未说明 FDR 是 reversal-only 而非基于 net_score。

**【Why it matters】** 一个药物若"逆转"与"模拟"同时很强，其 net_score≈0，却仍可能获得很低的 reversal-only FDR。即"通过 FDR 的命中数"（如 ligand 臂 24 个）度量的是"逆转显著性"，而稿件用以呈现"领先候选"的是 net_score。对 repositioning 解释而言，读者应知道 24 个 ligand 命中是"逆转显著"，未必是"net 逆转显著"。该点在当前领先候选（IFNγ ligand，net 与 fdr 均强）上不改变结论，但作为方法学透明度应在稿件中讲清。

**【Specific fix】** 在 Methods（line 62 附近）明确一句："BH-FDR 与 empirical p 基于 Fisher 合并的 reversal-only p_combined；候选排序/展示用 net_score = reversal − mimic。两者量纲不同，net 为正且 fdr 显著者方视为双向一致命中。"若可行，报告一个"net_score 显著（或 mimic 部分被扣除后）"的敏感性计数，以佐证 24 这个 ligand 命中数不是 reversal-only 假象。

---

## MANDATORY CHECK 2 — chemical 与 ligand 子库的分离，以及是否混淆

**【Problem】** 稿件在 Methods 与 Results 主体基本保持了"小分子 chemical 子库（n=1,069，FDR 命中=0）"与"细胞因子 ligand 子库（n=96，命中=24）"的分离，这是**正确的**；但**在讨论/结论中出现了将二者合并叙述的过度概括**（见 Item 9）。本项核验 chemical/ligand 分离本身是否成立。

**【Evidence (file:line / CSV)】**
- `layerD_immune_summary.json`：`libraries.chem.n_perturbagens = 1069`，`n_druglevel_fdr05 = 0`；`libraries.ligand.n_perturbagens = 96`，`n_druglevel_fdr05 = 24`。
- `layerD_whole_tx_summary.json`：chem `n_druglevel_fdr05 = 0`；ligand `n_druglevel_fdr05 = 0`（whole-tx 臂 ligand 也 0）。
- 稿件 Methods line 64："These two sub-libraries are analysed separately below and must not be conflated: the term 'drug-level FDR hits' refers to whichever sub-library is named"——与数据一致 ✓。
- 稿件 Results line 110 明确写了"the glucocorticoid class…"（chemical 臂）+ "six IFNγ perturbations…"（ligand 臂），二者未混淆 ✓。

**【Why it matters】** 分离本身是对的，应保留。问题在更高层表述（Item 9）。本项确认：底层数据与方法描述在 chemical vs ligand 层面**没有混淆**，chemical 臂确实 0 个药物水平 FDR 命中，ligand 臂 24 个。

**【Specific fix】** 维持分离措辞；仅在讨论/结论层补充 ligand 臂 24 命中（见 Item 9），避免"no compound cleared FDR"这一笼统否定。

---

## MANDATORY CHECK 3 — 糖皮质激素类的 p 值与百分位数校正诚实性

**【Problem】** 稿件将糖皮质激素类称为"marginal directional signal / class-level pass"，并称其**未达 a priori p<0.05 阈值**。经核验，这一**自我约束是诚实且准确的**；仅提示单侧 p 的脆弱性。

**【Evidence (file:line / CSV)】**
- `layerD_immune_control_class_tests.csv`：glucocorticoids `mean_percentile=58.46`，`null_mean=49.92`，`null_sd=4.87`，`z=1.751`，`perm_p_one_sided=0.0407`（≈稿件所述 0.041）；按双侧即 ≈0.080，与稿件 line 110 完全一致 ✓。
- `layerD_whole_tx_control_class_tests.csv`：glucocorticoids `perm_p_one_sided=0.0556`（稿件称 whole-tx 0.056）✓；成员：dexamethasone 95.1%、hydrocortisone valerate 92.0%、betamethasone acetate 90.9%、prednisolone hemisuccinate 87.2%——与稿件 line 110 所列 whole-tx 百分位数**逐字一致** ✓。
- immune-restricted 对应成员（同表）：dexamethasone acetate 76.6%、hydrocortisone 79.4%、betamethasone 27.4%、prednisolone 76.3%——与稿件 line 110 完全吻合 ✓；另注意 immune 臂中 hydrocortisone valerate 仅 5.1%、betamethasone 27.4%，说明类别均值由少数高位成员拉动。
- 稿件 line 110 明确"does not meet the a priori p < 0.05 threshold"且"directional, class-level signal rather than a definitive positive-control pass"——诚实 ✓。

**【Why it matters】** 这是少数作者主动把"边界显著"降级为"方向性信号"的案例，符合重定位诚实性要求。需提醒：单侧 p=0.041 紧贴 0.05，且 whole-tx 臂（p=0.056）单侧都不达阈值；该信号脆弱，不应被读者解读为"糖皮质激素保护 POD 的阳性对照通过"。

**【Specific fix】** 当前表述已可接受。可加半句强调单侧检验的边界性："one-sided p=0.041 is within noise of the 0.05 boundary and the whole-transcriptome arm (p=0.056) fails even one-sided, so the signal should not be read as a passed positive control."（稿件其实已隐含此意，仅建议显式化。）

---

## MANDATORY CHECK 4 — 他汀与 NSAID 作为失败预设阳性对照

**【Problem】** 核验确认：他汀与 NSAID 在两个臂均为**失败阳性对照**，且稿件"failed positive controls temper confidence"的框架**诚实**。

**【Evidence (file:line / CSV)】**
- `layerD_immune_control_class_tests.csv`：statins `perm_p_one_sided=0.853`；nsaids `perm_p_one_sided=0.784`。
- `layerD_whole_tx_control_class_tests.csv`：statins `perm_p_one_sided=0.798`；nsaids `perm_p_one_sided=0.931`。
- 稿件 line 110："statins (immune-restricted perm p = 0.85; whole-transcriptome p = 0.80) and NSAIDs (p = 0.78; 0.93)"——与 CSV 逐值一致 ✓，并称"failed positive controls and temper confidence… report the screen as hypothesis-generating, not as having fully passed its positive-control battery" ✓。

**【Why it matters】** 失败阳性对照被如实报告（而非 rescue 或隐匿），是重定位稿件可信度的关键加分项。两个最具临床相关性的抗炎/降脂类（他汀有保护性 RCT、NSAID 常作对照）均未逆转 signature，说明本屏的"敏感性"有限——作者据此下调结论强度，方向正确。

**【Specific fix】** 无需修正。可在 Limitations 再补一句：由于两个临床最相关的抗炎类（他汀、NSAID）均失败，即便糖皮质激素边际信号成立，也不能外推为"泛抗炎即有益"，与稿件讨论 line 131 的"specificity rather than a blanket anti-inflammatory effect"呼应即可。

---

## MANDATORY CHECK 5 — 6 个 IFNγ 扰动逆转 ISG 主导 signature，及 ISG 同义反复的可行动推断

**【Problem】** 6 个 IFNγ 扰动在 ligand 臂 FDR<0.05 逆转 signature 的**排序与最强/次强 FDR 可核验**；ISG 同义反复的警示**已存在且恰当**。但"可行动推断"（抑制 IFNγ/ I 型 IFN，如 JAK–STAT 阻断或 IFNAR 拮抗；"IFNγ 给药会加重 POD"）**本质上是对输入 signature 已有信息的重述**，LINCS 这一步提供的独立证据有限，应被明确降级为"假设生成/印证"，而非看似独立的新发现。

**【Evidence (file:line / CSV / web)】**
- `layerD_immune_summary.json: top25`：ifng-mcf7 `fdr=1.92e-4`（最强）、ifng-mdamb231 `fdr=3.36e-4`（次强）、ifng-bt20 `fdr=6.21e-4`、ifna-bt20 `fdr=7.79e-4`。稿件 line 110 称"strongest: ifng-mcf7 FDR = 1.9 × 10⁻⁴; second-strongest: ifng-mdamb231 FDR = 3.4 × 10⁻⁴" ✓；另称"weakest: ifng-hs578t FDR = 2.5 × 10⁻²"（6 个 IFNγ 跨 MCF7/MCF10A/HS578T/MDAMB231/BT20/SKBR3）。
- **核验范围说明**：允许的产物中未包含 ligand signatures CSV（仅有 chem 与 crispr_ko 的 signatures/genes CSV），故 6 个 IFNγ 中我直接从 immune summary 的 top25 核验到 3 个 ifng（mcf7、bt20、mdamb231）及 1 个 ifna（bt20，属 I 型 IFN 而非 IFNγ）；其余 3 个 ifng（mcf10a、hs578t、skbr3）及其 FDR 依赖稿件自身引述（最弱 ifng-hs578t=2.5e-2）与 ligand 臂 `n_druglevel_fdr05=24` 的容纳能力。该引述与 ligand 臂存在 24 个 FDR 显著扰动的内部一致性相符，故**接受但标注为部分核验**。
- 稿件 line 110 已写明"Because the disease signature is interferon-stimulated-gene dominated (IFNG, CXCL10, IRF7, IRF1, STAT1 among the 48 input up-genes), this reversal is partly tautological" ✓。
- 稿件 line 135 推断："inhibition of IFNγ / type-I-IFN signalling (e.g. JAK–STAT blockade or IFNAR antagonism)… IFNγ administration would be predicted to exacerbate rather than prevent POD" 并承认"limited by the tautology between the ISG-rich signature and the IFNγ perturbations"。

**【Why it matters】** 逻辑链本身**可防御**：若 POD 血液已呈活跃 IFNγ/ISG 程序（输入 signature 已含 IFNG、CXCL10、IRF7、IRF1、STAT1），则给予 IFNγ 在生物学上会强化而非逆转该程序——"IFNγ 加重 POD"是合理预测。但关键限定是：**这一预测几乎完全来自输入 signature 的 ISG 组成，而非 LINCS 的重定位步骤**。LINCS 在此只是"用 IFNγ 扰动再次确认 signature 是 ISG 驱动"，属同义反复式印证。因此把"可行动方向（JAK–STAT 阻断/IFNAR 拮抗）"呈现为 LINCS 的产出，会高估该步贡献。

**【Specific fix】**
1. 将 line 135 的"would be predicted to exacerbate"弱化为"is consistent with exacerbation"或"would be expected to exacerbate, pending experimental validation"，并显式说明该预测主要源自输入 signature 的 ISG 内容，LINCS 提供的是一致性印证而非独立证据。
2. 区分 IFNγ（II 型）与 ifna（I 型，top25 中 ifna-bt20 也显著）的受体/通路不同：IFNAR 拮抗针对 I 型 IFN，JAK–STAT 阻断则覆盖两者；建议把"IFNAR antagonism"标为针对 I 型、"IFNγ/IFNGR 阻断"为针对 II 型，避免把两类混为一谈。
3. 明确"24 个 ligand 命中里 IFN 家族占多少"，让读者看到该信号是否主要由干扰素类构成（若是，则进一步支持"ISG 主导"而非独立药物重定位发现）。

---

## MANDATORY CHECK 6 — KO 臂：1,109 vs 0 基因、BH-FDR 过宽、n_cell_lines 结构=1、非免疫顶命中

**【Problem】** 所有 KO 臂注意事项均**与数据一致且表述诚实**。仅有一处数字细节（"5208/5210"）与产物略有出入，且无关紧要。

**【Evidence (file:line / CSV)】**
- `layerD_immune_summary.json: crispr_ko.n_genelevel_fdr05 = 1109`；`layerD_whole_tx_summary.json: crispr_ko.n_genelevel_fdr05 = 0`——稿件 line 112"1,109 KO terms survived BH-FDR < 0.05 in the immune-restricted arm (vs 0 in the whole-transcriptome arm)" ✓。
- `layerD_immune_crispr_ko_genes.csv`：nupl2 `fdr=5.81e-8`、hla-dpa1 `2.21e-7`、itgb4 `4.24e-7`、pde4a `4.16e-6`、emb `4.16e-6`、plk5 `4.16e-6`、il4r `4.16e-6`——与稿件 line 112 所列 FDR 完全吻合 ✓；且 nupl2/itgb4/emb/plk5 在该 CSV 中 `n_cell_lines` 列=1、未标注免疫 plausibility，与"non-immune"判定一致 ✓。
- 稿件 line 112 称"the cell_line field is empty for 5208/5210 genes, so n_cell_lines is structurally 1"。我实测 `layerD_immune_crispr_ko_signatures.csv`：总数据行 5212，cell_line 为空的 5206 行（约 5206/5212，而非 5208/5210）。**差异由 GMT 条目数（5212）与基因数（5210）口径不同导致，属无关紧要的描述性小数；实质性结论"n_cell_lines 结构=1、无细胞系内重复结构可过滤"完全成立** ✓。
- 机制上 `layerD_real_lincs.py:133-141` 的 `parse_term` 对 KO 术语（仅"GENE"）回退为空 cell_line，`drug_level` 用 `g["cell_line"].nunique()`，空字符串计为 1 个唯一值——故"n_cell_lines 结构=1"是代码层面必然结果，理解正确 ✓。

**【Why it matters】** "BH-FDR 在微小重叠 Fisher 检验上过宽"（1109/5210≈21% 基因通过，而疾病 signature 仅 ~40 up / ~33 down 基因，生物学上不可能如此多真命中）的警示是恰当的方法学诚实。KO 臂被整体降级为"假设生成"，定位正确。

**【Specific fix】**
- "5208/5210"可改为"5206/5212 signatures（或 5208/5210 genes）"以与产物精确对齐；或简化为"essentially all KO terms carry an empty cell_line field, so n_cell_lines is structurally 1"——避免精确数字引发复核偏差。
- 可加一句说明 BH-FDR 过宽的另一来源：KO consensus 集合本身由单细胞系聚合、集合间非独立，多重检验相关结构未被 BH 假定（独立）满足，故 1109 应视为宽松上限。

---

## MANDATORY CHECK 7 — KO 阳性对照循环性（adora3/cd63/ltf/col18a1）

**【Problem】** KO 阳性对照包含 4 个锚点基因（adora3、cd63、ltf、col18a1）造成的**循环性声明准确**；各基因单独 KO 百分位与稿件一致；whole-tx z=−2.47 模拟失败**已正确陈述**。

**【Evidence (file:line / CSV / py)】**
- `layerD_real_lincs.py:71-74`：`KO_POSITIVE_CONTROLS = [..., "ltf", "cd63", "adora3", "col18a1"]`——4 个锚点确在对照列表 ✓。
- `layerD_immune_crispr_ko_genes.csv`：adora3 `percentile=25.28%`（fdr 0.428，net −1.11）、cd63 `17.10%`（fdr 0.946，net −1.64，强模拟）、ltf `41.55%`（fdr 0.0618）、col18a1 `44.17%`（fdr 0.557）——与稿件 line 146 所列 25.28% / 17.10% / 41.55% / 44.17% 完全一致 ✓。
- `layerD_immune_summary.json: positive_control_class_test.crispr_ko`：`z=0.427`，`perm_p_one_sided=0.342`——稿件称"immune-restricted z = +0.43, p = 0.34" ✓。
- `layerD_whole_tx_summary.json: positive_control_class_test.crispr_ko`：`z=-2.472`，`perm_p_one_sided=0.9934`——稿件 line 112"KO whole-tx z = −2.47 … p = 0.9934 (non-significant for … upper-tail)… post-hoc lower-tail P ≈ 0.0066 (two-sided ≈ 0.013)… in the wrong direction (mimicry), i.e. a failed positive control" ✓。

**【Why it matters】** 作者把"anchor 基因同时充当 KO 阳性对照"识别为循环性并建议"剔除 4 个锚点后重跑"，是比多数稿件更严格的处理。该框架正确，无需质疑。

**【Specific fix】** 仅一处可澄清：`positive_control_class_test` 的 `n_controls_found` 在两份 summary 中均为 **18**，而 `KO_POSITIVE_CONTROLS` 仅列 **12** 个基因——差值是 `control_set_enrichment` 用 `str.contains` 子串匹配导致（如 "stat3" 可能命中多个别名条目）。建议在方法或脚注说明该 18 的来源，避免读者误以为指定了 18 个对照。属 T3 清晰度问题。

---

## MANDATORY CHECK 8 — 对接（Layer E）：仅审计、无姿态；[20]/[23]/[24] 异质性与 ADORA3 PDB 数

**【Problem】** 对接仅做"实验结构可得性审计"，无 docking pose，该边界**诚实**；[20] 对 A3 拮抗方向的支撑、**三篇腺苷受体文献的亚型异质性声明、ADORA3=5 个 PDB、TMIGD3/COL13A1/SPATA13=0** 均**经核验准确**。

**【Evidence (file:line / CSV / web)】**
- `layerE_dockability_audit.csv` 与 `pdb_audit.json`：ADORA3 `n_pdb_entries=5`（8X16, 8X17, 9EBH, 9EBI, 9EHS）；TMIGD3/COL13A1/SPATA13 `=0`；COL18A1=9、CD63=2、LTF=10——与稿件 line 122 完全一致 ✓。
- 稿件 line 145："ADORA3 … A3 antagonism (counterintuitive, because A3 signalling counteracts the A2A-mediated inhibition of TLR responses in microglia [20])" + "[20] supports A3, [23] supports A2A, and [24] supports A1/A2A".
- **Web 核验（PubMed 19494284 / J Immunol 2009;182:7603）**：van der Putten 2009 摘要原文"uncover a role for A(3)R as dynamically regulated suppressors of A(2A)R-mediated inhibition of TLR-induced responses. This would suggest exploration of combinations of A(2A)R agonists and A(3)R antagonists to dampen microglial activation."——**明确支持 A3R 拮抗可抑制小胶质细胞活化**，稿件对 [20] 的引用方向正确 ✓。
- 异质性：检索确认 [23] Rebola 2011（A2A 调控神经炎症）、[24] Martí Navia 2020（A1 激动剂与 A2A 拮抗剂）——三篇确实分属 A3 / A2A / (A1+A2A) 不同亚型，稿件"不应被读作对 ADORA3 的统一背书"的提醒恰当 ✓。

**【Why it matters】** 对接层最大的诱惑是"既然有结构就暗示可成药"，稿件克制地仅做可行性审计并把无实验结构的 3 个锚点列为局限，方向正确。A3 拮抗的"反直觉"逻辑（A3 促进 TLR 响应，故拮抗 A3 抑制神经炎症）经 [20] 支撑成立，是可防御的假设。

**【Specific fix】**
- 建议在 Discussion 补充：ADORA3 在 PBMC 几乎不表达（Layer C），但其作为"脑/内皮"免疫调节靶点的定位需要一个**组织表达证据缺口**说明——即当前无 POD 脑/内皮转录证据直接支持 ADORA3 为枢纽， docking 可行性≠药理相关性。稿件 line 145 已提及"PBMC 非相关组织"，可再显式补一句"亦无 POD 脑组织 ADORA3 表达证据，故仍属假设"。
- TMIGD3 被描述为"transmembrane protein with no characterised ligand-binding pocket"——此判断超出本次核验范围（未在允许产物中），但作为 docking 优先级排序理由可接受；若可行引用结构注释来源更佳。

---

## MANDATORY CHECK 9 — "no compound cleared FDR at the drug level" 的准确性（**核心准确性缺陷，T2**）

**【Problem】** 稿件在 Discussion（line 135）与 Conclusions（line 153）多次称"no compound cleared FDR at the drug level"（无任何化合物在药物水平通过 FDR），但该陈述**与自身 ligand 臂数据直接矛盾**，且会误导读者认为整屏在药物水平全阴性。

**【Evidence (file:line / CSV)】**
- `layerD_immune_summary.json`：`libraries.ligand.n_druglevel_fdr05 = 24`（immune 臂 ligand 子库有 24 个 perturbagen 在 FDR<0.05 通过）。
- 同文件 `libraries.chem.n_druglevel_fdr05 = 0`（chemical 小分子臂确为 0）。
- 稿件 Results line 110 自己报告了"six IFNγ perturbations … reversed … at FDR < 0.05"——这 6 个正是 ligand 臂 24 命中里的成员。即稿件**既报告了 ligand FDR 命中，又在讨论中宣称"无化合物通过 FDR"**，内部不一致。
- `layerD_whole_tx_summary.json`：ligand `n_druglevel_fdr05 = 0`，chem `=0`——说明该矛盾仅存在于 immune-restricted 臂的 ligand 子库。

**【Why it matters】** 这是本次评审最重要的诚实性/准确性问题。若"compound"被严格读作"小分子药物"，则 chemical 臂 0 命中属实；但 ligand 臂的 24 个扰动（含 6 个 IFNγ 细胞因子）也是 Layer D 的产出，且稿件 Methods line 64 明言 ligand 是独立于 chemical 的子库、有自己的命中计数。笼统的"no compound cleared FDR at the drug level"会让读者（及编辑）误以为 Layer D 在药物水平无任何显著信号，低估了 ligand 臂产出的 ISG/干扰素证据，也与稿件自身 line 110 矛盾。更严重的是，它可能削弱稿件"诚实报告阴性"的立论基础——因为此处恰恰是把一个**阳性**子结果在高层叙述中抹去了。

**【Specific fix】** 将 line 135 与 line 153 的笼统否定改为精确表述，例如：
> "No **small-molecule (chemical-arm)** compound cleared FDR at the drug level (chem `n_druglevel_fdr05 = 0` in both signatures). The cytokine ligand arm, by contrast, produced 24 FDR-significant perturbations (including the 6 IFNγ entries), but these are endogenous signalling molecules rather than repurposable drugs and are interpreted as pathway-level—not compound-level—evidence; the virtual-knockout arm remains exploratory."
这样既保留了"chemical 臂全阴性"的事实，又不再与 ligand 臂 24 命中自相矛盾，且符合稿件一贯的"假设生成"基调。

---

## 额外观察（不属于 9 项强制核验，但影响重定位可信度）

**A. Results 小标题"chemical arm valid"（line 108）过度承诺（T2/T3）。**
正文随后即把糖皮质激素类降级为"marginal directional signal"、他汀/NSAID 标为失败、chemical 臂 0 药物水平命中。小标题"chemical arm valid, virtual-knockout arm hypothesis-generating"与正文张力过大。建议改为"chemical arm partially supportive / hypothesis-generating; virtual-knockout arm exploratory"，与正文一致。

**B. reversal-only FDR 与 net_score 排序的量纲不一致（见 Item 1，T3）。** 建议在 Methods 显式区分。

**C. 单侧 p=0.041 的边界性（见 Item 3，T3）。** 已诚实陈述，建议显式化其脆弱性。

**D. "nominal FDR"措辞（line 112）**。FDR 本身已是校正量，"nominal FDR"易歧义，建议改为"BH-adjusted FDR"或"top FDR values"。

**E. KO 对照 n_controls_found=18 vs 指定 12（见 Item 7，T3）。** 子串匹配导致，建议脚注说明。

**F. 重定位"可行动候选"的提炼不足。** 稿件正确地把领先候选定为"干扰素通路"与"糖皮质激素类"。但干扰素通路候选本质上是输入 signature 的回响（Item 5），糖皮质激素类是边界显著且临床证据混杂（dexamethasone meta 阴性 vs 单中心 RCT 阳性，稿件已引 [15]/[22]）。建议在摘要/结论更明确：本屏**未产出任何可直接进入临床验证的具体已获批药物候选**（chemical 臂 0 命中），"候选"是通路层面的（IFN 抑制、糖皮质激素类），而非分子层面的。这比当前"interferon-pathway / glucocorticoid-class pathways as the leading … candidates"更不夸大。

**G. 规模匹配零分布的实现细节（T3）。** `layerD_real_lincs.py:251-252`：`null_net = null_rev - null_rev.mean()`，即 net 的零分布被强制中心化；且 `z_vs_null` 基于 `reversal_score` 而非 `net_score`。这与 Item 1 同源，建议一句说明"z/emp_p 针对 reversal_score，net 仅用于排序展示"，以免方法透明性被质疑。

---

## § Stands up（稿件站得住脚之处，≥3）

1. **方法学局限的系统性诚实**：失败阳性对照（他汀 0.85/0.80、NSAID 0.78/0.93）、KO 全转录本方向性失败（z=−2.47 模拟）、KO 阳性对照循环性、BH-FDR 过宽、consensus 集合而非 CMap τ——均被主动报告，且逐值与产物一致。这在 LINCS 重定位稿件中罕见且可贵。
2. **chemical 与 ligand 子库在主体分析中保持分离**，底层 `n_druglevel_fdr05`（chem=0 / ligand=24）与方法叙述无混淆（高层结论除外，见 Item 9）。
3. **对接层严格限定为"可行性审计"**：ADORA3=5 PDB、TMIGD3/COL13A1/SPATA13=0、无实验结构者列为局限——边界清晰，无"有结构即暗示可成药"的越界。
4. **[20] 对 A3 拮抗方向的引用经 Web 核验准确**，且三篇腺苷受体文献的亚型异质性提醒恰当，未把异质子文献包装成对 ADORA3 的统一背书。
5. **糖皮质激素类的"边际方向性信号/未达 a priori 阈值"表述诚实**，百分位数（immune vs whole-tx）逐值与 CSV 吻合，未 rescue。
6. **所有定量主张可追溯**：本次核验的 9 项强制检查中，8 项的数值（FDR、百分位、z、PDB 计数、KO 命中数）均与 `analysis/results/` 产物逐字一致，未见数字捏造。

---

## § Questions for the authors（向作者提问）

1. 在 ligand 臂 `n_druglevel_fdr05 = 24` 中，干扰素家族（IFNγ/IFNα/IFNβ 等）占多少？若该信号主要由干扰素类构成，是否进一步支持"signature 是 ISG 主导"而非独立的药物重定位发现？
2. 24 个 ligand 命中里，剔除干扰素类后还剩哪些非冗余通路？是否有任何**非干扰素**的可行动（通路级）候选值得在讨论中点名？
3. KO 对照 `n_controls_found=18` 与 `KO_POSITIVE_CONTROLS` 指定的 12 个基因不一致，差异来源是否为 `str.contains` 子串匹配？能否报告剔除 4 个锚点基因后的非循环 KO 类别 z（如稿件 line 133/146 所承诺的"re-run excluding them"）？
4. reversal-only FDR 与 net_score 排序的量纲差异，是否改变 ligand 臂 24 命中或 chemical 臂 0 命中的结论？（预期不改变，但请确认。）
5. ADORA3 作为"脑/内皮"靶点的组织表达证据缺口：除 PBMC 低表达外，是否有任何 POD 相关脑/内皮转录或蛋白证据支持其作为枢纽？若无，是否同意在讨论中将其明确标为"纯假设级 docking 候选"？

---

## § What I actually checked（我实际核验的内容、网页检查与差异）

**读取并核验的稿件/产物文件：**
- `manuscript/POD_immune_hub_manuscript_v1.md`（全文 226 行）。
- `analysis/results/layerD/layerD_immune_summary.json`——chem/ligand/crispr_ko 的 n_perturbagens、n_druglevel_fdr05（chem=0、ligand=24、crispr_ko=1109）、top25（含 4 个 IFN ligand FDR）。
- `analysis/results/layerD/layerD_whole_tx_summary.json`——chem/ligand=0、crispr_ko=0；whole-tx 对照组（glucocorticoids p=0.0556、statins 0.798、nsaids 0.931）；KO 对照 z=−2.472。
- `analysis/results/layerD/layerD_immune_control_class_tests.csv`——糖皮质激素类 mean 58.46 / null 49.92 / z 1.751 / p 0.0407；成员百分位（dexamethasone acetate 76.6%、hydrocortisone 79.4%、betamethasone 27.4%、prednisolone 76.3%）；statins 0.853、nsaids 0.784。
- `analysis/results/layerD/layerD_whole_tx_control_class_tests.csv`——糖皮质激素类 p 0.0556；成员（dexamethasone 95.1%、hydrocortisone valerate 92.0%、betamethasone acetate 90.9%、prednisolone hemisuccinate 87.2%）；statins 0.798、nsaids 0.931。
- `analysis/results/layerD/layerD_immune_crispr_ko_genes.csv`——nupl2/hla-dpa1/itgb4/pde4a/emb/plk5/il4r 的 FDR；adora3 25.28%、cd63 17.10%、ltf 41.55%、col18a1 44.17%（含 net/fdr 一致）。
- `analysis/results/layerE/layerE_dockability_audit.csv` 与 `analysis/pdb_audit.json`——ADORA3=5、TMIGD3/COL13A1/SPATA13=0、COL18A1=9、CD63=2、LTF=10。
- `analysis/layerD_real_lincs.py`——combine_fisher χ²(4df)、BH on p_combined、KO_POSITIVE_CONTROLS 含 4 锚点、cell_line 空导致 n_cell_lines=1 的代码机制。

**网页核验（WebSearch/WebFetch）：**
- [20] van der Putten et al., *J Immunol* 2009;182:7603（PubMed 19494284）：确认"A(3)R suppresses A(2A)R-mediated inhibition of TLR responses"，并原文建议"A(3)R antagonists to dampen microglial activation"——稿件 A3 拮抗方向引用**准确**。
- 异质性：[23] Rebola 2011（A2A）、[24] Martí Navia 2020（A1 激动剂 + A2A 拮抗剂）确认分属不同亚型，稿件"三篇不应统一读作 ADORA3 背书"的提醒恰当。
- LINCS/Enrichr 方法学（Subramanian 2017 *Cell*；Enrichr Kuleshov 2016；LINCS workflow）：确认 consensus 基因集合 overlap-enrichment 是 LINCS/Enrichr 生态的标准用法，且存在全秩 CMap τ 的更高分辨率替代——稿件"consensus sets, not full-ranked CMap τ"的声明属实。

**差异/未决项（discrepancies）：**
- "no compound cleared FDR at the drug level"（line 135/153）与 ligand 臂 `n_druglevel_fdr05=24` 矛盾（详见 Item 9，T2）。
- "5208/5210 empty cell_line" 实测为 `layerD_immune_crispr_ko_signatures.csv` 中 5206/5212 空（口径差异，无关紧要，T3）。
- KO 对照 `n_controls_found=18` vs 指定 12（子串匹配，T3）。
- 6 个 IFNγ 中仅 3 个（mcf7/bt20/mdamb231）+1 个 ifna(bt20) 能从允许产物直接核验；其余 3 个 ifng 依赖稿件自引（部分核验，已标注）。

---

## 总体评语与层级（Overall verdict + T0–T3 tiers）

**总体评语：** 这是一份方法学上合理、且对局限异常诚实的 LINCS 重定位筛选稿。强制核验的 9 项中，8 项的数值与产物逐字一致，方法学框架（consensus 集合 + Fisher-one-sided-greater + Fisher 合并 χ²4df + BH-FDR + 3000 规模匹配零分布）符合 LINCS/Enrichr 生态标准用法，且作者明确声明未使用全秩 CMap τ。未发现数据捏造或致命方法学错误。主要缺陷集中在**两层表述**：(1) Discussion/Conclusions 的"no compound cleared FDR at the drug level"与 ligand 臂 24 命中自相矛盾（准确性/诚实性，T2）；(2) IFNγ 可行动推断本质是输入 signature 的回响，LINCS 贡献为同义印证，宜降级（T2）；(3) Results 小标题"chemical arm valid"过度承诺（T2/T3）。其余为清晰度/透明度问题（T3）。

**Tier 分级：**
- **T0（致命/必须撤回重做）：无。**
- **T1（重大/必须修回前解决）：无。**
- **T2（应修回前修正）：**
  - Item 9：将"no compound cleared FDR at the drug level"改为精确表述，区分 chemical 臂（0）与 ligand 臂（24，通路级而非化合物级）。
  - Item 5：把 IFNγ 可行动推断降级为"与输入 ISG signature 一致、LINCS 提供印证而非独立证据"，并区分 IFNγ(II 型) 与 IFNAR(I 型) 通路。
  - 额外 A：Results 小标题"chemical arm valid"改为与正文一致的"partially supportive / hypothesis-generating"。
- **T3（清晰度/透明度，建议修回时一并处理）：**
  - Item 1/G：显式说明 FDR 基于 reversal-only p_combined、排序用 net_score、z/emp_p 针对 reversal_score。
  - Item 3：显式化单侧 p=0.041 的边界脆弱性。
  - 额外 D："nominal FDR"→"BH-adjusted FDR"。
  - Item 7/E：说明 KO 对照 n_controls_found=18 的来源（子串匹配）。
  - Item 6：将"5208/5210"对齐为产物实际值或泛化为"essentially all"。

**推荐决定：** Minor revision（小修后接受）。无需要补实验或重跑核心分析；所有 T2 项均可通过文字精确化与层级下调解决，且不会削弱稿件"假设生成"的核心立场——事实上精确化后会**增强**其诚实性立论。
