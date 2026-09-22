# REVIEW_round8_2026-09-22 — 独立多专家评审合并报告

**对象稿件：** `manuscript/POD_immune_hub_manuscript_v1.md`（v1.1, Round-7 修订, 2026-09-22）
**目标期刊：** *Journal of Neuroinflammation* (JNI)
**评审模式：** 首次投稿姿态（fresh first-submission review），强制独立性
**评审日期：** 2026-09-22
**编辑（合并人）：** 小团 / 永新（作为编辑对最严重发现做亲自复核）

---

## 1. 独立性声明（机制 + 证据）

**机制。** 五名专家（C1 领域/临床、C2 设计/统计、C3 溯源/重算、C4 期刊/报告规范、C5 重定位药理/LINCS+对接）被分别派驻，各自只读取稿件正文、仓库产物 CSV/JSON/脚本以及公开网页，被**明令禁止**读取任何 `REVIEW_*.md`、`RESPONSE_*.md`、`REVISION_*.md`、`PROGRESS.md`、G1–G4 gate 结果、memory、以及本目录中其他专家的产出。每位专家被要求：把稿件当作首次投稿、每条论断必须来自自己读到的文本或源数据、凡能核实的数字必须亲自重算。

**证据独立性确实生效（来自各报告头部自陈）：**
- **C1** 明确写道："I did **not** open any `REVIEW_*.md` … `review_round8/` sibling files, `PROGRESS.md`, gate results, or other expert reports, per the review brief."
- **C2** 明确写道："I did **not** consult any prior review, response, or revision material."
- **C3 / C4 / C5** 均自陈为 "fresh first-submission review / 未读取任何前序评审材料"，并各自列出了实际读取的允许文件清单（不含任何前序 `REVIEW_*`）。

**独立性的诊断性信号：** 多名互不知情的专家从不同角度命中了同一类缺陷——(a) C1-F2 与 C3/C4 均指向 **[41] Evered 引用错误**；(b) C2-Check2 与 C3-§1b 均指向 **Layer A 组成校正表述的口径/归因问题**；(c) C2-Check8 与 C3 均独立发现 **Layer-C 最小 p 值算术错误**（虽 C2 给出 1/6、C3 未复核该项，但结论一致）。这说明独立性纪律有效，而非"各自发现互不相交的孤立问题"。

**编辑亲自复核（skill §5 强制）：** 在出报告前，编辑用**自己的脚本**直接读原始产物（不经任何审稿人中间结果）复核了以下"最严重发现"与"可能改变定级的发现"，并在表中以 **[编辑亲核]** 标注：
1. **[41] Evered 引用**：欧洲 PMC/`EXT_ID:30325806` 命中 Anesthesiology 2018;129(5):872–879，DOI 10.1097/ALN.0000000000002334 正确，作者列表不含 "Whittaker P / Pu Y"。稿件 30179813 指向一篇无关的短链氯化石蜡论文。**→ 编辑亲核确认 T1。**
2. **Layer-C 最小单侧 p**：对 n1=n2=2 做秩枚举，C(4,2)=6 种等可能分配，完全分离仅 1 种 → 精确单侧 p = **1/6 = 0.1667**，非稿件所写 1/3（差 2×）。**→ 编辑亲核确认 T2（算术错误）。**
3. **ligand 臂 `n_druglevel_fdr05`**：`layerD_immune_summary.json` 显示 chem=0 / ligand=**24**。稿件 "no compound cleared FDR at the drug level" 与自身 ligand 24 命中自相矛盾。**→ 编辑亲核确认 T2。**
4. **[6] Seki 2026 是否真实**：编辑 WebSearch 确认 Seki T et al. *Transl Psychiatry* 2026;16(1):349，**PMID 42143058 / DOI 10.1038/s41398-026-04067-6 真实存在**（与本项目 G1 门 E-utilities 核实一致）。C4 称 "EuropePMC hitCount 0" 系查询/索引滞后所致。**→ [6] 非缺陷，予以降级撤销（false alarm）。**
5. **免疫基因数 48/21/29 是否矛盾**：`bumphunter_adjusted_summary.csv` 中 immune_genes = 48/21/29 为权威值；三个 overlap CSV 各 79/24/36 **行**、每行 `n_immune_genes`=1（即每行重叠恰好 1 个免疫基因）。故 79/24/36 是 **DMR 行数**，48/21/29 是 **基因数（来自另一份汇总文件）**，稿件数字本身均正确，缺陷仅在 **line 50 的归因错置**与 **48/75/79 三文件不一致未解释**。**→ C3 "三 artefact 矛盾" 的"矛盾"措辞为 false alarm，降级为 T2 溯源澄清项。**

---

## 2.  verdict 表（逐专家 + 分布）

| 专家 | 角色 | verdict | 关键升级点 |
|---|---|---|---|
| **C1** | 领域/临床 + 文献 | **Major revision** | [41] 引用错误（T1）；缺神经炎症教条文献（T2）；"modifiable target"措辞过强（T2） |
| **C2** | 设计/统计 | **Major revision** | 组成校正"46–61%"非灵敏度区间（T1）；Set-1/Set-2 不可比对冲不足（T1）；Layer-C 1/3→1/6（T2） |
| **C3** | 溯源/重算 | **Minor revision** | 免疫基因数 line-50 归因（T2）；lymphoid −1.67 vs −1.60（T2）；其余数字基本复现 |
| **C4** | 期刊/报告规范 | **Revise（moderate）** | 缺强制 generative-AI 声明（T1）；伦理委员会未具名（T2）；[41] PMID（T2） |
| **C5** | 重定位药理 | **Minor revision** | "no compound cleared FDR" 与 ligand 24 自相矛盾（T2）；IFNγ 推断宜降级（T2） |

**分布：** 0 份 DESK-REJECT；2 份 Major/Revise（C1、C2）；3 份 Minor（C3、C5、C4 之"moderate"可归 Minor-to-Major 之间）。
**编辑合并 verdict：** **Major revision（修改后重投）**。依据 skill §7——C3/C5 的宽松 verdict 仅认证"溯源/重定位"单个层面，不覆盖设计层（C2）与领域层（C1）的真实问题；采用更严格定级并说明：宽松方认证的是"一层"，非"整稿"。

---

## 3. 交叉核实表（最高价值产物）

| # | 稿件位置 | 稿件声称 | 独立重算值 | 复核人 | 结论 |
|---|---|---|---|---|---|
| 1 | L102/106 | DMR 13,357 / 5,153 / 7,263；免疫重叠 79/24/36 | 一致 | C2,C3 | ✓ 匹配 |
| 2 | L50 | 免疫基因 48/21/29（称来自 overlap CSV） | 48/21/29 在 `bumphunter_adjusted_summary.csv` 权威；overlap CSV 仅含 79/24/36 **行**（每行 1 基因） | C3 | ✗ **归因错置**（数字对、出处错） |
| 3 | L106 | 中性粒 +3.38 pp, p=7.6e-4 | 根 `layerA_deconv_summary.csv`（centCAB100i_a12, CP_nnls_norm）复现 +3.384 pp / 7.6e-4；原始分数复算一致 | C3 **[编辑亲核读数]** | ✓ 可复现（仅被引文件路径需澄清，见 T3-10） |
| 4 | L106 | lymphoid −1.67 pp, p=0.003 | 同 config 复算 −1.599 pp（≈−1.60）；cover letter 写 −1.6 pp | C3 | ✗ **不一致**（稿件 vs cover letter） |
| 5 | L106 | monocyte −1.11 pp, p=1.2e-4 | 复算 −1.108 pp / 1.2e-4 | C3 | ✓ 匹配 |
| 6 | L96 | Layer-C 最小单侧 p = 1/3 | 精确枚举 n1=n2=2 → **1/6 = 0.1667** | C2 **[编辑亲核]** | ✗ **算术错误（2×）** |
| 7 | L106 | \|t\| 富集 0.968/0.903, p=2.5e-4 → 0.72/0.55 | 一致 | C2,C3 | ✓ 匹配 |
| 8 | L102,141 | Set-1 signed-t p=0.041；Set-2 0.296/0.098/0.141 | 两值均复现，但**为不同检验**（基因集与 t 定义均不同） | C2 | ⚠ **不可比**（非数值错） |
| 9 | L116 | GWAS 37 chr19 SNP；N=134,310；λ=1.0133；OR=1.86；off-chr19 p=0.277；全常染 p=0.081 | 全部复现 | C2,C3 | ✓ 匹配 |
| 10 | L110 | 糖皮质激素 0.041/0.056；双侧 0.080 | 单侧 0.0407/0.0556 复现；双侧应 = 2×0.0407 = **0.0814**（非 0.080） | C3 | ✗ **舍入小错** |
| 11 | L110 | 6 个 IFNγ 扰动 FDR<0.05 | ifng-mcf7 1.9e-4 / mdamb231 3.4e-4 / hs578t 2.5e-2 等均复现 | C3,C5 | ✓ 匹配 |
| 12 | L112 | KO 1,109（immune）vs 0（whole-tx）；顶基因 FDR 精确 | 复现 | C3,C5 | ✓ 匹配 |
| 13 | L135,153 | "no compound cleared FDR at the drug level" | ligand 臂 `n_druglevel_fdr05`=**24**（含 6 IFNγ）；chem=0 | C5 **[编辑亲核]** | ✗ **与自身数据自相矛盾** |
| 14 | L122 | 对接 PDB：ADORA3=5, COL18A1=9, CD63=2, LTF=10, TMIGD3/COL13A1/SPATA13=0 | 一致 | C3,C5 | ✓ 匹配 |
| 15 | L27,41 | [41] Evered 2018 PMID 30179813 | 正确 PMID=**30325806**；作者列表含伪造的 "Whittaker P, Pu Y" | C1 **[编辑亲核]** | ✗ **引用错误（T1）** |
| 16 | L43,48,100,122 | [6] Seki 2026 PMID 42143058 | 真实存在（WebSearch + G1 门核实）；C4 "查无" 系索引滞后 | C4(误) / **编辑亲核证真** | ✓ **非缺陷（撤销 C4 该项）** |
| 17 | L92 | Layer B 免疫基因集 MWU p=0.0079；one-sample t p=0.18 "lower bound" | 数值复现；但 "lower bound" 措辞不标准（p=0.18 比 0.0079 更不显著，方向倒置） | C2 | ⚠ **措辞误导** |

> 表中 "✗" = 稿件与重算不符；"⚠" = 数值可复现但表述/口径问题；"✓" = 复现一致。**所有 [编辑亲核] 行均由编辑用原始产物独立复核确认。**

---

## 4. 分级合并问题清单（去重，保留各专家精确位置与数字证据）

### Tier 0 — 结论颠覆性 / 撤稿级
**无。** 五名专家一致认为无致命缺陷。编辑确认。

### Tier 1 — 重大，接受前必须解决（分析需增补 / 表述需结构性重框）
- **T1-1｜[41] Evered 2018 引用错误（PMID + 作者列表）。** 位置 L41（参考文献 [41]）。证据：稿件 PMID 30179813 实为短链氯化石蜡无关论文；正确 PMID=30325806；作者列表含非作者 "Whittaker P, Pu Y" 且漏列 DeKosky/Rasmussen/Oh/Crosby/Berger/Eckenhoff。来源 C1-F2 + **[编辑亲核]**。这是稿件引用的"命名法权威"，错误会引发对全部参考文献的连锁质疑。**修复：** 整条替换为 C1-F2 给出的正确条目（含 PMID 30325806）。
- **T1-2｜缺失 BMC 强制 generative-AI 使用声明。** 位置：Methods（L35–84 无 AI 句）；cover letter 亦无。证据 C4-§4：BMC/SN 政策要求 LLM 使用须在 Methods（无 Methods 时在合适位置）记录。这是"return-to-author"级政策缺口，与科学内容无关但必须补。**修复：** 在 Methods 增加一句（选真实版本）："No generative-AI tools were used in the conception, analysis, or writing of this manuscript"（若仅用语言润色则写明 "AI-assisted copy editing for grammar/style only; no content was AI-generated"）。cover letter 同步。
- **T1-3｜Layer A "46–61% 组成可归因" 表述误作灵敏度区间。** 位置 L102/106/141。证据 C2-Check2：该区间两端来自 DMR 层**两个不同协变量规格**（仅 aNeu 61% vs aNeu+aEos+aBaso 46%），并非灵敏度包络；且免疫 DMR 占比在调整前后基本持平（0.59%→0.47–0.50%），说明组成校除掉的是**全局**术后甲基化偏移，免疫特异性结论实际由探针层 \|t\| 富集检验（p=2.5e-4→0.72）承载。B=0（无置换）使其为描述性上限。**修复：** 改写为"至少约 46%、至多约 61% 的**总体** DMR 负担可由组成解释（描述性上限，受去卷积测量误差影响）"；明确扁平免疫占比；建议探针层与 DMR 层统一用同一协变量集（推荐 3-lineage Neu+Lymph+Mono），或把两种协变量作为**显式分离情景**呈现，而非单一"46–61%"。
- **T1-4｜Set-1 与 Set-2 为不可比检验，对冲不足。** 位置 L102/141（Abstract 与 Results 首次出现 p=0.041 处）。证据 C2-Check3：Set-1（49 基因/737 探针，limma  moderated-t，未校正，p=0.041）与 Set-2（83 基因/1,195 探针，OLS-on-delta，校正后 0.296/0.098/0.141）在**基因集与 t 定义上均不同**，Set-2 既不能证实也不能证伪 Set-1。稿件"not-yet-composition-validated"对冲存在但易被低估。**修复：** 在 p=0.041 首次出现处加行内说明"Set-1，未校正，limma moderated-t，无该精确检验的组成校正对应项"；增一张四行 p 值总表（检验对象 / 基因集 / t 定义 / 是否校正 / p）使不可比性可见；Abstract 软化"directional probe-set test … p=0.041"为"unvalidated (… its only adjusted analogue, a different Set-2 OLS test, was null: 0.296/0.098/0.141)"。

### Tier 2 — 中等，应修正（多为措辞/归因/局部重框）
- **T2-1｜Layer-C 最小单侧 p "1/3" → "1/6≈0.167"。** 位置 L96。证据 C2-Check8 + **[编辑亲核]**：n1=n2=2 秩枚举得 1/6；n1=4,n2=3 下限 1/35≈0.029，观测最极端 p=0.23。定性结论（欠幂、描述性）不受影响，但地板值写错削弱了"读数即地板"的教学点。**修复：** 改为 "1/6 ≈ 0.167 (n1=n2=2)；library 级 (4v3) 下限 1/35≈0.029，观测最极端 p=0.23"；"read as a floor, not biological evidence" 框架保留。
- **T2-2｜lymphoid −1.67 pp 不可复现，与 cover letter −1.6 矛盾。** 位置 L106 vs cover_letter L11。证据 C3-§2b：同 config（centCAB100i_a12, CP_nnls_norm）复算 lymphoid 中位 Δ=−1.599 pp（≈−1.60）；cover letter 写 −1.6 pp；layerA_deconv_summary.csv 无 lymphoid 行使该 p=0.003 溯源最弱。**修复：** 改稿件 lymphoid 为 −1.6 pp 与 CP_nnls_norm 及 cover letter 对齐；或显式定义"total lymphoid"包含哪些细胞型并给出产生 −1.67 的 config，并补 lymphoid Wilcoxon p 到 summary 作派生行。
- **T2-3｜免疫基因数 line-50 归因错置 + 48/75/79 三文件不一致未解释。** 位置 L50。证据 C3-§1b + **[编辑亲核]**：48/21/29 权威（来自 `bumphunter_adjusted_summary.csv`，Set-2 探针→基因准则）；overlap CSV 各 79/24/36 **行**（每行 1 基因）→ 行数非基因数；`bumphunter_immune_gene_before_after.csv`=48 基因；`layerA_DMR_immune_gene_counts.csv`=75 行。**稿件数字本身均正确**，缺陷是 (a) line 50 把 48/21/29 错误归因于 overlap CSV，(b) 48/75/79 三文件差异无任何解释。**修复：** 删除 line 50 括号"(and 48/21/29 immune genes)"或改为精确表述（overlap CSV 含 79/24/36 免疫重叠 DMR，每行列其重叠基因；distinct 免疫基因数在 overlap CSV 为 79/24/36，而 `bumphunter_adjusted_summary.csv` 在更严格的 Set-2 探针→基因准则下报 48/21/29）；并解释 48/75/79 分歧来自不同探针→基因映射准则——正文取一种定义、脚注其余。
- **T2-4｜"no compound cleared FDR at the drug level" 与 ligand 臂 24 命中自相矛盾。** 位置 L135/153。证据 C5-Item9 + **[编辑亲核]**：`layerD_immune_summary.json` ligand `n_druglevel_fdr05`=24（含 6 IFNγ），chem=0。稿件既在 L110 报告 ligand 24，又在讨论中称"无化合物通过 FDR"。**修复：** 改为"无**小分子（chemical 臂）**化合物在药物水平通过 FDR（两 signature 均 chem=0）；cytokine ligand 臂则产生 24 个 FDR 显著扰动（含 6 IFNγ），但这些是内源信号分子而非可重定位药物，按通路级而非化合物级解读；虚拟敲除臂仍属探索性"。
- **T2-5｜IFNγ 可行动推断宜降级；区分 II 型(IFNγ/IFNGR) 与 I 型(IFNAR)。** 位置 L135。证据 C5-Item5：该预测几乎完全来自输入 signature 的 ISG 组成（IFNG/CXCL10/IRF7/IRF1/STAT1），LINCS 仅提供同义印证而非独立证据。**修复：** 将"would be predicted to exacerbate"弱化为"is consistent with exacerbation, pending experimental validation"，显式说明预测源自输入 ISG signature；区分 IFNγ(II 型)/IFNAR(I 型) 通路（IFNAR 拮抗针对 I 型，JAK–STAT 阻断覆盖两者）。
- **T2-6｜Results 小标题 "chemical arm valid" 过度承诺。** 位置 L108。证据 C5-A：正文随即将糖皮质激素降级为"marginal directional signal"、他汀/NSAID 标失败、chem 臂 0 药物水平命中。**修复：** 改为"chemical arm partially supportive / hypothesis-generating; virtual-knockout arm exploratory"。
- **T2-7｜缺失神经炎症教条文献（机制层）。** 位置 Background（L25–30 一带）。证据 C1-F8：Axis-2 所依"神经炎症/外周免疫激活"的奠基文献阙如——Cerejeira 2010 (*Acta Neuropathol* 119:737, PMID 20490615)、Cunningham & Maclullich 2013 (*Brain Behav Immun* 28:1, PMID 23088936)、Inouye 2014 (*Lancet* 383:911, PMID 24560607)、Terrando 围术期 HMGB1/DAMP 小鼠工作。**修复：** 在引入神经炎症处至少加 Cerejeira 2010、Cunningham 2013、Inouye 2014；在外周→中枢通路处加一篇 Terrando HMGB1。
- **T2-8｜"peripheral immune activation state axis" 应明确定义为整合构念。** 位置 Background/Conclusions。证据 C1-F10：该"轴"由三层弱且部分非独立信号拼成（方向性转录组 MWU p=0.0079 + 组成污染表观信号 + 同义 LINCS IFNγ），非单一测量端点。**修复：** 加一句"我们使用的'外周免疫激活状态轴'是一个工作性、整合性构念，跨越方向性血液转录组免疫上调、大体由白细胞组成关联的围术期甲基化偏移、以及 in-silico 免疫 signature 可逆性——非单一测量端点；各成分均属假设生成"。
- **T2-9｜"candidate, modifiable therapeutic target" 措辞强于所引证据。** 位置 L153 Conclusions。证据 C1-F11：唯一干预文献为混合至阴性（地塞米松 meta 阴性、瑞舒伐他汀阳性但 MoDUS 阴性、LINCS 无药物级 FDR 命中、KO 全转录本阳性对照失败）。**修复：** 改为"a candidate therapeutic target whose modifiability remains unproven and requires randomized interventional testing"。
- **T2-10｜[35] "therapeutically actionable" 略夸大。** 位置 L55 一带。证据 C1-F4：Cheng 2022 自身结论留"more clinical evidence urgently needed"。**修复：** 软化并加"although the cited review notes clinical confirmation remains incomplete"。
- **T2-11｜伦理声明缺原始 IRB 名称/编号或明确豁免依据。** 位置 L163 Declarations。证据 C4-§3a：BMC 要求命名原始伦理委员会及参考号（如适当）。**修复：** 具名 GSE330869 与 Armstrong GWAS（UKB/Bristol）的批准委员会及编号，或显式写"Not applicable — secondary analysis of public de-identified data; IRB approval waived per [institution]/[journal] policy"并给法律依据。
- **T2-12｜KO 全转录本 class test 的 post-hoc 下尾 p 须标 exploratory + 重述锚点循环性。** 位置 L112。证据 C2-Check9：pre-specified 上尾 p=0.99 正确报为失败阳性对照；但 post-hoc 下尾 p≈0.0066/双侧≈0.013 系看数据后挑选，须明确"探索性、非预设"，并在 z=−2.47 旁重述 KO class 含 4 个锚点基因故部分循环。**修复：** 标 "post-hoc, exploratory, not pre-specified"；在 z=−2.47 旁加循环性说明；保留"failed positive control / mimicry"框架。
- **T2-13｜bumphunter B=0 须每处标注 "(descriptive)"。** 位置 L102/106/140 等所有 "13,357 DMRs" 及组成削减%出现处。证据 C2-Check5：`layerA_bumphunter_adjusted.R` L92 确认 `B=0`（无置换 FDR）；组成削减%正派生自这些无控计数。**修复：** 在每一处 "13,357 DMRs" 与削减%后附 "(descriptive, B=0, no permutation FDR)"；加一句"the 46–61% reduction is therefore a descriptive proportion, not an inferential attribution"。
- **T2-14｜Layer B "one-sample t p=0.18 as a lower bound on the directional signal" 措辞误导。** 位置 L92。证据 C2-Check6：p 值不是效应量的界；p=0.18 比 p=0.0079 更不显著，称"lower bound"倒置直觉。**修复：** 改为"immune-gene t 分布相对背景上移（MWU p=0.0079），但平均单基因偏移不显著（one-sample t p=0.18），因个体免疫基因在 n=8 时效应小且异质——信号是 broad and directional 而非由少数大效应基因驱动"；删 "lower bound"。

### Tier 3 — 轻微/排版（建议一并处理）
- **T3-1** [3] Hsiao 缺 PMID **36827351**（C1-F3）。
- **T3-2** [15][33][39][40] 缺 PMID/DOI（C4-§5c）。
- **T3-3** 糖皮质激素双侧 p 0.080→**0.0814**（C3-§6a，折叠进 T2 数值，此处仅注记舍入）。
- **T3-4** 把 \|t\| 富集的 MWU p 加进 `adjust_summary.csv` 或显引 `adjust_report.txt`（C3-§3）。
- **T3-5** `pdb_audit.json` 与 `layerE_dockability_audit.csv` 的 TMIGD3/SPATA13 UniProt ID 不一致（C3-§7，非稿件缺陷，存档卫生）。
- **T3-6** 稿件 v1.1 vs Zenodo v1.0.0 版本对齐（C4-§3b）。
- **T3-7** Abstract 缩写最小化（C4-§2）。
- **T3-8** README/preregistration 中 "two-axis **model**" vs 稿件 "hypothesis" 术语统一（C4-§6）。
- **T3-9** 按 Editorial Manager 路线以 .docx/.tex 提交（C4-§7）。
- **T3-10** **澄清被引 `layerA_deconv_summary.csv` 路径**（根文件 vs `analysis/results/GSE330869/` 副本 Neu=NA 不同 schema）——C2 跨切面与 C3 均涉及；C3 已用根文件复现中性粒/单核，故数字可溯源，仅路径歧义需消。**降自 T2**（中性粒数字已证可复现）。
- **T3-11** 首次出现各 p 值时加 "signed-t (directional)" vs "\|t\| (magnitude)" 比较脚注（C2-Check4）。
- **T3-12** 在 Results 显列 GWAS 全常染 p=0.081（"weak, chr19-driven"）并加有效 N 调和句（C2-Check7）。
- **T3-13** "nominal FDR"→"BH-adjusted FDR"（C5-D）。
- **T3-14** KO 对照 `n_controls_found`=18 vs 指定 12（子串匹配）脚注说明（C5-E）。
- **T3-15** "5208/5210 empty cell_line" 对齐产物实际 5206/5212（C5-Item6）。
- **T3-16** Methods 显式区分 reversal-only FDR（基于 p_combined）与 net_score 排序、z/emp_p 针对 reversal_score（C5-Item1/G）。
- **T3-17** 显式化单侧 p=0.041 的边界脆弱性（C5-Item3）。
- **T3-18** [22] 方法注：delirium 用 MMSE + 精神科访谈（非 CAM-ICU）（C1-F7）。
- **T3-19** 确认 [23]/[24] 作者/年/刊字符串（C1-F6）。
- **T3-20** [41] F1 发病率 10–50% 与 POD/POCD/PND 区分成立，保留（C1-F1）。

---

## 5. 共识 / 互补 / 分歧

### 5.1 共识（多专家独立命中）
- **极度诚实的阴性/失败阳性对照报告**是稿件最大优点（C1、C2、C3、C5 均独立盛赞）：他汀/NSAID 失败、KO 全转录本方向失败、组成污染明示、GWAS 不可检验解离。
- **[41] Evered 引用错误**被 C1/C3/C4 共同指向（C2 未专门查文献但 C1 详证）→ 确证 T1。
- **Layer A 组成校正口径/归因问题**被 C2（covariate 混用）与 C3（line-50 归因）从两个角度共同确认。
- **JNI 期刊指标声明真实**（JCR 2025 JIF 11.5, Q1, SCIE, 全 OA，范围契合）被 C4 独立网页核实为真。

### 5.2 互补（各专家覆盖不同层，互不重叠地补强）
- C1 补**机制文献缺口**与**临床措辞过强**（modifiable target / therapeutically actionable）。
- C2 补**设计层不确定性量化**缺陷（组成区间非灵敏度、Set 不可比、Layer-C 算术、B=0 标注、KO post-hoc）。
- C3 补**溯源级数字复现**（除 line-50 归因与 lymphoid 外几乎全部复现）。
- C4 补**政策级缺口**（generative-AI 声明、伦理具名、版本对齐、figure 上传）。
- C5 补**重定位诚实性**（"no compound cleared FDR" 自相矛盾、IFNγ 同义反复、subheading 过度承诺）。

### 5.3 分歧与编辑裁定
- **C4 称 [6] Seki 2026 "查无此文"（T1）→ 编辑裁定：撤销（false alarm）。** 理由：编辑 WebSearch + 本项目 G1 门 E-utilities 核实均证 PMID 42143058 / DOI 10.1038/s41398-026-04067-6 真实存在；C4 的 EuropePMC hitCount 0 系 2026 新文索引滞后。稿件 [6] **无需修改**（仍建议作者投稿前再点一次 PubMed 确认，属 T3-verify 而非必改）。
- **C3 称 "48/21/29 与 overlap CSV 三 artefact 矛盾"（暗示 T1）→ 编辑裁定：降级为 T2。** 理由：编辑亲核显示 48/21/29 权威（来自 `bumphunter_adjusted_summary.csv`），79/24/36 是 overlap CSV 的 **DMR 行数**（每行恰好 1 基因），两者口径不同、均正确。稿件数字无事实错误；真实缺陷仅 (a) line-50 把基因数错误归因于 overlap CSV、(b) 48/75/79 三文件分歧未解释。故属 T2 溯源澄清，非 T1 事实矛盾。
- **C3 verdict "Minor revision" vs C1/C2 "Major revision" → 编辑采纳更严格定级（Major revision）。** 理由（skill §7）：C3 自限范围为"数字复现"，其宽松 verdict 认证的是"溯源层"而非整稿；设计层（C2 T1-3/T1-4）与领域层（C1 T1-1）的真实问题不因溯源层干净而消失。
- **Layer-C 1/3 算术**：仅 C2 显式复核并给出 1/6；C3 未复核该项但亦指 lymphoid 等其他不一致。编辑亲核确认 1/6，维持 T2-1。

---

## 6. 优先级 must-fix 清单（含 DESK-REJECT 旗标 + 增补分析 vs 仅改写 拆分）

**DESK-REJECT 旗标：无。** 所有发现均可在修改内解决，无撤稿级。

**优先级顺序（必改 → 应改 → 可选）：**

| 优先级 | 项 | 类型 | 动作 |
|---|---|---|---|
| **P0（接受前硬门槛）** | T1-1 [41] Evered | 仅改写（换参考文献条目） | 替换正确 PMID 30325806 + 正确作者 |
| **P0** | T1-2 generative-AI 声明 | 仅改写（加一句） | Methods + cover letter 加声明 |
| **P0** | T1-3 组成"46–61%"重框 | 改写为主 + 可选增补分析 | 改写为描述性上限；**建议**用统一 3-lineage 协变量重跑 DMR 层（非强制，但能夯实） |
| **P0** | T1-4 Set-1/Set-2 不可比对冲 | 仅改写 + 加表 | 加行内说明 + 四行 p 值总表 |
| **P1** | T2-4 "no compound cleared FDR" | 仅改写 | 区分 chem(0) vs ligand(24) |
| **P1** | T2-1 Layer-C 1/3→1/6 | 仅改写 | 改数字 |
| **P1** | T2-2 lymphoid −1.67→−1.60 | 仅改写 | 对齐 cover letter |
| **P1** | T2-3 line-50 归因 + 48/75/79 解释 | 仅改写 | 改归因 + 脚注定义分歧 |
| **P1** | T2-7 缺神经炎症教条文献 | 增补文献（非分析） | 加 3–4 篇引用 |
| **P1** | T2-11 伦理具名/豁免 | 仅改写 | 补声明 |
| **P1** | T2-12 KO post-hoc 标 exploratory | 仅改写 | 加标注 + 循环性重述 |
| **P1** | T2-13 B=0 每处标注 | 仅改写 | 加 "(descriptive, B=0)" |
| **P2** | T2-5 IFNγ 推断降级 | 仅改写 | 软化 + 通路区分 |
| **P2** | T2-6 subheading 过度承诺 | 仅改写 | 改小标题 |
| **P2** | T2-8 "state axis" 定义 | 仅改写 | 加一句 |
| **P2** | T2-9 "modifiable target" 软化 | 仅改写 | 改结论 |
| **P2** | T2-10 [35] 夸大 | 仅改写 | 加限定 |
| **P2** | T2-14 "lower bound" 误导 | 仅改写 | 改措辞 |
| **P3** | 全部 T3-1…T3-20 | 仅排版/补全 | 一并扫尾 |

**"必须增补分析" vs "必须改写" 拆分：** 本稿**无需补新实验**；唯一可选的增补分析是 T1-3 建议的"用统一 3-lineage 协变量重跑 DMR 层组成校正"（能直接回答 C2 的核心质疑，强烈建议但非硬门槛）。其余全部为**改写/补全引用/补声明**，可在单轮修订内完成。

---

## 7.  站得住脚之处（告诉作者**不要**改什么）

（自各专家 § Stands up 结转；编辑确认这些结论不受本轮任何 T1/T2 影响）
1. **DMR 总数与免疫重叠 DMR 数（13,357/5,153/7,263 与 79/24/36）**精确复现，组成削减叙述的骨架正确（C3-§1a）。
2. **去卷积中性粒 +3.38 pp (p=7.6e-4) 与单核 −1.11 pp (p=1.2e-4)** 从根 summary 与原始分数均精确复现——这是稿件的阳性对照，扎实（C3-§2a）。
3. **GWAS 层内部无瑕**：37 SNP 全 chr19、λ=1.0133、N=134,310、BETA/OR 符号逻辑自洽（OR=exp(0.621)=1.86/ε4 等位）（C2-Check7、C3-§5）。
4. **LINCS class 检验与 IFNγ/KO 数字**全部匹配产物（糖皮质激素 0.041、6 IFNγ FDR、KO 1,109 vs 0、顶基因 FDR 精确）（C3-§6、C5）。
5. **Layer B 转录组数字**（737/996、MWU 0.0079、mean t 0.253/−0.171、83 免疫基因）复现（C3-§8）。
6. **Set-1 signed-t (0.041) 与 \|t\| 富集 (2.5e-4) 正确区分为两个不同检验**，无混淆（C2-Check4）。
7. **JNI 指标声明真实且网页核实**（JCR 2025 JIF 11.5, Q1, SCIE, 全 OA，范围契合）（C4-§1）。
8. **Data availability 真实匹配**：Zenodo DOI 10.5281/zenodo.22896443 已发布、标题/作者/ORCID/MIT 一致（C4-§3b）。
9. **诚实的 Declarations**：单作者 B.M. 学位声明无虚构、IRB 豁免立场、无基金、无竞争利益（C4-§3c/d）。
10. **对接层严格限定为"可行性审计"**：ADORA3=5 PDB、TMIGD3/COL13A1/SPATA13=0、无结构者列局限，无"有结构即暗示可成药"越界（C5-§8）。
11. **腺苷受体异质性处理典范**：[20]/[23]/[24] 分属 A3/A2A/(A1+A2A)，稿件正确拒绝把异质子文献包装为对 ADORA3 的统一背书（C1-F6、C5-§8）。
12. **临床 RCT 限定精确**：MoDUS 正确排除于外科 POD 推论；地塞米松 null/瑞舒伐他汀 positive/Mardani positive 方向与 n 均正确（C1-F7）。

---

## 8.  推荐处理路径（A / B / C）

- **路径 A（推荐）：同文章类型（JNI Research article）修改后重投。** 所有 T1 项（T1-1 换文献、T1-2 加声明、T1-3 重框、T1-4 加表）与 P1/P2 改写可在**单轮修订**内完成；唯一可选增补是 T1-3 的统一协变量 DMR 重跑（非强制）。科学内容、阴性报告框架、GWAS 与去卷积阳性对照均无需动摇。**本稿适合 JNI，无需降型。**
- **路径 B（降型）：不适用。** 本稿多组学整合 + GWAS + LINCS 的方法学深度与 JNI 范围契合；设计层问题经改写可弥合，无降为 Brief Communication / Correspondence 的必要。
- **路径 C（仅改写、不补分析）：不足以单独成立。** 表面看多数项为改写，但 T1-3 若**只改写不补分析**，会留下"组成特异性结论仍由探针层 \|t\| 检验单方承载"的未决质疑；建议至少跑一次统一 3-lineage 协变量的 DMR 层校正以坐实。故 C 不能替代 A 中的可选增补。

**编辑给作者的直接建议（reframe, don't just criticise）：**
- **最站不住却非你之过的强论断 = DMR 层"46–61% 组成可归因"所暗示的"免疫特异性"结论。** 真相是免疫 DMR 占比在调整前后基本持平（0.59%→0.47–0.50%），所以"免疫特异性组成"实际由探针层 \|t\| 富集检验（p=2.5e-4→0.72）承载，而非 DMR 计数。这不是你的计算错误，而是描述性 B=0 DMR + 两层协变量口径错位导致的**表述性过伸**。
- **你被埋没的真正贡献 = (a) 围术期白细胞偏移作为可复现阳性对照（中性粒 +3.4 pp，多参考一致）、(b) APOE ε4 单态（GWAS 全 chr19、无 MR 可行、off-chr19 p=0.277）、(c) 一套罕见诚实的失败阳性对照/假设生成框架。** 把头条从"免疫甲基化是状态轴证据"换成"外周白细胞组成偏移是真实可复现的围术期信号；免疫表观信号为假设生成且大体可归因于细胞组成，非细胞内在表观重编程"，能把一篇脆弱信号稿转为稳定的结构/报告学发现稿。

---

## 9.  过程教训（自动化 gate 抓不到的，以及可扩展 gate 覆盖的设计层）

- **gate 全绿 ≠ 设计无瑕。** 本稿 Round-7 通过全部算术/溯源 gate（DMR 13,357 复现、数字溯源审计表齐备），但面板仍在设计层抓出：组成"区间"非灵敏度包络（T1-3）、Set-1/Set-2 不可比却并置（T1-4）、"lower bound" 误用 p 值（T2-14）、Layer-C 1/3 小样本精确 p 错（T2-1）。这些是 gates 不枚举的**设计层偏差**。
- **可新增的 gate 断言（扩展覆盖到设计层）：**
  1. **"声明↔源文件" 溯源映射 gate**：扫描稿件每个数字及其"见 X 文件"的归因，自动打开 X 文件验证该数字确实存在于 X（而非另一文件）。本 gate 会直接抓出 T2-3（line-50 把 48/21/29 归因于 overlap CSV，但该文件只含 79/24/36 行）与 T3-10（中性粒被引路径歧义）。
  2. **"相邻两个不同检验" 检测器**：当稿件在相邻句出现两个同对象但不同 t 定义/基因集的 p 值（如 signed-t 0.041 与 \|t\| 2.5e-4，或 Set-1 与 Set-2），强制要求显式"二者不可比"声明 + 对照表。本 gate 对应 T1-4。
  3. **小样本精确 p 校验器**：对所有 n≤4 的比较（如 Layer-C n1=n2=2），用枚举法重算精确 MWU 最小 p，对照稿件所写地板值。本 gate 会抓出 1/3→1/6（T2-1）。
  4. **"区间 vs 包络" 措辞 gate**：当稿件出现 "X–Y%" 形式且两端来自不同协变量/模型规格时，标记"是否误作灵敏度区间"。对应 T1-3。
- **独立性纪律的价值再次验证**：C4 的 [6] "查无" 与 C3 的 "48/21/29 矛盾" 两处 false alarm，均由编辑亲核（WebSearch / 读原始 CSV）及时撤销，避免作者浪费一轮去"修正"本不存在的缺陷——这正是 skill §5 "verify the worst claim yourself, do not forward it" 的意义。

---

*本报告由编辑合并五份独立专家报告（C1–C5）而成，所有 T1 严重项与可能改变定级的发现均经编辑用原始产物亲自复核并标注 [编辑亲核]。五份专家原文保留于 `review_round8/C{1-5}_*.md`，作为本轮修改为何改、改哪里的审计轨迹。*
