# HANDOFF：AML 骨髓“样本内部分工”主线（v1）

| 项 | 内容 |
|---|---|
| 文件名 | `HANDOFF_within_sample_v1.md` |
| 日期 | 2026-10-02 |
| 状态 | 方向和优先级已定；数值阈值是建议值，待冻结（见 §15） |
| 上游文件 | `AML_niche_CCC_OT_blueprint_zh.md`（旧主线）、`ARCHITECTURE.md`、`CODING_STANDARDS.md`、`AML_niche_CCC_dataset_inventory.xlsx`、`AML_ecosystem_research_dossier.md`（下文称“资料”） |
| 与旧主线的关系 | 旧主线（per-sample CCC 图 → FGW 对齐 → 评分）的 H1–H3 均为阴性，另行处理，不在本文件范围。本文件定义替换它的新主线 |
| 读法 | 只想开工：读 §0、§4、§16、§17。要写代码：再读对应 WP 的小节和附录 A。要写文章：读 §11、§13 |

标记说明：`✔` 表示本轮对话里检索核实过题目、期刊和年份；`○` 表示来自资料或既有知识，没有单独核实，使用前请核对。

---

## 0. 一页摘要

**要回答的问题**：同一份 AML 骨髓里，哪些细胞在提供信号，哪些细胞因此受益；突变细胞和残余正常细胞面对同一批信号时，反应差在哪里；这种分工在治疗之后还在不在。

**为什么换成这个问题**：旧主线有两步容易出问题。一是按细胞类型汇总配受体打分，容易被细胞比例带着走。二是拿不同样本互相比，容易被平台差异带着走。新主线把比较放到同一个样本内部，证据改看“接收细胞有没有真的起反应”。

**四个工作包（按优先顺序）**

| 顺序 | 工作包 | 一句话 | 数据是否够 |
|---|---|---|---|
| 1 | WP1 样本内部的分工 | 提供者和受益者是谁，再叠加突变或正常 | 最充足：只需诊断期样本 |
| 2 | WP2 治疗后仍在的多细胞程序 | 扣掉细胞比例后，多类细胞一起变的程序，缓解期是否仍异常 | 发现阶段充足，纵向部分偏少且以儿童为主 |
| 3 | WP3 代偿与共享瓶颈 | 两个上游输入是否可以互相替代 | 同 WP1，但要等 WP1、WP2 先成立 |
| 4 | WP4 验证模块 | 外部单细胞、bulk 队列、公开扰动数据、后续实验 | 只能支持方向 |

**开工第一件事**：做 §16 的三个可行性检验，不要先写分析主流程。

---

## 1. 已定的前提（本轮对话确认）

1. 目标是一条能出正向结果的新主线，可以替换旧主线，数据可以重选。
2. 核心分析只用 scRNA-seq 或 snRNA-seq，加样本临床元数据。不用空间转录组、单细胞表观多组学、自建 CRISPR 实验。
3. 不讲“分化停滞”，不讲“某两类细胞通讯异常”这类结论。
4. 计算分析是主体。实验只做验证，可做的类型：细胞系敲基因或加药、细胞因子或阻断抗体、共培养、特定突变的小鼠模型。
5. 公开 bulk 队列（TCGA-LAML、BeatAML）可以用于验证。
6. 突变类型不设限，按公开数据量来定。
7. 排序标准：成功率第一，其次看能拿到的数据数量和质量。
8. 文章卖点不预设，方法或发现都可以，哪个容易成做哪个。
9. 配受体通讯放在“解释原因”的位置，除非某个做法能证明自己不受细胞比例影响。
10. 有 GPU，随时可用。

---

## 2. 前情：旧主线留下的教训

以下是项目已有结论，细节以项目内的 findings 文档为准。

- 按细胞类型汇总的 CCC 距离、OT 距离，没有一种超过“7 类细胞比例”（`prop7`）这个基线。
- AML 与健康骨髓在网络拓扑上分不开，这不是样本量的问题。
- 同一患者前后的变化是真实的，但变化方向在患者之间不一致，所以“诊断 → MRD → 复发单调加深”这类检验做不成。
- 平台效应会明显污染跨数据集的梯度。
- inferCNV 会漏判静息的白血病干细胞，MRD 期的恶性细胞判定不可靠。
- 骨髓穿刺样本里基本没有基质细胞，“niche”实际指造血和免疫细胞。
- 分选或 MACS 富集样本的细胞比例反映的是操作，不能当作真实丰度。
- 现有的纵向配对是 11 对，都是诊断 → 复发：9 对来自 GSE227903，2 对来自 GSE201966。

---

## 3. 设计原则

- **P1 比较在样本内部完成**，或在同一患者的前后之间完成。跨样本的比较只出现在第二阶段，并且以患者为单位。
- **P2 通讯的证据看接收细胞的反应**。配体和受体都表达不算证据。
- **P3 统计单位是患者**。细胞数不能增加独立样本量，细胞层面的重采样只用来估计测量稳定性。
- **P4 每个结论先过三个对照**：细胞比例基线（`prop7`）、表达水平匹配的随机基因集、打乱配受体配对的假配对。
- **P5 先做小检验，再投入**。每个 WP 有通过门槛，不通过就换做法或停。
- **P6 程序先在单个研究内部发现，再跨研究重复**。不把不同平台的样本混在一起找规律。

---

## 4. 工作包来源对照

| 工作包 | 来自首轮方案 | 来自资料 | 主要改动 |
|---|---|---|---|
| WP1a 基因型层 | 场景 A（同一骨髓里的天然实验） | 场景 C 的“基因型版本” | 主检验改为“同一层级内，突变与野生型细胞的信号反应差”；用旁观细胞测环境信号强度 |
| WP1b 角色层 | ContactTracing 逻辑 | 场景 C（同状态细胞的不同生态角色） | 不用 scComm 的打分；“角色”拆成提供者（配体表达）和受益者（受体阳性细胞出现反应）；用患者间配体多寡做检验 |
| WP2 | 场景 C（扣掉比例后的协同程序） | 场景 A（治疗后持续的多细胞程序） | 加环境 RNA 控制；“持续存在”改成两个能证明的命题；预先固定的模块不超过 4 个 |
| WP3 | — | 场景 B（蛋白与代谢信号的代偿） | 不做跨样本通讯张量；改成样本内部的冗余检验；MEBOCOST、SigXTalk 只提供候选名单 |
| WP4 | bulk 与实验验证 | 场景 D（受干预约束的反事实预测） | CellOT 不当主线；公开扰动数据只用来核对方向 |
| 暂缓 | 场景 B（化疗当扰动）、场景 D（增长率） | 场景 B 原做法、场景 D 原做法 | 见 §9 |

---

## 5. WP1 样本内部的分工

WP1 分两层。WP1b 用全部诊断期患者，不依赖突变读数。WP1a 只用能判定基因型的患者子集。两层共用同一套细胞注释和信号打分。

### 5.1 假设

- **H-W1a**：同一份骨髓、同一层级状态内，突变细胞对炎症类信号的转录反应弱于残余的野生型细胞；环境信号越强，这个差距越大。
- **H-W1b**：恶性细胞内部存在分工。少数细胞高表达某些配体（提供者），带对应受体的细胞出现可测的反应（受益者）；提供者多的患者里，这种反应更强。非恶性细胞作为提供者的情况同样纳入检验。
- **H-W1c**：提供者和受益者在基因型上不对称（例如提供者多为突变细胞，被压制的多为野生型细胞）。

### 5.2 数据要求

| 要求 | WP1a | WP1b |
|---|---|---|
| 时间点 | 诊断期即可，不需要纵向 | 诊断期即可 |
| 细胞范围 | 未分选或 CD34 富集的骨髓；同一样本里要同时有突变和野生型干祖细胞 | 未分选骨髓（要有恶性细胞和多类非恶性细胞） |
| 原始 reads | 必须有（带 CB/UB 标签的 BAM） | 不需要 |
| 临床信息 | NPM1 状态（或能从 reads 里直接看到突变） | blast%、基因型、ELN 分层（有则用） |
| 患者数 | 粗估 15–25 位 NPM1 突变患者，需实测 | 项目清单内诊断期全细胞样本约 105 行，按患者去重后再计 |

候选数据集见 §10。WP1a 的首选是有原始 reads 的成人诊断期数据：GSE185381、GSE239721、GSE289435、GSE227903、GSE116256。Chen 2023 的数据用作阳性对照。

### 5.3 共用的预处理

1. **doublet 和环境 RNA 先处理干净**。受体阳性与阴性细胞的比较对这两样特别敏感：一个混入了别的细胞的液滴，会同时带上“受体”和“反应基因”，看起来像真的反应。项目里 doublet 预处理的状态（`01_preprocess` 的 per-sample QC）需要先确认。
2. **细胞注释沿用已有产出**：恶性标签用 `02_malignancy` 的共识结果；层级状态用 `03_hierarchy` 的 BoneMarrowMap 投影分箱；非恶性细胞用投影得到的类型。
3. **信号活性打分（每个细胞）**，用“反应特征”来估计细胞收到的信号，不用配体表达：
   - CytoSig：细胞因子反应特征矩阵加岭回归。
   - PROGENy 通路（可经 decoupleR）：TNFa、NFkB、JAK-STAT、TGFb、Hypoxia 等。
   - Hallmark 的 IFN-γ、IFN-α、TNFα/NFκB 反应基因集打分，作为最朴素的对照版本。
   - 三套结果方向一致才算稳。打分基因与下面“适应度”基因不得重叠。
4. **适应度代理（每个细胞）**：增殖评分（S 与 G2M 基因）、凋亡与 p53 通路评分、干性评分（项目已有的 LSC17、HSPC_core 等）。哪个是主读数见 §15。

### 5.4 WP1a 步骤（基因型层）

1. **逐细胞读 NPM1 位点**：按附录 A，从 BAM 里统计每个细胞在 NPM1 突变位点的突变 UMI 数和野生型 UMI 数。
2. **估计假阳性率**：T、NK、B 细胞理论上不带 NPM1 突变。它们里面的突变 UMI 比例，就是这个样本的环境污染率。
3. **后验判定**：用二项混合模型给每个细胞算“是突变细胞”的后验概率，分三类：突变、野生型、未定。未检出突变不等于野生型，未定的细胞不参与比较。细则见附录 A。
4. **非 NPM1 患者**：用 Numbat 的克隆后验做同样的事，比较 CNV 克隆与正常细胞。注意项目的语言约定：CNV 信号不是点突变，称 “CNV-based proxy”。
5. **只在同一层级分箱内比较**。已有报道里，HSPC 组分中的 NPM1 突变细胞多见于 LMPP、GMP，静息的干细胞多为健康或白血病前细胞。两类细胞都有的分箱才有可比性，具体是哪些分箱要看实测。这一步同时保证比的不是分化阶段。
6. **检验 1（主检验）：反应差**

   ```
   # per patient p, per hierarchy bin b (bins with >= N_MIN cells of each genotype)
   score_s(c) ~ mut(c) + log10(nUMI(c)) + pct_mito(c) + phase(c) + library(c)
   Delta[p,b,s] = coefficient of mut          # mutant − WT, for signal s
   Delta[p,s]   = inverse-variance weighted mean over bins
   ```

   跨患者检验 `Delta[.,s]` 是否偏离 0：符号翻转置换或 Wilcoxon，按信号做 BH 校正。逻辑是：同一样本同一状态的细胞，暴露在同一环境里，反应特征的差异就是反应能力的差异。要注意，这个差异也可能来自突变本身带来的基础表达差别，与外界信号无关。检验 2 就是用来区分这两种情况的：差距随环境信号增强而变大，才支持“反应能力不同”。
7. **检验 2：环境剂量关系**

   ```
   E[p,s] = mean response score of signal s in bystander cells of patient p
            (non-malignant T/NK and mature monocytes; cells not used in Delta)
   Delta[p,s] ~ E[p,s]        # within study; unit = patient; Spearman + linear fit
   ```

   用旁观细胞当“环境探针”，环境信号的测量和结果的测量来自不同细胞，避免循环论证。
8. **检验 3：代价是否不同（Tres 式交互）**

   ```
   fitness(c) ~ signal_s(c) * mut(c) + covariates      # per patient, per bin
   ```

   记录每位患者的交互项 t 值，跨患者汇总。问的是：同样的反应强度下，突变细胞的增殖或存活受损是否更小。
9. **耐受基因**：只在突变细胞里做 `fitness ~ signal_s × gene_g`，每个基因得到一个交互 t 值，跨患者取中位数排序。排在前面的是候选的“耐受基因”，交给 WP4 验证。

### 5.5 WP1b 步骤（角色层）

1. **配体可得性**（每位患者、每个配体 L）：

   ```
   A[L,p] = pseudobulk expression of L in a sender compartment of patient p
            (malignant cells by bin; each non-malignant type)
   ```

   主分析用“区室内的平均表达”，不乘以该区室的细胞占比；乘以占比的版本只做敏感性分析，因为细胞比例就是从这里混进来的。
2. **少数提供者结构**：每位患者、每个配体，算恶性细胞里检出该配体的比例和 Gini 系数。再问两件事：检出的不均匀程度是否超过同等平均表达下的随机期望（负二项零模型）；回归掉层级分箱、细胞周期、nUMI 之后，配体阳性细胞是否仍带有跨患者可重复的共表达程序。对照是表达水平匹配的随机基因。
3. **接收细胞的反应（第一阶段，每位患者内部）**：

   ```
   # receiver type T in patient p, receptor R
   Y_g(c) ~ CDR(c) + Rpos(c) [+ bin(c) + phase(c)]     # MAST hurdle model
   -> beta[g,p], se[g,p]
   ```

   `Rpos` 是受体是否检出，`CDR` 是细胞检出基因的比例。这一步只说明受体阳性和阴性细胞不同，还不能说是配体造成的。
4. **跨患者交互（第二阶段，单位是患者）**：

   ```
   beta[g,p] ~ A[L,p] + covariates(p)       # weights 1/se^2; within study
   -> b[g]; permutation p-value by shuffling A[L,·] across patients
   ```

   `b[g]` 显著的基因，就是“配体多的患者里才出现的受体依赖反应”，称为配体效应基因。每个（配体、受体、接收细胞类型）组合的统计量是配体效应基因的个数，和置换零分布比。这是 ContactTracing 的做法，只是把“实验条件”换成了“患者间配体的多寡”。原文的参考阈值：交互效应基因不少于 10 个（FDR < 0.25），配体差异 |log2FC| > 0.12（FDR < 0.05）。
5. **是否受益**：配体效应基因里有没有增殖、存活相关程序；`fitness ~ Rpos × A[L,p]` 的交互是否为正。
6. **叠加基因型**（WP1a 的子集）：提供者是突变还是野生型，受益者是突变还是野生型，列成 2×2 表，按患者汇总。
7. **跨研究**：每个研究内部先出结果，再做随机效应 meta 分析。

### 5.6 对照与反证

- **假配对**：把配体和受体的配对关系在表达水平相近的基因间打乱，重跑第二阶段。真实配对的显著数必须明显多于假配对。
- **混杂识别**：受体阳性与阴性细胞的差异如果在所有患者里都一样、和配体多寡无关，那是受体标记了某个亚群，不是配体效应。ContactTracing 原文也用这一点区分两者。
- **细胞比例**：把 `prop7` 和 blast% 放进第二阶段的协变量，结论不应消失。
- **污染**：配体效应基因如果主要是发送细胞的标志基因，按环境 RNA 或 doublet 处理，不算反应。
- **深度**：受体阳性细胞通常 UMI 更多。`CDR` 必须在模型里，并按 nUMI 分层重复一次。
- **WP1a 特有**：把基因型标签在同一分箱内随机打乱，`Delta` 应回到 0。野生型判定阈值放宽和收紧各跑一次，方向不变才算稳。

### 5.7 通过门槛（建议值，待冻结）

| 门槛 | 指标 | 建议值 | 不通过时 |
|---|---|---|---|
| GATE_W1.1 基因型可判定 | 每位患者在至少一个共有分箱里的高置信突变细胞数和野生型细胞数；T/NK/B 细胞的突变 UMI 比例 | 不少于 12 位患者各有 ≥30 个突变和 ≥30 个野生型细胞；污染率 ≤2% | WP1a 改用 CNV 克隆，或只保留 WP1b |
| GATE_W1.2 阳性对照 | 在 Chen 2023 数据上重复原文方向：残余正常干祖细胞的炎症特征高于健康骨髓里的对应细胞 | 方向一致 | 先查流程 |
| GATE_W1.3 检出力 | 通过预筛的（配体、受体、接收类型）组合数 | 至少 2 个研究里各有 ≥30 个组合；每个组合在 ≥10 位患者里受体阳性和阴性细胞各 ≥50 个 | 合并相近细胞类型，或只做恶性细胞作为接收方 |
| GATE_W1.4 假阳性 | 假配对下的显著组合比例 | 经验 FDR ≤ 0.2 | 收紧预筛，检查污染 |

### 5.8 交付物

- `npm1_genotype_percell.csv`：每个细胞的突变 UMI、野生型 UMI、后验概率、三分类标签。
- `signal_scores_percell`：三套信号活性打分和适应度打分。
- `w1a_delta_by_patient.csv`、`w1a_dose_response.csv`、`w1a_resilience_genes.csv`。
- `w1b_ligand_availability.csv`、`w1b_stage1_beta/`、`w1b_stage2_ligand_effects.csv`、`w1b_provider_structure.csv`。
- 每个门槛一份一页的结果说明（通过或不通过，附数字）。

### 5.9 已有工作与新意

- Chen 等 2023（Blood Cancer Discov）在 NPM1 突变 AML 里提出：炎症性 niche 重塑压制残余正常干细胞亚群，克隆细胞相对耐受。项目清单里这套数据是 6 例患者。它的判定规则是读到至少 1 条突变 read 算突变，至少 3 条野生型 reads 算可能的野生型。
- Jakobsen 等 2024（Cell Stem Cell）在克隆性造血里看到同样的现象：DNMT3A、TET2 突变的干细胞比同一样本里的野生型干细胞炎症反应更弱。
- CloneTracer（Beneyto-Calabuig 等 2023）在 19 例 AML 里说明静息干细胞多为健康或白血病前细胞。
- 所以“突变细胞更耐受炎症”这个概念不新。新的部分有四点：逐个信号定量；跨多个公开队列；用旁观细胞测环境剂量；加上“谁在提供信号”这一层。写作时要把这条边界讲清楚。

### 5.10 后续实验（只做验证）

- 细胞因子刺激突变与野生型干祖细胞，看反应差是否重现。
- 共培养加阻断抗体，针对 WP1b 排名靠前的配受体。
- 敲低候选耐受基因后，看细胞在炎症刺激下的增殖或存活。
- 突变小鼠的竞争移植加炎症刺激。

---

## 6. WP2 治疗后仍在的多细胞程序

### 6.1 假设

- **H-W2a**：扣掉细胞比例之后，AML 骨髓里存在跨细胞类型一起变的表达程序，并且与基因型或预后有关。
- **H-W2b**：缓解期或 MRD 期的骨髓，这些程序的得分仍高于健康供者。
- **H-W2c**：同一患者诊断期和缓解期的程序得分相关，缓解期得分高的患者更容易复发。

### 6.2 数据要求

| 用途 | 要求 | 候选 |
|---|---|---|
| 发现 | 多患者、未分选骨髓、诊断期；每个研究内部患者数尽量多 | GSE185381、GSE239721、GSE289435、GSE227903、GSE116256；AML scAtlas 里追溯到的原研究 |
| 缓解期对比 | 治疗后样本，保留多类非恶性细胞；同研究里有健康供者更好 | GSE227903 的 7 个 MRD 样本；GSE116256 的 19 个疗后样本（细胞数少）；Mumme 2023（GSE235923，儿童）；Lambo 2023（GSE235063，儿童） |
| 结局 | 复发或持续缓解的标注 | Mumme 2023 |
| 配对确认 | 诊断 → 复发 11 对 | GSE227903（9 对）、GSE201966（2 对） |
| 健康参照 | 健康骨髓 | 各 AML 研究内的健康供者、E-MTAB-11536 |

儿童和成人分开分析，骨髓和外周血分开分析。

### 6.3 步骤

1. **环境 RNA**：用未过滤的液滴矩阵估计每个样本的环境表达谱和污染比例（SoupX、DecontX 或 CellBender 任选其一并固定）。记录每个样本的污染比例，列出每个样本环境谱里排名靠前的基因，作为后面的黑名单。原始细胞多的样本，blast 的 RNA 会进入所有其他细胞，正好伪装成“多类细胞一起变”。
2. **分组**：恶性细胞按层级分箱；非恶性细胞按类型。每个样本每组不少于 30 个细胞才计入，某组要在七成以上样本里存在才保留（两个数都是建议值）。
3. **状态控制的 pseudobulk**：每个样本、每个组各自汇总成一个表达谱。这样“某种细胞变多”和“同一种细胞内部表达变了”被分开。这一步借自 GBM CARE 联盟的做法。
4. **研究内部发现程序**：scITD（供者 × 基因 × 细胞组的张量做 Tucker 分解）或 DIALOGUE（跨细胞类型的稀疏典型相关加多层模型）。因子个数用“对供者再抽样后的稳定性”来定。
5. **去混杂**：每个因子与 blast%、污染比例、细胞比例、nUMI、性别、年龄逐一检验。主要被 blast% 或污染比例解释的因子标记出来，不进入后续结论。因子的高载荷基因与环境黑名单的重叠不应超过随机水平。
6. **跨研究重复**：把一个研究里得到的载荷投到另一个研究，算得分，再看两边的载荷相关性和与表型的关联方向是否一致。
7. **是否超出细胞比例**：嵌套模型比较，`表型 ~ prop7` 对比 `表型 ~ prop7 + 程序得分`，按患者做交叉验证和 bootstrap。
8. **纵向**：
   - 缓解期或 MRD 期样本对比健康供者（单位是患者）。
   - 同一患者诊断期与缓解期得分的相关。
   - 11 对诊断 → 复发：只检验预先固定的模块，不超过 4 个。
   - Mumme 数据里，诱导结束时的得分对比复发与持续缓解两组。
9. **解释层**：某个程序里，A 类细胞的配体与 B 类细胞的程序部分相关时，回到 WP1b 的两阶段检验去确认受体依赖的反应。
10. **bulk 外推**：取程序里各细胞类型特异的基因组成签名，在 BeatAML、TCGA-LAML 里打分，看预后。bulk 分不开细胞类型，这一步只当旁证。

### 6.4 关于 11 对的检验力

11 对做符号翻转置换，共 2^11 = 2048 种排列。按最朴素的符号检验：11 对全部同向，双侧 p ≈ 0.001；10 对同向，p ≈ 0.012。做 Bonferroni 校正并容忍一位患者不一致，模块最多 4 个。模块必须在看这 11 对之前固定下来。

### 6.5 对照与反证

- 程序效应在控制 blast%、细胞比例、取样来源、处理流程之后消失，就不能说发现了独立的生态机制。
- 只有两个时间点，只能说“仍然异常”或“发生了转变”，不能说吸引子、迟滞、稳定反馈。
- “持续存在”不靠“诊断和缓解没有差别”来证明，那等于证明没变，11 对做不到。

### 6.6 通过门槛（建议值，待冻结）

| 门槛 | 指标 | 建议值 |
|---|---|---|
| GATE_W2.1 环境 RNA | 因子得分与污染比例的相关；高载荷基因与环境黑名单的重叠 | 校正后 \|ρ\| < 0.3；重叠不超过随机期望 |
| GATE_W2.2 超出比例 | 加入程序得分后的交叉验证增益 | 患者 bootstrap 的置信区间不含 0 |
| GATE_W2.3 可重复 | 跨研究的载荷相关 | 至少 2 个研究间为正且显著 |

### 6.7 已有工作与新意

- AML 里已有按细胞状态“共存模式”定义生态系统的工作（ACE，Brief Bioinform 2025，6 套数据、68 例患者），并在 bulk 队列里验证了预后意义。所以只讲“生态型”或“细胞共现”不新。
- 多发性骨髓瘤里也有按细胞亚群比例的 NMF 定义 TME 生态型的工作。
- 本工作包的角度只能是：扣掉比例之后的表达协同，以及它在缓解期是否仍异常。
- GBM CARE 的第二篇给出了参照写法：配对样本没有统一轨迹，但部分患者共享特定轨迹。

---

## 7. WP3 代偿与共享瓶颈

只有 WP1 或 WP2 已经给出一个站得住的接收细胞程序之后，才启动。

### 7.1 假设

同一个存活或耐受程序可以由多个上游输入维持。阻断其中一个，另一个可以顶上。这些输入在细胞内汇到同一批转录因子上，这批因子就是共享瓶颈。

### 7.2 做法

1. **候选输入名单**：从三处取，只当名单用，不用它们的打分。
   - 蛋白配受体：WP1b 第二阶段显著的组合。
   - 受体到转录因子到靶基因的汇聚：SigXTalk 的思路。
   - 代谢物与 sensor：MEBOCOST 的数据库。
2. **样本内部的冗余检验**（每位患者、每类接收细胞）：

   ```
   program(c) ~ R1pos(c) + R2pos(c) + R1pos(c):R2pos(c) + CDR(c) + bin(c)
   ```

   两个主效应都为正、交互项为负，说明两个输入对同一程序有饱和或互相替代的关系。每位患者得到一个交互估计，跨患者检验。
3. **只在 R1 阴性细胞里看 R2**：没有 R1 的细胞里，R2 的表达是否仍然伴随同一程序。这是代偿最直接的读法。
4. **共享瓶颈**：用转录因子活性打分（decoupleR 加 CollecTRI 或 DoRothEA），找同时与 R1、R2 的反应相关的因子。
5. **程序打分的基因与 R1、R2 及其已知下游不得重叠**，否则是循环使用同一组基因。

### 7.3 反证

- 关联如果只来自总体应激、同一患者的 blast 数或批次，拒绝机制解释。
- 代谢物是用酶的 mRNA 推出来的，只能叫候选，不能当作浓度或流量。
- 样本里基本没有基质细胞，来自基质的输入看不到。

### 7.4 验证

这一包最适合交给实验收尾：单独阻断 R1、单独阻断 R2、同时阻断，比较接收细胞程序和存活。你们做代谢实验有基础，代谢物一侧的候选优先留给实验。

---

## 8. WP4 验证模块

| 编号 | 内容 | 数据 | 能支持到什么程度 |
|---|---|---|---|
| V1 | 外部单细胞重复 | 发现阶段没用过的研究；AML scAtlas 追溯到的原研究；儿童数据集 | 结论是否跨队列成立 |
| V2 | bulk 队列 | BeatAML（预后、离体药敏）、TCGA-LAML（预后） | 签名与结局相关；bulk 分不开细胞类型 |
| V3 | 公开扰动数据 | CytoSig 数据库里的细胞因子处理表达谱；GSE146590（阿糖胞苷） | 配体效应基因的方向是否与真实处理一致 |
| V4 | 配对方向 | 11 对诊断 → 复发；MRD 样本 | 少量预先固定的检验 |
| V5 | 实验 | 见 §5.10、§7.4 | 因果方向 |

关于 V3 的两点限制：

- GSE146590 是未处理对照，加上阿糖胞苷处理后分选出的衰老样细胞。分选过的细胞不能代表处理后的全部细胞，供者也少，只能核对方向，不能训练或评估预测模型。
- 目前没有任何公开数据对应“阻断某条通讯”。这类结论只能等实验。

---

## 9. 暂缓的方向及原因

| 方向 | 原因 | 什么情况下重开 |
|---|---|---|
| 跨样本通讯张量（MEBOCOST、SigXTalk 打分拼成张量再分解；Tensor-cell2cell） | 分数是“发送细胞群均值 × 接收细胞群均值”，和旧主线同一结构，已知不超过细胞比例基线 | 不重开；问题本身改到 WP3 |
| CellOT 反事实预测当主线 | AML 没有合格的公开扰动训练数据；没有“阻断通讯”的数据 | 自己产生扰动数据之后 |
| 患者前后的“治疗效应”推断（LEMUR、CINEMA-OT、MrVI） | 患者前后是观察数据，药物作用、克隆选择、造血恢复混在一起，不能叫治疗效应 | 可作探索，统一叫“配对变化”；LEMUR 的设计写成 `~ patient + timepoint` |
| 把时间点当 ContactTracing 的条件 | 化疗直接改变接收细胞；原文的条件是肿瘤细胞内部的遗传操作，不直接作用于免疫细胞 | 只作探索，不进主结论 |
| TIGON 增长率 | 需要每位患者至少 3 个时间点，且细胞数能代表真实数量；现有数据不满足 | 有未分选的密集时间点数据之后 |
| 以细胞共现或生态型为主线（CoVarNet 式） | AML 已有同类工作（ACE）；比例受 blast 负荷和富集流程影响 | 只在 WP2 里作为基线对照 |
| scComm、scHyper 的打分 | scComm 把得分最高 5% 的细胞类型对当作真实事件来训练，标签来自自己的打分 | 不采用打分；“通讯角色”这个问题已并入 WP1b |
| 图神经网络、单细胞大模型（PINNACLE、Geneformer 等） | 2025 年 Nature Methods 的评测显示深度学习预测扰动效应未超过线性基线；AML 缺少可监督的标签 | 主线成立后再考虑 |
| 把治疗当工具变量（scIVCCC） | 化疗直接作用于多类接收细胞，排除限制不成立 | 不重开 |

---

## 10. 数据清单

样本级信息以项目里已经整理好的元数据表（`meta_v2.3`）为准。下表的行数来自 `AML_niche_CCC_dataset_inventory.xlsx`，是较早的版本，只用来判断量级。

### 10.1 项目清单内

| 数据集 | 研究 | AML 侧样本行 | 时间点 | 细胞范围 | 平台 | 原始 reads | 用于 | 注意 |
|---|---|---|---|---|---|---|---|---|
| GSE185381 | Lasry 2023 Nat Cancer | 39（另有健康 25） | 诊断，成人和儿童 | 未分选骨髓 | 10x 3′、5′，部分 CITE | SRA | WP1、WP2 发现 | 多数供者跨多个文库；部分文库同时含 AML 和健康供者，按 Patient_ID 拆分 |
| GSE239721 | IFNγ signaling in AML | 20 | 诊断 | 未分选 | 10x 5′ | 有 | WP1、WP2 发现 | 作者的患者编号对照表缺失 |
| GSE289435 | Zeng 2025 Blood Cancer Discov | 12（10 AML、2 急性红系白血病） | 诊断 | BMMC | 10x 5′ v1.1 | SRA | WP1、WP2 发现 | Numbat 补跑是否完成需确认 |
| GSE227903 | Ennis 2023 iScience | 诊断 10、MRD 7、复发 11 | 三个阶段 | 未分选 | 10x 3′ | 有 | WP1、WP2、V4 | 11 对中的 9 对来自这里 |
| GSE116256 | van Galen 2019 Cell | 诊断 16、疗后 19（另有健康） | D0 到 D171 | 骨髓单个核细胞 | Seq-Well | 有 | WP2 纵向；基因分型层只作外部核对 | 每样本细胞数少；同一天数可能对应相反的临床状态，要按 blast% 校正 |
| GSE201966 | 移植后复发的急性单核细胞白血病 | 初发 3、复发 3、缓解 1 | — | 未分选 | 10x | 未确认 | V4 | 11 对中的 2 对来自这里 |
| Zenodo 3345981 | Petti 2019 Nat Commun | 5（另有健康 4） | 诊断 | 未分选 | 10x 5′ | 原始 reads 在 dbGaP（受控），矩阵在 Zenodo | WP1a 参照 | 原文已给出逐细胞突变检出 |
| Mendeley gwjh3w6ztm.2 | Chen 2023 Blood Cancer Discov | 6 例 NPM1 突变 AML，各 2 个分选组分（另有健康 4 例） | 诊断 | CD34⁺ 组分；niche 加免疫组分 | 10x | 未确认 | GATE_W1.2 阳性对照 | 组分比例由分选决定，不能当真实比例 |
| GSE185991 | Naldini 2023 Nat Commun | 诊断 11、D14 11、D30 5、复发 4 | 四个节点 | 分选的原始细胞和祖细胞 | 10x 3′ | 有 | 恶性细胞程序的验证；WP1a 的 NPM1 读数可试 | 不含完整微环境；D14、D30 残存细胞很少 |
| GSE147989 | Riether 2020 Nat Med | 4 | 用药前后 | 外周血 CD34⁺ | 10x 3′ v2 | 有 | 暂不用 | 只有干祖细胞 |
| GSE207356 | Nicosia 2023 Cancer Cell | 3 | 用药前后 | 外周血 | 10x 3′ | SRA | 暂不用 | 单个患者 |
| GSE154109 | Bailur 2020 JCI Insight | 儿童 AML 8、B-ALL 7、健康 4 | 初诊 | 去除肿瘤细胞后的免疫细胞 | 10x | SRA | 对照 | 肿瘤细胞被人为去除 |
| E-MTAB-11536 | Domínguez Conde 2022 Science | 健康骨髓 | — | 免疫细胞 | 10x | 有 | 健康参照 | 少数供者贡献了大部分细胞，按供者均衡 |
| GSE253355 | Bandyopadhyay 2024 Cell | 健康骨髓 | — | 基质富集 | 10x 3′ v3.1 | SRA | 健康参照 | 组分比例由富集决定 |
| BoneMarrowMap | Zeng 2025 | 参考图谱 | — | — | — | — | 层级投影 | 不是队列 |

### 10.2 新增候选（都需要先做元数据整理）

| 数据 | 入口 | 内容 | 用于 | 注意 | 核实 |
|---|---|---|---|---|---|
| Mumme 2023 Nat Commun | GSE235923 | 儿童 AML 骨髓；诊断、诱导结束、复发；10x 3′ v3；带复发或持续缓解的结局 | WP2 缓解期对比和结局 | 儿童；配对数需查 | ✔ |
| Lambo 2023 Cancer Cell | GSE235063；EGAS00001007323 | 28 例儿童 AML；诊断、缓解、复发；骨髓或外周血 | WP2 | GEO 只有处理好的矩阵，原始 reads 受控；不是每位患者都有三个节点 | ✔ |
| AML scAtlas（Whittle 2025 eLife） | CELLxGENE 合集；AnnData 的 DOI 见资料 | 748,679 个细胞；159 例 AML、51 例健康；20 项研究 | 找研究、跨队列重复 | 汇集的是已有研究，要追溯到原研究并去重；与项目清单有重叠 | ✔（入口链接 ○） |
| 儿童 AML 诊断–复发配对 | GEO，登录号待核对 | 33 例患者的配对样本 | V1、V4 | 只在检索摘要里见到，细节未核实 | ○ |
| Li 等 2023 Leukemia | 登录号待查 | 难治或早期复发 AML，化疗过程中的纵向样本 | V1 | 数据是否公开未核实 | ○ |
| Duy 2021 Cancer Discov | GSE146590 | 患者来源 AML 细胞；未处理对照与阿糖胞苷处理后分选的衰老样细胞 | V3 | 分选过，只核对方向 | ✔ |
| Pei 2020 Cancer Discov | GSE143363 | 2 位患者诊断和复发配对的 CITE-seq | 不建议用 | 量太小；主题是单核样亚克隆 | ✔ |

### 10.3 外部资源

| 资源 | 用途 | 核实 |
|---|---|---|
| BeatAML（含离体药敏和生存） | V2 | ○ |
| TCGA-LAML | V2 | ○ |
| CytoSig 反应特征矩阵和细胞因子处理数据库 | 信号活性打分；V3 | ○ |
| PROGENy、CollecTRI（经 decoupleR） | 通路活性和转录因子活性 | ○ |
| 配受体数据库（CellChatDB、CellPhoneDB、CellTalkDB） | 候选配受体。ContactTracing 原文取 CellTalkDB 与 CellPhoneDB 的并集，并把复合体表达定义为各成员的最小值 | ✔ |
| MEBOCOST 的代谢物–sensor 数据库 | WP3 候选名单 | ✔ |

### 10.4 每个样本要补齐的字段

在 `meta_v2.3` 已有字段之外，本主线特别需要：是否有带 CB/UB 的 BAM；文库是 3′ 还是 5′；是否保留了未过滤的液滴矩阵；分选或富集流程；blast%；临床 NPM1、FLT3 状态；疗后样本距用药的天数和当时的缓解状态；复发或持续缓解的随访结局。未知的字段保留为未知。缓解、D14 残存、难治、MRD 阳性不能互换。

---

## 11. 文献清单

### 11.1 直接要实现或改写的方法

| 方法 | 文献 | DOI | 用在哪里，取它的哪一步 | 核实 |
|---|---|---|---|---|
| ContactTracing | Li 等，Non-cell-autonomous cancer progression from chromosomal instability，Nature 620:1080–1088（2023） | 10.1038/s41586-023-06464-z | WP1b。两步检验：同类细胞内受体阳性对阴性（MAST hurdle 模型的似然比检验）；再加“受体 × 条件”的交互项。不需要预先知道下游靶基因。原文用 CytoSig 数据库和配体处理实验做了验证 | ✔ |
| Tres | Zhang 等，A T cell resilience model associated with response to immunotherapy in multiple tumor types，Nat Med 28:1421–1431（2022） | 10.1038/s41591-022-01799-y | WP1a。先用 CytoSig 估计细胞受到的信号，再做“增殖 ~ 信号 × 基因”的交互检验，找削弱负相关的基因。原文用 168 个肿瘤、19 种癌 | ✔ |
| CytoSig | Jiang 等，Systematic investigation of cytokine signaling activity at the tissue and single-cell levels，Nat Methods（2021） | 待补 | 信号活性打分；V3 的参照数据库 | ○ |
| scITD | Mitchel 等，Coordinated, multicellular patterns of transcriptional variation that stratify patient cohorts are revealed by tensor decomposition，Nat Biotechnol 43:1192–1201（2025，2024 年 9 月在线） | 10.1038/s41587-024-02411-z | WP2。供者 × 基因 × 细胞类型的张量做 Tucker 分解；附带配体–受体分析 | ✔ |
| DIALOGUE | Jerby-Arnon 和 Regev，DIALOGUE maps multicellular programs in tissue from single-cell or spatial transcriptomics data，Nat Biotechnol 40（2022） | 10.1038/s41587-022-01288-0 | WP2 的另一种实现。多细胞程序（MCP）；代码 `github.com/livnatje/DIALOGUE` | ✔ |
| CARE 第一篇 | Nomura 等，The multilayered transcriptional architecture of glioblastoma ecosystems，Nat Genet 57:1155–1167（2025） | 10.1038/s41588-025-02167-5 | WP2 的“状态控制的 pseudobulk”：把细胞组成、细胞状态、状态内部的基线程序三层分开 | ✔ |
| CARE 第二篇 | Spitzer 等，Deciphering the longitudinal trajectories of glioblastoma ecosystems by integrative single-cell genomics，Nat Genet 57:1168–1178（2025） | 10.1038/s41588-025-02168-4 | 纵向配对的写法参照：全队列没有一致轨迹，部分患者共享特定轨迹。队列 59 例 | ✔ |

两篇 CARE 文章共用同一个队列，不能当作两次独立验证。

### 11.2 只提供候选名单或交叉核对

| 方法 | 文献 | DOI | 用法 | 核实 |
|---|---|---|---|---|
| SigXTalk | Dissecting crosstalk induced by cell-cell communication using single-cell transcriptomic data，Nat Commun 16:5970（2025） | 10.1038/s41467-025-61149-7 | WP3 候选：受体 → 信号因子或转录因子 → 靶基因的汇聚关系（超图学习） | ✔ |
| MEBOCOST | Zheng 等，MEBOCOST maps metabolite-mediated intercellular communications using single-cell RNA-seq，Nucleic Acids Res 53(12):gkaf569（2025） | 10.1093/nar/gkaf569 | WP3 候选：代谢物–sensor 数据库。它的通讯分数是发送细胞群酶均值乘以接收细胞群 sensor 均值，不采用 | ✔ |
| scTenifoldXct | A semi-supervised method for predicting cell-cell interactions and mapping cellular communication graphs，Cell Systems（2023） | 10.1016/j.cels.2023.01.004 | 可选：用胞内基因网络的一致性筛候选 | ○ |

### 11.3 探索用或暂缓

| 方法 | 文献 | DOI | 定位 | 核实 |
|---|---|---|---|---|
| LEMUR | Ahlmann-Eltze 和 Huber，Analysis of multi-condition single-cell data with latent embedding multivariate regression，Nat Genet 57:659–667（2025） | 10.1038/s41588-024-01996-0 | “配对变化”的探索；设计矩阵写患者加时间点；CPU 可跑 | ✔ |
| CINEMA-OT | Dong 等，Causal identification of single-cell experimental perturbation effects with CINEMA-OT，Nat Methods 20:1769–1779（2023） | 10.1038/s41592-023-02040-5 | 适合受控扰动数据；有处理细胞比例变化的重加权版本 | ✔ |
| MrVI | Boyeau 等，Deep generative modeling of sample-level heterogeneity in single-cell genomics，Nat Methods 22:2264–2274（2025） | 10.1038/s41592-025-02808-x | 探索：只出现在某些细胞亚群里的患者分层；需要 GPU | ✔ |
| CellOT | Learning single-cell perturbation responses using neural optimal transport，Nat Methods 20:1759–1768（2023） | 10.1038/s41592-023-01969-x | 暂缓，理由见 §9 | ○ |
| FlowSig | Almet 等，Inferring pattern-driving intercellular flows from single-cell and spatial transcriptomics，Nat Methods 21:1806–1817（2024） | 10.1038/s41592-024-02380-w | 暂缓。非空间数据需要对照加扰动的设计；输入含细胞通讯推断结果 | ✔ |
| TIGON | Sha 等，Reconstructing growth and dynamic trajectories from single-cell transcriptomics data，Nat Mach Intell（2024） | 10.1038/s42256-023-00763-w | 暂缓，数据不满足 | ✔ |
| CoVarNet | Shi 等，Cross-tissue multicellular coordination and its rewiring in cancer，Nature 643:529–538（2025） | 10.1038/s41586-025-09053-4 | 只作 WP2 的比例基线对照 | ✔ |
| scComm | Jin 等，scComm: a contrastive learning framework for deciphering cell–cell communications at single-cell resolution，Genome Biol 27:149（2026） | 10.1186/s13059-026-04043-9 | 不采用打分 | ✔ |
| scHyper | scHyper: reconstructing cell–cell communication through hypergraph neural networks，Brief Bioinform（2024） | 10.1093/bib/bbae436 | 不采用 | ○ |
| Tensor-cell2cell | Nat Commun（2022） | 10.1038/s41467-022-31369-2 | 不采用 | ○ |
| PINNACLE | Li 等，Contextual AI models for single-cell protein biology，Nat Methods 21:1546–1557（2024） | 10.1038/s41592-024-02341-3 | 不采用 | ✔ |
| scIVCCC | Brief Bioinform（2026），bbag139 | — | 不采用 | ○ |
| 扰动预测评测 | Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines，Nat Methods 22:1657–1661（2025） | 待补 | §9 的依据 | ✔（题目和出处） |

### 11.4 AML 及相关领域的已有工作（写作时必须交代区别）

| 文献 | 内容 | 与本主线的关系 | 核实 |
|---|---|---|---|
| Chen 等 2023，A Single-Cell Taxonomy Predicts Inflammatory Niche Remodeling to Drive Tissue Failure and Outcome in Human AML，Blood Cancer Discov（PMC10472197） | NPM1 突变 AML：炎症性 niche 重塑压制残余正常干细胞亚群，克隆细胞相对耐受；HSPC 组分里的 NPM1 突变细胞多见于 LMPP、GMP | WP1a 最接近的已有结论；GATE_W1.2 的阳性对照 | ✔ |
| Jakobsen 等 2024，Cell Stem Cell（10.1016/j.stem.2024.05.010） | 克隆性造血：突变干细胞相对同样本野生型干细胞炎症反应减弱 | WP1a 的概念先例 | ✔ |
| Beneyto-Calabuig 等 2023，Cell Stem Cell 30:706–721 | CloneTracer：19 例 AML 的克隆判定；静息干细胞多为健康或白血病前 | 基因型概率模型的参照 | ✔ |
| Petti 等 2019，Nat Commun 10:3660（10.1038/s41467-019-11591-1） | 10x 5′ 文库里读表达突变；两例样本 NPM1 位点有覆盖的细胞为 5,591/11,620 和 11,672/20,474；SNV 平均可在 22.7% 的细胞里检出 | WP1a 可行性依据 | ✔ |
| ACE（sciNMF），Brief Bioinform 26(1):bbaf028（2025）（10.1093/bib/bbaf028） | 6 套数据、68 例、256,352 个细胞；26 个细胞状态；按共存模式定义 AML 细胞生态系统并做预后模型 | WP2 的已有工作 | ✔ |
| Whittle 等 2025，eLife | AML scAtlas | 数据入口 | ✔ |
| Mumme 等 2023，Nat Commun 14:6209（10.1038/s41467-023-41994-0） | 儿童 AML 诊断、诱导结束、复发；复发组与持续缓解组的微环境比较 | WP2 数据和已有工作 | ✔ |
| Lambo 等 2023，Cancer Cell | 28 例儿童 AML 纵向图谱；复发时层级更原始 | WP2 数据 | ✔ |
| Duy 等 2021，Cancer Discov 11:1542–1561（10.1158/2159-8290.CD-20-1375） | 化疗后的衰老样细胞 | V3 数据 | ✔ |
| Pei 等 2020，Cancer Discov 10:536–551（10.1158/2159-8290.CD-19-0710） | VEN/AZA 耐药与单核样亚克隆 | 不建议用 | ✔ |
| Li 等 2023，Leukemia（10.1038/s41375-022-01789-6） | 难治或早期复发 AML 的化疗重编程 | 可能的外部数据 | ✔（题目和出处） |
| Lasry 等 2023，Nat Cancer；Ennis 等 2023，iScience；van Galen 等 2019，Cell；Naldini 等 2023，Nat Commun；Zeng 等 2025，Blood Cancer Discov | 项目清单内数据的原始文献 | 引用和区分 | ○ |

---

## 12. 统计与复现规则

- 外层按患者或按研究留出。同一患者的所有时间点放在同一分区。不能先用全部数据学出模块，再把结局预测称为独立验证。
- 每位患者、每个细胞类型做 pseudobulk 或低维汇总后再比较。不确定性用患者级 bootstrap 和配对置换来估计。
- 每个研究内部先估计效应，再比较或汇总。不同方案、年龄、取样时间和材料的样本不直接混成一个队列。
- 打分用的基因和被检验的基因尽量不重叠。每个结论配上基因集匹配、组成匹配、标签置换和简单方法基线。
- 没有原始采样比例时，不把富集样本的细胞频率当成真实骨髓丰度。
- 同一研究衍生出的图谱文件、GEO 文件、CELLxGENE 文件，同一位患者只算一次。
- 门槛和主读数在看结果之前写进 config 并提交，改动要留记录。
- 每个含随机性的步骤固定随机种子，种子取 `config_paths.R` 里的全局 `SEED`。

---

## 13. 结论的边界

可以说的：

- 同一样本、同一状态内，两类细胞的反应特征有差别，并且这种差别随环境信号强度变化。
- 某个配体多的患者里，带对应受体的细胞出现了特定的表达变化。
- 某个多细胞程序在缓解期仍高于健康水平，并与结局相关。

不能说的：

- RNA 推出的通讯不等于物理邻近、配体结合、代谢物流量。
- 同一分箱里的突变和野生型细胞“暴露相同”是假设。骨髓里如果两类细胞位置不同，这个假设会被破坏，文中要写明。
- 患者前后的差别不是治疗的因果效应。
- 删除模型里的一条边导致输出变化，只说明模型敏感，不说明真实干预的结果。
- 只有表达数据，不能得出“阻断某条通讯可以避免复发”。
- 临床报告的 NPM1 或 FLT3 状态是患者层面的标签，不是细胞层面的亚克隆标签。

最低可写的结果由三类证据定义：跨研究稳定；患者内有相符的纵向变化；独立的扰动数据或实验支持预测方向。缺第三类时，结论限定为可重复的关联和候选机制。

---

## 14. 项目约定与环境

- 对话用中文，代码注释用英文。写代码前一次只问一个澄清问题。
- 代码和脚本一律以可下载文件交付，不贴在对话里让人复制。
- 目录号等于阶段号，目录内脚本从 01 编号，命名 `NN_verb_noun.ext`，测试脚本用 `90_`、`91_` 前缀。
- 路径只在 `config/config_paths.{R,sh}` 里定义，用 `.here` 定位项目根；每个阶段一个 `config_<stage>.R`；业务脚本里不写路径字符串。
- 日志用 `message("[N] ...")`。每个产出步骤先查输出是否存在，存在就跳过。长任务用 `sbatch --requeue`。
- 全局 `SEED`：`CODING_STANDARDS.md` 写的是 `20260605L`。以 `config_paths.R` 里的实际值为准，开工时核对一次。
- 恶性标签只由 `02_malignancy/50_consensus_malignancy.R` 产出。CNV 信号称 “CNV-based proxy”，不称突变。NPM1 读数判定是另一列（`npm1_genotype`），不并入恶性共识，除非另行决定。
- 大对象放 `LARGE1`，表、图、脚本放 `FAST`（会被清理）。脚本要推到 GitHub。
- 计算环境：SLURM，分区 `gr10634c`，资源写法 `--rsc p=1:t=:c=:m=`；项目根 `/FAST/gr10634/gaozy/aml_niche_net`（脚本）和 `/LARGE1/gr10634/gaozy/aml_niche_net`（数据）；conda 环境 `/FAST/gr10634/gaozy/general_env`。
- 已有工具：Seurat v5、CellChat、GeneNMF、Symphony、data.table、here；Python 侧有 POT、pandas、numpy、scipy。本主线可能新增：MAST、decoupleR、CytoSig、pysam、SoupX 或 DecontX、scITD 或 DIALOGUE、lme4。

**建议的新目录**（不动 05–08，那是旧主线）：

```
scripts/
├── 11_genotype/         # NPM1 read-level genotype (the planned SNV evidence type)
├── 12_signal/           # per-cell signal activity and fitness scores
├── 13_within_sample/    # WP1a, WP1b
├── 14_mcp/              # WP2
├── 15_redundancy/       # WP3
└── 16_validation/       # WP4
```

**可直接复用的已有产出**：`01_preprocess` 的 per-sample QC 对象；`02_malignancy` 的逐细胞共识表和 STARsolo BAM；`03_hierarchy` 的逐细胞分箱和干性评分；`04_cnmf` 的 meta-program；`prop7` 基线和既有的配对检验框架。

旧蓝图的 L1（不做跨样本表达整合）在这里仍然成立：WP1 全在样本内部，WP2 用每个样本自己的 pseudobulk。LEMUR、MrVI 属于整合类方法，只放在探索项里。

---

## 15. 待冻结的决定

1. 从同一 scRNA 文库的原始 reads 里读突变，是否算在“只用 scRNA-seq”之内。建议算；van Galen 的靶向基因分型只作外部核对。
2. WP1a 的主检验。建议以“反应差”（检验 1）为主，适应度交互（检验 3）为次。
3. 适应度代理的主指标。建议增殖评分为主，凋亡和干性评分为辅。
4. 基因型后验的判定线。建议突变 ≥ 0.9，野生型 ≤ 0.1。
5. 各 GATE 的数值（§5.7、§6.6）。建议先跑 §16 的 T1–T5，再按实际数字定。
6. 儿童数据是否纳入。建议纳入，但与成人分开分析。
7. 是否申请受控数据（Lambo 的原始 reads、Petti 的 dbGaP 数据、CloneTracer 的 EGA 数据）。
8. 环境 RNA 用哪个工具，定下后不换。
9. 新目录编号（11–16）是否采用，蓝图文件是否同步更新。
10. WP2 预先固定的模块怎么选。建议取发现研究里可重复性最高的前 4 个。
11. WP2 用 scITD、DIALOGUE，还是自己写。
12. doublet 预处理状态的处理方式。这是 WP1 的前置条件。

---

## 16. 开工顺序

| 编号 | 任务 | 产出 | 对应门槛 |
|---|---|---|---|
| T0 | 读本文件；确认 §15 的第 1、12 项 | 两条决定 | — |
| T1 | 从 `meta_v2.3` 筛出满足 §5.2 的样本 | 每个数据集的患者数、文库类型、BAM 是否在、未过滤矩阵是否在 | — |
| T2 | 写 `11_genotype/01_count_npm1_site.py`（附录 A）；在一套 5′ 数据和一套 3′ 数据里各取 1–2 个 NPM1 突变样本试跑 | 位点覆盖率；T/NK/B 细胞的突变 UMI 比例；各分箱的突变与野生型细胞数 | GATE_W1.1 |
| T3 | 统计候选配体在恶性细胞里的检出比例、候选受体在各接收类型里的检出比例 | 能进入第二阶段的组合数 | GATE_W1.3 |
| T4 | 估计各样本的环境 RNA 污染比例 | 每样本一个数，加环境基因黑名单 | GATE_W2.1 的前置 |
| T5 | 在 Chen 2023 数据上重复原文方向 | 一页结果 | GATE_W1.2 |
| T6 | 按 T1–T5 的数字冻结 §15 的阈值，写 `config_genotype.R`、`config_signal.R`、`config_within.R` | 提交到仓库 | — |
| T7 | 写 WP1 主流程 | §5.8 的交付物 | GATE_W1.4 |
| T8 | WP1 有结果后启动 WP2；WP3 等 WP1、WP2 | — | — |

T1–T5 互不依赖，可以并行。任何一个门槛不通过，先停下来讨论，不往下写。

---

## 17. 新 session 的开场提示词

```
请先读 HANDOFF_within_sample_v1.md，以及 ARCHITECTURE.md 和 CODING_STANDARDS.md。
本次只做 §16 的 T0–T2，不写 WP1 主流程。
1. 先和我确认 §15 的第 1 项和第 12 项，一次只问一个问题。
2. 按 §5.2 的条件，告诉我你需要 meta_v2.3 的哪些列，我来提供或你来读。
3. 按附录 A 写 11_genotype/01_count_npm1_site.py，以可下载文件交付，遵守 CODING_STANDARDS.md。
   位点坐标和等位序列从比对所用的参考基因组和 GTF 里取，不要手写。
4. 给出在集群上试跑 1–2 个样本的 sbatch 脚本，输出覆盖率、T/NK/B 细胞的突变 UMI 比例、各分箱的细胞数。
不通过门槛时先停下来告诉我数字，不要自行换方案。
```

---

## 18. 文章骨架（暂定，随结果调整）

1. 设计图：为什么把比较放到样本内部；数据总览。
2. 突变细胞与野生型细胞的反应差；环境剂量关系；跨队列一致性。
3. 提供者与受益者：配体的不均匀表达；受体依赖的反应；与基因型的 2×2 关系。
4. 多细胞程序：扣掉比例之后的协同；缓解期是否仍异常；与结局的关系。
5. 代偿与共享瓶颈。
6. 验证：外部队列、bulk、公开扰动数据、实验。

---

## 附录 A：NPM1 位点的逐细胞读数与判定

### A.1 读数

1. **位点和等位序列从参考文件里取**。用比对时的参考基因组 FASTA 和 GTF 定位 NPM1 最后一个编码外显子里的插入位点，取两侧各约 40 bp。下面的写法只作提示，使用前核对：
   - 参考转录本 NM_002520；A 型 `c.860_863dupTCTG`，B 型 `c.863_864insCATG`，D 型 `c.863_864insCCTG`。A 型占大多数。
   - 另加一条通用规则：该位点处任意 4 bp 插入都记为突变候选，单独标记。
2. **取 reads**：用 pysam 取位点上下各约 100 bp 范围内的 reads，包括被软剪切的 reads。插入位点靠近 read 末端时，比对软件常把它剪掉，所以不能只看 CIGAR。
3. **按 k-mer 分类**：对每条 read 的原始序列，找跨过插入点的等位特异 k-mer（k 取 20–24）。跨过插入点且两侧各有不少于 8 bp 的 read 才算有信息。分为突变、野生型、无信息三类。
4. **按 UMI 合并**：同一个 `CB` 加 `UB` 下的 reads 取多数；出现冲突的 UMI 丢弃。
5. **输出**：每个细胞的 `n_mut_umi`、`n_wt_umi`。

3′ 文库能否读到这个位点要实测。已发表的数字来自 5′ 文库。

### A.2 判定

记 `n = n_mut_umi + n_wt_umi`，`m = n_mut_umi`。每个样本估三个量：

- `theta`：突变细胞里突变 UMI 的期望比例。杂合突变且两个等位表达相近时约为 0.5，用数据估计。
- `eps`：野生型细胞里出现突变 UMI 的比例，来自环境 RNA 和错误。用 T、NK、B 细胞估计。
- `pi_b`：层级分箱 b 里突变细胞的先验比例，用 EM 估计。

```
P(mut | m, n) = pi_b * Binom(m; n, theta)
                / ( pi_b * Binom(m; n, theta) + (1 - pi_b) * Binom(m; n, eps) )

label = "mut"  if P >= 0.9
        "wt"   if P <= 0.1
        "undetermined" otherwise
```

### A.3 为什么不能用“读到 3 条野生型就算野生型”

取 `theta = 0.5`、`eps = 0.01`，细胞里一条突变 UMI 都没有：

| 分箱里突变细胞的先验比例 | 野生型 UMI 数 | 仍是突变细胞的后验概率 |
|---|---|---|
| 0.9 | 3 | 0.54 |
| 0.9 | 7 | 0.07 |
| 0.5 | 3 | 0.11 |
| 0.5 | 4 | 0.06 |

在突变细胞占多数的分箱里，3 条野生型 UMI 远远不够，要 7 条左右。所以判定必须带上分箱的先验，不能用一个固定条数。

### A.4 核对

- 健康供者样本里突变 UMI 应接近 0。
- T、NK、B 细胞里的突变 UMI 比例就是污染率，逐样本报告。
- van Galen 数据有独立的靶向基因分型，可用来核对读数判定的一致性。
- Petti 数据原文给出了逐细胞的突变检出，可对比覆盖率量级。

### A.5 非 NPM1 患者

用 Numbat 的克隆后验代替，流程相同。比较的是 CNV 克隆与正常细胞，措辞用 “CNV-based proxy”。Numbat 给出 `no_CNV_detected` 的样本记为缺失，不记为正常。

