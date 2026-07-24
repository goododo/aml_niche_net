# 白血病干细胞向骨髓微环境通讯的保守依赖轴及其在治疗压力下的拓扑重塑

### —— 一个跨异质性 AML 单细胞队列的最优传输图对齐发现框架（Project Design Blueprint v1.0）

> **文档定位**：本稿是**可执行的科学方案蓝图 + 面向 NCS 级别论文的 project design 底稿 + 个人/组内 pipeline 归档**。主体按"**先有发现、方法服务于发现**"的逻辑组织；摘要/意义/创新性力求凝练有冲击力，而 Methods / Validation / Expected Outcomes 一律写成**可证伪的假设检验**，不使用"结果已出"的口吻。
>
> 全文夹带三类 **co-designer 批注**，用于直面 NCS 审稿与自我执行的卡点： **【R#｜审稿风险】**（reviewer 最可能攻击的点）、**【D#｜设计决策】**（已锁定的取舍）、**【M#｜对原框架的修订】**（我建议改动你原始设计的地方）。
>
> 配套文件：`AML_niche_CCC_dataset_inventory.xlsx`（339 个样本逐条分层，研究×角色汇总，图例）。本稿正文只写**类别**，逐数据集明细见该 Excel 与附录 A 精简表。
>
> 数学公式用规范 LaTeX；正文引用用 `(作者, 年份)`，文末按编号列出带 DOI/PMID 的完整文献。本版为中文主稿，设计冻结后再生成严格双语镜像版。

------------------------------------------------------------------------

## 0. 中心发现命题（Central Discovery Claim）

> **中文**：在遗传背景高度异质的 AML 中，**LSC 样恶性状态会收敛到一组保守、可跨平台与跨队列复现的微环境通讯依赖轴**；疾病信号编码在「恶性细胞 → 骨髓微环境」**通讯网络的拓扑结构**之中，而非任何单一 ligand–receptor（L–R）对；该拓扑相对健康造血的偏离程度，沿**治疗压力轴（诊断 → 残留/MRD → 复发）单调加深**。
>
> **EN**：Across genetically heterogeneous AML, LSC-like malignant states converge onto a small, reproducible set of conserved microenvironment-communication dependency axes; the disease signal is encoded in the **topology** of the malignant→microenvironment communication network rather than in any single ligand–receptor pair; and the magnitude of this topological deviation from healthy hematopoiesis increases **monotonically along the treatment-pressure axis (diagnosis → residual/MRD → relapse).**

这一命题把两个直觉合并成一个更强、可证伪的原则：(i) *收敛性*——不同遗传亚型的 LSC 样细胞利用一组**共同**的 niche 依赖轴；(ii) *拓扑性*——进展/复发不是单个 L–R 改变，而是通讯网络的**系统性重排**。方法（最优传输 / Fused Gromov–Wasserstein 图对齐）只是把"拓扑偏离"变成**可度量、可比较、可证伪**的量。

> **【D1｜"niche" 的操作性重定义】** 经数据盘点（见 §5），公开 AML 抽吸数据几乎不含结构性血管–成骨 niche。本项目把 **"microenvironment / niche" 主定义为造血-免疫与旁分泌微环境**（T、NK、单核/巨噬、B、DC、巨核、粒系），结构性基质（MSC/内皮/成骨/周细胞）降级为 **cohort 层辅助证据**。全文措辞用 "BM microenvironment"，避免 Frenette/Morrison 式结构 niche 的过度声称。这是消解头号审稿攻击的前提（见 **R1**）。

------------------------------------------------------------------------

## 1. 摘要（Abstract）

急性髓系白血病（AML）以克隆、转录组与微环境的多重异质性著称；骨髓微环境通过 L–R 通讯为白血病干细胞（LSC）提供庇护，驱动治疗抵抗与复发 (Ennis et al., 2023; Chen et al., 2023)。然而，现有 cell–cell communication（CCC）比较框架大多停留在 **cell-type 级**的静态描述，既未把通讯拓扑锚定到 **恶性细胞状态层级**与**恶性置信度**，也无法区分"跨队列真正保守的通讯结构"与"平台/注释/抽样造成的假象"。本研究提出一个三层发现框架：(L1) 在**每个队列内部、不做跨样本表达整合**地定义恶性 meta-program 并量化其跨研究复现性；(L2) 构建**锚定于恶性层级的患者特异 CCC 图**；(L3) 用**最优传输 / Fused Gromov–Wasserstein（FGW）图对齐**在拓扑空间比较患者图与"健康/AML/治疗压力"共识重心（barycenter），输出可证伪的患者级拓扑偏离评分。我们把公开 AML 队列的异质性**当作对稳健性的生物学压力测试**：只有跨平台、跨中心、跨患者构成反复出现的通讯模块，才被认定为核心 AML–微环境依赖。框架显式分层，使每层都有可独立成文的产出，从而对更高、技术更难的层做风险对冲。我们将检验三个可证伪假设：LSC 样状态的通讯依赖轴在亚型间**收敛**（H1）；网络**拓扑**比任意单一 L–R 或细胞构成更能判别 AML/健康（H2）；该拓扑偏离沿**治疗压力轴单调加深**（H3）。计算发现为主线，最后提出少量可行的空间/功能/扰动实验以检验候选轴的因果性。

------------------------------------------------------------------------

## 2. 科学问题与中心假设

**核心科学问题。** 在恶性细胞状态层级上，LSC 样 / blast 样状态如何重组其与骨髓微环境的通讯？哪些通讯模块在异质 AML 间**可复现地保守**（而非平台/注释/抽样伪影）？患者特异的 CCC 网络**拓扑偏离**能否定量刻画疾病状态与治疗压力下的演变？

三个**可证伪假设**（每个都给出"被证伪"的明确判据，见 §9 评估）：

- **H1（收敛性 / Convergence）**：存在一组 LSC 样 → 微环境通讯轴，在 ≥ 60% 独立研究中复现，且在多个遗传亚型（NPM1、TP53、KMT2A 重排、单核细胞型等）中均富集，显著超过置换零分布。 *证伪*：保守轴数量不超过零分布，或保守性可被平台/批次完全解释。
- **H2（拓扑性 / Topology \> single pair）**：网络拓扑偏离（HDS/ATS，经 FGW/GW 距离）在严格跨队列验证下，判别 AML vs 健康的能力**优于**最佳单一 L–R 模型与"仅细胞构成"模型。 *证伪*：单 L–R 或构成模型达到与拓扑模型相当的判别力（即拓扑无增量信息）。
- **H3（治疗压力单调性 / Monotonic intensification）**：在配对患者内，拓扑偏离沿 **Dx \< MRD/残留 \< 复发**单调增加。 *证伪*：无单调趋势，或趋势可被 blast 比例单独解释。

> **【M1｜把"复发"从主 claim 降级为验证终点】** 数据盘点显示真正可用于全骨髓 CCC 的"诊断–复发"配对仅 \~12 例、且高度集中在 1–2 个研究（见 §5），无法支撑"跨异质队列保守的复发原则"这一主 claim——否则与"跨队列保守"卖点自相矛盾。故**主时间轴改为"治疗压力"（Dx→MRD/post-tx，功效更足）**，复发（RLS）作为**探索性验证终点**。H3 因此以治疗压力为主、复发为极端点。

------------------------------------------------------------------------

## 3. 意义与创新性

- **从"用了复杂方法"转向"发现一条普适原则"**：本框架的卖点不是"我做了 FGW"，而是"AML 的微环境依赖在分子多样性之上**收敛**、且疾病编码在**网络拓扑**而非单点"。方法是为这条原则服务的必要工具。
- **把异质性从"噪声"转为"证据"**：跨平台/中心/构成仍复现的模块 = 核心依赖；这把公开数据的最大短板（异质性）变成稳健性论证（见 (Gavish et al., 2023) 在泛癌 meta-program 上的同类逻辑）。
- **拓扑级疾病编码的可证伪量化**：HDS/ATS/RLS/SCS 是一组可被临床 metadata 检验的患者级读数，把"网络重排"从定性叙述变成可统计的量。
- **方法学新意（相对先行工作）**：(Nagai et al., 2025) 的 scACCorDiON 在 **cell-type 级**有向图上做 OT 重心；本框架首次把图节点**锚定到参考投影的恶性细胞状态层级 + meta-program 状态 + sender/receiver 角色**，并用 **FGW**（同时编码拓扑与节点特征）而非纯 Wasserstein，使"保守边/涌现边"具备生物学可解释性。
- **分层可发表性**：L1（可复现 MP 词汇表）、L2（层级锚定 CCC 图）、L3（拓扑偏离评分）任一层都能独立成文，降低高风险层的项目风险。

------------------------------------------------------------------------

## 4. 设计哲学与方法学定位

把*表达层的程序发现*与*图层的拓扑比较***解耦**，是全框架的轴心：

1.  **不在 CCC 推断之前做跨样本表达整合。** CCC 推断（LIANA+ / CellChat）仍在**样本层、基于表达**的 L–R 信号上进行 (Dimitrov et al., 2024; Jin et al., 2021)。我们移除的只是"先把所有细胞投进共享表达空间"这一最易产生 L–R 假阳性的步骤。
2.  **图对齐发生在 per-sample CCC 推断之后**，比较的是**拓扑结构**而非协调表达值——这恰是 OT/FGW 的用武之地 (Vayer et al., 2020)。
3.  **OT/FGW 不替代 CCC 推断内部的批次校正**；它替代的是"跨样本表达整合"这一前置步骤。这一点必须在论文中讲清，否则会被 (R3/R6) 攻击。

> **【D2｜平台不变性由 rank-based 距离承担】** 数据涵盖 Seq-Well、10x 3′/5′ 多版本、BD Rhapsody CITE、MutaSeq/Smart-seq（见 §5）。强度→距离主用 **rank-percentile**（§Phase 5），对深度/平台不变；per-sample 推断 + 不整合表达；并在所有评分-表型关联中把**平台与 study 作为协变量/随机效应**。

------------------------------------------------------------------------

## 5. 数据基础与样本分层（焊入真实清单）

完整 339 样本逐条分层见 `AML_niche_CCC_dataset_inventory.xlsx`（sheet `02_Per_sample_inventory`），研究×角色汇总见 `01_Study_x_Role_summary`，角色定义见 `00_Legend`。下表为**类别级**精简（附录 A 给研究级摘要）。

**五类角色（blueprint role）：**

1.  **Reference-scaffold（投影参考）**：BoneMarrowMap 图谱 (Zeng et al., 2025)（≈263k 健康 BM 细胞，55 态）用于层级投影；CloneTracer/Triana CITE 健康参考作表面标志辅助 (Beneyto-Calabuig et al., 2023; Triana et al., 2021)。**不作发现/验证队列**。
2.  **Healthy-control（健康重心 B_healthy）**：`E-MTAB-11536` (Domínguez Conde et al., 2022)、van Galen 健康 BM (van Galen et al., 2019)、Petti 健康 (Petti et al., 2019)、Lasry 对照 (Lasry et al., 2023)、`GSE253355` (Bandyopadhyay et al., 2024)。
3.  **Discovery-AML（横断面，主队列；B_AML + L1 + L2）**：`GSE239721`（20 例，5′）、Petti AML（5 例，含突变/核型/恶性计数）、`GSE289435` (Zeng2025_BCD，\~12 例 BMMC)、Lasry AML（大量，免疫 receiver 背景丰富但**恶性逐细胞标签弱**）、Chen 2023（6 例，含 niche-immune 富集库）。
4.  **Treatment-axis（主时间轴，Dx→MRD/post-tx）**：van Galen Dx→疗后随访（全 MNC，强 L2）、`GSE185991` (Naldini et al., 2023)（Dx→D14/D30，**CD34/CD117 sorted**→强 L1、弱 L2）、Ennis 的 MRD 时点 (Ennis et al., 2023)、Riether CD34 LSC 治疗 (Riether et al., 2020)。
5.  **Validation: relapse（探索性 RLS）**：`GSE227903` (Ennis et al., 2023)（\~9 例 Dx/(MRD)/R，全 MNC，最适合）、`GSE201966`（移植后单核细胞 AML，Primary→Relapse，\~3 例）。
6.  **Auxiliary-stroma（cohort 层辅助，不进 per-patient 图）**：`GSE253355`（健康股骨头、富集基质）+ Chen 2023 `*_Niche_Immune` 富集库（AML，NPM1，仅诊断期）。
7.  **Comparator / Excluded**：Bailur B-ALL 与肿瘤去除的儿童 AML 免疫 (`GSE154109`) 作特异性对照；PDX、细胞系、annotation/genotype/nanopore 层排除出人 CCC 建模。

**数据现实对设计的三条硬约束：**

- **结构性 niche 几乎只存在于 `GSE253355`（仅健康）+ Chen 2023 niche-immune（仅 AML 诊断期）**。→ 锁定 **D1**：(a) 造血-免疫为主、(b) 结构基质为辅。
- **`malignant_cell_count` 几乎全空、`bam_available` 几乎全 No（仅 FASTQ via SRA）**。→ "≥200 恶性细胞"阈值只能在自建恶性注释后判定；多数样本需**从 FASTQ 重比对建 BAM**（见 **M2**）。
- **平台高度混杂**。→ rank 距离 + per-sample 推断 + 平台协变量（**D2**）。

> **【R1｜审稿："你说 niche 其实没测到结构 niche"】** 应对：术语全程用 "BM microenvironment"；结构基质仅在 `GSE253355`/Chen niche-immune 上做 cohort 层辅助分析并显式标注"受供体数据集构成混杂"；主 claim 不依赖结构 niche。 **【R9｜审稿：sorted/tumor-depleted 数据会扭曲 CCC】** 应对：在分层表里把 **L1-only**（CloneTracer、Riether、Naldini、Chen CD34 库——缺 receiver 区室，仅供恶性程序发现）与 **L2-capable**（含 sender+receiver 的全 MNC 库）严格区分；per-patient CCC 图**只**来自 L2-capable 库。

------------------------------------------------------------------------

## 6. 概念架构：三层框架

| 层级 | 问题 | 主要产出 | 失败时的回退 |
|----|----|----|----|
| **L1 恶性状态-程序层** | 哪些恶性 meta-program 跨队列复现？ | 带 global/conditional 复现性评分（SRRS）、HSI、healthy-match 评分的 MP 词汇表 | 回退到仅 hierarchy-bin 图（不依赖 MP） |
| **L2 CCC 图层** | 恶性层级状态如何与微环境区室通讯？ | 锚定 hierarchy-bin × MP 状态 × 角色的 per-sample CCC 图 | 用稳健的 Main 图（仅 bin×角色） |
| **L3 拓扑偏离层** | 患者拓扑是否偏离健康/AML/治疗压力共识？ | 健康/AML barycenter、emergent edge 集合、患者级 HDS/ATS/RLS/SCS | 用更简单的图距离 baseline 仍可成文 |

------------------------------------------------------------------------

## 7. 方法与分析计划（Phase 1–8）

> 全节写成**待执行的假设检验流程**。所有阈值为预注册的起始值，将在 robustness（Phase 8）中做敏感性分析；非"已得结果"。

### Phase 1 — 数据策划与共识恶性注释

**study 级划分**：按**研究（而非样本）**分层随机为 **Discovery（≈70% 研究）/ Validation（≈30% 研究）**，按遗传亚型分层；全部纵向/复发数据（Ennis、GSE201966、van Galen 疗后、Naldini）保留至 Validation/treatment 轴。

**per-sample QC**：最低 UMI、线粒体% 阈值、doublet 移除（Scrublet, (Wolock et al., 2019)）；样本层阈值 **总细胞 ≥ 500、恶性细胞 ≥ 200**；未达标进 `P2-exclude`。

**恶性注释（三方共识）**：在 BAM/FASTQ 可用样本上取 **inferCNV + Numbat + VarTrix** 共识。inferCNV 以表达相对参考的偏离推断大尺度 CNV (Patel et al., 2014; Tirosh et al., 2016)；Numbat 整合表达、等位基因比与群体单倍型定相（allele counts 由 cellsnp-lite 生成）推断等位基因特异 CNV 并重建克隆系统发生 (Gao et al., 2023; Huang & Huang, 2021)；VarTrix 在已知突变样本（如 Petti、Riether）上提供 SNV 级证据 (10x Genomics, VarTrix)。**判定规则**：≥ 2/3 证据一致即记为 malignant；无已知 SNV 的样本用 CNV+allelic 双证据。

> **【M2｜BAM 几乎不可得 → 流程前置一步重比对】** 由于绝大多数样本只有 FASTQ（via SRA），Phase 1 需先 **STARsolo 重比对 → 带 CB/UB 的 BAM**，再跑 cellsnp-lite → Numbat 与 VarTrix。这与现有单细胞 CNV/SNV 栈一致，但应在计划与算力预算中显式列出（每样本一次重比对）。 **【M3｜恶性标签弱/缺失的数据集的处置】** Lasry（无逐细胞标签、置信低）**不**用于定义 malignant sender，仅作**免疫 receiver 背景**与 B_AML 构成（标注 caution）；肿瘤去除的 Bailur 仅作 comparator。这避免把"未分辨恶性"误当 sender 污染 SCS。 **【R2｜审稿：同一 hierarchy-bin 内如何分恶性 vs 正常？】** 应对：恶性身份由 **CNV/SNV 共识**层叠加在层级标签之上（hierarchy 来自投影、恶性来自基因组证据），二者正交；报告每 bin 的 malignant fraction 与置信度。

### Phase 2 — 层级锚定的细胞状态分配（定义节点词汇）

经**参考投影**到 BoneMarrowMap 健康图谱 (Zeng et al., 2025) 分配 **5–8 个 hierarchy bin**（HSC 样、LMPP/GMP 样、早幼/blast 样、单核样、红系、淋巴系），用 Symphony (Kang et al., 2021) 或 scArches (Lotfollahi et al., 2022)，**避免队列内重整合**。辅助注释：Palantir pseudotime 与熵 (Setty et al., 2019)；多种 stemness 评分（LSC17 (Ng et al., 2016)、van Galen HSC 样签名 (van Galen et al., 2019)）。**LSC 样细胞操作性定义**：位于 HSC/MPP/LMPP bin 且 stemness 评分高的**恶性**细胞。pseudotime 仅作 **bin 内排序**，不作主要分化轨迹。仅可视化整合（scVI/Harmony）允许、但不进下游建模。

> **【R7｜审稿：pseudotime 不可靠】** 应对：pseudotime 仅用于 bin 内精排序；并用 Palantir + 扩散映射多法一致性作敏感性；主结论不依赖单一轨迹。

### Phase 3 — per-sample meta-program 提取与复现性筛选

对每样本恶性细胞跑 **GeneNMF / cNMF**（$K = 4\text{–}9$），**不跨样本整合表达** (Kotliar et al., 2019)；跨样本聚合为 meta-program（MP），方法学沿用泛癌共识程序 (Gavish et al., 2023)。

**双层 SRRS（Study-level Reproducibility Score）**： - **Global SRRS**：MP 在多少比例的**独立研究**中被复现（程序 top-基因集合 Jaccard $\ge \tau$，$\tau$ 预设 0.2，敏感性 0.15–0.30）： $$\mathrm{SRRS}_{\mathrm{global}}(\mathrm{MP}) = \frac{1}{|\mathcal{S}|}\sum_{s\in\mathcal{S}} \mathbb{1}\!\left[\max_{p\in s}\ \mathrm{Jaccard}(\mathrm{MP}, p) \ge \tau\right]$$ 其中 $\mathcal{S}$ 为独立研究集合，$p$ 为研究 $s$ 内的程序。 - **保留判据**：$\mathrm{SRRS}_{\mathrm{global}} \ge 0.6$。

> **【M4｜Conditional SRRS 的阈值必须下调】** 你原框架要求"每个亚型分层 ≥ 5 个独立研究"。按真实清单，除 **NPM1**（Chen、Naldini/GSE185991 等）外，多数亚型（TP53、KMT2A、单核型）**凑不到 5 个独立研究**，该硬门槛不可满足。**修订**：(i) 以 Global SRRS 作硬门槛；(ii) **conditional 复现性仅在 ≥ 3 个独立研究的分层内计算**，否则标注"single-study, 仅注释不设门"；(iii) 用**证据加权**而非硬阈值： $$w_{\mathrm{cond}}(\mathrm{MP}, c) = \frac{n^{\mathrm{recovered}}_{c}}{n^{\mathrm{total}}_{c}}\cdot \log\!\left(1 + n^{\mathrm{studies}}_{c}\right)$$ 每个保留 MP 标注 HSI（hierarchy specificity index）、healthy-match 评分（区分"被挪用"vs"新涌现"程序）、所属分层。

### Phase 4 — CCC 图构建

两种互补图： - **Main 图（稳健、全样本）**：节点 $=$ Hierarchy-bin $\times$ Role（sender/receiver），共 **10–16 节点**；MP 活性作节点特征。 - **Secondary 图（更丰富、按样本可行性）**：节点 $=$ Hierarchy-bin $\times$ 主导 MP 状态 $\times$ Role，共 **30–50 节点**；检验 AML 特异程序态能否细化拓扑。每样本 **\< 15–20 细胞**的节点合并或标 missing。

每样本节点层 pseudobulk 表达 → 用 **LIANA+（aggregated rank consensus）** 推断 L–R 强度矩阵 $S$ (Dimitrov et al., 2024)，以 **CellChat** 作敏感性对照 (Jin et al., 2021; Jin et al., 2024)，必要时用 **NicheNet** 做下游靶基因佐证 (Browaeys et al., 2020)。

> **【R3｜审稿：CCC 假阳性（深度/构成驱动）】** 应对：边进入 barycenter 前要求 **LIANA+ 与 CellChat 双法一致**；rank 距离已抑制深度膨胀；并构造**构成置换零分布**（打乱细胞标签重算 S）以剔除纯丰度驱动的边。 **【R9 续】** per-sample 图只来自 L2-capable 库（含 sender+receiver）。

### Phase 5 — 强度→距离变换

把 $S$ 转为图内距离矩阵 $C$： - **主（rank-based，平台/深度不变）**： $$C_{ij} = 1 - \mathrm{rank}_{\mathrm{pct}}\!\left(S_{ij}\right)$$ - **副（log 归一化，再 min–max 到** $[0,1]$）： $$C_{ij} = -\log\!\left(\tilde{S}_{ij} + \epsilon\right)$$

由 sender/receiver 节点拆分保证对称性；节点质量 $p_i$ 正比于该节点细胞数。两种 $C$ 定义下结果稳定性作内置稳健性检验。

### Phase 6 — Fused Gromov–Wasserstein 对齐与重心

主对齐用 **FGW**（同时编码图拓扑与节点特征：MP 活性向量、细胞丰度、hierarchy 位置）(Vayer et al., 2020)：
$$\mathrm{FGW}_{\alpha}(\mu,\nu) = \min_{T\in\Pi(p,q)}\left[(1-\alpha)\sum_{i,j} D(x_i,y_j)\,T_{ij} + \alpha\sum_{i,j,k,l}\big|C_1[i,k]-C_2[j,l]\big|^2 T_{ij}T_{kl}\right]$$
$\alpha$ 在 $\{0,0.25,0.5,0.75,1\}$ 扫描，主设 $\alpha=0.5$；$\alpha=0$（纯特征 OT）、$\alpha=1$（纯 GW）为端点。节点缺失时用 **unbalanced FGW**。健康/AML barycenter 在 Main 图（固定节点词汇）上算，用 GW/FGW 重心算法 (Peyré et al., 2016) + Sinkhorn 正则 (Cuturi, 2013)。FGW 与更简单距离（Frobenius、cosine、纯 Wasserstein on edge mass）做 benchmark，证明保守边排序稳健、且 FGW 更可解释。

> **【R4｜审稿：barycenter 在强异质下不可解释/被少数样本主导】** 应对：固定节点词汇 + unbalanced FGW + 样本级 bootstrap CI（Phase 8）+ 置换零分布 + 与简单距离 benchmark；并报告每条 barycenter 边的样本支持度（多少样本贡献质量）。

### Phase 7 — Emergent edge 识别与患者级评分

emergent edge 三法互证： - **逐元素差（共享节点词汇下为主）**：$\Delta C[i,j] = C_{\mathrm{AML}}[i,j] - C_{\mathrm{Healthy}}[i,j]$，保留 top-$k$。 - **耦合质量偏移**：算两重心间 FGW 耦合 $T^{*}$，标记非对角线质量偏移最大的边为"重塑边"。 - **置换零分布**：随机打乱 Healthy/AML 标签重算重心（$n \ge 1000$），检验真实 emergent 边显著性。

**患者级偏离评分**（GW/FGW 距离）： $$\mathrm{HDS} = \mathrm{GW}\!\left(G_{\text{patient}}, B_{\text{healthy}}\right), \qquad \mathrm{ATS} = \mathrm{GW}\!\left(G_{\text{patient}}, B_{\text{AML}}\right)$$ $$\mathrm{RLS} = \mathrm{GW}\!\left(G_{\text{patient}}, B_{\text{Dx}}\right) - \mathrm{GW}\!\left(G_{\text{patient}}, B_{\text{relapse}}\right) \quad(\text{仅 Validation 配对队列})$$ $$\mathrm{SCS} = \sum_{(i,j)\in \mathcal{E}_{\mathrm{LSC}\to\mathrm{env}}} T^{*}_{ij}\, w_{ij} \quad(\text{限定 LSC 样 sender}\to\text{免疫/旁分泌 receiver 边的质量加权耦合})$$

> **【M5｜SCS 改义】** 因 niche 重定义为免疫/旁分泌（D1），SCS 由"干性-niche"改为 **"干性-微环境通讯评分"**：LSC 样 sender → 免疫/旁分泌 receiver 边的质量加权耦合，名实相符。 **【R5｜审稿：拓扑评分只是 blast 比例的代理】** 应对：所有评分-表型关联中**回归掉 blast/恶性比例与细胞构成**（偏相关 / 含 blast-fraction 协变量的混合模型），并构造**构成匹配的零分布**；H2/H3 的判别力须在"构成匹配"后仍成立。

### Phase 8 — 稳健性与交叉验证

- **拓扑扰动**：随机边丢弃 30% + 节点遮蔽 20%，100 次 bootstrap；barycenter 高质量边须 **≥ 70% 恢复率**。
- **样本级 bootstrap**：80% 样本子集重算重心（$n=100$），输出边级 CI。
- **外部验证**：Validation 重心 vs Discovery 重心的 FGW 距离；emergent 边独立复现检验。
- **平台/批次稳健性**：把平台与 study 作随机效应；rank vs log 两种 $C$ 一致性。
- **空间交叉检验（可选）**：对 top emergent L–R 边，用公开骨髓空间数据（10x Xenium/Visium）验证对应状态的物理邻近。

> **【R6｜审稿：你测的是批次不是生物学】** 应对：study 级 train/val 划分 + 跨队列复现 + 平台协变量 + rank 不变性 + 置换零分布，五重防线。 **【R8｜审稿：cNMF 在低细胞数不稳】** 应对：$K$ 跨 4–9 取共识 + 恶性细胞下限 ≥ 200 + 回退到仅 hierarchy 的 Main 图（不依赖 MP）。

------------------------------------------------------------------------

## 9. 预期产出与可证伪评估

> 一律写成"将检验 / 判据"，不写"将证明"。

1.  **AML 特异、状态锚定的 MP 词汇表**（带 Global SRRS、conditional 证据权重、HSI、healthy-match）。 *评估*：H1 —— 保守轴数量与跨亚型富集须显著超过置换零分布；被证伪若可由平台解释。
2.  **健康/AML 共识 barycenter + emergent edge 集合**（经扰动测试验证为候选不变边）。 *评估*：emergent 边须在置换零分布下显著、且在外部队列复现（Phase 8）。
3.  **患者级拓扑偏离评分（HDS/ATS/RLS/SCS）**及其与亚型/治疗反应/结局的关联（metadata 允许时）。 *评估*：H2 —— 拓扑模型 vs 单 L–R / 仅构成模型的判别力比较（嵌套交叉验证、构成匹配后）；H3 —— 配对患者内 **Dx\<MRD\<R** 单调性（Jonckheere–Terpstra 趋势检验），并在回归掉 blast 比例后仍成立。

------------------------------------------------------------------------

## 10. 验证策略（计算发现为主线 + 可选实验）

主线为**计算发现**；最后提出少量**可行**实验以检验候选轴因果性（有条件与时间则 A+B+C，但**不承诺具体实验必须完成**）：

- **A. 组织/空间层面**：IHC / IF / RNAscope，验证 top L–R 轴对应细胞的**邻近与表达定位**（与 Phase 8 空间交叉检验衔接）。
- **B. 细胞功能层面**：共培养、colony assay、流式，检验候选微环境轴是否影响 **LSC 样表型**（分化、克隆形成、存活）。
- **C. 机制扰动层面**：CRISPRi / siRNA / 阻断抗体 ± **perturbation scRNA-seq**，检验该轴是否**驱动**网络重排（扰动后重算 CCC 图与 HDS/ATS，看拓扑是否回移）。

> *可证伪表述*：若敲低/阻断候选轴后 LSC 样表型与网络拓扑偏离均无显著变化，则该轴被判为"相关但非驱动"，从中核模块剔除。

------------------------------------------------------------------------

## 11. 风险与应对（汇总）

| 风险 | 应对 |
|----|----|
| 程序发现不稳（低 $n$） | $K$ 共识、恶性细胞下限、回退 hierarchy-only 图（R8） |
| CCC 假阳性 | LIANA+/CellChat 双法一致 + 构成置换零分布（R3） |
| 稀疏亚型分层 | conditional 改为 ≥3 研究 + 证据加权（M4） |
| 重心不可识别/被主导 | unbalanced FGW + bootstrap CI + 置换 + 简单距离 benchmark（R4） |
| 拓扑=blast 比例混杂 | 回归掉构成 + 构成匹配零分布（R5） |
| 批次 vs 生物学 | study 级划分 + 跨队列复现 + 平台协变量（R6） |
| 复发数据薄 | 复发降级为探索性 RLS；主轴用治疗压力（M1） |
| 结构 niche 缺失 | niche 重定义 + 结构基质仅 cohort 层辅助（D1, R1） |
| BAM 缺失 | FASTQ→STARsolo 重比对前置（M2） |
| 分层失败保底 | L1/L2/L3 各自可独立成文 |

------------------------------------------------------------------------

## 12. 时间线（示意）

| 时段 | 里程碑 | 层 |
|----|----|----|
| 第 1–6 月 | 策划、QC、FASTQ→BAM 重比对、三方恶性注释、层级投影分箱 | Phase 1–2（L1 基础） |
| 第 6–14 月 | per-sample cNMF、MP 聚合、SRRS 筛选；MP 词汇表 v1 | Phase 3（L1） |
| 第 12–22 月 | CCC 图、强度→距离、FGW barycenter、emergent 边 | Phase 4–7（L2–L3） |
| 第 20–30 月 | 患者级评分、治疗压力轴 H3、外部验证、（可选）空间交叉、少量功能/扰动实验 | Phase 7–8 + 验证 |

------------------------------------------------------------------------

## 附录 A — 研究级数据集分层（精简）

> 完整逐样本明细见 `AML_niche_CCC_dataset_inventory.xlsx`。下表只列研究级角色与关键限制。

| 研究 / Accession | 角色 | 时间轴 | niche 捕获 | 关键限制 |
|----|----|----|----|----|
| Zeng 2025 BoneMarrowMap | Reference-scaffold | 参考 | 造血图谱 | 仅投影参考，非队列 |
| Domínguez Conde 2022 `E-MTAB-11536` | Healthy-control | 参考 | 免疫为主 | 仅免疫，无恶性 |
| Petti 2019 `Zenodo 3345981` | Healthy + Discovery-AML(突变锚定) | 横断面 | 无基质 | 有突变/计数；5 例 AML |
| `GSE239721`（IFNγ-AML） | Discovery-AML | 横断面 | 造血/免疫 | 20 例诊断；无纵向/基质 |
| Chen 2023 `gwjh3w6ztm` | Discovery-AML + Auxiliary-stroma | 横断面 | niche+免疫富集 | 6 例 NPM1；仅诊断；唯一 AML 富 niche |
| Ennis 2023 `GSE227903` | Validation: 治疗轴+复发 | 治疗压力/复发 | 造血/免疫(基质有限) | \~9 例 Dx/(MRD)/R 全 MNC，最适合 L2 |
| van Galen 2019 `GSE116256` | Discovery + 治疗轴 | 治疗压力 | 造血(基质有限) | Seq-Well；Dx→疗后随访 |
| Beneyto-Calabuig 2023 CloneTracer | Reference + Discovery(克隆/HSPC) | 参考/横断面 | sorted-CD34 | CD34 only → L1，无 receiver |
| Naldini 2023 `GSE185991` | Treatment-axis | 治疗压力(少量 REL) | sorted-blast/progenitor | CD34/CD117 sorted → 强 L1 弱 L2 |
| `GSE201966`（单核细胞 AML 复发） | Validation: 复发 | 复发 | BMMC 无基质 | Primary→Relapse \~3 例 |
| Riether 2020 `GSE147989` | Treatment-axis(CD34 PB) | 治疗压力 | sorted-CD34 | LSC 治疗反应；无全骨髓 CCC |
| Bandyopadhyay 2024 `GSE253355` | Reference + Auxiliary-stroma | 参考 | **基质富集** | 仅健康；结构 niche 参考 |
| Lasry 2023 `GSE185381` | Discovery-AML(免疫 receiver) | 横断面 | 免疫为主 | 大队列；恶性逐细胞标签弱 |
| Zeng 2025 `GSE289435` | Discovery-AML(BMMC) | 横断面 | 无基质 | \~12 例；PDX 排除 |
| Nicosia 2023 `GSE207356` | Exploratory | 治疗压力 | 血液无基质 | 单患者；仅生成假设 |
| Bailur 2020 `GSE154109` | Comparator | — | 肿瘤去除/免疫 | B-ALL + 儿童 AML，对照用 |

------------------------------------------------------------------------

## 参考文献（正文 Author–Year；下列为编号完整版，附 DOI/PMID 或数据 accession）

1.  Zeng AGX, Iacobucci I, Shah S, …, Mullighan CG, Dick JE. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of aberrant AML differentiation (BoneMarrowMap). *Blood Cancer Discov*. 2025;6(4):307–324. <doi:10.1158/2643-3230.BCD-24-0342>.
2.  van Galen P, Hovestadt V, Wadsworth MH II, …, Bernstein BE. Single-cell RNA-seq reveals AML hierarchies relevant to disease progression and immunity. *Cell*. 2019;176(6):1265–1281.e24. <doi:10.1016/j.cell.2019.01.031>. <PMID:30827681>. [数据：GSE116256]
3.  Ng SWK, Mitchell A, Kennedy JA, …, Wang JCY. A 17-gene stemness score for rapid determination of risk in acute leukaemia (LSC17). *Nature*. 2016;540:433–437. <doi:10.1038/nature20598>. <PMID:27926740>.
4.  Gavish A, Tyler M, Greenwald AC, …, Tirosh I. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. *Nature*. 2023;618(7965):598–606. <doi:10.1038/s41586-023-06130-4>.
5.  Kotliar D, Veres A, Nagy MA, …, Sabeti PC. Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq (cNMF). *eLife*. 2019;8:e43803. <doi:10.7554/eLife.43803>.
6.  Setty M, Kiseliovas V, Levine J, Gayoso A, Mazutis L, Pe'er D. Characterization of cell fate probabilities in single-cell data with Palantir. *Nat Biotechnol*. 2019;37(4):451–460. <doi:10.1038/s41587-019-0068-4>. <PMID:30899105>.
7.  Jin S, Guerrero-Juarez CF, Zhang L, …, Nie Q. Inference and analysis of cell–cell communication using CellChat. *Nat Commun*. 2021;12:1088. <doi:10.1038/s41467-021-21246-9>.
8.  Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell–cell communication from single-cell and spatially resolved transcriptomics (CellChat v2). *Nat Protoc*. 2024. <doi:10.1038/s41596-024-01045-4>.
9.  Dimitrov D, Schäfer PSL, Farr E, …, Saez-Rodriguez J. LIANA+ provides an all-in-one framework for cell–cell communication inference. *Nat Cell Biol*. 2024;26(9):1613–1622. <doi:10.1038/s41556-024-01469-w>.
10. Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nat Methods*. 2020;17:159–162. <doi:10.1038/s41592-019-0667-5>.
11. Gao T, Soldatov R, Sarkar H, …, Kharchenko PV. Haplotype-aware analysis of somatic copy number variations from single-cell transcriptomes (Numbat). *Nat Biotechnol*. 2023;41:417–426. <doi:10.1038/s41587-022-01468-y>.
12. Patel AP, Tirosh I, Trombetta JJ, …, Bernstein BE. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. *Science*. 2014;344(6190):1396–1401. <doi:10.1126/science.1254257>.
13. Tirosh I, Izar B, Prakadan SM, …, Garraway LA. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq (CNV-from-expression basis of inferCNV). *Science*. 2016;352(6282):189–196. <doi:10.1126/science.aad0501>.
14. Huang X, Huang Y. Cellsnp-lite: an efficient tool for genotyping single cells. *Bioinformatics*. 2021;37(23):4569–4571. <doi:10.1093/bioinformatics/btab358>.
15. 10x Genomics. VarTrix: extracting single-cell variant information from 10x data (software). <https://github.com/10XGenomics/vartrix>.
16. Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Syst*. 2019;8(4):281–291.e9. <doi:10.1016/j.cels.2018.11.005>.
17. Kang JB, Nathan A, Weinand K, …, Raychaudhuri S. Efficient and precise single-cell reference atlas mapping with Symphony. *Nat Commun*. 2021;12:5890. <doi:10.1038/s41467-021-25957-x>.
18. Lotfollahi M, Naghipourfar M, Luecken MD, …, Theis FJ. Mapping single-cell data to reference atlases by transfer learning (scArches). *Nat Biotechnol*. 2022;40:121–130. <doi:10.1038/s41587-021-01001-7>.
19. Nagai JS, Maié T, Schaub MT, Costa IG. scACCorDiON: a clustering approach for explainable patient-level cell–cell communication graph analysis. *Bioinformatics*. 2025;41(5):btaf288. <doi:10.1093/bioinformatics/btaf288>.
20. Vayer T, Chapel L, Flamary R, Tavenard R, Courty N. Fused Gromov-Wasserstein distance for structured objects. *Algorithms*. 2020;13(9):212. <doi:10.3390/a13090212>.
21. Peyré G, Cuturi M, Solomon J. Gromov-Wasserstein averaging of kernel and distance matrices. *Proc. 33rd Int. Conf. Machine Learning (ICML)*, PMLR. 2016;48:2664–2672.
22. Cuturi M. Sinkhorn distances: lightspeed computation of optimal transport. *Adv. Neural Inf. Process. Syst. (NeurIPS)*. 2013;26:2292–2300.
23. Kunisaki Y, Bruns I, Scheiermann C, …, Frenette PS. Arteriolar niches maintain haematopoietic stem cell quiescence. *Nature*. 2013;502(7473):637–643. <doi:10.1038/nature12612>. <PMID:24107994>.
24. Ennis S, Conforte A, O'Reilly E, …, Szegezdi E. Cell-cell interactome of the hematopoietic niche and its changes in acute myeloid leukemia. *iScience*. 2023;26(6):106943. <doi:10.1016/j.isci.2023.106943>. <PMID:37332612>. [数据：GSE227903]
25. Chen L, Pronk E, van Dijk C, …, Raaijmakers MHGP. A single-cell taxonomy predicts inflammatory niche remodeling to drive tissue failure and outcome in human AML. *Blood Cancer Discov*. 2023;4(5):394–417. [数据：Mendeley 10.17632/gwjh3w6ztm.2]
26. Naldini MM, Casirati G, Barcella M, …, Genovese P. Longitudinal single-cell profiling of chemotherapy response in acute myeloid leukemia. *Nat Commun*. 2023;14:1285. <doi:10.1038/s41467-023-36969-0>. [数据：GSE185991]
27. Beneyto-Calabuig S, Merbach AK, Kniffka JA, …, Velten L. Clonally resolved single-cell multi-omics identifies routes of cellular differentiation in acute myeloid leukemia (CloneTracer). *Cell Stem Cell*. 2023;30(5):706–721. <doi:10.1016/j.stem.2023.04.001>.
28. Petti AA, Williams SR, Miller CA, …, Ley TJ. A general approach for detecting expressed mutations in AML cells using single-cell RNA sequencing. *Nat Commun*. 2019;10:3660. <doi:10.1038/s41467-019-11591-1>. [数据：Zenodo 10.5281/zenodo.3345981; phs000159]
29. Lasry A, Nadorp B, Fornerod M, …, Aifantis I. An inflammatory state remodels the immune microenvironment and improves risk stratification in acute myeloid leukemia. *Nat Cancer*. 2023;4(1):27–42. <doi:10.1038/s43018-022-00480-0>. [数据：GSE185381]
30. Bandyopadhyay S, Duffy MP, Ahn KJ, …, Tan K. Mapping the cellular biogeography of human bone marrow niches using single-cell transcriptomics and proteomic imaging. *Cell*. 2024;187(12):3120–3140.e29. <doi:10.1016/j.cell.2024.04.013>. [数据：GSE253355]
31. Domínguez Conde C, Xu C, Jarvis LB, …, Teichmann SA. Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*. 2022;376(6594):eabl5197. <doi:10.1126/science.abl5197>. [数据：E-MTAB-11536]
32. Triana S, Vonficht D, Jopp-Saile L, …, Velten L. Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states. *Nat Immunol*. 2021;22(12):1577–1589. <doi:10.1038/s41590-021-01059-0>.
33. Riether C, Pabst T, Höpner S, …, Ochsenbein AF. Targeting CD70 with cusatuzumab eliminates acute myeloid leukemia stem cells in patients (单细胞 LSC 治疗反应). *Nat Med*. 2020;26(9):1459–1467. <doi:10.1038/s41591-020-0910-8>. [数据：GSE147989]
34. [数据集，无独立期刊论文] Single-cell transcriptome landscape of relapse fingerprint in acute monocytic leukemia post-transplantation. GEO accession **GSE201966**.
35. [数据集] Nicosia L, et al. scRNA-seq of CCS1477 (inobrodib) therapy time course in leukaemia. GEO accession **GSE207356**.
36. [数据集] Bailur JK, et al. Single-cell profiling of pediatric leukemia bone marrow T cells. *JCI Insight*. 2020. GEO accession **GSE154109**.
37. [数据集] IFNγ signaling characterization in AML, single-cell. GEO accession **GSE239721**.

> 说明：第 34–37 项以**期刊（如有）+ 数据 accession** 作硬标识；其 DOI 未经我逐条核验，故不臆造，引用时以 accession 为准、需要时可补全 DOI。其余 1–33 项 DOI/PMID 均已核验。
