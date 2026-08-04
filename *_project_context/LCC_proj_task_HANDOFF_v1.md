# HANDOFF — TP53-adjacent PI side题 (PVRL4/Nectin-4 与骨髓纤维化)

> **怎么用这份文件 / How to use this**
> 这是 `aml_niche_net` 主线之外的一个**机会性分支**——老板实验室两个 TP53 AML 相关题目，
> 评估能不能用现有单细胞数据掺一脚。新开对话时把这份文件发给 Claude 即可，不需要一并带
> `HANDOFF_phase2.md`（除非要交叉引用主线 Phase 1b/2 的产出路径，那部分本文件已摘录了必要的）。
> Working style: Chinese chat, English code comments/paths, co-designer mode。

---

## 0. 背景与老板的两个题目

老板实验室的两个意外发现，想知道能不能用 `aml_niche_net` 现有的 13 个数据集帮上忙：

**题目 A（PVRL4/Nectin-4）**：老板团队在自己的 bulk RNA-seq + cell line 里发现 PVRL4（Nectin-4）
在血液里本应低表达，但部分病人**病理切片 IHC（蛋白层级）染色很强**——意外发现。想知道：
TP53 突变的 AML 相较其他 AML，PVRL4 表达是否更高？

**题目 B（TP53 → 骨髓纤维化）**：老板团队发现带 TP53 突变的 AML/MDS 病人，形成骨髓纤维化的
比率更高。想知道：这类病人是否有更高的纤维化相关基因表达，或与间质细胞有更明显的纤维化 interaction？

---

## 1. 结论（截至本文件撰写时仍然成立，无新数据推翻）

### 题目 A — PVRL4：**可以做，但要换个问法**

**硬伤**：13 个数据集里，被标注 TP53 突变的样本目前只有 **1 个**（Petti 809653，complex
karyotype + del(17q)）。有逐样本突变注释的数据集只有 4 个（GSE185991、Chen2023、Petti2019、
GSE147989），加起来 TP53 也只此一例；有核型信息的样本总共 12 个。**「TP53 AML vs 其他 AML」
的直接分组比较，现有公开数据做不出来（n=1）。**

**但两条路能救活**：

- **(a) 用 CNV 流程自造 TP53 代理组**：TP53 突变高度富集 17p 缺失 / complex karyotype。
  Phase 1b 的 inferCNV + Numbat 本来就在逐细胞调 CNV，17p13.1（TP53 位点）丢失是能直接读出来的。
  可以从恶性注释结果里定义一个「17p-loss / complex karyotype」群，作为 TP53-aberrant 代理组，
  再比 PVRL4 表达。**Caveat**：17p loss ≠ TP53 突变，是富集代理，不是等价，写论文/报告时必须
  明确这一点。
- **(b) 单细胞的真正增值点：PVRL4 的细胞来源**。老板的 bulk RNA-seq + IHC 回答不了「是 blast
  表达、还是某个 blast 亚群、还是基质/免疫细胞」——但单细胞能。若能证明 PVRL4 富集在 LSC 样 /
  恶性细胞而非旁观细胞，对一个已知 ADC 靶点（Nectin-4 是 enfortumab vedotin 的靶）有直接临床
  意义，比单纯"TP53 组更高"更有价值。

**必须最先做的事（5 分钟体检，任何下游分析之前）**：PVRL4 血液本底表达低，10x 3′ dropout 会很
狠——先确认它在现有对象里到底捕没捕到、哪些细胞有非零计数。若几乎全 0，后面的分析都无意义，
必须趁早知道。

### 题目 B — 纤维化：**现有数据基本做不了，建议不接或明确告知老板数据缺口**

死结是蓝图 D1 决策点名的问题：**公开 AML 抽吸数据几乎不含结构性间质 niche**。

- 真正 stroma-enriched 的样本只有 `GSE253355` 的 12 个——**全是健康骨髓**。
- 有 AML 间质的只有 Chen2023 的 niche-immune 库（14 个，**全是 NPM1、全是诊断期**）。
- 纤维化基因（COL1A1/COL3A1/FN1 等）由间质/成纤维细胞产生，blast 不产——没捕到 MSC 就看不到
  这些基因，也算不出与间质的纤维化 interaction。
- **TP53 ∩ 间质捕获 = 0**：唯一有 AML 间质的 Chen2023 是 NPM1，一个 TP53 都没有。

抽吸物天然冲掉间质，这题需要骨芯活检 / 带骨片样本或空间数据 + TP53 状态，公开数据里基本没有
成规模的。勉强做只能在 Chen2023 的健康/NPM1 间质里看纤维化基因，**完全碰不到 TP53 这个核心
假设**。建议直接跟老板说清楚：要做这题得用他自己的活检/空间数据，不能用公开 scRNA 抽吸数据
硬凑，否则是过度声称（over-claiming）。

---

## 2. 如果要参与，题目 A 的具体做法

### 2.1 定位：机会性支线，不阻塞主线

这是 `aml_niche_net` 主线（Phase 1b 恶性注释 → Phase 2 hierarchy bin → Phase 3+ CCC/FGW）之外
的**旁支**。不占用主线关键路径的算力/时间预算，可以在主线 array job 挂后台跑的间隙推进。
Step 0（体检）现在就能做，不依赖任何主线进度；Step 2/3（TP53 代理组）依赖 Phase 1b 恶性注释
的产出（进行中）；细胞来源问题（1(b)）最好等 Phase 2 hierarchy bin 做完再叠加，但也可以先用
粗粒度的 malignant/non-malignant 二分做初版。

### 2.2 Step 0 — PVRL4/Nectin 家族检出体检（现在就能做，不依赖任何流程进度）

直接在 QC RDS（`/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/01_per_sample_qc/<dataset>/<sample>.rds`）
上跑一个轻量脚本：

- 对每个数据集抽几个有代表性的样本（不用全跑），检查 `PVRL4` 基因是否在 counts matrix 的
  rownames 里（不同参考基因组/注释版本，symbol 可能是 `PVRL4` 或 `NECTIN4`，两个都要试）。
- 计算：全局非零细胞比例、pseudobulk 平均表达、按（粗粒度）细胞类型/cluster 拆分的非零比例。
- 顺手看一眼 Nectin 家族其他成员（PVRL1/NECTIN1, PVRL2/NECTIN2, PVRL3/NECTIN3）作为背景对照，
  确认不是基因名/参考版本导致的系统性丢失。
- 输出一张小表：dataset × sample → detection rate，一眼看出这题还能不能做。

若检出率普遍趋近 0，题目 A 基本判死刑，直接如实回报给老板；若有可观测的非零信号（哪怕只在
少数样本/少数细胞类型），再往下走 Step 1/2。

### 2.3 Step 1 — PVRL4 的细胞来源（1(b)，价值更高的部分）

依赖：Phase 1b 恶性标签（进行中，见下方"依赖主线的产出接口"）；理想情况下再叠加 Phase 2
hierarchy bin（未开始）区分恶性细胞属于哪个分化层级。

- 先用粗粒度二分（恶性 vs 非恶性，来自 Phase 1b c50 consensus）：PVRL4+ 细胞是否显著富集在
  恶性细胞里？
- 若 Phase 2 hierarchy bin 就绪，再看 PVRL4+ 的恶性细胞集中在哪个 bin（LSC 样/HSC-MPP-like
  还是更分化的 blast 样）。
- 这一步的产出即使 Step 2（TP53 代理组）做不出显著性，本身也是独立可报的发现。

### 2.4 Step 2 — 用 CNV 流程自造 TP53 代理组（17p-loss / complex karyotype）

依赖：Phase 1b 的 Numbat / inferCNV 逐细胞输出（进行中）。

**待与 Claude 在新对话里敲定的设计点（D# 风格，一次一个问题）**：

- **代理组的判定粒度**：按**样本**（该样本恶性细胞群体是否携带 17p 缺失克隆）还是按**细胞**
  （单细胞级别的 17p 状态，允许同一样本内异质）？前者噪声小、样本量小；后者反过来。建议先按
  样本（更稳），细胞级作为敏感性分析。
- **17p 缺失的读出来源**：
  - Numbat：`clone_post`/`segs_consensus_*.tsv`/`bulk_clones_final.tsv.gz` 里应该有按染色体
    臂/区段的拷贝数状态（此前验证 6323_R 时见过 `GT_opt` 形如 `7d,11e,19a,2b` 这类基因型字符串，
    需要确认其中是否覆盖 17p、以及具体格式怎么解析出"17p 缺失"这个判定）。
  - inferCNV：HMM 状态矩阵（`run.final.infercnv_obj`）按基因组区段给出拷贝数状态，17p13.1
    区域的基因子集可以直接读取。
  - 这两个源头目前**都还没针对"提取 17p 状态"这个具体查询验证过格式**，是新对话要做的第一件
    技术活。
- **complex karyotype 的操作性定义**：如果只用 17p 太窄（假阴性多），可以放宽到"多条染色体臂
  异常的克隆"（呼应蓝图里 TP53 常见的 complex karyotype 富集），但要先定一个可复现的量化阈值
  （比如"≥3 条染色体臂有 CNV 事件"），而不是主观判断。
- **代理组 vs 真实 TP53 突变样本（n=1, Petti 809653）的一致性检查**：如果流程正确，这唯一一个
  真阳性样本应该被代理组规则命中——这是最起码的 sanity check，做之前先把这一步跑通。

### 2.5 Step 3 — 比较 PVRL4 表达：代理-TP53 组 vs 其他

- 需要控制混杂变量：数据集/平台（study 作协变量，呼应蓝图 D2 的平台不变性设计）、恶性细胞比例、
  样本层 UMI 深度。
- 由于代理组是弱信号（CNV 富集代理，不是突变本身），统计功效会比真实分组低很多——预期效应量
  上要打折扣，报告时必须清楚写明这是 proxy-based 而非 genotype-confirmed 的比较。
- 输出建议：pseudobulk 层面的组间比较（更稳） + 单细胞层面的分布可视化（辅助、非主要判据）。

---

## 3. 数据与路径速查（沿用主线的存储约定）

| 内容 | 路径 |
|---|---|
| 数据集清单（含 key_mutations/cytogenetics/subtype 列） | `/mnt/project/AML_niche_CCC_dataset_inventory.xlsx`（sheet `02_Per_sample_inventory`）|
| 每样本 QC RDS（Step 0 直接读这个） | `/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/01_per_sample_qc/<dataset>/<sample>.rds` |
| Phase 1b 恶性共识标签（c50 输出，Step 1/3 要用） | `/FAST/gr10634/gaozy/aml_niche_net/results/tables/03_malignancy/<dataset>/` |
| Numbat 逐样本输出（Step 2 的候选来源之一） | `/LARGE1/gr10634/gaozy/aml_niche_net/03_cnv_snv/numbat/<dataset>/<sample>/numbat/` |
| inferCNV 逐样本输出（Step 2 的候选来源之二） | `/LARGE1/gr10634/gaozy/aml_niche_net/05_cnv_snv/infercnv/<dataset>/<sample>/` |
| inferCNV burden CSV（只有全基因组总分，无区段级信息，Step 2 用不上，仅供参考） | `/LARGE1/gr10634/gaozy/aml_niche_net/05_cnv_snv/infercnv_burden/<dataset>/<sample>_infercnv_burden.csv` |

已知的唯一真阳性 TP53 样本：**Petti 2019，样本 809653**（`key_mutations` = "TP53 E286G, CEBPA
R142fs, NRAS G12D, FLT3-ITD negative"；`cytogenetics` = complex karyotype with del(5q),-7,+8,
del(12p),**del(17q)**,t(1-15),der(22)t(1-22)）。

Chen2023 niche-immune 富集库（题目 B 唯一有 AML 间质捕获的来源，但全 NPM1、全诊断期，对题目 B
的 TP53 假设无法验证，仅作背景参考）。

---

## 4. 依赖主线的产出接口（供交叉引用，非本文件重点）

- Phase 1b 恶性共识标签：schema `cell, malignant(0/1), score, method, sample` + per-sample
  `evidence_tier` (A/B/C)，见主线 `c50_consensus_malignancy.R` 输出。GSE227903 目前是 Tier A
  锚点数据集（进行中）。
- Phase 2（未开始）：BoneMarrowMap 投影 + hierarchy bin，输出 `cell, hierarchy_bin,
  projection_confidence` + per-sample × per-bin 恶性比例汇总表。Step 1 的精细版分析要等这个。

如果这次新对话里发现需要主线的实时进度（比如 Phase 1b 具体跑到哪个样本了），把
`HANDOFF_phase2.md` 一并带过去；如果只是做题目 A 的探索性分析，本文件应该够用。

---

## 5. 下一步建议（新对话开场可以直接问）

1. 先做 Step 0 体检脚本（PVRL4/NECTIN4 + 家族其他成员的检出率表），这个不依赖任何其他进度，
   可以立刻出结果。
2. 体检结果如果可行，再决定是先做 Step 1（细胞来源，价值更高但依赖 Phase 1b/2）还是先做
   Step 2（TP53 代理组定义，技术上要先摸清 Numbat/inferCNV 的区段级输出格式）。
3. 题目 B：建议这次新对话里先把"如何跟老板措辞说明数据缺口"过一遍，而不是硬做分析。
