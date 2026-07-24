# aml_niche_net — Blueprint v1.1 Patch
 
> **文档定位**：本文件是对 `AML_niche_CCC_OT_blueprint_zh.md` (v1.0) 的**增量修订补丁**，不是替代品。
> 每一节标注 `[REPLACE §x]` / `[INSERT AFTER §x]` / `[AMEND §x]`，可直接套用到主稿。
> 生成日期：2026-07-22。触发原因：2026 年 4–7 月发表的 4 篇关键工作 + 1 篇直接竞争预印本。
>
> **本版新增编号**：`M6–M9`（设计修订）、`R10–R12`（审稿风险）、`D3–D4`（设计决策）。
> v1.0 的 D1/D2、M1–M5、R1–R9 全部保留有效，除下文明确标注修订者外不做改动。
>
> **执行影响**：本补丁的全部改动位于 Phase 5–8（`06_distance` / `07_fgw` / `08_scoring`）与文字层，
> **不需要重跑 CellChat**。159 样本 `05_ccc` array 可照常提交，唯一前置动作见 §12 QC-0。
 
---
 
## 0. 文献核实状态表（引用前必读）
 
本补丁引用的新文献，按核实程度分三级。**Tier-1 已逐段读过正文**，可直接作为承重引用；**Tier-2 仅核实标题/摘要/DOI**；**Tier-3 未独立核实，引用前必须自行确认**。
 
| # | 文献 | DOI / ID | 级别 | 用途 |
|---|---|---|---|---|
| N1 | Wang C, Pan Y, …, Li L, Kang X. *Nat Cell Biol* 2026;28(5):890–902 | 10.1038/s41556-026-01939-3 | **Tier-1** | R10 / R12 / M9 |
| N2 | Ochi Y, …, Lehmann S, Ogawa S (eCHROMA). *Nature* 2026, online 07-08 | 10.1038/s41586-026-10703-4 | **Tier-1** | R11 / H1 分层 / D4 |
| N3 | Struyf N, …, Erkers T. *bioRxiv* 2026-07-01 (preprint) | 10.64898/2026.07.01.735780 | **Tier-1** | §3 竞争定位 / H1 佐证 |
| N4 | Li Z, …, Yu J, Caligiuri MA. *Nat Commun* 2026 | 10.1038/s41467-026-68582-2 | Tier-2 | M7 功能轴锚定 |
| N5 | CD81 supports LSC function. *Signal Transduct Target Ther* 2026 | 10.1038/s41392-026-02697-2 | Tier-2 | M8 节点特征 |
| N6 | STCase. *Nat Comput Sci* 2025 | 10.1038/s43588-025-00809-6 | Tier-2 | Related work |
| N7 | SCILD. *Commun Biol* 2026 | 10.1038/s42003-025-09413-w | Tier-2 | Related work |
| N8 | COMMOT. *Nat Methods* 2023;20:218–228 | 10.1038/s41592-022-01728-4 | Tier-2 | Related work |
| N9 | miR-126 vascular niche ABM. *npj Syst Biol Appl* 2026 | 10.1038/s41540-026-00675-6 | **Tier-3** | Discussion only |
| N10 | CellRefiner. *Nat Commun* 2026 | 10.1038/s41467-026-70090-2 | **Tier-3** | Limitation only |
| N11 | SUCNR1/succinate. *Nat Commun* 2026 | (accession 待补) | **Tier-3** | Limitation only |
| N12 | Extramedullary/lung AML niche. *Nat Immunol* 2026 (N&V) | 10.1038/s41590-026-02598-0 | **Tier-3** | 不引用为承重证据 |
 
> **待核对数据点**：N5 报告的 HR。你手头记的是 RFS 2.44 / OS 2.31；早期版本（*Blood* 2023 摘要）记录为 **RFS HR 2.01 (1.16–3.48)、OS HR 2.38 (1.32–4.29)**，且限定 **non-CBF AML**。写进正文前请核对 STTT 终稿数字与队列限定。
 
---
 
## 1. `[AMEND §0]` 中心发现命题 — 术语与范围收窄
 
**保留**原命题三要素（收敛性 / 拓扑性 / 治疗压力轴），但作两处收窄：
 
**(a) 全文 niche 术语统一为 `aspirate-accessible bone marrow microenvironment`。**
v1.0 的 D1 已把 niche 重定义为造血-免疫+旁分泌。N1 证明真正决定 LSC 命运的轴之一是**纵向解剖位置**（metaphysis vs central marrow），其关键细胞是 N-cadherin⁺ MSC，位于骨内膜——这一区室在骨髓抽吸样本中**结构性不可及**。因此不只是"结构基质降级为辅助证据"，而是要明说**主分析的采样面本身是有边界的**。详见 R10。
 
**(b) "单调加深"从命题正文降为待检验的三个竞争模型之一。** 见 M6。
 
修订后的命题（中文正文替换段）：
 
> 在遗传与表观背景高度异质的 AML 中，LSC 样恶性状态会收敛到一组保守、可跨平台与跨队列复现的**抽吸可及微环境**通讯依赖轴；疾病信号编码在「恶性细胞 → 微环境」通讯网络的**拓扑结构与偏离方向**之中，而非任何单一 ligand–receptor 对；该拓扑相对健康造血的偏离，沿治疗压力轴（诊断 → 残留/MRD → 复发）**按可区分的轨迹模型演化**，其中复发是否沿诊断期同一方向加深，本身即为可证伪的检验对象。
 
---
 
## 2. `[REPLACE §2]` 三个可证伪假设 — 全部改写
 
### H1（收敛性）— 分层从"遗传亚型"扩展为三重分层
 
**v1.0 版本**：保守轴在 ≥60% 独立研究中复现，且在多个**遗传亚型**中富集。
 
**v1.1 版本**：
 
> 存在一组 LSC 样 → 微环境通讯轴，在 ≥60% 独立研究中复现，且在**三种正交分层方案**下均富集：
> (i) **遗传亚型**（NPM1 / TP53 / KMT2A-r / 单核细胞型）；
> (ii) **表观亚组**（eCHROMA 16 个 ATAC subgroup，经 RNA-based ClaNC 分类器指派，见 D4）；
> (iii) **转录状态类**（本项目 04_cnmf 产出的 MP 词汇表聚类）。
> 显著性以置换零分布为基准。
>
> **证伪**：保守轴数量不超过零分布；或保守性在任一分层方案下消失；或可被平台/study 完全解释。
 
**修订理由（M4 的正式升级）**：v1.0 的 M4 因为多数遗传亚型凑不到 ≥3 独立研究而不得不下调 conditional SRRS 门槛，这是一个功效妥协。N2 提供了两个出路：
 
1. **ATAC 亚组现在可从 RNA 指派**（见 D4），分层维度从 4 个稀疏的遗传亚型扩展到 16 个（合并后可用 6–8 个）表观亚组，且分层变量与遗传变量正交。
2. N2 明确证明**遗传分层本身就是错的分层**：除 A/B/C/I 四个融合基因定义的亚组外，穷尽的决策树分析无法用已知 driver 及其组合唯一定义 ATAC 亚组；一个表观亚组可映射到多个 WHO/ICC 遗传亚型，反之亦然。这把"我们为什么不只按突变分层"从辩解变成了引用。
**同时新增一条 H1 的外部支持（写入 §3）**：N2 在 36 个覆盖全部 16 个 ATAC 亚组的样本、281,167 个细胞的 scRNA/ATAC 多组学分析中发现，高 LSC 分数在**所有**亚组中都一致落在分化轨迹早期，即 LSC 样细胞存在于每一种表观背景下，仅丰度不同。**这是"LSC 样状态跨异质背景收敛到共同层级位置"的 1,563 例队列级独立验证**，是 H1 的前提条件。本项目的贡献因此定位为：*层级位置的收敛已确立，本研究检验通讯依赖是否同样收敛*。
 
> ⚠️ **循环论证防火墙**：MP 状态类（分层 iii）可用于 CCC 拓扑分析的分层，**不得**用作 SRRS 自身的分层变量。SRRS 的分层只用 (i) 遗传亚型 与 (ii) ATAC 亚组。
 
---
 
### H2（拓扑性）— 增加"方向"维度
 
**v1.1 版本**：
 
> 网络拓扑偏离在严格跨队列验证下，判别 AML vs 健康的能力优于：(a) 最佳单一 L–R 模型、(b) 仅细胞构成模型、(c) 仅节点特征模型（α=0 端点）。
> **且**：偏离的**方向**（见 M6 定义的 $\cos\theta$ 与正交残差 $\mathrm{HDS}_\perp$）携带超出**标量距离**的增量判别信息。
>
> **证伪**：单 L–R / 构成 / 纯特征模型达到相当判别力（拓扑无增量）；或方向分量在构成匹配后不携带增量信息（则退回标量 HDS，FGW 相对 GW 无优势）。
 
**新增 (c) 的理由**：M8 大幅扩容节点特征后，必须证明 FGW 的融合项不是唯一贡献者，否则审稿人会问"你其实在做特征分类，图结构是装饰"。α=0 端点扫描本来就在 Phase 6 计划里，现在把它提升为 H2 的正式对照臂。
 
---
 
### H3（治疗压力轴）— 从单调性改为三模型竞争检验（**M6**）
 
> **【M6｜H3 从"单调加深"改为三轨迹模型竞争检验】**
>
> **动机（数据内因）**：本项目 03_hierarchy 已观察到 V 形恶性 LSC 架构——Dx 高 primitive-bin 富集 → MRD 谷 → Relapse 重建。若预注册单调递增而结果为 V 形，将被判为事后合理化。
>
> **动机（文献外因）**：N3 用 20-gene state 签名追踪 van Galen 纵向数据时观察到 **MLS / TIM 状态在 D10–25 早期一过性升高后回落**——非单调动力学；同时观察到诊断期 state 构成高度异质，而**复发样本构成彼此趋同**。这两条同时指向"MRD 是第三态"与"复发是收敛态"。
 
**H3 v1.1（预注册三个竞争模型）**：
 
| 模型 | 陈述 | 主检验量 |
|---|---|---|
| **A 单调强化** | $\mathrm{HDS}_{\mathrm{Dx}} < \mathrm{HDS}_{\mathrm{MRD}} < \mathrm{HDS}_{\mathrm{Rel}}$ | Jonckheere–Terpstra 趋势检验 |
| **B 反弹 / V 形** | $\mathrm{HDS}_{\mathrm{MRD}} < \mathrm{HDS}_{\mathrm{Dx}}$ 且 $\mathrm{HDS}_{\mathrm{MRD}} < \mathrm{HDS}_{\mathrm{Rel}}$；MRD 由再生-炎症拓扑主导，表现为 $\mathrm{HDS}_\perp$ 在 MRD 显著升高（第三方向） | 二次项 LRT + $\mathrm{HDS}_\perp$ 组间比较 |
| **C 迟滞 / hysteresis** | $\mathrm{HDS}_{\mathrm{Rel}} > \mathrm{HDS}_{\mathrm{Dx}}$ 但**方向不同**：$\cos\theta_{\mathrm{Rel}} < \cos\theta_{\mathrm{Dx}}$，复发进入新的拓扑区域而非沿 Dx 轴外推 | $\cos\theta$ 配对比较 |
 
**模型选择**：以 patient 为随机效应、study 与 platform 为随机效应、`frac_malignant` 为固定协变量的混合效应模型拟合三者，按 AIC/BIC 与嵌套 LRT 比较，**预注册比较方案，不做事后择优**。
 
**证伪**：三个模型均不优于零模型（时点无信息）；或最优模型的效应在回归掉 blast/恶性比例后消失（退回 R5）。
 
**H3 的外部佐证（可写入正文）**：N5 报告的配对诊断–复发样本中，CD81⁺ LSC 比例由 33% 升至 63%——复发期 LSC 的**受体状态**发生了改变，而非仅仅数量变化。这与模型 C（方向改变）一致，且直接连到 M8 的 receptor competence 节点特征。
 
---
 
## 3. `[REPLACE §3]` 意义与创新性 — 必须加入竞争定位
 
v1.0 的 §3 把先行工作定位为 scACCorDiON（cell-type 级 OT 重心）。**这已经不够了。** N3 使用了与本项目**完全相同的技术栈**：BoneMarrowMap + Symphony 投影、LSC17 + Eppert 干性签名、**CellChat v2.2.0**、van Galen (GSE116256) 治疗轴、Pei 2023 复发对照。审稿人一定会问"和 Struyf et al. 有什么不同"。
 
**新增段落（§3 末尾）**：
 
> **相对最近先行工作的定位。** 近期已有工作将 AML 恶性细胞状态与药敏、可溶性蛋白组及微环境信号关联（N3），并确立了"转录状态而非免疫表型定义 LSC 功能"这一方向。本框架与之在三个层面正交：
>
> 1. **样本规模与队列结构**：先行工作以 6 例诊断期样本定义状态词汇表；本框架在 13 个独立研究、159 个 L2-capable 样本上建立词汇表，并以 study 级 train/validation 划分检验其跨队列复现性。前者无法区分"AML 的状态"与"这 6 个病人的状态"。
> 2. **患者间比较的形式化**：先行工作将多例细胞合并后运行一次 CCC 推断，输出"哪个状态是 signaling hub"——这是**队列平均的描述**。本框架构建 per-patient 图并用 FGW 在拓扑空间对齐，输出**患者级、可与临床 metadata 关联的偏离评分**。跨患者拓扑比较在该文献中完全空白。
> 3. **恶性身份的基因组锚定**：先行工作仅依赖参考投影与聚类定义细胞状态，未做逐细胞恶性判定。其据以建立核心论述的两个"stromal" 状态（vascular-like / mesenchymal-like）来自 BMMC 抽吸物——该制备中基本不含真正的基质细胞，这两个状态更可能是**表达了异常基质样程序的恶性细胞**（该文亦观察到 "VLS 尽管具有基质转录身份，却呈现 progenitor 样互作谱" 这一内部矛盾）。
>
> 这构成一个可推广的方法学论点：**在缺乏逐细胞恶性标签的情况下，恶性细胞的异常转录程序会被误判为基质区室，并据此建立错误的 niche 因果叙事。** 本框架以 inferCNV/Numbat 共识证据将恶性身份（基因组来源）与层级身份（投影来源）正交分离，从设计上排除这一失效模式。
 
**另需在 §3 增补的一句（针对 N1 的两个阴性结果）**：
 
> 共表达并不蕴含功能。N1 在同一模型中证明：N-cadherin⁺ MSC 同时共表达 *Cxcl12*、*Gas6*、*Angpt1*、*Kitl*、*Tgfb1*、*Vcam1*，但条件敲除 *Kitl*（SCF）对 LSC 维持、AML 植入与分布**无任何影响**；敲除 Nestin⁺ MSC 的 *Cxcl12* 同样无影响。任何基于配体–受体共表达的 CCC 推断都会把 KITL–KIT 排入 top edge，而它在该系统中功能无关。**这说明单边强度排序不可靠，而跨队列复现性与拓扑收敛才是可用的过滤器**——正是本框架的核心设计前提。
 
> **【R13｜审稿：你的框架被 Struyf 抢先了】** 应对：以上三点写进 §3 正文而非 Discussion；并在 §12 时间线中把"L2 层跨队列状态词汇表 + per-patient 图"列为可独立成文的第一篇（不等 L3），以规模与形式化为卖点。
 
---
 
## 4. `[AMEND §5]` 数据基础 — 新增第四条硬约束
 
v1.0 §5 列了三条硬约束（结构 niche 缺失 / 恶性计数缺失 / 平台混杂）。**新增第四条**：
 
> **4. 抽吸样本对骨髓纵向解剖轴是坍缩的，且对 metaphysis 区室系统性欠采样。**
> → 主 claim 限定为 aspirate-accessible 区室；MRD 时点的恶性细胞估计需双重折扣（见 R10）。
 
### `[INSERT AFTER R9]` 新增三条审稿风险
 
> **【R10｜审稿：真正决定 LSC 命运的是空间位置，而你把它平均掉了】**
>
> **攻击形态**：N1 证明 LSC 命运由骨髓**纵向位置**决定——静息 LSC 优先定位于 metaphysis，被动员到 central marrow 后进入循环、失去干性并凋亡；PM/DM 来源 L-GMP 的 tumour-initiating cell 频率为 1/76.6 与 1/102，CM 来源为 1/394（PM vs CM *P* = 1.8×10⁻⁴）。你全部数据是骨髓抽吸，既打散空间结构，又对骨内膜/metaphysis 区室系统性欠采样。你的 CCC 图恰好在你声称最重要的轴上做了平均。
>
> **这是有效攻击，必须正面承认**，不能靠术语回避。三层应对：
> 1. **范围声明**：主 claim 全程限定为 "aspirate-accessible BM microenvironment"；N-cad⁺ MSC / GPC3 / 骨内膜区室在设计上不可及，明确写入 Limitation。
> 2. **把它转化为 MRD 解释的护栏**：本项目已独立观察到「MRD 期 HSC_MPP bin 中 stemness 最高的细胞未被 inferCNV 判为恶性」。现在有**两个正交机制**共同解释这一盲区——(i) 静息 LSC 转录输出低，CNV 检出率下降；(ii) 静息 LSC 优先位于 metaphysis，在抽吸物中本就欠采样。**双重不可见性**把 MRD 谷从"方法学 caveat"升级为有机制支撑的预期偏倚，直接支撑 M6 的模型 B。
> 3. **量化敏感性**：按时点报告「stemness 高分位但未被判为恶性」细胞的比例，作为采样+检出盲区的直接指标；该指标随时点的变化必须与 HDS 轨迹一并呈现。
>
> **不做的事**：不声称能从抽吸 scRNA-seq 恢复空间信息；不把 CXCL12–CXCR4 边的强度解释为 retention 的定量代理（见 R12）。
 
> **【R11｜审稿：AML 的真实分层是表观的，你按突变分层等于分错了层】**
>
> **攻击形态**：N2 在 1,563 例中定义 16 个 ATAC 亚组，具有独立于 ELN 的预后价值与特异药敏；除 A/B/C/I 外，穷尽决策树无法用已知 driver 组合唯一定义这些亚组。你按 NPM1/TP53/KMT2A 分层声称"跨亚型收敛"，可能只是在错误分层上做平均。
>
> **应对**：见 H1 三重分层 + D4。核心是**这不再是无法回应的攻击**——ATAC 亚组可从 RNA 指派。同时把 N2 的 LSC 跨亚组收敛结果引为 H1 前提的正面支持。
 
> **【R12｜审稿：CellChat 看不到配体的翻译后加工，你的 CXCL12 边是反的】**
>
> **攻击形态**：N1 的机制是 **CXCL12–DPP4–GPC3** 三元轴：AML 细胞高表达的 DPP4 截短并失活 CXCL12，N-cad⁺ MSC 上的 GPC3 又抑制 DPP4。决定性证据是——*Dpp4* 敲除后 **CXCR4 表达未变、体外迁移能力未变**，改变的只有梯度本身。**配体表达不变、受体表达不变、功能完全逆转。** 任何基于 L–R mRNA 共表达的推断在 CM 区室会把 CXCL12–CXCR4 打成高分，而那里恰恰是 CXCL12 正在被灭活的地方。
>
> **应对（可执行，成本近零）**：新增 **modifier-enzyme 节点特征层**（见附录 C 第 5 块），把配体加工/降解酶与其内源抑制物作为节点特征进入 FGW 融合项，并在所有 top emergent edge 的解释中做 modifier-aware 标注。措辞纪律：CCC 推断输出称 **"communication propensity"**，不称 "signaling activity"。
>
> **附带收益**：这同时部分吸收了代谢物–受体轴的问题（N11）——腺苷轴 `ENTPD1 / NT5E → ADORA2A/2B` 是唯一在 AML 中既重要、又能用外切酶表达近似捕捉的代谢信号轴。
 
---
 
## 5. `[AMEND Phase 4]` CCC 图构建 — 新增功能轴标注（**M7**）
 
> **【M7｜边的功能轴分解：把 HDS 从标量变成向量】**
>
> **动机**：N4 证明人源 ILC1 通过 TNF 阻断 LSC 向白血病祖细胞及 M2 巨噬样支持性细胞分化，通过 IFN-γ 部分抑制后续分化步骤——即一条 CCC 边可以控制 LSC 的**去向**而非生死。survival 与 differentiation-control 是正交的功能轴。v1.0 的 HDS/ATS 是单一标量，把这两类混在一起，既不可解释也无法定位偏离来源。
>
> **方案**：CellChat 的 `subsetCommunication()` 输出已含 `pathway_name`。建立一张**预注册的** pathway → functional axis 映射表（附录 B），冻结于版本控制，然后用最优耦合 $T^*$ 按轴聚合，得到分轴偏离评分：
>
> $$\mathrm{HDS}^{(a)} = \sum_{(i,j)\in\mathcal{E}_a} T^{*}_{ij}\,\big|C_{\text{patient}}[i,j] - C_{B_{\text{healthy}}}[i,j]\big|,\qquad a \in \mathcal{A}$$
>
> 其中 $\mathcal{A}$ = {Survival, Retention, Differentiation, ImmuneEvasion, Inflammatory}，$\mathcal{E}_a$ 为归属轴 $a$ 的边集。标量 HDS 保留为全轴之和，用于与简单基线做 benchmark。
>
> **实现代价**：post-hoc join，**不需要重跑 CellChat**。
>
> **必须写入的防御**：
> - 功能轴标注是**先验注释、非结论**；IFN-γ / TNF / TGF-β 在 AML 中方向双重（本队列的 GSE239721 本身即 IFN-γ AML 研究，GSE185381/Lasry 为炎症状态研究）。
> - 双重归属 pathway 采用 **primary + secondary 双列**（附录 B），主分析用 primary-only，敏感性分析用 0.5/0.5 权重拆分，报告两者一致性。
> - Therapy-adaptation **不是先验轴**：它由 MRD/Relapse 时点的 emergent edge 事后定义，标注为 derived，不参与 H1 检验。
 
> **【D3｜节点词汇表冻结】** 7-bin 词汇表（`HSC_MPP / LMPP_GMP / Mono_DC / Erythroid / Megakaryocyte / T_NK / B_Plasma`）**不因本补丁改变**。ILC1（N4）在 BMMC 中占比约 0.1–0.5% 且与 NK/T 难以分离，**本框架不承诺 ILC1 节点**；其影响以 `T_NK` 节点的 ILC1-signature 特征分 + IFN-II/TNF 轴边的形式进入模型。任何扩节点的提议都会触发 159 样本重跑，本阶段一律拒绝。
 
---
 
## 6. `[AMEND Phase 5–6]` 节点特征扩容与 FGW 融合项（**M8**）
 
> **【M8｜节点特征从 3 维扩容为分块向量：让 FGW 的 "F" 真正起作用】**
>
> **动机（方法学软肋）**：v1.0 的节点特征为 `(frac_malignant, mean_stemness, n_cells)` 三维。FGW 的融合项
> $$(1-\alpha)\sum_{i,j} D(x_i,y_j)\,T_{ij}$$
> 在三维特征下几乎不承载信息——实质上等价于纯 GW，α 扫描会近乎平坦，且**答不上"为什么用 FGW 而不是 GW"**。
>
> **动机（生物学）**：N5 证明 CD81（tetraspanin，参与膜蛋白组织与黏附）支持 LSC 功能，且配对样本中 CD81⁺ LSC 比例 Dx 33% → Relapse 63%。这把 LSC intrinsic state、niche attachment 与临床复发连成一条链，说明**受体端的"接收能力"是一个独立于配体供给的维度**。
>
> **方案**：节点特征扩为 6 个语义块（完整清单见附录 C）：
> B1 结构/构成 · B2 干性与周期 · B3 receptor competence · B4 黏附/tetraspanin · B5 免疫逃逸配体能力 · B6 modifier enzymes（R12）· B7 receiver response state（M9）。
>
> **两条不可省的实现纪律**：
> 1. **全部特征在样本内取 rank-percentile，不用归一化表达值。** 否则会把平台效应从边权重（已由 D2 的 rank 距离压制）重新灌进特征项，等于绕过 D2。
> 2. **特征距离按块等权聚合**，不按基因数等权：
> $$D(x_i, y_j) = \frac{1}{|\mathcal{B}|}\sum_{b\in\mathcal{B}} d_b\!\left(x_i^{(b)}, y_j^{(b)}\right)$$
> 否则基因数最多的 B3/B6 会因维数而支配整个融合项。
>
> **对 H2 的连带影响**：α=0（纯特征 OT）端点现在是 H2 的正式对照臂之一（见 §2）。若 α=0 达到与 α=0.5 相当的判别力，则拓扑无增量信息，H2 被证伪——这是扩容特征后必须付出的诚实代价。
 
---
 
## 7. `[REPLACE Phase 7]` 患者级评分 — 加入方向分量（**M6** 数学部分）
 
v1.0 的 HDS/ATS 是 GW 距离，**是标量，丢失了方向**。"复发离健康更远"与"复发往一个新方向走"是两个完全不同的生物学命题，而只有后者能体现 FGW 相对简单距离的优势。
 
**定义（新增）**。设健康重心 $B_H$、AML 重心 $B_A$，患者图 $G_p$。用 $G_p$ 与 $B_H$ 的最优耦合 $T^*$ 把 $G_p$ 的代价矩阵拉回到 $B_H$ 的固定节点词汇上，得到对齐后的 $\tilde{C}_p$。取上三角向量化：
 
$$\mathbf{d}_p = \mathrm{vec}_{\triangle}\!\left(\tilde{C}_p - C_{B_H}\right), \qquad
\mathbf{u} = \frac{\mathrm{vec}_{\triangle}\!\left(C_{B_A} - C_{B_H}\right)}{\left\lVert \mathrm{vec}_{\triangle}\!\left(C_{B_A} - C_{B_H}\right)\right\rVert}$$
 
于是患者级读数由 1 个标量扩为 4 个：
 
| 量 | 定义 | 含义 |
|---|---|---|
| $\mathrm{HDS}$ | $\lVert\mathbf{d}_p\rVert$（等价于原 GW 距离） | 偏离**幅度** |
| $\mathrm{HDS}_\parallel$ | $\langle \mathbf{d}_p, \mathbf{u}\rangle$ | 沿 healthy→AML 主轴的投影 |
| $\mathrm{HDS}_\perp$ | $\sqrt{\lVert\mathbf{d}_p\rVert^2 - \mathrm{HDS}_\parallel^2}$ | **正交残差**：走到主轴之外的程度 |
| $\cos\theta_p$ | $\mathrm{HDS}_\parallel / \lVert\mathbf{d}_p\rVert$ | 偏离**方向**与主轴的一致性 |
 
$\mathrm{ATS}$、$\mathrm{RLS}$ 同法分解。$\mathrm{SCS}$（M5 定义的干性-微环境通讯评分）保留，并按 M7 分轴报告。
 
**这组量直接映射到 M6 的三模型**：模型 B 的判据是 MRD 的 $\mathrm{HDS}_\perp$ 升高；模型 C 的判据是 $\cos\theta_{\mathrm{Rel}} < \cos\theta_{\mathrm{Dx}}$。**没有方向分量，B 和 C 无法与 A 区分。**
 
---
 
## 8. `[INSERT INTO Phase 8]` 接收端响应验证（**M9**）
 
> **【M9｜receiver-side response validation：给 CCC 推断一个生物学证伪面】**
>
> **动机**：v1.0 对 CCC 假阳性的防线（R3）全是统计的——双法一致、构成置换零分布。N1 提供了一个**生物学**防线的材料。
>
> **材料**：N1 中两个独立扰动（`Dpp4`-KO 与 `N-cad; Cxcl12`⁻/⁻）导致恶性细胞转录组**收敛**：PCA 上两者与各自对照分离并彼此聚拢；GSEA 一致显示 IL6–JAK–STAT3、MAPK、NF-κB 负富集，代谢过程正富集；蛋白层 p-STAT3 / p-ERK1/2 / p-p38 / p-p65 一致下降。即"**niche 支持被剥夺**"具有可复现的转录指纹。RNA-seq 公开于 SRA `SRP323430`（BioProject `PRJNA1077712`）。
>
> **构建**：`NSD score`（Niche-Support-Deprived signature）= 两组比较的 DEG 交集，鼠→人 ortholog 映射，方向：↓IL6-JAK-STAT3 / ↓MAPK / ↓NF-κB / ↓stemness / ↑OXPHOS / ↑cell cycle。
>
> **可证伪预测**：在每个样本内，`HSC_MPP` 恶性组分的 NSD score 应与该节点接收的 **Retention + Survival 轴入边权重之和负相关**（跨样本 Spearman，以 platform/study 为随机效应）。
>
> **证伪与后果**：若无相关，说明本框架推断出的 retention/survival 边不反映接收细胞的实际信号状态 → **边级生物学解释被削弱**（须降格为 propensity 描述），但**拓扑级跨患者比较仍可成立**（H2/H3 不受此检验直接影响）。这是一个分层的、诚实的证伪设计。
>
> **必须标注的限制**：小鼠来源签名 + ortholog 映射；MLL-AF9 / AML-ETO9a 模型对应人类 KMT2A-r 与 t(8;21)，外推到其他亚型需谨慎。因此 NSD 检验作为**支持性证据**呈现，不作为 H1–H3 的主判据。
 
### `[AMEND Phase 8]` 新增两项稳健性检验
 
- **S9 — 功能轴归属稳健性**：primary-only vs. primary/secondary 0.5 拆分，两套映射下分轴 HDS 的排序一致性（Spearman ρ ≥ 0.8 为通过）。
- **S10 — 特征块消融**：逐块移除 B2–B7，检验 H2 的判别力是否由单一特征块驱动。若移除 B3（receptor competence）后判别力崩塌，说明结论是特征驱动而非拓扑驱动，须在正文明说。
---
 
## 9. `[AMEND §9]` 预期产出与评估判据 — 更新
 
| 产出 | v1.1 变更 | 评估 / 证伪判据 |
|---|---|---|
| 1. MP 词汇表 | 增加 ATAC 亚组标注 | H1：三重分层（遗传 / ATAC / MP 类）下均显著超置换零分布 |
| 2. 健康 / AML barycenter + emergent edges | 增加功能轴标注与 modifier-aware 注释 | 置换零分布显著 + 外部队列复现；S9 通过 |
| 3. 患者级评分 | **由 4 个标量扩为 4×(1+5) 组读数**：{HDS, HDS∥, HDS⊥, cosθ} × {全轴 + 5 功能轴} | H2：拓扑模型 vs 单 L–R / 仅构成 / α=0 三个基线，构成匹配后仍优；S10 通过 |
| 4. **治疗压力轨迹模型**（新） | A/B/C 三模型竞争拟合结果 | H3：最优模型显著优于零模型，且在回归掉 `frac_malignant` 后保持 |
| 5. **NSD 一致性检验**（新） | receiver-side 验证结果 | 支持性证据；不通过则边级解释降格，不影响 H2/H3 |
 
---
 
## 10. `[AMEND §11]` 风险表增补
 
| 风险 | 应对 |
|---|---|
| **空间坍缩（抽吸不含 metaphysis / N-cad⁺ MSC）** | 范围限定 aspirate-accessible；转化为 MRD 双重不可见性护栏；量化 stemness-high-not-malignant 比例（R10） |
| **表观分层错位（16 ATAC 亚组）** | RNA-based ClaNC 指派 + H1 三重分层；引 N2 的 LSC 跨亚组收敛为正面支持（R11, D4） |
| **配体翻译后加工不可见（DPP4/GPC3）** | modifier-enzyme 节点特征块 + 措辞降级为 communication propensity（R12） |
| **被 Struyf et al. 抢先** | §3 三点正交定位；L2 层先行成文（R13） |
| **功能轴归属主观** | 预注册映射表 + primary/secondary 双列 + S9 稳健性（M7） |
| **特征扩容导致"拓扑是装饰"** | α=0 端点作为 H2 正式对照臂 + S10 块消融（M8） |
| **配体竞争未建模** | CellChat 逐边独立计算，不建模多 receiver 对同一配体池的竞争；列为 Limitation，可选 collective-OT 敏感性分析（N8） |
 
---
 
## 11. 新增决策与附录
 
> **【D4｜表观亚组以 RNA 指派，不追加 ATAC 实验】**
>
> N2 提供了基于表达的 ATAC 亚组预测模型（ClaNC 最近质心分类器，**每亚组 15 个标记基因**），并在 TCGA / BeatAML / Leucegene 等 **4 个外部队列共 1,079 例** 上验证；另提供一个 **30 基因高危亚组预测模型**。基因列表见该文 Supplementary Table 11。
>
> **本项目做法**：对每个样本的**恶性组分 pseudobulk** 应用 ClaNC 模型，得到 ATAC 亚组标签 + 高危评分，作为**样本级协变量与分层变量**（不是节点特征）。
>
> **必须做的 sanity check**：模型训练于诊断期高 blast 的 **bulk** RNA-seq，投到 scRNA-seq 恶性 pseudobulk 存在 domain shift。以 Petti2019（有已知核型/突变）与 Chen2023（NPM1）验证：由融合基因唯一定义的亚组 A(`PML::RARA`) / B(`RUNX1::RUNX1T1`) / C(`CBFB::MYH11`) 与 I(`CEBPA`-bZIP) 是否被正确指派。**若这四个可验证亚组指派失败，D4 整体作废，退回仅遗传+MP 双重分层。**
>
> **不做**：pySCENIC regulon 活性作为节点特征列为**条件性追加**（159 样本成本高），不进关键路径。
 
> **【D5｜空间验证限定为 imaging-based 单细胞分辨率】**
>
> 若 Phase 8 执行可选空间交叉验证，**只用 Xenium / CosMx / MERFISH 等成像法单细胞空间数据，不用 Visium 等 spot-level 数据**。理由：spot-level 需要反卷积或单细胞空间重建（N10 处理的正是这一层），邻接关系的不确定性会直接污染 CCC 结论，引入一个本可避免的失效环节。此决策现在锁定，避免后期临时选型。
 
---
 
### 附录 B — CCC 功能轴映射表（M7，预注册）
 
**冻结规则**：本表在查看任何 barycenter/emergent edge 结果**之前**冻结，存为 `scripts/05_ccc/ccc_functional_axis_map.tsv`，进版本控制。任何后续修改必须留 git 记录并在正文说明。
 
**先验轴（5 条，参与 H1–H3）**
 
| 轴 | 生物学定义 | 主要 CellChat pathway（primary） | secondary |
|---|---|---|---|
| **Survival** | 抗凋亡、增殖与代谢支持 | `CSF`, `KIT`, `FLT3`, `IGF`, `HGF`, `EGF`, `FGF`, `PDGF`, `VEGF`, `ANGPTL`, `GAS`, `PROS`, `MK`, `PTN`, `VISFATIN`, `MIF` | `IL6`, `SPP1` |
| **Retention** | 黏附、滞留与静息维持 | `CXCL`, `SELL`, `SELPLG`, `SELE`, `VCAM`, `ICAM`, `ITGAL-ITGB2`, `FN1`, `COLLAGEN`, `LAMININ`, `CDH`, `CADM`, `ESAM`, `JAM`, `PECAM1`, `ANGPT`, `THBS`, `TENASCIN`, `VTN` | `TGFb`, `SPP1` |
| **Differentiation** | 分化推动或阻断、状态转换控制 | `IFN-II`, `TNF`, `NOTCH`, `WNT`, `ncWNT`, `BMP`, `ACTIVIN`, `GDF`, `EPO`, `THPO`, `OSM`, `LIFR` | `IL6`, `TGFb` |
| **ImmuneEvasion** | 抑制免疫杀伤与识别 | `GALECTIN`, `PD-L1`, `SIRP`, `CD226`, `PVR`, `NECTIN`, `TIGIT`, `MHC-I`, `CD48`, `CD70`, `CD80`, `CD86`, `ICOS`, `CD39`, `VISTA`, `BTLA`, `LAIR1` | `TGFb`, `APP` |
| **Inflammatory** | 炎症/应激重塑（非特异免疫抑制） | `CCL`, `CX3C`, `IL1`, `IL17`, `IFN-I`, `ANNEXIN`, `COMPLEMENT`, `GRN`, `RESISTIN`, `SAA` | `IL6`, `TNF` |
 
**衍生轴（1 条，事后定义，不参与 H1）**
 
| 轴 | 定义方式 |
|---|---|
| **TherapyAdaptation** | 在 MRD / Relapse 时点相对 Dx barycenter 显著涌现、且在置换零分布下显著的边集合。标注为 `derived`，只用于描述与假设生成。 |
 
**未分配（`UNASSIGNED`）**：其余 pathway 标注为未分配，进入全轴 HDS 但不进入任何分轴 HDS。分配率必须在正文报告（若 UNASSIGNED 占总耦合质量 > 30%，说明轴体系覆盖不足，需扩表）。
 
> ⚠️ **上表的 pathway 名称必须先与你实际安装的 CellChatDB (v2.2.0.9001) 校验**，不同版本命名有出入。校验一行：
> ```r
> setdiff(map$pathway_name, unique(CellChatDB.human$interaction$pathway_name))
> ```
> 非空即为拼写/版本不符，必须逐条修正后再冻结。
 
---
 
### 附录 C — 节点特征向量规范（M8 + R12 + M9）
 
**通用规则**
1. 除 B1 外，所有特征在**样本内**取 rank-percentile（`frank(x)/.N`），不用归一化表达值 —— 保护 D2 平台不变性。
2. 恶性/正常组分**分别**计算，作为同一节点的两组子特征（节点不拆分，遵循既有决策）。
3. 特征距离按**块等权**聚合，不按基因数等权（见 M8）。
4. 某块基因在样本中全部缺失 → 该块记 `NA`，由 unbalanced FGW 处理，不填 0。
| 块 | 内容 | 依据 |
|---|---|---|
| **B1 结构/构成** | `frac_malignant`, `log(n_cells)`, `frac_of_sample` | v1.0 保留 |
| **B2 干性与周期** | LSC17, van Galen HSC-like, Zeng2025 LSC score, G0/quiescence score, G2M score | v1.0 扩展；与 N2 同一套签名 |
| **B3 receptor competence** | `CXCR4, KIT, FLT3, IL3RA, CSF2RA, CSF2RB, MPL, TEK, AXL, MERTK, NOTCH1, NOTCH2, TGFBR1, TGFBR2, IFNGR1, IFNGR2, TNFRSF1A, TNFRSF1B, IL6R, IL6ST, ADORA2A, ADORA2B` | N5（CD81 逻辑推广）；H3 模型 C |
| **B4 黏附 / tetraspanin** | `CD81, CD82, CD9, ITGA4, ITGB1, ITGAL, ITGB2, CD44, SELPLG, ADGRG1, VLA4 复合体成员` | **N5 核心** |
| **B5 免疫逃逸配体能力** | `CD274, PVR, NECTIN2, LGALS9, CD47, CD24, CD70, CD86, HLA-E` | N3（TIGIT 轴）；N4 |
| **B6 modifier enzymes** | `DPP4, GPC3, MMP2, MMP9, MMP14, ADAM10, ADAM17, HPSE, CTSG, ELANE, PLAU, PLAUR, TFPI, ENTPD1, NT5E` | **N1 核心（R12）** |
| **B7 receiver response state** | NSD score（M9）, HALLMARK_IL6_JAK_STAT3, HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_KRAS/MAPK, HALLMARK_OXPHOS | **N1（M9）** |
 
**样本级（非节点级）协变量**：ATAC subgroup label + 30-gene 高危评分（D4）；platform；study；`frac_malignant`（全样本）。
 
> **B5/B6 与 TP53 侧项目的连接**：此前测得 `PVRL4/NECTIN4` 检出率 ~0–0.95% 且无恶性富集。**注意 NECTIN4 与 NECTIN2 (PVRL2) 不是同一分子**——N3 报告的 AML 免疫抑制轴是 `NECTIN2–TIGIT` 与 `PVR–TIGIT`，二者在 BMMC 中表达远高于 NECTIN4。侧项目的阴性结果**不排除** NECTIN2/PVR 轴，建议在剩余可行分析中补测。
 
---
 
## 12. 实施清单
 
### QC-0 —— **array 提交前必须跑掉的一行检查**
 
技术基因移除 pattern 与 CellChatDB 基因集的交集。N3 报告的 `MT-RNR2 (Humanin) – FPR2/FPR3` 边中，`MT-RNR2` 命中你 `CODING_STANDARDS §8` 的线粒体 pattern `^MT[-.]`（核编码旁系 `MTRNR2L*` 不命中）。虽然该 pattern 目前只声明用于 04_cnmf，仍需确认 05_ccc 输入路径未继承：
 
```r
db_genes <- unique(c(CellChatDB.human$interaction$ligand,
                     CellChatDB.human$interaction$receptor))
db_genes <- unique(unlist(strsplit(db_genes, "_")))
grep("^MT[-.]|^RP[SL]|^HB[^P]", db_genes, value = TRUE)
```
非空 → 记录被移除的 CellChatDB 成员，在 Limitation 中声明；并确认 05_ccc 的输入对象未做该过滤。
 
### 与 159 样本 array 的关系
 
| 项 | 是否阻塞 array | 落点 |
|---|---|---|
| QC-0 gene pattern 检查 | **是**（一行，先跑） | 05_ccc |
| M7 功能轴映射表 | 否（post-hoc join） | 05_ccc 输出后 |
| M8 节点特征提取 | 否（独立于 CellChat） | 06_distance 输入 |
| M9 NSD 签名构建 | 否（外部数据） | 09_robustness |
| D4 ATAC 亚组指派 | 否（bulk pseudobulk） | 样本级 metadata |
| R10/R11/R12 文字 | 否 | 主稿 |
| 节点词汇表 | — | **冻结，不改**（D3） |
 
### 新增脚本占位（按 CODING_STANDARDS 命名）
 
```
scripts/05_ccc/ccc_functional_axis_map.tsv          # M7 预注册映射表（冻结）
scripts/05_ccc/06_annotate_functional_axis.R        # M7 post-hoc join
scripts/06_distance/03_build_node_features.R        # M8 分块特征提取
scripts/06_distance/node_feature_blocks.tsv         # M8 基因清单（B3–B6）
scripts/07_fgw/04_deviation_direction.R             # M6 方向分量 HDS∥/HDS⊥/cosθ
scripts/08_scoring/05_trajectory_models.R           # M6 三模型竞争拟合
scripts/09_robustness/03_nsd_receiver_validation.R  # M9
scripts/09_robustness/04_atac_subgroup_assign.R     # D4 ClaNC 指派 + sanity check
```
 
---
 
## 13. 参考文献增补（接续 v1.0 编号 38–）
 
38. Wang C, Pan Y, Dong R, …, Li L, Kang X. Longitudinal localization of leukaemic stem cells between the metaphysis and central marrow governs their behaviour. *Nat Cell Biol*. 2026;28(5):890–902. doi:10.1038/s41556-026-01939-3. [数据：BioProject PRJNA1077712；RNA-seq SRP323430]
39. Ochi Y, Liew-Littorin M, Nannya Y, …, Lehmann S, Ogawa S. Chromatin landscape and epigenetic heterogeneity of acute myeloid leukaemia (eCHROMA). *Nature*. 2026 (online 8 July). doi:10.1038/s41586-026-10703-4.
40. Li Z, et al. Human type-1 innate lymphoid cells control leukemia stem cell differentiation and limit acute myeloid leukemia development. *Nat Commun*. 2026. doi:10.1038/s41467-026-68582-2.
41. Surface CD81 supports leukemia stem cell function and reveals a therapeutic vulnerability in acute myeloid leukemia. *Signal Transduct Target Ther*. 2026. doi:10.1038/s41392-026-02697-2. 〔HR 数值与队列限定待核对终稿〕
42. Struyf N, Hartmanis L, Rico Pizarro L, …, Kallioniemi O, Erkers T. Transcriptionally defined AML cell states associate with treatment response and microenvironmental remodeling. *bioRxiv*. 2026. doi:10.64898/2026.07.01.735780. 〔**预印本，未经同行评审**；作为竞争定位对象引用，不作承重证据〕
43. Cang Z, et al. Screening cell–cell communication in spatial transcriptomics via collective optimal transport (COMMOT). *Nat Methods*. 2023;20:218–228. doi:10.1038/s41592-022-01728-4.
44. Interpretable niche-based cell–cell communication inference using multi-view graph neural networks (STCase). *Nat Comput Sci*. 2025. doi:10.1038/s43588-025-00809-6.
45. Advancing spatial cellular communication inference with ligand diffusion and transport model (SCILD). *Commun Biol*. 2026. doi:10.1038/s42003-025-09413-w.
46. 〔Tier-3，引用前自行核实〕miR-126 vascular niche agent-based model. *npj Syst Biol Appl*. 2026. doi:10.1038/s41540-026-00675-6.
47. 〔Tier-3〕CellRefiner. *Nat Commun*. 2026. doi:10.1038/s41467-026-70090-2.
48. 〔Tier-3〕SUCNR1/succinate limits haematopoiesis and AML progression. *Nat Commun*. 2026. accession 待补。
49. 〔Tier-3，News & Views〕Extramedullary AML niche in lung. *Nat Immunol*. 2026. doi:10.1038/s41590-026-02598-0.
---
 
## 14. v1.0 → v1.1 变更摘要（一页版）
 
| 编号 | 类型 | 一句话 |
|---|---|---|
| M6 | 修订 | H3 从"单调加深"改为 A/B/C 三轨迹模型竞争；HDS 增加方向分量 HDS∥/HDS⊥/cosθ |
| M7 | 修订 | HDS 按 5 条先验功能轴分解，映射表预注册冻结 |
| M8 | 修订 | 节点特征从 3 维扩为 7 个语义块；样本内 rank-percentile；块等权距离 |
| M9 | 新增 | 用 N1 的"niche 支持剥夺"签名做 receiver-side 响应验证 |
| R10 | 新增风险 | 空间坍缩：抽吸不含 metaphysis / N-cad⁺ MSC；转化为 MRD 双重不可见性护栏 |
| R11 | 新增风险 | 16 个 ATAC 表观亚组挑战遗传分层；以 RNA 指派回应 |
| R12 | 新增风险 | 配体翻译后加工（DPP4/GPC3）对表达法不可见；modifier-enzyme 特征块 |
| R13 | 新增风险 | Struyf et al. 技术栈重合；以规模/形式化/基因组锚定三点定位 |
| D3 | 决策 | 7-bin 节点词汇表冻结，拒绝加 ILC1 节点 |
| D4 | 决策 | ATAC 亚组以 RNA-based ClaNC 指派，不追加实验；四个融合亚组作 sanity check |
| D5 | 决策 | 空间验证限定 imaging-based 单细胞，排除 spot-level |
| H1 | 假设改写 | 分层扩为 遗传 × 表观 × 转录状态 三重 |
| H2 | 假设改写 | 增加 α=0（纯特征）对照臂与方向分量增量检验 |
| H3 | 假设改写 | 见 M6 |
 