# Candidate revision — M10: co-option lens at the communication layer (L3)

> **状态：CANDIDATE（草稿，待 review）。** 本文件**不修改** blueprint 本体，也不改 v1.1 patch；它是一条待评估的补丁提案，计划在**组会之后 + Phase-1 标注/元数据修复之后**再决定是否并入 `AML_niche_CCC_OT_blueprint_zh.md`。在此之前，H1 的正式结论仍以现有 emergent-edge 检验为准。
>
> 建议合入位置：`[AMEND H1]`（v1.1 patch §2）+ `[INSERT Phase 7–8]`（新增模块 **M10**）+ `[AMEND §9]` 判据。
> 模块号 **M10** 当前未被占用（v1.1 patch 中的 "M10" 是 ADAM10 子串，非模块）。
>
> 提出日期：2026-07-26。先导证据脚本：`scripts/10_figures/cooption_redistribution.py`（已跑，见 §5）。

---

## 0. 动机 — 我们测的是"涌现"，没测"挪用"

现有 H1 检验（Phase 7 `emergent_edges`）算的是每条边的 **差异**：

```
dC(edge) = C_AML(edge) − C_healthy(edge)   →  置换检验
```

这是一个**差异检验**：只有当一条通讯边的强度在 AML 与健康之间**改变**时，它才会被标记为显著。结果为 0/49 边显著、5 个 block 全部 q = 0.93、pre-registered 的 primitive→immune 轴 p = 0.86 —— 一个干净的 null。

问题在于：**差异检验对"被挪用（co-opted）的依赖"是结构性失明的。** 如果一条通路在健康造血里**本来就存在**，AML 只是**换个用法**在用它（同等强度、但改变发送方/接收方/网络情境），那么 `dC ≈ 0`，emergent 检验必然给出 null —— 而这恰恰就是我们观察到的结果。换句话说，**当前的 null 没有排除 co-option；它是 co-option 会自然产生的结果。**

**H1 自身逻辑里还有一处错位需要点明**：H1 假设的文字是"存在一组**保守、可复现**的 LSC→niche 轴（convergence）"，但代码把它操作成了"**AML 与健康不同**的边（emergence）"。一条保守的依赖轴按定义就应当在健康里也有，`dC` 自然接近 0 —— **假设问的是"保守"，检验测的是"差异"，二者不是一回事。**

**这个概念 blueprint 其实已经有了，只是放错了层。** L1（meta-program）层已有 healthy-match 评分，用来区分 **"被挪用（co-opted）" vs "新涌现（emergent）"** 的基因程序（见 blueprint Phase 3）。M10 的核心就是：**把这把"co-opted vs emergent"的尺子，从基因程序层（L1）搬到通讯边层（L3）。**

---

## 1. `[AMEND H1]` — 给 H1 增加一条 co-option 检验臂

v1.1 patch 的 H1 证伪条款目前只覆盖"涌现/差异"。建议扩写为：

> H1 的通讯依赖检验分为两臂：
> **(a) Emergent 臂（现状）**：保守轴在 AML 中相对健康**显著增强**（`dC` 检验）。
> **(b) Co-option 臂（新增 M10）**：一条在 AML 与健康中**同等普遍存在**（prevalence 保守）的通路，其**通讯组成（sender/receiver composition）在 AML 中显著向 LSC/primitive 隔室重新分配**。
>
> **证伪（补充）**：若一条候选轴 (a) 不涌现、(b) 组成也不移位、且在跨 study 对比中无法与 study/platform 效应区分，则该轴不被计为保守依赖。**单靠 emergent 臂的 null 不足以否证 H1。**

理由：co-option 是肿瘤利用正常发育/稳态通路的常见方式；把它排除在检验之外，会让"生物学上最可能为真"的那类信号被系统性漏掉。

---

## 2. `[INSERT Phase 7–8]` M10 — L3 co-option lens（正式定义）

对每个样本 `s`、每条通路（以 ligand `L` 为单位），在其显著边（`pval < 0.05`）上定义三个**通路内、样本内归一化**的组成量（scale-free，不受 cell-count / edge-count 混淆）：

```
total(s,L)             = Σ prob
primitive_sender_share = Σ(prob | sender ∈ PRIMITIVE) / total     # 谁在发这条通路
immune_receiver_share  = Σ(prob | receiver ∈ IMMUNE)  / total     # 谁在收
prim_to_immune_share   = Σ(prob | sender ∈ PRIMITIVE & receiver ∈ IMMUNE) / total   # SCS 轴
```

其中 `PRIMITIVE = {HSC_MPP, LMPP_GMP}`（LSC 样/干-祖隔室），`IMMUNE = {T_NK, B_Plasma, Mono_DC}`。

**两步判据：**

1. **CONSERVED?（保守性）** — 通路在 AML 与健康中的 **prevalence**（≥1 条显著边的样本比例）都要高（预注册门槛，建议两组均 ≥ 60%）。这是"非涌现/本来就有"的证据。
2. **SHIFTED?（移位/再分配）** — 组成量的分布在 AML vs 健康之间**显著不同**（Mann-Whitney U，two-sided），且**方向指向 primitive/LSC 隔室**（`primitive_sender_share` 或 `prim_to_immune_share` 在 AML 中更高）。

**跨 study 混淆的处理（必须）**：健康与 AML 样本来自不同 study，pooled 检验被 study 混淆。因此 SHIFTED 检验**必须**同时在"同时含健康与 AML 的数据集"内部报告（当前为 `GSE116256, GSE185381, Petti2019`）。**pooled 显著 + within-study 同向**才计为候选；单靠 pooled 不足。

**零模型**："组成移位"要区别于噪声。预注册零模型二选一（或都做）：
- **打乱组内标签**：在同一 study 内打乱 healthy/AML 标签，重算组成量差异，得经验 p；
- **打乱 sender 归属**：在样本内保持每条通路的边集与总强度不变、打乱 sender_bin，检验 `primitive_sender_share` 的观测值是否超出零分布。

---

## 3. 与已有模块的关系（避免重复造轮子）

| 已有/相邻模块 | 与 M10 的关系 |
|---|---|
| **L1 healthy-match（Phase 3）** | M10 = 把同一"co-opted vs emergent"判据从**程序层**移到**通讯边层**。可共用 healthy-match 评分公式。 |
| **M8 方向分量（HDS⊥ / cosθ）** | M8 捕捉"边强度不变但网络**情境/邻居**变了"的 co-option（拓扑式挪用）。M10 捕捉"**发送/接收方**变了"的 co-option（组成式挪用）。二者互补，不重叠。 |
| **M9 NSD 接收端验证** | M10 只能找**候选**（组成移位）；"是否为**依赖**"需 M9 的接收端响应签名或湿实验来钉死。**co-option = M10 找候选 → M9 验依赖，两步。** |
| **H2 拓扑检验** | H2 测 bin-level 图的拓扑；M10 测通路的发送/接收组成。M10 的**最强版本**（见下）超出 H2 范围。 |
| **Secondary graph（bin × 恶性态 × role，blueprint 有但未建）** | **M10 的最强版本**：把节点拆成"恶性 vs 正常"再算通讯，直接问"一条健康里也有的边，在 AML 里是否**不成比例地由恶性细胞在发**"。这**依赖可靠标注** —— 因此 Phase-1 标注修复 + secondary graph 不只救 H2/H3，它**解锁了一个 H1 能出信号的检验**。这也部分调和了"标注对 H1 重要"的直觉（机制不是 FPR 污染 dC，而是 co-option 需要恶性-正常分解）。 |

---

## 4. `[AMEND §9]` 预注册判据 — 一条"确认的 co-opted 轴"需要

1. prevalence 在 AML 与健康两组均 ≥ 60%（保守）；**且**
2. 组成量向 primitive/LSC 的移位在 pooled 检验显著（多重比较校正后）；**且**
3. 该移位在 within-study 对比中**同向**（功效允许时亦显著）；**且**
4. 通过 §2 的零模型；**且**
5. （升级为"依赖"时）M9 接收端响应签名支持。

任一不满足 → 记为"候选/趋势"，不记为确认。

---

## 5. 先导证据（pilot，已跑）

脚本 `scripts/10_figures/cooption_redistribution.py`，输入现有 148 样本 CCC tensor + `fgw_input_index.csv` 的 healthy flag，对 top-3 配体（LGALS9 / MIF / CD99）跑了 §2 的两步判据。输出见 `results/tables/05_ccc/cooption/`。

**CONSERVED？—— 稳健成立。** 三个主导配体在 AML 与健康中都无处不在：

| 通路 | prevalence AML | prevalence healthy |
|---|---|---|
| MIF | 98% | 93% |
| Galectin-9 (LGALS9) | 80% | 89% |
| CD99 | 96% | 96% |

→ 这些是**保守、非涌现**信号；H1 的 emergent 差异检验结构上就看不见它们。这是对整个 co-option 论点最干净、最不受混淆的支持。

**SHIFTED？—— pooled 强、within-study 为趋势，且方向一致。** 组成量（中位数，Mann-Whitney）：

| 通路 · 组成量 | AML | healthy | p (pooled) | p (within-study) | 方向 |
|---|---|---|---|---|---|
| **Galectin-9** · primitive_sender_share | 0.47 | 0.19 | 8.6e-4 | 0.205 | ↑ AML |
| **Galectin-9** · prim_to_immune_share | 0.26 | 0.11 | 2.0e-3 | **0.071** | ↑ AML |
| **CD99** · primitive_sender_share | 0.62 | 0.18 | 2.9e-5 | 0.138 | ↑ AML |
| **CD99** · immune_receiver_share | 0.35 | 0.76 | 3.3e-5 | 0.066 | ↓ AML |
| **MIF** · primitive_sender_share | 0.48 | 0.47 | 0.478 | 0.246 | ≈（阴性对照） |

（within-study 池 = GSE116256 / GSE185381 / Petti2019；healthy n = 14–16，AML n = 31–45。）

**图**（`results/figures/13_cooption/cooption_sender_composition.png`，脚本 `scripts/10_figures/fig_cooption_composition.py`）：三个保守通路的**发送隔室组成**堆叠条，AML vs healthy。隔室定义：**Primitive** = HSC_MPP + LMPP_GMP（LSC 样/干-祖）；**Immune** = Mono_DC + T_NK + B_Plasma；**Other** = Erythroid + Megakaryocyte。条上标注的是 primitive 发送份额。

> ⚠️ **看图须知（如实标注）**：此图展示的是**队列平均份额（pooled mean）**，仅供直观；它**不是**显著性声明。Galectin-9、CD99 的 primitive 份额在 AML 中明显更高、MIF 几乎不变——但 pooled 差异**受 study 混淆**，within-study 对比中只是**趋势**（Galectin-9 primitive→immune p=0.071、CD99 immune-receiver p=0.066，均未达 0.05）。显著性与混淆判断以上表为准，勿以此图的视觉差直接下结论。

**解读（诚实版）：**
- **保守性**：铁证 —— 主导信号在健康里同样普遍。
- **再分配**：**Galectin-9 和 CD99** 在 AML 中显著地更多由 primitive/LSC 隔室发出（pooled 极显著），方向与 co-option 假设一致；within-study 对比中**同向但仅为趋势**（Galectin-9 primitive→immune p=0.071，CD99 immune-receiver p=0.066），因为只有 3 个数据集、healthy 端欠功效。pooled 的强度部分来自 study 混淆，不能单独采信。
- **阴性对照**：**MIF** 保守但**几乎不再分配** —— 说明这个再分配信号不是所有配体通用的假象。

**Galectin-9 是当前最强的单分子 co-option 候选**：保守（两组都在）、pooled 极显著地向 LSC→immune 轴重新分配、within-study 同向。它也正是 AML 里经典的免疫抑制配体，生物学上自洽。

**局限**（写入论文时必须带）：(a) 组成移位 ≠ 依赖，依赖需 M9/湿实验；(b) pooled 受 study 混淆，within-study 欠功效；(c) 仅 CellChat，需 LIANA+/NicheNet 复现；(d) bin-level，未做恶性-正常分解（M10 最强版本，待 secondary graph）。

---

## 6. 脚本占位（按 CODING_STANDARDS 命名）

```
scripts/10_figures/cooption_redistribution.py        # M10 先导：组成移位两步判据（已存在）
scripts/10_figures/fig_cooption_composition.py       # M10 先导：发送隔室组成堆叠图（已存在）
scripts/09_robustness/05_cooption_null.py            # M10 零模型（§2，待建）
scripts/05_ccc/07_secondary_graph_malignant.R        # M10 最强版：bin×恶性态×role 通讯（依赖 Phase-1 标注，待建）
```

输出表（已生成）：
```
results/tables/05_ccc/cooption/pathway_summary.csv        # 保守? 移位? 每通路
results/tables/05_ccc/cooption/sender_composition.csv     # 通路×组×sender_bin 平均份额（供堆叠图）
results/tables/05_ccc/cooption/prevalence_top_ligands.csv # top 配体在 AML/healthy 的 prevalence
results/tables/05_ccc/cooption/persample_shares.csv       # 逐样本原始份额（透明留痕）
```

---

## 7. 一句话总结

M10 把 H1 的故事从"**null → 大概没信号**"改写为"**我们只测了涌现；更可能的生物学（挪用保守通路）此前没测**"。先导证据表明：主导 CCC 信号确为保守，且 Galectin-9/CD99 显示出向 LSC→immune 轴的方向一致的组成移位 —— 一个 emergent 差异检验按定义看不见的信号。这条修订同时给"修标注 + 建 secondary graph"补上了一个直接的科学动机。
