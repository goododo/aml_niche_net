# aml_niche_net — Blueprint v1.2 Patch

> **状态:正式补丁(组会后授权)。** 本补丁在 v1.0 基线 + v1.1 补丁之上叠加,合入锚点用 `[AMEND §X]` / `[INSERT]` 标注。它把 `CANDIDATE_M10_cooption_lens.md` 从候选**升为正式模块 M10**,并纳入组会 PI(Daichi Inoue)提出的重构:**把打击目标从"配体"上移到"配体在 LSC 内汇聚的共同下游 hub",并把单靶失败解释为"冗余/代偿"。**
>
> 提出日期:2026-07-27。新模块号 **M11 / M12** 此前未占用。先导证据脚本:`scripts/10_figures/cooption_redistribution.py`(已含 sender/receiver/autocrine),文献支撑:`CCC_TOP_LIGANDS_LITERATURE.md`。

---

## 0. 本次升级的动因

1. **组会反馈(Daichi)**:Gal-9 / CD99 / MIF 三个分子单独看都被发表烂了,且 TIM-3/sabatolimab 临床失败 —— **单打没有 novelty**。真正的问题是:**这三条不同的受体输入(TIM-3 / CD99 / CD74)是否在 LSC 内汇聚到同一个可打的下游节点,以及为什么单靶封锁会被代偿。**
2. **数据已经指向这个方向**:M10 的接收侧探针发现 CD99、Gal-9 在 AML 里塌缩成 **autocrine LSC↔LSC 环**,MIF 保持广泛(见 §5)。这正好落进 Daichi 的"同一 LSC vs 代偿性 niche"二分。
3. 因此把 M10 升为正式模块,并新增 **M11(下游汇聚)** 与 **M12(同一 LSC vs 代偿)**。

---

## 1. `[AMEND §0 / §3]` 中心命题与创新性 — 收窄到"汇聚 hub + 代偿"

**v1.2 中心命题(取代 v1.1 的"被挪用的保守轴"表述):**

> 在 AML/LSC 中,一组**保守的**通讯配体(以 Gal-9 / CD99 / MIF 为代表)不再是各自独立的旁分泌信号,而是**通过不同受体(TIM-3 / CD99 / CD74)汇聚到 LSC 内的一个共同下游调控节点(regulon/hub)**;这种汇聚在**治疗压力下**强化,并通过**多轴冗余**使单靶封锁被代偿。

**创新性定位(写入 §3):**
- 单个配体/受体都已发表,且 TIM-3/sabatolimab 已临床失败 —— 我们**不**主打任何单分子。
- **未被发表的是:** (i) 这些信号**汇聚到同一下游 hub**;(ii) 它们在 AML 里塌缩成**同一 LSC 的自分泌环**;(iii) **冗余/代偿**是单靶失败的机制解释。
- **把 sabatolimab 的失败从"坏消息"反转为论点核心**:单靶被代偿 → 必须打汇聚节点或组合封锁。

---

## 2. `[AMEND H1]` 收敛性 — 从"通讯轴收敛"扩展到"下游 hub 收敛",并加入 autocrine 与二分

H1(收敛性)现在分**两个层级**检验:

> **(L3 边层级 — M10):** 一组保守配体在 AML 中把**发送方与接收方同时移向 LSC/primitive 隔室**,塌缩成 **autocrine LSC↔LSC 环**(相对健康显著)。
> **(下游层级 — M11):** 上述多受体输入在 LSC 内**汇聚到同一个 regulon / 靶基因程序**(相对置换零分布显著)。
>
> **预注册分叉(M12):** 汇聚是 **(A) 同一 LSC 亚群**(单细胞共表达多受体、单一 hub)还是 **(B) 不同 niche 细胞的代偿性支持**(冗余多路)。**两个结果都可发表**:A → "可打的汇聚节点";B → "冗余 niche 支持,需组合封锁"。
>
> **证伪:** 无共同下游 regulon 超过零分布;或 sender/receiver/autocrine 移位在 within-study / 多方法下不成立;或 in-silico 敲除单轴不被其他轴代偿(即无冗余)。

---

## 3. `[INSERT Phase 7–8]` M10(正式)— L3 co-option / convergence lens

对每个样本、每条通路(ligand),在显著边(pval<0.05)上定义**通路内、样本内归一化**的 scale-free 组成量:

```
primitive_sender_share    = Σ(prob | sender ∈ PRIM) / total        # 谁发
primitive_receiver_share  = Σ(prob | receiver ∈ PRIM) / total      # 谁收
primitive_autocrine_share = Σ(prob | sender ∈ PRIM & receiver ∈ PRIM) / total   # LSC->LSC 环 (M10 headline)
prim_to_immune_share      = Σ(prob | sender ∈ PRIM & receiver ∈ IMMUNE) / total # SCS 轴
immune_receiver_share     = Σ(prob | receiver ∈ IMMUNE) / total
```
`PRIM = {HSC_MPP, LMPP_GMP}`,`IMMUNE = {Mono_DC, T_NK, B_Plasma}`。

**判据:** CONSERVED(两组 prevalence 均 ≥60%)→ SHIFTED(组成量 AML vs healthy 的 Mann-Whitney,方向指向 LSC)→ CONVERGE(autocrine 份额升高)。**pooled 显著 + within-study 同向**才计候选;并配打乱标签/打乱 sender 的零模型。脚本 `cooption_redistribution.py`(已存在,输出 `results/tables/05_ccc/cooption/`)。

---

## 4. `[INSERT 新模块 M11]` 下游汇聚(SCENIC + NicheNet)

**问题:** TIM-3 / CD99 / CD74 三条受体输入在 LSC 内是否汇聚到**同一个下游调控节点**?

- **SCENIC**:在 LSC(恶性 primitive)中推断 regulon,检验三条轴的下游是否落在**同一 regulon**;候选 hub 有先验 —— β-catenin/TCF、NF-κB、Src-family/ERK 效应子、干性 TF(HOXA9/MEIS1)。
- **NicheNet ligand→target**:Gal-9 / CD99 / MIF 是否预测**重叠的靶基因程序**(target 层的汇聚)。
- **LIANA+ = 跨方法对照,不作结论**(Daichi 明确)。
- **机制先验(非瞎捞):** 三条已发表下游本就重叠于 Src-family/ERK 与 β-catenin/NF-κB(见 `CCC_TOP_LIGANDS_LITERATURE.md` §1–3)。
- **判据:** 三轴共享的下游 regulon/target 富集显著超过置换零分布,且在 LSC 中比非 LSC 更强。

---

## 5. `[INSERT 新模块 M12]` 同一 LSC vs 代偿性 niche

三条互补检验,回答 Daichi 的 A/B 分叉:

1. **单细胞受体共表达:** 单个恶性 LSC 是否**同时**表达 HAVCR2(TIM-3)+ CD99 + CD74?(同一细胞 = 支持 A)。可在现有 scRNA 上直接做。
2. **恶性特异 secondary graph**(bin × 恶性态 × role):修标注后(VarTrix/Numbat 用于有 BAM 的样本,表达型 consensus 用于其余),检验 autocrine 环是否**由恶性细胞驱动**。
3. **in-silico 组合敲除:** 在 FGW/拓扑上**单边 vs 组合边**敲除,看去掉一条轴是否被其他轴**代偿**(冗余 = 支持 B,并预演湿实验的组合封锁)。

---

## 6. `[AMEND R10 / Phase 8]` 空间共定位 — 升为必需验证

`in silico` 只生成假设(Daichi)。新增必需项:**空间转录组 / 多重成像**证明这些 hub(LSC 与发送隔室)在**耐药骨髓**中**物理共定位**,而非仅在解卷积网络里相邻。

---

## 7. `[AMEND §9]` 实验验证设计 — 组合封锁,而非单抑制剂

单抑制剂(如 ISO-1 单用)**无法**验证拓扑重排/冗余假设。预注册 **ex vivo 共培养 + 组合封锁**(如 **anti-CD99 + MIF 抑制剂** ± anti-TIM-3),读出**代偿性通路激活**(某轴被封后其他轴/下游 hub 是否上调)。

---

## 8. 先导证据(已在手,来自 M10)

**发送 / 接收 / autocrine 组成量(中位数,pooled Mann-Whitney;within-study 见括注):**

| 通路 | primitive 发送 (AML/健康) | primitive 接收 | **autocrine LSC↔LSC** | pooled p (autocrine) |
|---|---|---|---|---|
| **CD99** | 0.62 / 0.18 | 0.60 / 0.06 | **0.34 / 0.00** | **3.6e-5** |
| **Galectin-9** | 0.47 / 0.19 | 0.34 / 0.24 | **0.16 / 0.045** | **2.9e-3** |
| **MIF** | 0.48 / 0.47 | 0.26 / 0.18 | 0.13 / 0.10 | 0.044(弱,阴性对照) |

**读法:** CD99 在 AML 里几乎凭空长出一个自分泌 LSC 环(0→0.34),Gal-9 次之(与 Kikushige 2015 的 TIM-3/Gal-9 autocrine 环吻合),MIF 基本不动(广泛 niche 因子)。**这把三个分子分到 Daichi 的两边:CD99/Gal-9 = 候选"同一 LSC 汇聚";MIF = 候选"代偿性 niche"。**

**必须带的 caveat(诚实):** (a) **within-study 全部未达 p<0.05**(autocrine bothDS:CD99 0.15、Gal-9 0.57、MIF 0.29 —— 3 数据集欠功效);(b) 仅 CellChat;(c) **bin 层级、非单细胞亚群**(M12 才是单细胞版);(d) pooled 受 study 混淆;(e) **治疗压力轴数据最薄**(见 §9)。因此这是**假设生成级**证据,须经 M11/M12 + §9 硬化。

---

## 9. 优先级与测序(诚实的数据现实)

- **诊断优先,压力其后。** "治疗压力下的汇聚"恰是我们纵向数据最薄处(H3 blocker:仅 GSE227903 有干净三联、6323_MRD 一患者占 62% MRD 细胞、M6 方向分数未建、配对患者 ~6)。因此**先在 diagnosis(121 个 AML 样本)钉死汇聚/autocrine,再向耐药/压力语境延伸**,不把整个故事压在薄弱纵向数据上。
- **可靠性门槛链(承 `CCC_TOP_LIGANDS_LITERATURE.md` §5):** 生物学可信度(✅ 文献已复现 split)→ present-node 混淆检查 → LIANA+/NicheNet 复现 → 零模型 + bootstrap → within-study / SRRS → 恶性特异 secondary graph。过 1–5 的边才交实验组。

---

## 10. 脚本占位(按 CODING_STANDARDS 命名)

```
scripts/10_figures/cooption_redistribution.py       # M10:sender/receiver/autocrine(已存在)
scripts/10_figures/fig_cooption_convergence.py      # M10:sender/receiver/autocrine 三联图(待建)
scripts/09_robustness/05_cooption_null.py           # M10 零模型(打乱标签/sender,待建)
scripts/11_convergence/01_scenic_regulons.R         # M11:LSC regulon 推断 + 三轴汇聚(待建)
scripts/11_convergence/02_nichenet_ligand_target.R  # M11:ligand→target 重叠(待建)
scripts/11_convergence/03_liana_crosscheck.R        # 跨方法对照(待建)
scripts/05_ccc/07_secondary_graph_malignant.R       # M12:bin×恶性态×role(依赖标注修复,待建)
scripts/05_ccc/08_receptor_coexpression.R           # M12:单细胞 TIM-3/CD99/CD74 共表达(待建)
scripts/09_robustness/06_edge_knockout_compensation.py  # M12:in-silico 单/组合敲除(待建)
```

---

## 11. v1.1 → v1.2 变更摘要(一页)

1. **中心命题上移**:配体 → **配体在 LSC 汇聚的共同下游 hub**;单靶失败 → **冗余/代偿**(§1)。
2. **H1 收敛性**扩为两层(L3 边收敛 + 下游 regulon 收敛)+ autocrine 环 + 同一 LSC/代偿分叉(§2)。
3. **M10 升为正式模块**,新增 receiver + autocrine 度量(§3;先导证据 §8)。
4. **新增 M11**(SCENIC/NicheNet 下游汇聚,LIANA+ 仅对照)与 **M12**(同一 LSC vs 代偿:单细胞共表达 + 恶性 secondary graph + in-silico 组合敲除)(§4–5)。
5. **空间共定位升为必需验证**(§6);**实验改为组合封锁**(§7)。
6. **诚实排期**:诊断优先、压力其后;within-study 尚未显著,证据为假设生成级(§8–9)。
