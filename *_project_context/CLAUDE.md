# AML Niche Network — 项目说明

13 个公共 AML / 健康骨髓 scRNA-seq 数据集的整合分析管线：
**摄入 → QC → 恶性细胞标注 → 造血层级投影(7 bin) → CellChat 细胞通讯图 → FGW 图对齐 → 统计检验。**

- 项目根：`/FAST/gr10634/gaozy/aml_niche_net`（脚本 + 表格 + 图，可清理层）
- 大对象：`/LARGE1/gr10634/gaozy/aml_niche_net`（Seurat 对象、CNV、cNMF，**785G**）
- 所有路径与全局 `SEED=491638` 由 `scripts/config/config_paths.R` 单点定义；
  集群端镜像 `config_paths.sh`。**两份必须同步修改。**
- 无 git，无正式文档。每个脚本头部的注释块是唯一的规格说明，改动前先读它。
  （例外：`.claude/settings.local.json` 是本仓库的 agent 权限配置；
  `00_project/tmp/mp_bin_check.txt` 与根目录的 `recon_*` 文件是一次性核查记录，不是规格。）

---

## 既定前提（不要重新论证）

### 1. 恶性标注以 inferCNV 为主线；Numbat 只覆盖少数数据集

**决策**：Numbat 需要 FASTQ/BAM 才能跑，而队列里**只有两三套数据满足条件**。
它因此天然只能覆盖一小部分样本，**不能作为全队列的证据来源**——
这是数据可得性的限制，不是 Numbat 本身不可靠。
恶性标注的主线是 inferCNV（唯一能覆盖全队列的臂）：

```
20_refnorm_identify.R → 44_infercnv_run_one.R (+45 array) → 41_infercnv_to_percell.R → 50_consensus_malignancy.R
```

`50_consensus_malignancy.R` 是**唯一的恶性标注器**（多臂 union 投票 + 证据等级）。

**但磁盘与代码的实际状态与这条决策并不完全一致，接手前必须知道：**

- **130 个样本中 125 个是 `arms=infercnv` / `tier=C_single`（单证据）。
  另外 5 个 GSE227903 样本的标签里已经含有 Numbat 证据**：
  `1886_Dg`、`1886_R`、`3853_Dg`、`6323_Dg` 为 `B_multi_partial`，`6323_R` 为 `A_concordant`。
  这 5 个样本的 `conflict_frac` 高达 0.14–0.76。
- **Numbat 臂在代码里没有关闭**：`run_consensus_all.sh:50` 只要发现 numbat percell 文件
  就会自动追加 `--numbat`；`50_consensus_malignancy.R:137` 把它作为 `allele_call` 送进 union 投票。
- **bamchain / Numbat 目前仍在跑**（截至 2026-07-24 仍有 SLURM array `bamchain`
  处于 R/PD 状态，`logs/bamchain.*` 覆盖 07-16 ~ 07-24）。
  GSE289435 已产出 5 个 `__numbat_percell.csv`（07-16 ~ 07-23），
  **但从未合入 consensus** —— 该数据集 12 个样本的标签仍全是 07-14 的 `arms=infercnv`。

> **已决（2026-07-24）**：bamchain array **不终止**，让它跑完——多几个样本有等位基因
> 证据是净收益。但 Numbat 不作为全队列策略。
>
> ⚠️ **仍待决**：GSE289435 那 5 个已产出的 Numbat 结果何时合入。
> 重跑 `50_consensus` 会改变该数据集标签并使全部下游（03–08）失效，
> 因此必须与 metadata 修复后的重跑**合并成一次**，不要单独触发。

**提升标签质量的下一步（已定，组会后执行）**：跑 `42_exprcnv_run.R`（CopyKAT/SCEVAN）。
它直接读 per-sample QC 的 `.rds`，**不需要 FASTQ/BAM**，因此是唯一能覆盖全队列的第二条
独立证据臂——可把大批样本从 `tier=C_single` 抬到多证据，直接针对 HSC_MPP 的高假阳性率。
这是 blueprint Phase 1「≥2/3 三方共识」在当前数据条件下最现实的近似。

其他臂的状态：
- `42_exprcnv_run.R`（CopyKAT/SCEVAN）：从未运行。
- `43_author_to_percell.R`：**跑过**（`/LARGE1/.../03_cnv_snv/author/` 下有 15 个
  `__author_percell.csv`，Chen2023 ×1 + GSE239721 ×14，06-19/06-25），
  但从未通过 `--author` 传给 `50_consensus`（日志里 130 个样本全是 `[author] absent`）。

恶性标签绝大多数是单方法结论，做下游解释时要按"单证据"对待。已有两个负控：
`96_malignancy_fpr_healthy.R`（健康供者假阳性率）、`94_diagnose_label_quality.R`（测序深度混杂）。
`51_malignancy_plan.py` 可扫描每个样本当前可达的 evidence tier。

### 2. niche 已重定义为造血-免疫微环境（不含基质）

**队列里有两套数据带基质成分，但它们无法构成对照**：

| 数据集 | 角色（`config_qc.R` 的 ROLE_TABLE 原文） | 基质细胞数 |
|---|---|---|
| Chen2023 | `Discovery-AML+AuxStroma` | ~6,284（19 样本） |
| GSE253355 | `Reference-scaffold+AuxStroma` | ~25,370（12 样本） |

**Chen2023 是唯一带基质的 AML 数据集；GSE253355 虽然基质更多，但是纯健康供者参考
（Bandyopadhyay 2024，`Role : healthy reference + auxiliary stroma`，全部 baseline 无时间点）。**
也就是说跨队列的「AML vs 健康」基质-造血通讯对比只有单臂，做不了。
因此项目已把 "niche" **重定义为造血-免疫微环境**，基质细胞不进入分析。

这条前提已经编码在代码里，不是口头约定：

- `scripts/03_hierarchy/bmm_bin_map.tsv`：`Stromal → in_ccc_graph = FALSE`
  （注释原文：`non-hematopoietic; NOT a CCC-graph node (B-layer add-on, Chen2023 later)`）
- `scripts/config/config_ccc.R:20`：`CCC_NODES` 固定为 7 个造血 bin —
  `HSC_MPP, LMPP_GMP, Mono_DC, Erythroid, Megakaryocyte, T_NK, B_Plasma`

FGW 的 7×7 图结构、barycenter、08_scoring 的 49 条边检验全部建立在这 7 个节点上。
**增删节点会使 06/07/08 三个阶段的全部产物失效。**

`scripts/判断样本MSC含量.r` 是当初做出这个判断的那次筛查脚本（一次性，非管线步骤）。

---

## 目录说明

### `scripts/` — 全部代码（205 个文件）

活跃代码约 103 个；其余 102 个在 `以前06_hierarchy/` 归档里（其中 85 个是旧 SLURM 日志）。

| 目录 | 内容 |
|---|---|
| `config/` | 10 个配置库，只被 `source()`，从不单独运行。`config_paths.R`/`.sh` 是路径唯一真相源 |
| `00_ingest/` | 13 个 `ingest_<数据集>.R` + 共享读取器 `00_common_readers.R`（**真正写文件的地方**）+ MASTER 汇总 |
| `01_preprocess/` | 角色清单 → 70/30 研究级拆分 → `03_per_sample_qc.R`（核心，SLURM array 13 任务）→ 落地核查 |
| `02_malignancy/` | **24 个文件**，最复杂的一层。审计/规划工具：`51_malignancy_plan.py`（扫描各样本可达 tier 并生成 `_generated/run_consensus.sh`）、`91_check_infercnv.R`、`92_check_numbat.R`、`93_check_bams.sh` |
| `03_hierarchy/` | BoneMarrowMap 投影 → per-bin 恶性率 → 分布漂移检验 → 干性打分。严格线性 01→02→03，04 挂在 01 上 |
| `04_cnmf/` | GeneNMF 元程序（分支 A）。`mp_labels.tsv` 是手工标注表 |
| `05_ccc/` | CellChat 通讯图（分支 B，主线）。`02_run_cellchat.R` 走 SLURM array |
| `06_distance/` | 单脚本：LR 张量 → 有向 7×7 边权 → 样本内 rank 距离 `C = 1 - rank` |
| `07_fgw/` | FGW 图对齐（R 建输入 → Python 算 barycenter/HDS/ATS → H3 复发检验） |
| `08_scoring/` | **11 个文件**（9 个 `.py` + 2 个 `.sbatch`），统计推断层，**当前活跃前沿** |
| `10_figures/` | 1 个共享 config（`f00_fig_config.R`）+ 5 个出图脚本（含 `f03_*_test.R` 试验版），**已写但从未跑过** |
| `以前06_hierarchy/` | **归档旧版**（"以前" = previous），已被 `03_hierarchy/` + `04_cnmf/` 取代。不要修改，不要放进执行顺序 |

顶层两个杂项：`00.0_Job_submit.sh`（SLURM/tmux 命令备忘，**不是驱动脚本**）、
`判断样本MSC含量.r`（一次性 MSC 筛查，见上文前提 2）。

### 散落的 `recon_*` 探针（一次性诊断，非管线步骤）

`00_ingest/recon_GSE185381_structure.R`、`05_ccc/recon_ccc_probe.R`、`05_ccc/recon_node_feasibility.R`、
`07_fgw/recon_fgw_probe.py`、`07_fgw/recon_fgw_directed_probe.py`。

**它们的输出写在项目根目录而不是 `results/tables/`**：`recon_node_counts.csv`(100K)、
`recon_ccc_probe.txt`、`recon_fgw_probe.txt`、`recon_fgw_directed_probe.txt`、
`recon_GSE185381_{donor_qc,donor_span,library_mix}.csv`。
不要把它们当阶段产物解释，也不要放进执行顺序。

### `results/tables/` — 各阶段产物

`00_ingest/` `01_preprocess/` `02_malignancy/` `03_hierarchy/` `04_cnmf/`
`05_ccc/` `06_distance/` `07_fgw/` `08_scoring/` `09_robustness/`(空) `side_tp53_pvrl4/`

管线当前头部是 `results/tables/08_scoring/`。
事后负控产物：`02_malignancy/{sample_quality_vs_malignancy,malignancy_fpr_healthy,malignancy_fpr_by_bin}.csv`、
`00_ingest/timepoint_rule0_dryrun.csv`。

**注意：`03_hierarchy` 的 per-cell 产物不在这里** —— `__bmm_percell.csv` /
`__stemness_percell.csv` 写在 `/LARGE1/.../02_seurat_objects/03_bmm_projected/`（见下）。

### `results/figures/`

只有 `01_preprocess/`（13 张 png）和 `04_cnmf/`（12 个文件 = 5 组 pdf+png ＋ 2 张 png）有内容。
`00_ingest/` 和 `02_malignancy/` 是**空目录**（占位，不代表出图脚本跑过）。

### 日志

`logs/`（344 个, 11M）+ `00_project/logs/`（26 个, 3.7M）= 370 个 / 14M。
另有 85 个旧日志在 `scripts/以前06_hierarchy/logs/`（归档，忽略）。

---

## 不要读取的大数据目录

**以下路径一律不要用 Read 打开，也不要 `cat` 其中的文件。**
需要了解内容时，用 `ls` / `du` / `find -printf` 看文件名、大小、mtime，或读脚本头部的 schema 说明。

### 绝对禁止（单文件即可撑爆上下文）

| 路径 | 大小 | 说明 |
|---|---|---|
| `/LARGE1/gr10634/gaozy/aml_niche_net/` | **785G** | 整棵树。下面逐项穷举 |
| `└─ 05_cnv_snv/` | 374G | inferCNV 运行目录 + burden |
| `└─ 03_cnv_snv/` | 365G | Numbat / cellsnp / STARsolo 工作区（含 BAM）；另有 `author/` 15 个 percell csv |
| `└─ 00_raw/` | 23G | 公共数据原始下载（GEO/Zenodo/ArrayExpress） |
| `└─ 01_processed_counts/rds/` | 8.8G | 13 个合并 Seurat 对象 |
| `└─ 06_cnmf/` | 6.7G | 旧版 cNMF 输入/运行（legacy） |
| `└─ 02_seurat_objects/` | 4.7G | 220 个 per-sample QC 对象 + **350 个纯文本 CSV**（见下方警告） |
| `└─ 05_ccc_graphs/` | 2.6G | 148 个 CellChat 对象 |
| `└─ reference/` | 200M | SingleR/inferCNV 外部参考。**`gencode_GRCh38_gene_order.txt` 是 1.2M 纯文本，扩展名清单挡不住它，同样不要读** |
| `└─ 04_cnmf/` | 309M | NMF 结果 `nmf_res.rds` |
| `└─ objects/`、`└─ 99_archive/` | 空 | 列出以示穷举 |
| `/FAST/.../SRR32323369/` | **30G** | `MLL_16703.bam` — 遗留下载，不是管线产物 |
| `/FAST/.../SRR32323368/` | 2.4G | `MLL_17746.bam.tmp` — **未完成的下载**（有 .lock） |

### 按扩展名跳过（无论在哪个目录）

`.rds` `.RDS` `.h5ad` `.h5` `.mtx` `.bam` `.bam.tmp` `.npz`

> ⚠️ **扩展名清单只挡二进制。以下是纯文本但同样不能读：**
> - `/LARGE1/.../02_seurat_objects/03_bmm_projected/` 下 **350 个 `.csv`**
>   （`__bmm_percell.csv` / `__stemness_percell.csv`，单个最大 **2.3M**）
> - `/LARGE1/.../reference/gencode_GRCh38_gene_order.txt`（1.2M，约 6 万行）

### FAST 侧需要谨慎的（不是二进制，但很大）

| 路径 | 大小 | 处理方式 |
|---|---|---|
| `results/tables/00_ingest/*_qc_percell.csv.gz` | 合计 90M | 单个最大 **49M**（GSE185381）。只读 `_qc_summary.csv`，不要解压 percell |
| `results/tables/02_malignancy/` 的 9 个子目录 | 49M / 265 个文件 | 单个最大 1.4M。不要遍历 per-sample，读根目录 4 个汇总表即可 |
| `results/tables/05_ccc/ccc_edge_distance.csv` | **728K** | 已知问题 5 的重复副本 |
| `results/tables/05_ccc/ccc_node_features.csv` | 104K | 需要时用 `head`/`awk` 取样 |
| `results/tables/05_ccc/tensors/` | 148 个 LR 张量 | 只读单个样本，不要 glob 全部 |
| `results/tables/06_distance/edge_distance.csv` | **728K** | 只能取样 |
| `results/tables/07_fgw/fgw_edges_long.csv` | **447K** | 只能取样 |
| `results/tables/07_fgw/fgw_nodes_long.csv` | 135K | 只能取样 |
| `logs/` | 344 个 / 11M | 按前缀分组统计；单个日志最多读头部 30 行 |
| `00_project/logs/` | 26 个 / 3.7M | **单文件最大 521K**（`01_qc_3622821_*.err`）。绝不整读，只用 `tail -30` / `grep` |

### 安全的读取入口（都是小表，可直接读）

```
results/tables/{01_preprocess,03_hierarchy,04_cnmf,08_scoring}/*.csv
results/tables/00_ingest/*_qc_summary.csv          # 不含 *_qc_percell.csv.gz
results/tables/02_malignancy/*.csv                 # 仅根目录 4 个表，全部 <11K，不递归
results/tables/06_distance/edge_qc.csv             # 9K
results/tables/07_fgw/{fgw_input_index,patient_scores,paired_rls_scs,rls_grouped}.csv
results/tables/05_ccc/ccc_edge_qc.csv              # 9K
```

---

## 执行顺序

**没有任何 SLURM `--dependency`/`afterok`。阶段之间完全靠手工按顺序执行。**
真正串联多脚本的驱动只有两个：`run_consensus_all.sh`（41→50）和
`10_prealign_one_sample.sh`（下载→STARsolo→numbat pileup_and_phase→30_numbat_run.R，
由 `11_submit_array.slurm:74` 提交）。另有 `51_malignancy_plan.py` 生成（但不执行）一份 41/42/43→50 的计划。

```
[1] 00_ingest/ingest_<13数据集>.R → ingest_MASTER_summary.R ;（旁支审计）95_test_timepoint_rule0.R
[2] 01_preprocess: 01 → 02 → 04 ; (预检 91/90) → 03_per_sample_qc.R [sbatch array] → 05
        ↓ per-sample QC .rds  ← 全项目中心枢纽，7 个下游阶段都读它
[3] 02_malignancy: 01_gene_order ; 20_refnorm → 90_preflight → 44+45[array] → 91_check ↺
                   → run_consensus_all.sh (41 → 50) → 60_rollup
                   （旁支）11+10[array] → 30_numbat  ← 仍在跑，见前提 1
                   （事后负控，可独立重跑）94_diagnose_label_quality / 96_malignancy_fpr_healthy
[4] 03_hierarchy:  01_bmm_project → 02_per_bin → 03_shift ; 01 → 04_stemness
[5A] 04_cnmf:  01(normal) → 02(normal) → 01(malignant) → 02(malignant) → 03 → 04
[5B] 05_ccc:   01_define → 02_run_cellchat[array] → 03_node_features
[6] 06_distance/01  →  [7] 07_fgw: 01(R) → 02(py) → 03 → 04
[8] 08_scoring: 01 → 02 → 03 → 04 → 05 → 06 → 07/08/09   ← 当前前沿
```

### 容易踩的坑

- **`02_malignancy` 的编号不是执行顺序**：`20` 必须先于 `40/44`；
  `30_numbat_run.R` 从不直接调用（被 `10_prealign_one_sample.sh:272` 内部调用）；
  `02/03/04/10/11/30` 是并行的 FASTQ 旁支，不是主线的一步。
- **`04_cnmf/02` 有自依赖**：normal 那次运行的输出是 malignant 那次的输入，顺序不能反。
- **`08_scoring` 的编号是叙事顺序，不是数据依赖**：除 `01→05` 外，
  `02/03/04/06/07` 只读 `07_fgw` 的三个 CSV；**`08/09` 还额外读 `06_distance/edge_distance.csv`**
  （只用 nodes/index，不读 edges）。这些脚本彼此可独立并行重跑。
- **`07_fgw/02_fgw_align.py` 不 source `config_fgw.R`**，而是内联重声明同一批常量
  （注释写着 "keep in sync"）。改 FGW 参数要改两处。
- `config_paths.R` 的 `PROJ_OBJ_DIR` 指向 `04_bmm_projected` 是**死常量**；
  真实路径是 `config_hierarchy.R` 里的 `03_bmm_projected`。

---

## 已知问题（截至 2026-07-24，尚未处理）

1. **`results/tables/02_malignancy/ref_norm_summary.csv` 在磁盘上不存在**，
   而 `40/44/91` 三个脚本 `stopifnot` 它，`90_preflight` 则检查后 `quit(status=1)`
   并提示 "run 20_refnorm_identify.R first"。
   路由表本身可以用 `20_refnorm_identify.R` 重建（输入齐全：220 个 QC rds + 106 个
   `ref_norm_cells.txt`），但 QC 对象在 07-15 已刷新，**重建出的路由表不再等于
   07-14 标签所依据的那一份 —— 丢失的是原次运行的可复现性。**
2. **时间戳与阶段编号矛盾**：`02_malignancy` 的 per-sample 共识产物是 07-14
   （rollup `ALL_consensus_summary.csv` 是 07-17，负控 CSV 是 07-22），
   而 `00_ingest`/`01_preprocess` 是 07-15 —— ingest 与 QC 在恶性标注**之后**重跑过，
   磁盘上的恶性标签对应的是旧一批 QC 细胞。`90_preflight_infercnv.R` 的注释确认这是已知情况。
3. **GSE289435 有 5 个新 Numbat 结果未合入**（07-16 ~ 07-23），
   但 `50_consensus` 从未重跑：该数据集 12 个样本的 `__consensus_percell.csv` 仍全是 07-14 18:12。
   因此**下游（03/04/05/06/07/08）与磁盘上的标签是一致的**；
   不一致的是「已产出但未合入」的 Numbat 证据。见前提 1 的待决项。
4. **队列不完整**：13 个数据集都有 ingest + QC（220 个 QC 对象），
   但只有 9 个数据集 / 130 个样本有恶性共识（缺 E-MTAB-11536、GSE147989、GSE185991、GSE253355）。
   覆盖数对照：220 个 `__bmm_percell` / 130 个 `__stemness_percell` / 148 个 CCC 张量。
   **stemness = 130 不是缺口**：`04_stemness_score.R` 取 QC ∩ consensus ∩ projection 的交集，
   与那 130 个有 consensus 的样本逐一对应。CCC 的 148 则不以 consensus 为门控，是另一套筛选。
5. **`results/tables/05_ccc/ccc_edge_distance.csv` 与 `ccc_edge_qc.csv` 是两份重复副本**，
   与 `06_distance/edge_distance.csv`、`edge_qc.csv` 字节完全相同（md5 一致），
   mtime 反而早 19 分钟，全树无脚本引用这两个文件名。可安全删除。
6. **`scripts/10_figures/` 从未运行**：不 source 项目 config，硬编码的是旧版路径
   （`results/tables/01_qc/`、`03_malignancy/`），当前会在 `stopifnot` 处直接失败。
   另外 `03_qc_report__ALL.csv` 也不存在（只有 13 个 per-dataset 版本）。
7. **`results/tables/09_robustness/` 是空目录，且没有对应的脚本目录**。
8. **`results/tables/side_tp53_pvrl4/` 是孤儿产物**，全树 grep 不到生产它的脚本。

---

## 工作约定

- 改配置常量时先确认它是否在别处被重复声明（`config_paths.R`/`.sh` 成对；
  `config_fgw.R` 与 `02_fgw_align.py` 成对；`config_cnmf.R` 复制了 `HIER_PROJ_DIR`）。
- 新增/修改脚本时保持现有的头部注释块风格：目的、INPUT、OUTPUT、调用方式、WHY。
  这是本项目唯一的文档。
- 重跑任何阶段前，先看该阶段脚本是否有 `--force` 标志与"输出已存在则跳过"的逻辑，
  多数阶段是 resume-safe 的。
- **重跑 `50_consensus_malignancy.R` 会使全部下游失效**（03/04/05/06/07/08），
  且会改变现有标签。动它之前先读前提 1 的待决项。
