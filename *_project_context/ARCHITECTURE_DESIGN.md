# aml_niche_net — 最终架构设计 (Architecture Design)
 
> 全局重构蓝图。定义:目录树 + 编号方案(对齐 blueprint Phase)+ config 层级 + 恶性标签统一。
> 这是**方向文档**,审定后再产出"逐层执行手册 + 治理文档"。
> 原则:目录编号 = blueprint 的科学 Phase = 项目叙事。一个路径真相,一套风格,一套恶性标签。
 
---
 
## 1. 核心设计原则
 
1. **目录 = 科学阶段**:脚本目录编号直接对应 blueprint 的 Phase 1–8,读目录树即读项目逻辑。
2. **单一路径真相**:所有路径根定义在一个地方(R + bash 两个镜像),其余全部派生。
3. **分层 config**:`config_paths` (只有根) + 每阶段专属 config (源 paths + 加自己的常量)。
4. **一套恶性标签**:只有 c50 (union_mode 合并 d35 洞见);废弃 d35。
5. **一套风格**:ALL_CAPS 常量、`message("[N] ...")` 日志、`make_/read_/load_` helper、`----` section、`=`-form CLI、共享 utils、SEED=20260605。以现 `d00_config.R` 为模板。
6. **生成物与源码分离**:自动生成的执行脚本(如 c45_run_consensus.sh)不进源码目录。
---
 
## 2. Blueprint Phase → 目录编号映射
 
blueprint 定义 Phase 1–8。当前脚本目录 `01_qc / 02_preprocess / 05_cnv_snv / 06_hierarchy / 10_figures` 与之错位。目标映射:
 
| blueprint Phase | 科学内容 | 当前位置 | **新目录** |
|---|---|---|---|
| (前置) 数据摄入 | 原始→per-dataset merged RDS | `01_qc` | `00_ingest` |
| Phase 1 (part A) | per-sample QC + doublet + 角色/split | `02_preprocess` | `01_preprocess` |
| Phase 1 (part B) | 重比对→共识恶性注释 (CNV/SNV) | `05_cnv_snv` | `02_malignancy` |
| Phase 2 | 层级锚定分箱 (BMM 投影 + stemness) | `06_hierarchy` (d) | `03_hierarchy` |
| Phase 3 | per-sample cNMF + MP 聚合 + SRRS | `06_hierarchy` (e) | `04_cnmf` |
| Phase 4 | CCC 图构建 | (未来) | `05_ccc` |
| Phase 5 | 强度→距离 | (未来) | `06_distance` |
| Phase 6 | FGW 对齐 + barycenter | (未来) | `07_fgw` |
| Phase 7 | emergent edge + 患者评分 | (未来) | `08_scoring` |
| Phase 8 | 稳健性 / QC 对照 | `06_hierarchy` (f) | `09_robustness` |
| (贯穿) | 图表 | `10_figures` | `10_figures` |
 
> **决策点 A**:是否接受"数据摄入=00_ingest,Phase 1 拆成 01_preprocess + 02_malignancy"?
> 理由:blueprint 的 Phase 1 同时含"QC 策划"和"共识恶性注释"两大块,体量都很大,拆成两个目录比挤在一个更清晰。QC(01)是恶性注释(02)的上游,顺序自然。
> 备选:严格按 blueprint 保持 Phase 1 = 单目录 `01_phase1_qc_malignancy`(不拆)。我**倾向拆**(01/02),因为两块独立成文、独立运行。
 
---
 
## 3. 目标目录树
 
```
aml_niche_net/
├── scripts/
│   ├── config/                      # ← 分层 config (单一路径真相)
│   │   ├── config_paths.R           # 所有路径根 + 派生 (R 侧)
│   │   ├── config_paths.sh          # 同上镜像 (bash 侧, 供 02_malignancy 的 .sh)
│   │   ├── config_qc.R              # source paths + QC 阈值/gene sets
│   │   ├── config_malignancy.R      # source paths + inferCNV/Numbat/consensus 参数
│   │   ├── config_hierarchy.R       # source paths + 投影/stemness/bin 常量
│   │   ├── config_cnmf.R            # source paths + cNMF/SRRS 常量
│   │   └── utils.R                  # 共享 helper: get_counts / fwrite_safe / core16 / ...
│   │
│   ├── 00_ingest/                   # (旧 01_qc) 原始 → merged RDS
│   │   ├── 00a_common_readers.R     # (旧 00_common_qc.R 的 readers/make_seurat)
│   │   ├── ingest_<DATASET>.R       # 13 个 (旧 01_qc_<DATASET>.R)
│   │   └── ingest_MASTER_summary.R
│   │
│   ├── 01_preprocess/               # (旧 02_preprocess) Phase 1A
│   │   ├── 01_dataset_roles.R
│   │   ├── 02_study_split.R
│   │   ├── 03_per_sample_qc.R
│   │   ├── 04_export_role_split_xlsx.R
│   │   ├── 05_qc_landing_check.R    # header 修正为 05
│   │   ├── submit_03_per_sample_qc.sbatch
│   │   └── tests/                   # 90/91 挪这里 (或保持前缀, 见决策 B)
│   │       ├── 90_dryrun_chen2023.R
│   │       └── 91_smoke_test_all.R
│   │
│   ├── 02_malignancy/               # (旧 05_cnv_snv) Phase 1B
│   │   ├── c00_make_gene_order.R    # (旧 c00b)
│   │   ├── c01_fetch_srr_table.sh
│   │   ├── c02_build_samples_sheet.py   # (旧 c01b)
│   │   ├── c03_vet_read_structure.sh    # (旧 c02)
│   │   ├── c10_run_one_sample.sh    # 前置链: 下载→STARsolo→cellsnp→Numbat
│   │   ├── c11_submit_array.slurm   # (旧 c20)
│   │   ├── c20_ref_norm_identify.R  # (旧 c05) autologous normal ref
│   │   ├── c21_ref_norm_diagnose.R  # (旧 c05b)
│   │   ├── c30_run_numbat.R         # allele-CNV (含 --max_entropy)
│   │   ├── c40_run_infercnv.R       # expr-CNV
│   │   ├── c41_infercnv_to_percell.R
│   │   ├── c42_run_expr_cnv.R       # (旧 c44) copykat/scevan
│   │   ├── c43_author_to_percell.R
│   │   ├── c44_run_infercnv_one.R   # (旧 c46)
│   │   ├── c45_submit_infercnv.slurm # (旧 c47)
│   │   ├── c48_vartrix.R            # [PLANNED, 未激活] SNV 臂占位
│   │   ├── c50_consensus_malignancy.R  # ★ 唯一恶性标签 (+ --union_mode)
│   │   ├── c51_malignancy_plan.py   # (旧 c45_malignancy_plan.py) 生成计划
│   │   └── _generated/              # ← 自动生成的执行脚本 (不 commit, gitignore)
│   │       └── run_consensus.sh     # (旧 c45_run_consensus.sh)
│   │
│   ├── 03_hierarchy/                # (旧 06_hierarchy/d) Phase 2
│   │   ├── d10_project_bmm.R
│   │   ├── d15_derive_mapping_thresholds.R
│   │   ├── d20_assign_bins.R
│   │   ├── d25_stemness.R
│   │   ├── d30_per_bin_malignant.R
│   │   ├── d40_rollup.R             # (d35 已废弃并入 c50)
│   │   ├── bmm_bin_map.tsv
│   │   ├── dataset_platform.tsv
│   │   └── stemness_signatures.tsv
│   │
│   ├── 04_cnmf/                     # (旧 06_hierarchy/e) Phase 3
│   │   ├── e05_export_cnmf_input.R  # ← 改读 c50 输出 (03_malignancy)
│   │   ├── e10_run_cnmf.py
│   │   ├── e10_array.sbatch
│   │   └── e20_aggregate_programs.R
│   │
│   ├── 09_robustness/               # (旧 06_hierarchy/f) Phase 8 (QC 对照)
│   │   ├── f05_platform_deviation_control.R
│   │   └── f07_validate_numbat_rescue.R
│   │
│   ├── 10_figures/
│   │
│   └── docs/                        # ← 治理文档
│       ├── ARCHITECTURE.md
│       ├── CODING_STANDARDS.md
│       ├── PIPELINE.md              # 跑一个新数据集的完整顺序
│       └── NEW_DATASET_RUNBOOK.md   # c10 前置链 + 3 个 USER-CONFIRM
│
├── results/tables/                  # 输出按新阶段编号
│   ├── 00_ingest/                   # (旧 01_qc 的 percell/summary)
│   ├── 01_preprocess/               # (旧 01_qc 里混装的 manifest/split/qc_report)
│   ├── 02_malignancy/               # c50 consensus (唯一恶性标签)
│   ├── 03_hierarchy/
│   ├── 04_cnmf/
│   └── 09_robustness/
│
└── (LARGE1) 大对象目录同理按阶段编号
```
 
> **决策点 B**:测试脚本 90/91 —— 你之前说"靠命名区分,不挪目录"。上图我暂放 `tests/`,但可改回**同目录 + 90/91 前缀**。你定:(i) `tests/` 子目录,还是 (ii) 同目录留 90/91 前缀?
 
> **决策点 C**:`02_malignancy` 里 c 系列**是否重编号**?它跨 4 种语言、有前置链依赖。选项:(i) 全重编号成连续 c00–c51(上图),(ii) 保留现有 c 编号只补文档。我倾向 (i) 连续化,但 c 系列重编号会动很多 .sh/.slurm 内部引用,风险较高——**这是最该谨慎的一批**。
 
---
 
## 4. 分层 config 设计
 
### 4.1 单一路径真相 (两个镜像文件)
 
**`config/config_paths.R`** (R 侧):
```r
# 唯一路径根。改这里,全项目生效。
FAST_DIR   <- "/FAST/gr10634/gaozy/aml_niche_net"
LARGE1_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net"
REF_DIR    <- "/LARGE1/gr10634/gaozy/reference"
ENV_PREFIX <- "/FAST/gr10634/gaozy/general_env"
SEED       <- 20260605L
# 派生的阶段目录 (输出)
DIR_INGEST      <- file.path(FAST_DIR, "results/tables/00_ingest")
DIR_PREPROCESS  <- file.path(FAST_DIR, "results/tables/01_preprocess")
DIR_MALIGNANCY  <- file.path(FAST_DIR, "results/tables/02_malignancy")
DIR_HIERARCHY   <- file.path(FAST_DIR, "results/tables/03_hierarchy")
DIR_CNMF        <- file.path(FAST_DIR, "results/tables/04_cnmf")
# 大对象
OBJ_QC_DIR      <- file.path(LARGE1_DIR, "02_seurat_objects/01_per_sample_qc")
...
```
 
**`config/config_paths.sh`** (bash 侧, 供 02_malignancy 的 .sh/.slurm):
```bash
export FAST_DIR="/FAST/gr10634/gaozy/aml_niche_net"
export LARGE1_DIR="/LARGE1/gr10634/gaozy/aml_niche_net"
export REF_DIR="/LARGE1/gr10634/gaozy/reference"
export ENV_PREFIX="/FAST/gr10634/gaozy/general_env"
# ... 与 .R 镜像一致
```
> 两个文件手工保持同步(内容少、变动少)。或用一个 `.env` 两边 source——但 R source .env 需解析,反而复杂。**倾向两镜像 + 一句注释"改一处必改另一处"**。
 
### 4.2 阶段 config (source paths + 加常量)
 
每个阶段一个,风格统一为 `d00_config.R` 那样(ALL_CAPS + 决策注释 + 路径语义):
- `config_qc.R` — QC 阈值 (MAD/doublet/min_cells) + gene sets (MITO_SHORT 等)
- `config_malignancy.R` — inferCNV/Numbat 参数 + consensus tier 词表 + ref-norm 参数 (合并旧 c00_refnorm_config.R)
- `config_hierarchy.R` — 投影/stemness/bin/mapping-QC (旧 d00 的 Phase 2 部分)
- `config_cnmf.R` — cNMF/SRRS (旧 d00 的 Phase 3 部分)
### 4.3 共享 utils
 
`config/utils.R` 收纳现在重复/分散的 helper:
- `get_counts()` (v5/v4 兼容, 现在 d00 + c05 各有一版) → **统一一版**
- `fwrite_safe()` (现 d00)
- `core16()` (现 d30/d35/e05/f07/c50 **五处**, 两种实现!) → **统一一版, 加长度校验**
- `list_qc_samples()` / `list_datasets()`
- logging: 若保留时间戳需求,做 `msg <- function(...) message(sprintf("[%s] ...", ...))` 统一包装
---
 
## 5. 恶性标签统一 (c50 收编 d35)
 
### 5.1 目标:一个脚本、一套输出、一套 tier
 
废弃 `d35_infercnv_label.R`;`c50_consensus_malignancy.R` 成为唯一恶性标签生产者。
 
### 5.2 c50 要吸收的 d35 洞见
 
1. **`--union_mode` 选项**:任一 valid 证据类型阳性即恶性(对应 d35 并集/盲区互补)。默认仍是 type-level majority (≥2/3);union_mode 用于当前"两证据、盲区互补"场景。
2. **Numbat 双 schema 修复** [🔴 bug]:c50 现在对 degraded Numbat (`no_CNV_detected`) 当 `malignant=0`(误判"正常")。**必须改成 NA**(该样本 Numbat 无意见,不投票),否则污染共识。这是 d35 已正确处理、c50 待修的 bug。
3. **tier 词表统一**:合并两套 tier 命名,首字母对齐 `TIER2CONF (A=high/B=medium/C=low)`:
   ```
   A_concordant   : ≥2 证据类型 且 一致  (旧 c50 A_allelic_or_snv ∪ d35 A_cnv_concordant)
   B_multi_partial: ≥2 证据类型 但部分一致 (旧 B_expr_cnv_multi ∪ d35 B_cnv_union)
   C_single       : 单一证据类型          (旧 C_single_type ∪ d35 C_cnv_single)
   ```
4. **SNV 臂 [PLANNED]**:c50 保留 `--vartrix` 接口,注释 "planned; most datasets lack VarTrix; currently expr+allele two-evidence"。tier 逻辑不因 SNV 缺失报错。c48 占位脚本。
### 5.3 下游改动
 
- `e05_export_cnmf_input.R`:`--malignancy_dir` 默认从 `MALIGNANCY_INFERCNV_DIR` (d35) 改为 `DIR_MALIGNANCY` (c50)。
- `d30_per_bin_malignant.R`:`--malignancy_dir` 同上;它消费的 schema (`cell,malignant,confidence` + summary `evidence_tier`) c50 已兼容。
- `d00` 里 `MALIGNANCY_INFERCNV_DIR / TIER_CNV_* / PILOT_EVIDENCE_TIER` 删除(迁入 config_malignancy.R 的统一 tier)。
### 5.4 等价验证 (你已同意直接重跑)
 
废弃 d35 前,拿 pilot 的 1216_Dg / 3853_Dg / 6323_Dg(覆盖 Numbat 有效/inferCNV 盲区/Numbat 盲区三种情况)跑 c50 `--union_mode`,确认恶性标签与 d35 输出逐样本一致(malignant_frac、tier)。一致 → 废弃 d35 + 全量重跑 e05/e10/e20。
 
---
 
## 6. 逐层重构执行顺序 (你定的数据流顺序)
 
按 **manifest → config → 00_ingest → 01_preprocess → 02_malignancy → 03_hierarchy → 04_cnmf → 09_robustness**,每层改完立即验证(重跑该层、确认输出正确)再下一层。详细步骤在后续"逐层执行手册"。
 
| 批次 | 内容 | 验证方式 | 动数据? |
|---|---|---|---|
| 0 | 建 config/ (paths + utils + 阶段 config) + docs 骨架 | source 无报错 | 否 |
| 1 | 00_ingest: 重命名 + 改 source + Timepoint 词表 + uid_patient 上移 + QC min/max | 重跑 1 个数据集, 比对 RDS | 迁移 RDS |
| 2 | 01_preprocess: 改 source + seed 统一 + 文件名/header + 输出目录分离 | 重跑 manifest+split+QC | 迁移表 |
| 3 | 02_malignancy: 重编号 + c50 收编 d35 + Numbat schema 修复 + config 合并 | 等价验证(§5.4) | 迁移 CNV 产物(备份) |
| 4 | 03_hierarchy: 改 source + 删 d35 + core16 用 utils | 重跑 d10–d40 | 迁移表 |
| 5 | 04_cnmf: e05 改读 c50 + 改 source | 重跑 e05–e20, 比对 151 MP | 迁移 cNMF 产物 |
| 6 | 09_robustness + 10_figures + 文档定稿 | 重跑 f05/f07 | 否 |
 
> 每批"迁移数据"= mv 到新路径 + 校验数量/大小,**Numbat/inferCNV 结果保留备份**(最贵)。
 
---
 
## 7. 待你审定的决策点汇总
 
- **决策 A**:Phase 1 拆成 `00_ingest + 01_preprocess + 02_malignancy` 三目录(我倾向拆)?还是保持更少目录?
- **决策 B**:测试脚本 90/91 → `tests/` 子目录,还是同目录留前缀?
- **决策 C**:02_malignancy 的 c 系列是否连续重编号(动 .sh/.slurm 引用,风险较高)?还是保留现编号+补文档?
- **决策 D**:config_paths 用"R+sh 双镜像"?确认可接受手工同步?
- **决策 E**:总编号风格 —— 目录用 `00_/01_/.../10_`,目录内脚本保留字母前缀(ingest 无前缀 / preprocess 数字 / malignancy=c / hierarchy=d / cnmf=e / robustness=f)?还是目录内也统一?
审定这 5 个决策 + 目录树方向,我再产出**逐层执行手册 + CODING_STANDARDS.md + ARCHITECTURE.md 定稿**。