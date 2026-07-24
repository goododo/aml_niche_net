# aml_niche_net — 逐层重构执行手册 (Refactor Execution Handbook)
 
> 配套 `ARCHITECTURE_DESIGN.md`。这里是**可执行步骤 + 完整重命名映射表**。
> 决策已锁定:A(Phase1 拆三目录)、B(90/91 同目录留前缀)、C(c 系列重编号)、
> D(config R+sh 双镜像)、E1(每目录内从 01 起,废字母前缀)。
> 执行顺序:config → 00_ingest → 01_preprocess → 02_malignancy → 03_hierarchy → 04_cnmf → 09_robustness。
> 每批:重命名 + 改 source + 迁移数据(校验) + 重跑 + 验证,通过后再下一批。
 
---
 
## ★ 完整重命名映射表 (旧 → 新)
 
> 按此表 `mv`。**注意:改名后必须同步改脚本内部的 `source()` / 调用路径 / header**(见各批"改 source"步骤)。
> 表内"内部引用"列标注该文件是否被别的脚本硬编码引用(改名后要一起改)。
 
### 批 0 — config/ (新建,非重命名)
| 新文件 | 来源 |
|---|---|
| `config/config_paths.R` | 新建 (抽取自 d00 的路径根 + 各 config 的路径) |
| `config/config_paths.sh` | 新建 (镜像 c00_cluster_config.sh 的 export) |
| `config/config_qc.R` | 抽取自旧 02_preprocess/00_config.R (QC 阈值/doublet/gene sets) |
| `config/config_malignancy.R` | 合并旧 c00_refnorm_config.R + c50 tier 词表 + inferCNV/Numbat 参数 |
| `config/config_hierarchy.R` | 抽取自旧 d00_config.R (Phase 2 部分) |
| `config/config_cnmf.R` | 抽取自旧 d00_config.R (Phase 3 部分) |
| `config/utils.R` | 收纳 get_counts/fwrite_safe/core16(统一)/list_qc_samples/list_datasets |
 
### 批 1 — 00_ingest/ (旧 01_qc/)
| 旧 | 新 |
|---|---|
| `01_qc/00_common_qc.R` | `00_ingest/00_common_readers.R` |
| `01_qc/01_qc_Chen2023.R` | `00_ingest/ingest_Chen2023.R` |
| `01_qc/01_qc_E-MTAB-11536.R` | `00_ingest/ingest_E-MTAB-11536.R` |
| `01_qc/01_qc_GSE116256.R` | `00_ingest/ingest_GSE116256.R` |
| `01_qc/01_qc_GSE147989.R` | `00_ingest/ingest_GSE147989.R` |
| `01_qc/01_qc_GSE185381.R` | `00_ingest/ingest_GSE185381.R` |
| `01_qc/01_qc_GSE185991.R` | `00_ingest/ingest_GSE185991.R` |
| `01_qc/01_qc_GSE201966.R` | `00_ingest/ingest_GSE201966.R` |
| `01_qc/01_qc_GSE207356.R` | `00_ingest/ingest_GSE207356.R` |
| `01_qc/01_qc_GSE227903.R` | `00_ingest/ingest_GSE227903.R` |
| `01_qc/01_qc_GSE239721.R` | `00_ingest/ingest_GSE239721.R` |
| `01_qc/01_qc_GSE253355.R` | `00_ingest/ingest_GSE253355.R` |
| `01_qc/01_qc_GSE289435.R` | `00_ingest/ingest_GSE289435.R` |
| `01_qc/01_qc_Petti2019.R` | `00_ingest/ingest_Petti2019.R` |
| `01_qc/01_qc_MASTER_summary.R` | `00_ingest/ingest_MASTER_summary.R` |
 
> **内部引用**:每个 ingest_*.R 都 `source(".../01_qc/00_common_qc.R")` → 改为 source 新的 `00_common_readers.R` + `config/config_paths.R` + `config/config_qc.R` + `config/utils.R`。
 
### 批 2 — 01_preprocess/ (旧 02_preprocess/)
| 旧 | 新 |
|---|---|
| `02_preprocess/00_config.R` | 拆解 → `config/config_qc.R` + `config/config_paths.R` (不再单独存在) |
| `02_preprocess/01_dataset_roles.R` | `01_preprocess/01_dataset_roles.R` |
| `02_preprocess/02_study_split.R` | `01_preprocess/02_study_split.R` |
| `02_preprocess/03_per_sample_qc.R` | `01_preprocess/03_per_sample_qc.R` |
| `02_preprocess/04_export_role_split_xlsx.R` | `01_preprocess/04_export_role_split_xlsx.R` |
| `02_preprocess/05_qc_landing_check.R` | `01_preprocess/05_qc_landing_check.R` (header 从 "04_" 改回 "05_") |
| `02_preprocess/90_dryrun_chen2023.R` | `01_preprocess/90_dryrun_chen2023.R` (前缀保留, 决策 B) |
| `02_preprocess/91_smoke_test_all.R` | `01_preprocess/91_smoke_test_all.R` |
| `02_preprocess/submit_03_per_sample_qc.sbatch` | `01_preprocess/submit_03_per_sample_qc.sbatch` |
 
> 编号内部已连续(01–05),无需改号。主要改动:source 换成 config/ 分层;`.log()` → `message("[N]")`;seed 引用统一常量;输出目录从 `results/tables/01_qc/` 改为 `results/tables/01_preprocess/`。
 
### 批 3 — 02_malignancy/ (旧 05_cnv_snv/, c 系列重编号)
| 旧 | 新 | 内部引用 (改名后要改) |
|---|---|---|
| `c00_cluster_config.sh` | → `config/config_paths.sh` (合并) | 被 c10 source |
| `c00_refnorm_config.R` | → `config/config_malignancy.R` (合并) | 被 c05 source |
| `c00b_make_gene_order.R` | `02_malignancy/01_make_gene_order.R` | — |
| `c01_fetch_srr_table.sh` | `02_malignancy/02_fetch_srr_table.sh` | — |
| `c01b_build_samples_sheet.py` | `02_malignancy/03_build_samples_sheet.py` | — |
| `c02_vet_read_structure.sh` | `02_malignancy/04_vet_read_structure.sh` | — |
| `c10_run_one_sample.sh` | `02_malignancy/10_run_one_sample.sh` | source config_paths.sh; 调 c30→新号 |
| `c20_submit_array.slurm` | `02_malignancy/11_submit_array.slurm` | 调 10_run_one_sample.sh |
| `c05_ref_norm_identify.R` | `02_malignancy/20_ref_norm_identify.R` | source config_malignancy.R |
| `c05b_ref_norm_diagnose.R` | `02_malignancy/21_ref_norm_diagnose.R` | — |
| `c30_run_numbat.R` | `02_malignancy/30_run_numbat.R` | 被 10_run_one_sample.sh 调 |
| `c40_run_infercnv.R` | `02_malignancy/40_run_infercnv.R` | 被 44/45 调 |
| `c41_infercnv_to_percell.R` | `02_malignancy/41_infercnv_to_percell.R` | 被 51 生成的脚本调 |
| `c43_author_to_percell.R` | `02_malignancy/43_author_to_percell.R` | 被 51 调 |
| `c44_run_expr_cnv.R` | `02_malignancy/42_run_expr_cnv.R` | 被 51 调 (copykat/scevan) |
| `c46_run_infercnv_one.R` | `02_malignancy/44_run_infercnv_one.R` | 被 45 调 |
| `c47_submit_infercnv.slurm` | `02_malignancy/45_submit_infercnv.slurm` | 调 44 |
| (新) VarTrix 占位 | `02_malignancy/48_vartrix.R` [PLANNED] | — |
| `c50_consensus_malignancy.R` | `02_malignancy/50_consensus_malignancy.R` (+ union_mode + schema 修复) | 被 51 生成的脚本调 |
| `c45_malignancy_plan.py` | `02_malignancy/51_malignancy_plan.py` | 生成 _generated/run_consensus.sh |
| `c45_run_consensus.sh` | `02_malignancy/_generated/run_consensus.sh` (生成物, gitignore) | — |
 
> **这是最易错的一批**:10/11 调 30;44/45 互调;51 生成的脚本调 41/43/42/50。改号后**逐一核对内部引用**。建议改完先 `grep -rn "c30\|c40\|c41\|c43\|c44\|c46\|c50" 02_malignancy/` 确认没有残留旧号。
 
### 批 4 — 03_hierarchy/ (旧 06_hierarchy/d, 废字母前缀)
| 旧 | 新 |
|---|---|
| `d00_config.R` | → `config/config_hierarchy.R` (Phase2 部分) + `config/config_cnmf.R` (Phase3 部分) |
| `d10_project_bmm.R` | `03_hierarchy/01_project_bmm.R` |
| `d15_derive_mapping_thresholds.R` | `03_hierarchy/02_derive_mapping_thresholds.R` |
| `d20_assign_bins.R` | `03_hierarchy/03_assign_bins.R` |
| `d25_stemness.R` | `03_hierarchy/04_stemness.R` |
| `d30_per_bin_malignant.R` | `03_hierarchy/05_per_bin_malignant.R` (--malignancy_dir → c50 输出) |
| `d35_infercnv_label.R` | **废弃** (并入 02_malignancy/50) |
| `d40_rollup.R` | `03_hierarchy/06_rollup.R` |
| `bmm_bin_map.tsv` | `03_hierarchy/bmm_bin_map.tsv` |
| `dataset_platform.tsv` | `03_hierarchy/dataset_platform.tsv` |
| `stemness_signatures.tsv` | `03_hierarchy/stemness_signatures.tsv` |
 
### 批 5 — 04_cnmf/ (旧 06_hierarchy/e, 废字母前缀)
| 旧 | 新 |
|---|---|
| `e05_export_cnmf_input.R` | `04_cnmf/01_export_cnmf_input.R` (读 c50 输出 DIR_MALIGNANCY) |
| `e10_run_cnmf.py` | `04_cnmf/02_run_cnmf.py` |
| `e10_array.sbatch` | `04_cnmf/02_run_cnmf.sbatch` |
| `e20_aggregate_programs.R` | `04_cnmf/03_aggregate_programs.R` |
 
### 批 6 — 09_robustness/ (旧 06_hierarchy/f, 废字母前缀)
| 旧 | 新 |
|---|---|
| `f05_platform_deviation_control.R` | `09_robustness/01_platform_deviation_control.R` |
| `f07_validate_numbat_rescue.R` | `09_robustness/02_validate_numbat_rescue.R` |
| `c45_malignancy_plan.tsv` | (生成物,不进源码;放 _generated/ 或 results/) |
 
---
 
## 逐批执行步骤
 
### 批 0 — config/ (不动数据,先夯地基)
 
1. 建 `scripts/config/`。
2. 写 `config_paths.R` + `config_paths.sh`(双镜像,内容见 ARCHITECTURE_DESIGN §4.1)。顶部注释:"改一处必改另一处"。
3. 写 `utils.R`:把 d00 的 `get_counts/fwrite_safe`、统一版 `core16`(加长度校验)、`list_qc_samples/list_datasets` 收进来。
4. 拆 d00 → `config_hierarchy.R` + `config_cnmf.R`;拆旧 02 的 00_config.R → `config_qc.R`;合并 c00_refnorm → `config_malignancy.R`。
5. **验证**:`Rscript -e 'source("config/config_paths.R"); source("config/config_qc.R"); source("config/utils.R"); cat("OK\n")'` 无报错;`bash -c 'source config/config_paths.sh; echo $FAST_DIR'` 正确。
6. 产出 docs/ 骨架(4 个 .md 空文件占位)。
### 批 1 — 00_ingest/
 
1. 按映射表 `mv` 15 个文件到 `00_ingest/`。
2. 每个 ingest_*.R:`source` 改为 `config/config_paths.R + config_qc.R + utils.R + 00_common_readers.R`;删除各自的 `SCRIPTS_DIR + dir.create` 样板(改由 config 处理)。
3. `00_common_readers.R`:make_seurat 里加 `seu$uid_patient <- paste0(dataset, ":", patient)`(C3 上移);加 Disease_state 默认值(C4);add_qc_metrics 加三 pattern 全空的 warning(A1);normalize_symbols 加 anyDuplicated 检查(A2);summarise_qc 加 min/max(反馈 #4)。
4. Timepoint 受控词表(C1):在 config 定义 `CANONICAL_TIMEPOINTS`,各 ingest 脚本映射到它。
5. **迁移数据**:旧 `results/tables/01_qc/*_qc_{percell,summary}.csv` → `results/tables/00_ingest/`;旧 RDS 路径不变(已在 LARGE1/01_processed_counts)。
6. **验证**:重跑 1 个数据集(如 Petti2019,小),确认新 RDS 的 meta 有 uid_patient、Timepoint 是受控值、QC 表有 min/max。
### 批 2 — 01_preprocess/
 
1. `mv` 9 个文件。
2. source 换 config 分层;`.log()` → `message("[N]")`(或 utils 里的时间戳包装);删 `list_datasets` 里 `setdiff(rds_exclude_files)` 那行(#2);seed 用统一常量。
3. 05_qc_landing_check.R header "04_"→"05_"(B1);输出目录 `01_qc/`→`01_preprocess/`(B2)。
4. **迁移**:manifest/split/qc_report 表 → `results/tables/01_preprocess/`。
5. **验证**:重跑 01_dataset_roles + 02_study_split + 03(1 数据集),确认 manifest/split 正确、QC 报告落新目录。
### 批 3 — 02_malignancy/ (最谨慎)
 
1. `mv` + 重编号(见表);`c00_*` 并入 config;`c45_run_consensus.sh` → `_generated/`。
2. **c50 收编 d35**:加 `--union_mode`;修 Numbat degraded schema(→NA 非 0);统一 tier 词表(A_concordant/B_multi_partial/C_single);SNV 臂标 planned。
3. 改所有内部引用(10→30、44↔45、51→41/43/42/50);`grep` 核对无残留旧号。
4. config_paths.sh 的 3 个 USER-CONFIRM(STAR_INDEX/NB_SNP_VCF/STAR_TMP_BASE)写进 NEW_DATASET_RUNBOOK.md(A3)。
5. **迁移(备份!)**:Numbat/inferCNV 产物 mv 到新路径,**旧目录保留为备份**(最贵,20h/样本);校验数量/大小。
6. **等价验证**:1216_Dg/3853_Dg/6323_Dg 跑 `50_consensus --union_mode`,比对旧 d35 输出(malignant_frac、tier)。一致 → 废弃 d35。
7. **验证**:全量重跑 c50(经 51 生成的脚本),确认恶性标签与 pilot 一致。
### 批 4 — 03_hierarchy/
 
1. `mv` + 废 d 前缀(见表);删 d35。
2. source 换 config_hierarchy;core16 用 utils;d30 的 `--malignancy_dir` 默认 → `DIR_MALIGNANCY`(c50)。
3. **迁移**:percell/perbin/thresholds 表 → `results/tables/03_hierarchy/`。
4. **验证**:重跑 01_project_bmm → 06_rollup(1 数据集),确认 bin/stemness/malignant 标签正确,rollup 总表对得上。
### 批 5 — 04_cnmf/
 
1. `mv` + 废 e 前缀。
2. `01_export_cnmf_input.R`:读 c50 输出(不再是 d35);source config_cnmf;e10.py 顶部常量注释同步。
3. **迁移**:cNMF input/runs/programs → 新路径。
4. **验证**:重跑 e05→e10→e20,确认 **1599 programs / 151 MP / 9 跨数据集** 复现(这是 pilot 金标准,数字必须对得上)。
### 批 6 — 09_robustness/ + 收尾
 
1. `mv` f05/f07 → 09_robustness/(废 f 前缀)。
2. source 换 config;core16 用 utils。
3. **验证**:重跑 f05(平台对照,确认 T_NK baseline 数字对得上)、f07。
4. 定稿 4 份 docs;`10_figures` 检查 source。
---
 
## 验证总清单 (每批的"金标准"数字)
 
| 批 | 金标准(重构后必须复现) |
|---|---|
| 1 | Petti2019 RDS 有 uid_patient + 受控 Timepoint |
| 2 | 13 数据集 manifest;70/30 study split |
| 3 | 1216_Dg 恶性率、tier 与 d35 一致(等价验证) |
| 4 | 1216_Dg bin 分布(HSC_MPP 72%/T_NK 0.8%)、stemness 方向 |
| 5 | **1599 programs / 151 MP / 9 跨数据集 MP**(MP040 LSC-like 等) |
| 6 | f05: T_NK baseline 5.15%(E-MTAB)→8.03%(GSE227903) |
 
> 任一批的金标准对不上 → 停,排查,不进下一批。这是"逐步夯实"的安全锁。
 