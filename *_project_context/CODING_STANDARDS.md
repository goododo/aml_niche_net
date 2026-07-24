# aml_niche_net — 编码规范 (CODING_STANDARDS)
 
> 本文件定死项目代码风格。**任何新对话中,读此文件即可产出与现有 codebase 一致的代码。**
> 规范提炼自四阶段 review 后收敛的最佳实践;模板脚本 = 旧 `d00_config.R`(风格最佳)。
 
---
 
## 1. 目录与编号
 
- 目录号 = blueprint Phase:`00_ingest / 01_preprocess / 02_malignancy / 03_hierarchy / 04_cnmf / 05_ccc / 06_distance / 07_fgw / 08_scoring / 09_robustness / 10_figures`。
- **目录内脚本从 01 起编号**(双层:目录号=Phase,脚本号=该 Phase 内步骤)。**不用字母前缀。**
- 脚本命名:`NN_verb_noun.ext`,如 `03_assign_bins.R`、`10_run_one_sample.sh`。
- 测试脚本:`90_/91_` 前缀,与主流程同目录(便于对照)。
- 自动生成的脚本:进 `_generated/`,gitignore,**不进源码**。
---
 
## 2. 配置 (分层)
 
- **单一路径真相**:所有路径根只在 `config/config_paths.{R,sh}` 定义;双镜像手工同步(顶部注明"改一处必改另一处")。
- **路径定位用 `.here` 锚点(零硬编码)**:项目根放一个空文件 `.here`;`config_paths.R` 里 `library(here)` + `here::here()` 向上查找 `.here` 锁定项目根,派生 `CONFIG_DIR`。脚本从任何工作目录运行都正确。
- 每阶段一个 `config_<stage>.R`,内部先 source `config_paths.R`(经 here),再定义本阶段常量。
- 业务脚本开头统一(**用 here,不硬编码 CONFIG_DIR**):
  ```r
  suppressPackageStartupMessages({ library(here) })
  source(here::here("scripts", "config", "config_paths.R"))   # 定义 CONFIG_DIR/SEED/路径
  source(here::here("scripts", "config", "config_<stage>.R")) # 本阶段常量
  source(here::here("scripts", "config", "utils.R"))          # 共享 helper
  ```
- **禁止**:在业务脚本里裸写路径字符串;硬编码 CONFIG_DIR(用 here);跨阶段 source 别的阶段的 config(全部经 config/)。
---
 
## 3. 命名
 
- **路径常量**:`ALL_CAPS`,如 `FAST_DIR`、`DIR_MALIGNANCY`、`QC_RDS_DIR`。
- **可调参数常量**:`ALL_CAPS`,如 `SEED`、`MIN_PROJ_PROB`、`CNMF_K_RANGE`。
- **helper 函数**:`make_xxx()`(构造)、`read_xxx()`(读盘)、`load_xxx()`(读+处理)、动词开头。
- **不用** `cfg$list` 风格(旧 02_preprocess 的模式,已废弃);统一裸 ALL_CAPS 常量。
- 元数据列名约定:`uid_patient = dataset:Patient_ID`(命名空间键,防碰撞);`orig.ident = Patient_ID`;作者注释前缀 `anno_`。
---
 
## 4. 日志 (统一)
 
- **唯一风格**:`message("[N] ...")`,N 是步骤号,与脚本 section 对应。
- 需要时间戳时,用 utils 里的包装:`msg("[3] ...")` → `[HH:MM:SS] [3] ...`。
- **禁止**:`cat()` 做日志(cat 留给数据输出);`.log()`(旧 02 的自定义,已废弃);`print()` 做进度(留给最终结果表)。
- 日志走 stderr(message 默认),不污染 stdout 数据流。
---
 
## 5. 脚本结构
 
- **Header 块**(每个脚本开头):
  ```r
  # NN_verb_noun.R ----
  # <一句话作用>
  # INPUT  : <路径/格式>
  # OUTPUT : <路径/格式>
  # Usage  : Rscript NN_verb_noun.R --arg=val
  # <关键决策/坑的注释,引用 blueprint 决策号如 [D2]/[M4]>
  ```
- **Section 标题**:以 ` ----` 结尾(RStudio Outline 兼容),如 `## -- Step 1. load ----`。
- **决策注释**:关键取舍处标 blueprint 决策号(`# [DECISION Q2a] ...`),保证可追溯。
- CLI 参数:`=`-form(`--dataset=GSE227903`),用 optparse。
- **resume-skip**:每个产出文件的步骤先检查输出是否存在,存在则 skip(幂等、可续跑)。
- 长流程用 `sbatch --requeue`,不用 login 节点 tmux。
---
 
## 6. 可复现性
 
- **全局 SEED = 20260605L**(唯一;不允许其他 seed,除非注释说明特殊理由并经确认)。
- 每个含随机性的步骤 `set.seed(SEED)`;并行用 `future.seed = TRUE` 并注释说明。
---
 
## 7. 数据栈与约定
 
- R:`data.table`(I/O + 处理)、`Seurat`/`SeuratObject`、`ggplot2` + `patchwork`(图)、`igraph`(图/Leiden)。
- `fwrite_safe()`(utils)写表:自动建父目录(FAST 是 purgeable scratch)。
- `get_counts()`(utils):v5 layer / v4 slot 双兼容,不直接调 GetAssayData。
- `core16()`(utils):barcode 16bp core 提取,**唯一实现 + 长度校验**(不在各脚本重定义)。
- 存储分区:`LARGE1` 放大对象(BAM/CNV/大 RDS),`FAST` 放表/图/脚本(可清除)。
---
 
## 8. 基因 pattern (约定)
 
- 线粒体:`^MT[-.]`(+ MITO_SHORT 兜底剥前缀数据集)
- 核糖体:`^RP[SL]`
- 血红蛋白:`^HB[^P]`
- cNMF 技术基因移除(保守):mito/ribo/HB + `^MIR[0-9] ^SNOR ^SCARNA ^RNU[0-9] ^RN7S`
- **保留**:lncRNA(AC/AL/LINC/AS)、IEG(FOS/JUN)、cell-cycle —— 由 SRRS 复现性过滤,不手删(避免 cherry-pick)。
---
 
## 9. 恶性标签 (唯一来源)
 
- **只有 `02_malignancy/50_consensus_malignancy.R`** 生产恶性标签。d35 已废弃。
- 三证据类型:expression-CNV / allele-CNV / SNV;`--union_mode` 用于盲区互补。
- 语言纪律:inferCNV/Numbat 检测的是 **CNV/LOH 信号,非点突变**。结果一律称 "malignancy proxy / CNV-based",**不称 "TP53 mutation group"** 等。
- Numbat degraded schema(`no_CNV_detected`)= NA(无意见),**不是 malignant=0**。
- Numbat 迭代:最高号 `segs_consensus_{N}.tsv` 为准。
- tier 词表:`A_concordant / B_multi_partial / C_single`(首字母对齐 TIER2CONF A=high/B=medium/C=low)。
---
 
## 10. 语言
 
- **对话用中文,代码注释用英文。**
- co-designer 模式:主动指出 reviewer 攻击点、执行卡点;写代码前一次一个澄清问题。
---
 
## 附:一个合规脚本骨架
 
```r
# 03_assign_bins.R ----
# Assign each cell to one of 8 hierarchy bins from the BMM projection.
# INPUT  : PROJ_OBJ_DIR/<ds>/<sample>.rds  (d 01 projection output)
# OUTPUT : PERCELL_DIR/<ds>/<sample>_percell.csv
# Usage  : Rscript 03_assign_bins.R --dataset=GSE227903
# [DECISION] bins from broad(24)->8 map in bmm_bin_map.tsv (author-defined, not invented).
 
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))
 
opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character")
)))
 
message("[1] listing QC samples for ", opt$dataset)
samples <- list_qc_samples(opt$dataset)
 
for (i in seq_len(nrow(samples))) {
  out <- file.path(PERCELL_DIR, opt$dataset, paste0(samples$sample[i], "_percell.csv"))
  if (file.exists(out)) { message("    [skip] ", samples$sample[i]); next }   # resume-skip
  set.seed(SEED)
  # ... work ...
  fwrite_safe(res, out)
  message("    [write] ", out)
}
message("[done] ", opt$dataset)
```
 