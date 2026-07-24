# aml_niche_net — 架构参考 (ARCHITECTURE)
 
> 定稿版(决策全部锁定)。DESIGN 文档里的开放问题已在此确定。
> 这是项目的**正式结构参考**:目录树 + 数据流 + config 层级 + 恶性标签。
 
---
 
## 1. 项目一句话
 
跨异质性 AML 单细胞队列,检验 LSC 样恶性状态是否**收敛**到保守的骨髓微环境通讯依赖轴,疾病信号是否编码在通讯网络**拓扑**而非单一 L–R,以及拓扑偏离是否沿**治疗压力轴单调加深**。计算发现为主线,方法(cNMF meta-program → per-sample CCC 图 → FGW 拓扑对齐)服务于发现。整合无关(integration-free):先建图,后对齐拓扑。
 
---
 
## 2. 目录树 (定稿)
 
```
scripts/
├── config/                    # 分层配置,单一路径真相
│   ├── config_paths.R / .sh   # 路径根 (双镜像)
│   ├── config_qc.R
│   ├── config_malignancy.R
│   ├── config_hierarchy.R
│   ├── config_cnmf.R
│   └── utils.R                # get_counts/fwrite_safe/core16/list_qc_samples
├── 00_ingest/                 # 原始 → per-dataset merged RDS
├── 01_preprocess/             # Phase 1A: per-sample QC + doublet + 角色/split
├── 02_malignancy/             # Phase 1B: 重比对 → 三证据共识恶性注释
├── 03_hierarchy/              # Phase 2: BMM 投影 + 8-bin + stemness
├── 04_cnmf/                   # Phase 3: per-sample cNMF + MP 聚合 + SRRS
├── 05_ccc/         (未来)     # Phase 4: LIANA+/CellChat per-sample CCC 图
├── 06_distance/    (未来)     # Phase 5: 强度 → 距离 (rank-based)
├── 07_fgw/         (未来)     # Phase 6: FGW 对齐 + barycenter
├── 08_scoring/     (未来)     # Phase 7: emergent edge + HDS/ATS/RLS/SCS
├── 09_robustness/             # Phase 8: 稳健性 + 平台对照 QC
├── 10_figures/
└── docs/
    ├── ARCHITECTURE.md        # 本文件
    ├── CODING_STANDARDS.md
    ├── PIPELINE.md            # 跑一个新数据集的完整顺序
    └── NEW_DATASET_RUNBOOK.md # c10 前置链 + 3 USER-CONFIRM
```
 
---
 
## 3. 数据流 (每 Phase 输入 → 输出)
 
```
[原始 GEO/Zenodo/ENA]
   │  00_ingest: 各数据集 reader → make_seurat (uid_patient, 受控 Timepoint, QC 指标)
   ▼
LARGE1/01_processed_counts/rds/<dataset>.rds   (merged, 未过滤, 每数据集一个)
   │  01_preprocess: 角色/split manifest;per-sample MAD QC + 双 doublet 共识
   ▼
LARGE1/02_seurat_objects/01_per_sample_qc/<ds>/<sample>.rds   (PASS 样本, per-sample)
   │  02_malignancy: FASTQ→STARsolo→cellsnp→Numbat;inferCNV;c50 三证据共识
   ▼
FAST/results/tables/02_malignancy/<ds>/<sample>__consensus_{percell,summary}.csv
   │  03_hierarchy: BMM 投影 → 8-bin;stemness;per-bin 恶性叠加 → rollup
   ▼
FAST/results/tables/03_hierarchy/{percell,perbin}/...  + rollup 总表
   │  04_cnmf: 恶性/fallback 细胞 → cNMF(K=4-9) → 程序 → 跨样本聚合 → MP + SRRS
   ▼
FAST/results/tables/04_cnmf/programs/...  + MP 词汇表 (151 MP, 9 跨数据集)
   │  [未来] 05_ccc → 06_distance → 07_fgw → 08_scoring
   ▼
患者级拓扑偏离评分 (HDS/ATS/RLS/SCS)
```
 
**贯穿**:09_robustness 对各阶段做 QC 对照(如 f05 平台偏离控制,每个新数据集入组时跑)。
 
---
 
## 4. config 层级
 
- `config_paths.{R,sh}`:唯一路径根(FAST/LARGE1/REF/ENV)+ SEED + 所有派生阶段目录。R/bash 双镜像手工同步。
- 阶段 config 各 source paths + 加本阶段常量:
  - `config_qc.R`:MAD/doublet/min_cells + gene sets
  - `config_malignancy.R`:inferCNV/Numbat 参数 + ref-norm + consensus tier 词表
  - `config_hierarchy.R`:投影/stemness/8-bin/mapping-QC 阈值
  - `config_cnmf.R`:cNMF/SRRS 参数 + 技术基因 pattern
- `utils.R`:共享 helper(get_counts v5/v4、fwrite_safe、core16 唯一实现、list_qc_samples）。
---
 
## 5. 恶性标签 (唯一)
 
- 唯一生产者:`02_malignancy/50_consensus_malignancy.R`(d35 已废弃并入)。
- 三证据类型:expression-CNV(inferCNV/copykat/scevan/author)、allele-CNV(Numbat)、SNV(VarTrix,[planned])。
- 投票:type-level majority(≥2/3);`--union_mode` = 任一 valid 类型阳性即恶性(盲区互补,当前两证据场景用)。
- tier:`A_concordant / B_multi_partial / C_single` → confidence `high/medium/low`。
- Numbat 双 schema:full(compartment_opt)可用;degraded(no_CNV_detected)= NA(无意见,不投票)。
- **语言纪律**:CNV/LOH 信号,非点突变;称 "proxy",不称 "mutation group"。
---
 
## 6. 关键设计决策 (blueprint 对应)
 
| 决策 | 内容 |
|---|---|
| integration-free (locked L1) | 先建图,后对齐拓扑;不做跨样本表达整合 |
| D1 niche 重定义 | "BM microenvironment" = 造血-免疫+旁分泌;结构基质仅 cohort 辅助 |
| D2 平台不变性 | rank-based 距离 + per-sample 推断 + 平台协变量 |
| M1 复发降级 | 主时间轴 = 治疗压力(Dx→MRD);复发为探索性验证终点 |
| M2 BAM 前置 | FASTQ→STARsolo 重比对 → Numbat/VarTrix |
| M4 conditional SRRS | ≥3 研究才算 conditional;否则证据加权 |
| 恶性标签统一 | c50 唯一;union_mode 收编 d35;degraded Numbat→NA |
 
---
 
## 7. 数据集角色 (13 活跃)
 
- Reference-scaffold:BoneMarrowMap/Zeng 2025(投影参考,非队列)
- Healthy-control:E-MTAB-11536、GSE253355(+ 各 AML 研究内健康供体)
- Discovery-AML:GSE239721、Petti2019、GSE289435、GSE185381(Lasry)、Chen2023
- Treatment-axis:GSE116256(van Galen)、GSE185991(L1-only)、GSE147989(L1-only)
- Validation-relapse:GSE227903(主 L2 验证)、GSE201966
- Exploratory:GSE207356
- 平台高度混杂(10x 3'/5'、Seq-Well、CITE)→ 是稳健性压力测试,不是缺陷。
---
 
## 8. 当前进度 (重构基线)
 
- Phase 1–3 pilot 已跑通(2 数据集:GSE227903 + GSE239721):
  - 恶性注释(inferCNV∪Numbat,当前 d35;重构后 → c50)
  - 8-bin 投影 + 4 签名 stemness + prob-based mapping QC
  - cNMF pilot:**1599 programs / 41 样本 → 151 MP,9 跨数据集**(含 MP040 LSC-like)
- QC 对照:f05 证明"投影偏离梯度"主要是平台效应(T_NK baseline 5.15%→8.03%)。
- 重构进行中(本文档定义的目标结构);之后再扩数据集。
 