# 加载必要的包
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)

# ====================================================================
# 1. 定义文件路径与目标基因字典
# ====================================================================
base_dir <- "/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/01_per_sample_qc"

# 将基因与细胞类型映射，方便后续在表格中查看
gene_dict <- list(
  MSC        = c("LEPR", "CXCL12", "KITLG", "PDGFRA", "NT5E", "THY1", "NES"),
  EC         = c("PECAM1", "CDH5", "EMCN", "VWF"),
  Osteoblast = c("BGLAP", "SP7", "RUNX2"),
  Pericyte   = c("RGS5", "ACTA2", "MCAM")
)

# 展平字典生成基因参考数据框
gene_ref <- stack(gene_dict) %>% 
  rename(Gene = values, Cell_Type = ind) %>% 
  mutate(Gene = as.character(Gene))

target_genes <- gene_ref$Gene

# ====================================================================
# 2. 递归获取所有 rds 文件的完整路径
# ====================================================================
# recursive = TRUE 会自动穿透 Chen2023, GSE116256 等所有子目录
rds_files <- list.files(path = base_dir, pattern = "\\.rds$", 
                        full.names = TRUE, recursive = TRUE)

cat(sprintf("✅ 共检索到 %d 个 .rds 文件。开始批量处理...\n", length(rds_files)))

# ====================================================================
# 3. 定义单个 rds 文件的处理函数
# ====================================================================
# ====================================================================
# 3. (升级版) 定义单个 rds 文件的处理函数
# ====================================================================
process_single_rds <- function(file_path) {
  dataset_name <- basename(dirname(file_path))
  sample_name  <- str_remove(basename(file_path), "\\.rds$")
  
  res <- tryCatch({
    obj <- readRDS(file_path)
    
    # 兼容处理 Assay
    if("RNA" %in% Assays(obj)) {
      DefaultAssay(obj) <- "RNA"
    }
    
    # 获取表达矩阵 (如果 data 为空，尝试获取 counts)
    expr_mat <- GetAssayData(obj, layer = "data")
    if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
       expr_mat <- GetAssayData(obj, layer = "counts")
    }
    
    total_cells <- ncol(expr_mat)
    obj_genes <- rownames(expr_mat)
    
    # ---------------------------------------------------------
    # 核心修复：大小写不敏感匹配 (兼容人类 LEPR 和小鼠 Lepr)
    # ---------------------------------------------------------
    # match 寻找 target_genes 在 obj_genes 中的位置
    match_idx <- match(toupper(target_genes), toupper(obj_genes))
    
    # 剔除找不到的 NA 值
    valid_target_genes <- target_genes[!is.na(match_idx)]
    valid_obj_genes    <- obj_genes[na.omit(match_idx)]
    
    if (length(valid_obj_genes) == 0) {
      # 如果连忽略大小写都匹配不到，大概率行名是 Ensembl ID 或被破坏
      cat(sprintf("⚠️ [跳过] %-20s (Dataset: %s): 矩阵中未找到任何目标基因！请检查行名是否为 ENSG... ID。\n", 
                  sample_name, dataset_name))
      rm(obj, expr_mat)
      gc()
      return(NULL)
    }
    
    # 提取匹配成功的基因矩阵
    sub_mat <- expr_mat[valid_obj_genes, , drop = FALSE]
    pos_counts <- rowSums(sub_mat > 0)
    
    # 清理内存
    rm(obj, expr_mat, sub_mat)
    gc()
    
    # 将提取出来的小鼠/人类基因名，统一映射回你规定的全大写规范名，以便后续合并表格
    names(pos_counts) <- valid_target_genes
    
    data.frame(
      Dataset        = dataset_name,
      Sample         = sample_name,
      Gene           = names(pos_counts),
      Total_Cells    = total_cells,
      Positive_Cells = as.numeric(pos_counts)
    ) %>%
      mutate(Percentage = (Positive_Cells / Total_Cells) * 100)
      
  }, error = function(e) {
    # 如果是因为 Seurat 对象损坏或版本极度不兼容报错，打印出来
    cat(sprintf("❌ [报错] %-20s (Dataset: %s): %s\n", sample_name, dataset_name, e$message))
    return(NULL)
  })
  
  return(res)
}

# ====================================================================
# 4. 循环处理所有文件并合并结果
# ====================================================================
all_results_list <- list()

for (i in seq_along(rds_files)) {
  cat(sprintf("[%d/%d] 正在处理: %s\n", i, length(rds_files), basename(rds_files[i])))
  
  df_res <- process_single_rds(rds_files[i])
  if (!is.null(df_res)) {
    all_results_list[[i]] <- df_res
  }
}

# 合并所有样本的结果，并与细胞类型匹配
final_long_df <- bind_rows(all_results_list) %>%
  left_join(gene_ref, by = "Gene") %>%
  # 调整列顺序，看起来更舒服
  select(Dataset, Sample, Cell_Type, Gene, Total_Cells, Positive_Cells, Percentage)

cat("✅ 所有文件处理完成！\n")

# ====================================================================
# 5. 生成宽表 (方便人工在 Excel 中查看) 及保存
# ====================================================================
# 创建一个按 样本 x 基因 呈现阳性比例(Percentage)的宽表
final_wide_df <- final_long_df %>%
  select(Dataset, Sample, Total_Cells, Gene, Percentage) %>%
  pivot_wider(
    names_from = Gene, 
    values_from = Percentage,
    values_fill = NA # 若某样本未测到该基因填NA
  )

# 保存为 CSV 文件
write.csv(final_long_df, file = "/FAST/gr10634/gaozy/tmp/Niche_Markers_Summary_Long.csv", row.names = FALSE)
write.csv(final_wide_df, file = "/FAST/gr10634/gaozy/tmp/Niche_Markers_Summary_Wide.csv", row.names = FALSE)

cat("💾 结果已保存为 Niche_Markers_Summary_Long.csv 和 Niche_Markers_Summary_Wide.csv\n")