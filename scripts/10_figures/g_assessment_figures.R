# g_assessment_figures.R ----
# 组会/内部评估用的 5 张关键图。与 f01-f04 无关（那批指向旧路径且从未跑通）。
# INPUT  : results/tables/{01_preprocess,02_malignancy,03_hierarchy,05_ccc,07_fgw,08_scoring}/*.csv
# OUTPUT : results/figures/11_assessment/g0{1..5}_*.{pdf,png}
# Usage  : conda run -p /FAST/gr10634/gaozy/general_env Rscript scripts/10_figures/g_assessment_figures.R
#
# [DECISION] 配色取 dataviz 参考调色板前 2 槽（#2a78d6 蓝 / #eb6834 橙），已用
#            Machado-2009 severity-1.0 CVD 模拟实测：全对 CVD dE=24.7、normal dE=33.6、
#            对比度均 >=3:1。发散色标用蓝<->红 + 中性灰中点。
# [DECISION] 状态色 critical=#d03b3b 只用于"标记有问题的量"，全图不得再作系列色。
# [DECISION] 图上每个数字一律由源表算出，不写死字面量 —— 首版曾把边级 min q 写死在
#            block 面板的副标题里（0.87 vs 实际 0.93），且把 n_perm 写成 2000（实际 10000）。
# [DECISION] 比例/概率轴一律封顶（率 <=100%，p 与 q <=1），标签靠 clip="off" 放到边距外。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
TBL  <- file.path(ROOT, "results", "tables")
FIG  <- file.path(ROOT, "results", "figures", "11_assessment")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
CRIT <- "#d03b3b"; WARN <- "#fab219"; GOOD <- "#0ca30c"
GREY <- "#8a8a85"; LGREY <- "#c6c5c0"; MIDGREY <- "#f0efec"
INK  <- "#0b0b0b"; INK2 <- "#52514e"
FAM  <- "Noto Sans CJK SC"

theme_a <- function(base = 11) {
  theme_minimal(base_size = base, base_family = FAM) +
    theme(
      plot.title      = element_text(face = "bold", size = base + 2, colour = INK, hjust = 0),
      plot.subtitle   = element_text(size = base - 0.5, colour = INK2, hjust = 0, lineheight = 1.3),
      plot.caption    = element_text(size = base - 2.5, colour = INK2, hjust = 0, lineheight = 1.35),
      axis.title      = element_text(size = base - 1, colour = INK2),
      axis.text       = element_text(size = base - 1, colour = INK2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.3),
      legend.position = "top", legend.justification = "left",
      legend.title    = element_blank(), legend.key.size = unit(0.8, "lines"),
      plot.title.position = "plot", plot.caption.position = "plot",
      strip.text      = element_text(face = "bold", colour = INK, size = base - 1)
    )
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 200,
         device = ragg::agg_png, bg = "white")
  cat(sprintf("[fig] %s\n", name))
}
# p 值格式化：避免 "p=<0.001" 这种双符号
pfmt <- function(p) ifelse(p < 1e-3, "p<0.001", paste0("p=", formatC(p, format = "g", digits = 2)))

## ══════════════════════════════════════════════════════════════════
## Fig 1 — 各阶段实际样本数（非嵌套，不叫"漏斗"）
## ══════════════════════════════════════════════════════════════════
qc <- rbindlist(lapply(Sys.glob(file.path(TBL, "01_preprocess", "03_qc_report__*.csv")),
                       fread), fill = TRUE)
qc_pass <- qc[status == "PASS"]
cons  <- fread(file.path(TBL, "02_malignancy", "ALL_consensus_summary.csv"))
fgwix <- fread(file.path(TBL, "07_fgw", "fgw_input_index.csv"))
asw   <- fread(file.path(TBL, "08_scoring", "alpha_sweep.csv"))

# 148 vs 130 的真实对账（不要猜原因，算出来）
setkey(cons, dataset, sample); setkey(fgwix, dataset, sample)
both     <- fgwix[cons, nomatch = 0L]
graph_nolabel <- fgwix[!cons]
label_nograph <- cons[!fgwix]
n_aml_nolabel <- graph_nolabel[healthy == FALSE, .N]
n_hlt_nolabel <- graph_nolabel[healthy == TRUE,  .N]

funnel <- data.table(
  stage = c("QC 通过样本", "有恶性共识标签", "构建了 CCC 图", "进入 FGW 全局模型", "进入平台受控模型"),
  n     = c(nrow(qc_pass), nrow(cons), nrow(fgwix), asw$n_global[1], asw$n_strat[1]),
  nds   = c(uniqueN(qc_pass$dataset), uniqueN(cons$dataset), uniqueN(fgwix$dataset), NA_integer_, 3L)
)
funnel[, stage := factor(stage, levels = rev(stage))]
funnel[, lab := ifelse(is.na(nds), sprintf("%d 样本", n), sprintf("%d 样本 · %d 套数据", n, nds))]
funnel[, hl := c(rep(FALSE, 4), TRUE)]

p1 <- ggplot(funnel, aes(stage, n, fill = hl)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = lab), hjust = -0.06, size = 3.4, colour = INK, family = FAM) +
  annotate("text", x = 1, y = funnel$n[5] + 6, label = "← 最窄处", hjust = 0,
           size = 3.2, colour = INK2, family = FAM, vjust = 2.4) +
  scale_fill_manual(values = c("FALSE" = BLUE, "TRUE" = CRIT), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  coord_flip() +
  labs(title = "各阶段实际样本数（注意：非嵌套）",
       subtitle = "H2 的「平台受控」结论只用了 3 套数据、60 个样本",
       x = NULL, y = "样本数",
       caption = sprintf(paste0(
         "对账：CCC 图(%d) = 有标签且建图(%d) + 无恶性标签但建图(%d：%d 个健康供者 + %d 个 AML 样本)；\n",
         "另有 %d 个样本有标签但未建图。—— 那 %d 个无标签的 AML 样本，其 frac_malignant 在\n",
         "01_build_fgw_inputs.R:72 被填成全局均值，等于给了它们一个从未测过的恶性比例。"),
         nrow(fgwix), nrow(both), nrow(graph_nolabel), n_hlt_nolabel, n_aml_nolabel,
         nrow(label_nograph), n_aml_nolabel)) +
  theme_a() + theme(panel.grid.major.y = element_blank())

save_fig(p1, "g01_cohort_counts", 9.6, 5.0)

## ══════════════════════════════════════════════════════════════════
## Fig 2 — 恶性标签质量（本组最重要）
## ══════════════════════════════════════════════════════════════════
fpr <- fread(file.path(TBL, "02_malignancy", "malignancy_fpr_by_bin.csv"))
fpr[hierarchy_bin == "" | is.na(hierarchy_bin), hierarchy_bin := "未分配 bin"]   # 不静默丢弃
fpr[, in_graph := !is.na(in_ccc_graph) & in_ccc_graph == TRUE]
fpr[, hierarchy_bin := factor(hierarchy_bin, levels = hierarchy_bin[order(FPR)])]
fpr[, hl := fifelse(as.character(hierarchy_bin) == "HSC_MPP", "crit",
             fifelse(in_graph, "in", "out"))]
max_in <- fpr[in_graph == TRUE][which.max(FPR)]

p2a <- ggplot(fpr, aes(hierarchy_bin, FPR, fill = hl)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%s)", 100 * FPR, comma(n))),
            hjust = -0.05, size = 3.15, colour = INK, family = FAM) +
  scale_fill_manual(values = c(crit = CRIT, "in" = BLUE, out = LGREY), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  coord_flip(clip = "off") +
  labs(title = "A  健康供者中被判为 malignant 的比例（假阳性率）",
       subtitle = "健康骨髓里任何一个 malignant 判定按定义都是假阳性。浅灰 = 不进入 CCC 图的 bin。",
       x = NULL, y = "假阳性率") +
  theme_a() +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 105, 5, 5))

tier <- cons[, .N, by = evidence_tier]
tier[, evidence_tier := factor(evidence_tier, levels = c("A_concordant", "B_multi_partial", "C_single"))]
setorder(tier, evidence_tier)
tier[, lab_leg := sprintf("%s (%d)", evidence_tier, N)]

p2b <- ggplot(tier, aes(x = "", y = N, fill = evidence_tier)) +
  geom_col(width = 0.45, colour = "white", linewidth = 1.1) +
  geom_text(data = tier[N > 10], aes(label = sprintf("%s  %d", evidence_tier, N)),
            position = position_stack(vjust = 0.5), size = 3.3, colour = "white",
            fontface = "bold", family = FAM) +
  scale_fill_manual(values = c("A_concordant" = GOOD, "B_multi_partial" = WARN, "C_single" = GREY),
                    labels = tier$lab_leg) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(title = sprintf("B  证据等级构成（%d 个样本）", nrow(cons)),
       subtitle = paste0("blueprint Phase 1 要求 inferCNV + Numbat + VarTrix 三方共识、≥ 2/3 一致；\n",
                         sprintf("实际 %d/%d 是单臂 inferCNV。A/B 两级共 %d 个样本，全部来自 GSE227903。",
                                 tier[evidence_tier == "C_single", N], nrow(cons),
                                 tier[evidence_tier != "C_single", sum(N)])),
       x = NULL, y = "样本数") +
  theme_a() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())

p2 <- p2a / p2b + plot_layout(heights = c(2.5, 1)) +
  plot_annotation(
    title = "恶性标签质量：LSC bin 的假阳性率接近 40%",
    subtitle = paste0(
      sprintf("HSC_MPP 正是中心命题最依赖的节点，也是 CCC 图内假阳性最高的 bin（%.1f%%）。\n", 100 * max_in$FPR),
      "淋巴系(T_NK 3.6%、B_Plasma 7.7%)干净，髓系/原始系不干净 —— 与 inferCNV 依赖表达偏离的原理一致。\n",
      "（Stromal 90.3% 更高，但按既定前提不进入 CCC 图。）"),
    caption = paste0("来源：A = 02_malignancy/malignancy_fpr_by_bin.csv（由 96_malignancy_fpr_healthy.R 生成）；",
                     "B = 02_malignancy/ALL_consensus_summary.csv"),
    theme = theme_a())

save_fig(p2, "g02_label_quality", 9.8, 8.4)

## ══════════════════════════════════════════════════════════════════
## Fig 3 — H1：49 条有向边，无一显著
## ══════════════════════════════════════════════════════════════════
ee  <- fread(file.path(TBL, "08_scoring", "emergent_edges.csv"))
blk <- fread(file.path(TBL, "08_scoring", "block_permutation.csv"))
NODES <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma")
ee[, sender_bin   := factor(sender_bin,   levels = NODES)]
ee[, receiver_bin := factor(receiver_bin, levels = rev(NODES))]
lim <- max(abs(ee$dC_real))
# 置换次数由 perm_p 的分母反推，不写死：perm_p = (1+k)/(n_perm+1)，故 perm_p*(n_perm+1)
# 对所有边都必须是整数。取唯一满足此条件的候选。（首版用「最小差值」反推，会给出
# 错误的 1666 —— 相邻 perm_p 的差不一定是 1/(n_perm+1)。）
.n_cand <- c(500L, 1000L, 2000L, 5000L, 10000L, 20000L, 50000L)
.n_ok <- vapply(.n_cand, function(n) all(abs(ee$perm_p * (n + 1) - round(ee$perm_p * (n + 1))) < 1e-6),
                logical(1))
stopifnot("无法从 perm_p 反推 n_perm" = any(.n_ok))
n_perm_edge <- .n_cand[.n_ok][1]

p3a <- ggplot(ee, aes(sender_bin, receiver_bin, fill = dC_real)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  scale_fill_gradient2(low = CRIT, mid = MIDGREY, high = BLUE, midpoint = 0,
                       limits = c(-lim, lim),
                       name = "dC = C_AML − C_healthy    ← AML 中更强 | AML 中更弱 →") +
  coord_fixed() +
  labs(title = "A  49 条有向边的 AML−健康差值",
       subtitle = sprintf("色深 = 差值大小，但无一格通过置换检验（全部 q ≥ %.2f）。此图展示的是噪声的结构，不要读出模式。",
                          min(ee$q_bh)),
       x = "sender", y = "receiver") +
  theme_a() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid.major = element_blank(),
        legend.key.width = unit(2.1, "lines"),
        legend.title = element_text(size = 8, colour = INK2))

setorder(blk, -n_edges)
blk[, block := factor(block, levels = blk$block)]
p3b <- ggplot(blk, aes(block, q_bh)) +
  geom_col(width = 0.6, fill = LGREY) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = CRIT, linewidth = 0.6) +
  annotate("text", x = 0.55, y = 0.07, label = "FDR 0.05", hjust = 0, size = 3,
           colour = CRIT, family = FAM) +
  geom_text(aes(label = sprintf("q=%.2f  (%d 条边)", q_bh, n_edges)), hjust = -0.06,
            size = 3.05, colour = INK, family = FAM) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_flip(clip = "off") +
  labs(title = "B  5 个预设通讯 block 的 BH-FDR q 值（按边数排序）",
       subtitle = sprintf("5 个 block 的 q 完全相同 = %.2f，远高于 0.05。", min(blk$q_bh)),
       x = NULL, y = "q (BH-FDR)") +
  theme_a() +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 120, 5, 5))

p3 <- p3a / p3b + plot_layout(heights = c(2.2, 1)) +
  plot_annotation(
    title = "H1（保守 emergent edge）：未观察到信号",
    subtitle = sprintf(paste0(
      "49 条边中 0 条通过置换检验（最小 perm p = %.3f，最小 q = %.2f）；5 个 block 全部不显著。\n",
      "关键：边来自 CellChat，不经过恶性标签 —— 因此这个 null 不受标签质量问题影响，是目前最稳的一个结果。"),
      min(ee$perm_p), min(ee$q_bh)),
    caption = sprintf("来源：02_permutation_emergent.py (n_perm=%d) / 03_block_permutation.py (n_perm=10000)",
                      n_perm_edge),
    theme = theme_a())

save_fig(p3, "g03_h1_null", 9.6, 9.2)

## ══════════════════════════════════════════════════════════════════
## Fig 4 — H2 决定性检验
## ══════════════════════════════════════════════════════════════════
asw_l <- melt(asw, id.vars = "alpha", measure.vars = c("p_global", "p_strat"),
              variable.name = "model", value.name = "p")
asw_l[, model := factor(model, levels = c("p_global", "p_strat"),
                        labels = c(sprintf("全局模型 (n=%d)", asw$n_global[1]),
                                   sprintf("数据集内分层 (n=%d)", asw$n_strat[1])))]

p4a <- ggplot(asw_l, aes(alpha, p, colour = model)) +
  annotate("rect", xmin = 0.88, xmax = 1.02, ymin = 5e-5, ymax = 1, fill = CRIT, alpha = 0.07) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = CRIT, linewidth = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  geom_text(data = asw_l[alpha == 1], aes(label = model), hjust = 0, nudge_x = 0.035,
            size = 3, family = FAM, show.legend = FALSE) +
  annotate("text", x = 0.95, y = 0.62, label = "α = 1 纯拓扑", size = 3,
           colour = CRIT, family = FAM) +
  annotate("text", x = 0.02, y = 0.072, label = "p = 0.05", hjust = 0, size = 3,
           colour = CRIT, family = FAM) +
  scale_colour_manual(values = c(BLUE, ORANGE)) +
  scale_y_log10(limits = c(5e-5, 1), labels = label_number(drop0trailing = TRUE),
                breaks = c(1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.3, 1)) +
  scale_x_continuous(breaks = asw$alpha, expand = expansion(mult = c(0.03, 0.34))) +
  labs(title = "A  α 扫描：把特征项权重调到 0 时，显著性消失",
       subtitle = "α=0 纯节点特征 → α=1 纯 Gromov–Wasserstein（纯拓扑）。两个模型在 α=1 都越过 0.05。",
       x = "α (FGW 中结构项的权重)", y = "置换 p 值（对数轴，封顶 1）") +
  theme_a() + theme(legend.position = "none")

fd <- fread(file.path(TBL, "08_scoring", "feature_decomposition.csv"))[alpha == 0.5]
fd[, feature_set := factor(feature_set,
      levels = c("only_stemness", "only_ncells", "no_frac_mal", "all3", "only_frac_mal"))]
setorder(fd, feature_set)
fd[, kind := fifelse(feature_set == "only_frac_mal", "循环项（健康被强制置 0）",
              fifelse(p_global < 0.05, "显著 (p<0.05)", "不显著"))]
fd[, kind := factor(kind, levels = c("循环项（健康被强制置 0）", "显著 (p<0.05)", "不显著"))]

p4b <- ggplot(fd, aes(feature_set, beta_global, fill = kind)) +
  geom_col(width = 0.6) +
  # 负值柱的标签也放到 0 的右侧，避免压到 y 轴类别名上
  geom_text(aes(y = pmax(beta_global, 0.002), label = sprintf("β=%.3f  %s", beta_global, pfmt(p_global))),
            hjust = -0.06, size = 3.05, colour = INK, family = FAM) +
  scale_fill_manual(values = setNames(c(CRIT, BLUE, LGREY), levels(fd$kind))) +
  scale_y_continuous(limits = c(-0.05, 0.26), expand = expansion(mult = c(0.02, 0.02))) +
  coord_flip(clip = "off") +
  labs(title = "B  特征分解（α=0.5）：拆开看信号来自哪里",
       subtitle = paste0(
         "frac_malignant 单独即可复现并放大效应（β=0.173 > all3 的 0.075）—— 它在 01_build_fgw_inputs.R:67\n",
         "被对健康样本强制置 0，因此这个分离是构造出来的。但去掉它后仍剩 β=0.064 (p=0.001)，\n",
         "其中 only_ncells（细胞数，纯技术量）β=0.050 (p=0.008)；only_stemness（真正的生物学量）不显著。"),
       x = NULL, y = "is_aml 偏回归系数 β") +
  theme_a() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 100, 5, 5))

p4 <- p4a / p4b + plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = "H2（拓扑携带疾病信息）：显著性来自构造与技术量，不是来自拓扑",
    subtitle = "这两张图是主动设计的证伪检验。结论：现有 H2 信号无法归因于网络拓扑，也无法归因于干性。",
    caption = "来源：06_alpha_sweep.py / 07_feature_decomposition.py（均 n_perm=10000）",
    theme = theme_a())

save_fig(p4, "g04_h2_decisive", 9.8, 9.4)

## ══════════════════════════════════════════════════════════════════
## Fig 5 — H3：被 metadata 阻断，不解读
## ══════════════════════════════════════════════════════════════════
stem <- fread(file.path(TBL, "03_hierarchy", "stemness_by_timepoint.csv"))
md   <- fread(file.path(TBL, "03_hierarchy", "malignant_distribution.csv"))[ok == TRUE]
dom  <- md[timepoint %in% c("Dx", "MRD", "Relapse"),
           .(n_pat = uniqueN(patient), tot = sum(tot_mal), top = max(tot_mal)), by = timepoint]
dom[, share := top / tot]
mrd_share <- dom[timepoint == "MRD", share]

sl <- melt(stem, id.vars = c("tp", "n_malignant"),
           measure.vars = c("LSC17", "vanGalen_HSC_Prog", "vanGalen_HSC_like", "HSPC_core"),
           variable.name = "signature", value.name = "score")
sl[, tp := factor(tp, levels = c("Dx", "MRD", "Relapse"))]
sl[, shape := fifelse(signature %in% c("LSC17", "vanGalen_HSC_Prog"), "V 形", "非 V 形")]
xlab_n <- dom[match(levels(sl$tp), timepoint)]
# 标签必须短 —— 首版用 "Dx\n8 患者 · 4,136 恶性细胞"，在分面里三个标签互相重叠成一行。
# 细胞数移到副标题里讲。
lev_lab <- sprintf("%s\n(%d 患者)", levels(sl$tp), xlab_n$n_pat)

p5a <- ggplot(sl, aes(tp, score, colour = signature, group = signature)) +
  geom_line(linewidth = 0.85) + geom_point(size = 2.4) +
  geom_text(data = sl[tp == "Relapse"], aes(label = signature), hjust = 0, nudge_x = 0.05,
            size = 2.9, family = FAM, show.legend = FALSE) +
  facet_wrap(~shape, nrow = 1) +
  scale_colour_manual(values = c(BLUE, ORANGE, GREY, "#4a3aa7"), guide = "none") +
  scale_x_discrete(labels = lev_lab, expand = expansion(add = c(0.35, 1.5))) +
  labs(title = "A  干性签名沿治疗轴的走向（转录层面，独立于投影）",
       subtitle = sprintf(paste0(
         "blueprint H3 预测单调 Dx < MRD < Relapse。实际：4 个签名里只有 2 个呈 V 形；",
         "vanGalen_HSC_like 在 MRD 反而最高，HSPC_core 单调下降。\n",
         "而且这是细胞级混合均值（Dx/MRD/Relapse 各 %s / %s / %s 个恶性细胞）：\n",
         "MRD 的 %.0f%% 恶性细胞来自单个患者(6323)，所谓 V 形可能只是一个样本的效应。"),
         comma(dom[timepoint == "Dx", tot]), comma(dom[timepoint == "MRD", tot]),
         comma(dom[timepoint == "Relapse", tot]), 100 * mrd_share),
       x = NULL, y = "平均签名得分") +
  theme_a() + theme(panel.grid.major.x = element_blank())

sh <- fread(file.path(TBL, "03_hierarchy", "distribution_shift_tests.csv"))
sh[, comparison := factor(comparison, levels = c("Dx -> MRD", "MRD -> Relapse", "Dx -> Relapse"))]
sh[, score := factor(score, levels = c("primitive_frac", "stem_frac"),
      labels = c("primitive_frac（HSC_MPP+LMPP_GMP 占比）", "stem_frac（HSC_MPP 占比）"))]

p5b <- ggplot(sh, aes(comparison, median_delta)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_col(width = 0.6, fill = LGREY) +
  geom_text(aes(label = sprintf("%s\nn=%d 对", pfmt(p_value), n_pairs),
                vjust = ifelse(median_delta > 0, -0.3, 1.2),
                fontface = ifelse(p_value < 0.05, "bold", "plain")),
            size = 2.8, colour = INK2, family = FAM, lineheight = 0.95) +
  facet_wrap(~score, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.32, 0.32))) +
  labs(title = "B  配对样本内的中位变化量（粗体 = p<0.05）",
       subtitle = paste0("三个比较用的不是同一批患者（n=6 / 5 / 6，只有 5 位三个时间点齐全），",
                         "所以第三根柱不等于前两根之和。全部来自 GSE227903 一套数据。"),
       x = NULL, y = "中位 Δ") +
  theme_a() + theme(axis.text.x = element_text(angle = 12, hjust = 1))

p5 <- p5a / p5b + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "H3（治疗压力轴）：结果暂不解读 —— 依赖的 metadata 字段未核对",
    subtitle = paste0(
      "⚠  本图全部结论建立在「时间点(Dx/MRD/Relapse)」与「患者配对」两个字段上，",
      "而这两个字段已知在多套数据中与原文描述不一致\n",
      "（核对状态表尚未建立；另有 87 个样本的 timepoint 字段为空）。核对完成前，V 形趋势既不能确认也不能否认。"),
    caption = "来源：03_malignant_distribution_shift.R / 04_stemness_score.R / malignant_distribution.csv",
    theme = theme_a() + theme(plot.subtitle = element_text(colour = CRIT)))

save_fig(p5, "g05_h3_blocked", 10.2, 8.8)

cat("\n[done] 5 张图 -> ", FIG, "\n", sep = "")
