
<!-- TOC START -->
## 📑 目录

- [1] 叶克强教授AD研究及创新药科普 — 2026-05-16
- [2] 斑块相关微胶质细胞GRN空间分析流程 — 2026-05-17

<!-- TOC END -->



---
## [1] 叶克强教授AD研究及创新药科普
> 【AD_Pathology_Notes】来源: 2026-05-16-xhs-draft-yekeqiang-纯文本版.md | 合并日期: 2026-05-16
---
🌟 阿兹海默症的"破局者"！深理工叶克强教授全科普 🧠
今天必须安利一位神经科学领域的"真·大佬"——叶克强教授！他不仅发文章发到手软，更是真正把科研变成救命药的人。看完你就知道什么叫做科学家该有的样子 👇
✦ ✦ ✦
🌍 学术履历有多硬核？
20 年埃默里大学（Emory）执教经验，却毅然全职回国。现在是深圳理工大学生命健康学院讲席教授 + 中科院深圳先进院研究员。几十年只做一件事：攻破阿尔茨海默病、帕金森病等神经退行性疾病。学术天花板，就是这么顶 🏆
✦ ✦ ✦
🔬 两大颠覆性发现，每一个都够发 Cell
1️⃣ 锁定致病"元凶" —— Tau N368 片段
Tau 蛋白是 AD 研究的关键，叶教授团队发现天冬酰胺内肽酶（AEP）剪切产生的 Tau N368 片段，既是早期诊断的金标志物，更是药物靶点的"靶心" 🎯
2️⃣ 肠脑轴的奇妙证据 🦠➡️🧠
团队在《Cell》子刊证实：Tau 和 α-突触核蛋白的病理竟能从肠道"顺藤摸瓜"传到大脑！从此神经退行性疾病多了一个全新视角
✦ ✦ ✦
📚 代表性文献速览
▎Nature Medicine（2014）
⭐ 发现 δ-secretase 剪切 Tau N368 → 神经缠结
▎Nature（2022）
阻断 FSH 改善 AD 小鼠认知，解释女性高发机制
▎Nat Struct Mol Biol（2022）
DOPEGAL 修饰 Tau 促聚集，揭示蓝斑核易损之谜
▎Cell（2023）
全球首个 α-突触核蛋白 PET 示踪剂
▎Nat Commun ×4（2015–23）
δ-secretase 抑制剂、C/EBPβ 通路、FSH+ApoE4
📌 核心学术主线：C/EBPβ → δ-secretase → Tau N368，一条从机制到药物的完整链条
✦ ✦ ✦
💊 不只是发论文，他做出了真药
教授 + 创业者双重身份，2015 年创办 博芮健制药，落户深圳脑创中心。核心管线 BrAD-R13 片 是全球首款进入临床的口服 TrkB 受体激动剂，精准激活脑内保护通路，从源头减少致病蛋白。I 期临床已顺利收官，给全球千万 AD 患者带来真实可触的希望 🌅
✦ ✦ ✦
👨‍🔬 为什么他值得被更多人知道？
做基础研究的人很多，但敢跨越"死亡之谷"、把论文变成药片的人太少。从上海到深圳，叶教授一路追寻最好的转化生态，用硬核成果证明了中国科学家的担当 💪
✦ ✦ ✦
让我们一起期待创新药早日上市，造福千万患者 ❤️
✦ ✦ ✦
#阿尔茨海默症 #神经科学 #科研日常 #药学生物 #深圳理工大学 #中科院深圳先进院 #创新药 #Tau蛋白 #脑科学 #科学家故事


---
## [2] 斑块相关微胶质细胞GRN空间分析流程
> 【AD_Pathology_Notes】来源: 2026-05-17-plaque_associated_GRN_clean.ipynb | 合并日期: 2026-05-17
---
[markdown]
# Plaque-Associated Microglia: Gene Expression & GRN Spatial Analysis
---
**Analysis Pipeline**
1. Environment & Data Loading
2. Part A: Epigenetic gene expression vs plaque distance
3. Part B: GRN (regulon) activity vs plaque distance — linear regression + visualization

[markdown]
## 1. Environment & Data Loading

[code]
# ==========================================
# 1.1 Libraries (all at once)
# ==========================================
library(Signac)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(mgcv)
library(stats)
library(cowplot)
library(patchwork)
library(ggpubr)

print("Libraries loaded.")

[code]
# ==========================================
# 1.2 Read spatial metadata & clean barcode
# ==========================================
spatial <- read.csv('/home1/yzhang/Project_SN_Multiome/Analysis/Spatial_integration/sn_plaqueD_metadata.csv')
spatial$Clean_Barcode <- gsub("@.*$", "", spatial$X)
spatial <- spatial %>%
  select(X, Clean_Barcode, sample_id, min_center_dist) %>%
  mutate(dist_um = min_center_dist / 2)

print(paste("Spatial:", nrow(spatial), "spots"))
head(spatial, 3)

[code]
# ==========================================
# 1.3 Read multiome Seurat object & extract meta
# ==========================================
combined_full <- readRDS('/public/home/yzhang/zy_data/Project_AD_st/Multiome/Raw_data_0423_full/QC_0512/Combine/Combine_RNA_ATAC_anno_rnafromatac_260219.rds')
DefaultAssay(combined_full) <- "RNA"

meta <- combined_full@meta.data %>%
  select(cell_id, genotype, month, brain, class_id_label, Final_subclass_name)

print(paste("Combined:", ncol(combined_full), "cells"))
combined_full

[code]
# ==========================================
# 1.4 Read GRN data: AUC matrix, DGRN, RSS
# ==========================================
# AUC matrix — keep only +/+ regulons
AUC_gene <- read.csv('/home1/yzhang/Project_SN_Multiome/Analysis/Scenic+/Version_2_mouse/gene_AUC_matrix.csv', check.names = FALSE)
cell_col_name <- colnames(AUC_gene)[1]
plus_cols <- colnames(AUC_gene)[grepl("\\+/\\+", colnames(AUC_gene))]
AUC_gene <- AUC_gene[, c(cell_col_name, plus_cols)]

# Fix cell barcodes: "barcode-sample" → "sample_barcode"
AUC_gene[[cell_col_name]] <- gsub("^([^-]+)-(.*)$", "\\2_\\1", AUC_gene[[cell_col_name]])

# Differential GRN & RSS
DGRN <- read.csv('/home1/yzhang/Project_SN_Multiome/Analysis/Scenic+/Version_2_mouse/Differential_GRN_0515.csv')
RSS <- read.csv('/home1/yzhang/Project_SN_Multiome/Analysis/Scenic+/Version_2_mouse/RSS/RSS_class_id_label.csv', row.names = 1)
RSS <- as.data.frame(t(RSS))

print(paste("AUC matrix:", ncol(AUC_gene), "columns | DGRN:", nrow(DGRN), "rows | RSS:", ncol(RSS), "cell types"))

[code]
# ==========================================
# 1.5 Merge spatial + AUC + meta
# ==========================================
spatial <- spatial %>%
  inner_join(AUC_gene, by = c("Clean_Barcode" = "Cell")) %>%
  inner_join(meta, by = c("Clean_Barcode" = "cell_id"))

print(paste("Merged spatial:", nrow(spatial), "spots ×", ncol(spatial), "columns"))

[markdown]
## 2. Part A: Gene Expression vs Plaque Distance
Epigenetic regulators (Tet2, Asxl1, Kmt2d, Atrx, Cbl) in microglia

[code]
# ==========================================
# 2.1 Subset microglia & extract expression
# ==========================================
mg <- subset(combined_full, Final_subclass_name == '334 Microglia NN')

target_genes <- c("Tet2", "Asxl1", "Kmt2d", "Atrx", "Cbl")
expr_data <- FetchData(mg, vars = target_genes)
expr_data$cell_id <- rownames(expr_data)

print(paste("Microglia:", ncol(mg), "cells"))

[code]
# ==========================================
# 2.2 Binned mean expression vs plaque distance
#     (replaces original cells 15-22 with cleaner floor() approach)
# ==========================================
bin_size <- 25
max_dist <- 200

df_binned <- spatial %>%
  mutate(dist_um = min_center_dist / 2) %>%
  filter(dist_um <= max_dist) %>%
  inner_join(expr_data, by = c("Clean_Barcode" = "cell_id")) %>%
  mutate(bin_center = floor(dist_um / bin_size) * bin_size + bin_size / 2) %>%
  pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Expression") %>%
  group_by(bin_center, brain, Gene) %>%
  summarise(
    Mean_Expr = mean(Expression, na.rm = TRUE),
    SE_Expr   = sd(Expression, na.rm = TRUE) / sqrt(n()),
    Count     = n(),
    .groups   = "drop"
  ) %>%
  filter(Count >= 3)

head(df_binned)

[code]
# ==========================================
# 2.3 Plot: gene expression gradient
# ==========================================
ggplot(df_binned, aes(x = bin_center, y = Mean_Expr, color = Gene)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean_Expr - SE_Expr, ymax = Mean_Expr + SE_Expr),
                width = 5, alpha = 0.5) +
  facet_wrap(~ brain, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title    = "Epigenetic regulator expression vs plaque distance",
    subtitle = paste0("Bin size: ", bin_size, " µm | microglia only"),
    x        = "Distance to plaque center (µm)",
    y        = "Mean normalized expression"
  ) +
  theme(
    strip.background = element_rect(fill = "#f2f2f2", color = "black"),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "bottom"
  )

[markdown]
## 3. Part B: GRN Activity vs Plaque Distance
Linear regression + heatmap + single-TF line plot

[code]
# ==========================================
# 3.1 Filter & build TF mapping
# ==========================================
spatial_sub <- spatial %>% filter(min_center_dist < 200)

grn_cols <- grep("\\+/\\+", colnames(spatial_sub), value = TRUE)

# Pure TF name → complex GRN name dictionary
pure_tfs  <- gsub("_direct_\\+/\\+.*$", "", grn_cols)
tf2grn    <- setNames(grn_cols, pure_tfs)

cell_types <- unique(na.omit(spatial_sub$class_id_label))
months     <- unique(na.omit(spatial_sub$month))

print(paste("Cell types:", length(cell_types), "| Months:", paste(months, collapse = ", "),
            "| GRN columns:", length(grn_cols)))

[code]
# ==========================================
# 3.2 Precompute valid GRNs per cell type
#     (DGRN significant + RSS Top20 union)
# ==========================================
sig_diff_grn <- DGRN %>% filter(p_val_adj < 0.05 & abs(log2FC) > 0.1)

valid_grns_by_ctype <- lapply(cell_types, function(ctype) {
  # a. Differentially significant GRNs
  diff_grns <- sig_diff_grn %>%
    filter(CellType == ctype) %>%
    pull(Regulon) %>% unique()
  if (length(diff_grns) > 0 && !any(grepl("\\+/\\+", diff_grns))) {
    diff_grns <- na.omit(unname(tf2grn[diff_grns]))
  }
  
  # b. RSS Top 20
  if (ctype %in% colnames(RSS)) {
    top20 <- na.omit(unname(tf2grn[
      rownames(RSS)[order(RSS[[ctype]], decreasing = TRUE)[1:20]]
    ]))
  } else {
    top20 <- character(0)
  }
  
  intersect(unique(c(diff_grns, top20)), colnames(spatial_sub))
})
names(valid_grns_by_ctype) <- cell_types

print("GRN counts per cell type:")
sapply(valid_grns_by_ctype, length)

[code]
# ==========================================
# 3.3 Linear regression: GRN ~ distance
#     Triple loop: cell_type → month → GRN
# ==========================================
results_list <- list()

for (ctype in cell_types) {
  grns_to_test <- valid_grns_by_ctype[[ctype]]
  if (length(grns_to_test) == 0) next
  
  for (m in months) {
    sub_data <- spatial_sub %>% filter(class_id_label == ctype, month == m)
    if (nrow(sub_data) < 30) next
    
    for (grn in grns_to_test) {
      tmp_df <- na.omit(data.frame(
        activity = sub_data[[grn]],
        dist     = sub_data$min_center_dist
      ))
      if (nrow(tmp_df) < 30) next
      
      fit <- lm(activity ~ dist, data = tmp_df)
      s <- summary(fit)$coefficients
      
      if ("dist" %in% rownames(s)) {
        results_list[[paste(ctype, m, grn, sep = "_")]] <- data.frame(
          CellType     = ctype,
          Month        = m,
          GRN          = grn,
          Beta_Distance = s["dist", "Estimate"],
          Pval         = s["dist", "Pr(>|t|)"],
          N_Cells      = nrow(tmp_df),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Aggregate & FDR correction
grn_stats <- bind_rows(results_list)
if (nrow(grn_stats) > 0) {
  grn_stats$FDR <- p.adjust(grn_stats$Pval, method = "BH")
  grn_stats <- grn_stats %>% arrange(FDR)
}

print(paste("Total tests:", nrow(grn_stats), "| Sig (FDR<0.05):", sum(grn_stats$FDR < 0.05)))
print(head(grn_stats, 10))

[code]
# ==========================================
# 3.4 Heatmap: GRN activity Z-score by distance bin
#     Glial cells → per brain region; Neurons → merged
# ==========================================
dist_breaks <- c(seq(0, 200, by = 40), Inf)
bin_labels  <- c(paste0(seq(0, 160, by = 40), "-", seq(40, 200, by = 40)), ">=200")

for (ctype in names(valid_grns_by_ctype)) {
  ctype_stats <- grn_stats %>% filter(CellType == ctype & Pval < 0.05)
  sig_grns <- unique(ctype_stats$GRN)
  if (length(sig_grns) == 0) next
  
  is_glial <- !grepl("Glut|Gaba", ctype, ignore.case = TRUE)
  regions  <- if (is_glial) unique(na.omit(spatial$brain[spatial$class_id_label == ctype])) else "All"
  
  for (reg in regions) {
    pd <- spatial %>% filter(class_id_label == ctype)
    if (reg != "All") pd <- pd %>% filter(brain == reg)
    if (nrow(pd) == 0) next
    
    pd <- pd %>%
      select(month, min_center_dist, all_of(sig_grns)) %>%
      mutate(Dist_Bin = cut(min_center_dist, breaks = dist_breaks,
                            labels = bin_labels, include.lowest = TRUE)) %>%
      drop_na(Dist_Bin) %>%
      pivot_longer(all_of(sig_grns), names_to = "GRN", values_to = "Activity") %>%
      group_by(month, Dist_Bin, GRN) %>%
      summarise(Mean_Activity = mean(Activity, na.rm = TRUE), .groups = "drop") %>%
      group_by(GRN, month) %>%
      mutate(Z_Score = as.numeric(scale(Mean_Activity))) %>%
      ungroup() %>%
      left_join(ctype_stats %>% select(Month, GRN, Pval, Beta_Distance),
                by = c("month" = "Month", "GRN" = "GRN")) %>%
      mutate(
        Sig_Star = case_when(
          is.na(Pval) ~ "", Pval < 0.001 ~ "***",
          Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""
        ),
        Label = ifelse(Sig_Star != "", paste0(Sig_Star, "\n", signif(Beta_Distance, 2)), "")
      )
    
    first_bin <- levels(pd$Dist_Bin)[1]
    ann <- pd %>% filter(Dist_Bin == first_bin)
    
    title_str <- if (reg == "All") paste("GRN Activity in", ctype)
                else paste("GRN Activity in", ctype, "(", reg, ")")
    
    p <- ggplot(pd, aes(x = Dist_Bin, y = GRN, fill = Z_Score)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_text(data = ann, aes(label = Label),
                color = "black", size = 2.5, fontface = "bold", hjust = 0, nudge_x = -0.4) +
      scale_fill_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0) +
      scale_y_discrete(labels = function(x) gsub("_direct_\\+/\\+.*", "", x)) +
      facet_wrap(~ month, ncol = length(unique(pd$month))) +
      labs(title = title_str, x = "Distance (µm)", y = "GRN", fill = "Z-Score") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid   = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA)
      )
    print(p)
  }
}

[code]
# ==========================================
# 3.5 Single-TF line plot example: Spi1 in Immune cells
#     Facet labels show β ± significance per month
# ==========================================
target_ctype <- "34 Immune"
target_tf    <- "Ebf1"

ctype_stats <- grn_stats %>% filter(CellType == target_ctype)
sig_grns <- ctype_stats %>%
  filter(grepl(paste0("^", target_tf, "_"), GRN)) %>%
  pull(GRN) %>% unique()

if (length(sig_grns) > 0) {
  plot_data <- spatial %>%
    filter(class_id_label == target_ctype) %>%
    select(month, min_center_dist, all_of(sig_grns)) %>%
    mutate(Dist_Bin = cut(min_center_dist, breaks = dist_breaks,
                          labels = bin_labels, include.lowest = TRUE)) %>%
    drop_na(Dist_Bin) %>%
    pivot_longer(all_of(sig_grns), names_to = "GRN", values_to = "Activity") %>%
    group_by(month, Dist_Bin, GRN) %>%
    summarise(Mean_Activity = mean(Activity, na.rm = TRUE), .groups = "drop") %>%
    group_by(GRN, month) %>%
    mutate(Z_Score = as.numeric(scale(Mean_Activity))) %>%
    ungroup() %>%
    left_join(ctype_stats %>% select(Month, GRN, Pval, Beta_Distance),
              by = c("month" = "Month", "GRN" = "GRN")) %>%
    mutate(
      Clean_GRN = gsub("_direct_\\+/\\+.*$", "", GRN),
      Sig_Star  = case_when(
        is.na(Pval) ~ "", Pval < 0.001 ~ "***",
        Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ "ns"
      ),
      Facet_Label = paste0(month, "M | β=", signif(Beta_Distance, 2), Sig_Star)
    )
  
  ggplot(plot_data, aes(x = Dist_Bin, y = Z_Score,
                        color = Clean_GRN, group = Clean_GRN)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ Facet_Label, ncol = length(unique(plot_data$month))) +
    scale_color_manual(values = "#d73027") +
    labs(
      title = paste("GRN Trajectory:", target_tf, "in", target_ctype),
      x = "Distance to plaque (µm)", y = "Activity (Z-Score)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.background = element_rect(fill = "#E5E7EB", color = "black"),
      legend.position  = "bottom"
    )
}
