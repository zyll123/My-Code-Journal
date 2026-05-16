

# ==========================================
# 【Seurat_R_Pipeline】来源: 2026-05-16-Differential_GRN_clean.ipynb
# 合并日期: 2026-05-16
# ==========================================
[markdown]
# Differential GRN Analysis: AD vs WT
---
**Goal:** For each cell type × month, identify regulons with significantly different activity between AD and WT (Mann-Whitney U test), filtered to regulons with ≥10% activation rate in at least one group.

**Input:**
- `combined`: Seurat object with metadata (genotype, month, class_id_label)
- `activity`: Binarized TF activity matrix (0/1 per cell)
- `AUC_gene`: Regulon AUC activity matrix (continuous scores)

[markdown]
## 1. Libraries & Data Loading

[code]
# ==========================================
# 1.1 Libraries
# ==========================================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

print("Libraries loaded.")

[code]
# ==========================================
# 1.2 Read all input data
# ==========================================
# Seurat object (metadata only needed)
combined <- readRDS('/public/home/yzhang/zy_data/Project_AD_st/Multiome/Raw_data_0423_full/QC_0512/Combine/Combine_RNA_ATAC_anno_rnafromatac_260219.rds')

# Binarized TF activity matrix (0/1)
activity <- read.csv('/home1/yzhang/Project_SN_Multiome/Analysis/Scenic+/Version_2_mouse/binarize/b_matrix(1).csv')

# Regulon AUC matrix — keep only +/+ columns
AUC_gene <- read.csv('/home1/yzhang/Project_SN_Multiome/Analysis/Scenic+/Version_2_mouse/gene_AUC_matrix.csv', check.names = FALSE)
cell_col_name <- colnames(AUC_gene)[1]
plus_cols <- colnames(AUC_gene)[grepl("\\+/\\+", colnames(AUC_gene))]
AUC_gene <- AUC_gene[, c(cell_col_name, plus_cols)]

print(paste("AUC columns (after +/+ filter):", ncol(AUC_gene)))

[markdown]
## 2. Align Cell IDs Across Datasets
Convert `barcode-sample` → `sample_barcode` format to match.

[code]
# ==========================================
# 2.1 Fix cell ID format & align
# ==========================================
fix_cell_ids <- function(ids) {
  gsub("^([^-]+)-(.*)$", "\\2_\\1", ids)
}

# Extract & fix metadata
meta <- combined@meta.data
meta$cell_id <- rownames(meta)

# Fix activity and AUC_gene IDs
activity$X <- fix_cell_ids(activity$X)
rownames(activity) <- activity$X

AUC_gene[[cell_col_name]] <- fix_cell_ids(AUC_gene[[cell_col_name]])
rownames(AUC_gene) <- AUC_gene[[cell_col_name]]

# Intersect
common_cells <- Reduce(intersect, list(meta$cell_id, rownames(activity), rownames(AUC_gene)))
print(paste("Common cells:", length(common_cells)))

meta <- meta[common_cells, ]
activity <- activity[common_cells, ]
AUC_gene <- AUC_gene[common_cells, ]

[markdown]
## 3. Build TF → Regulon Mapping

[code]
# ==========================================
# 3.1 Map TF names to AUC regulon column names
#     e.g. "Ahctf1" → "Ahctf1_direct_+/+_(39g)"
# ==========================================
auc_cols <- setdiff(colnames(AUC_gene), cell_col_name)

regulon_mapping <- data.frame(
  AUC_name = auc_cols,
  TF_name  = sapply(strsplit(auc_cols, "_"), `[`, 1),
  stringsAsFactors = FALSE
) %>%
  filter(TF_name %in% colnames(activity))

print(paste("Mapped regulons:", nrow(regulon_mapping)))

[markdown]
## 3.5 Redefine Cell Types for Differential Analysis
Neurons (class_id_label ending in Glut/GABA) remain unchanged. Non-neurons are refined:
- **OPC-Oligo** → split into **OPC** and **Oligo** based on `Final_subclass_name`
- **Immune** → split into **Micro** (Microglia) and **Other Immune** based on `Final_subclass_name`
- Set `split_by_brain <- TRUE` to further split non-neuronal types by brain region (`brain` column)

[code]
# ==========================================
# 3.5 Redefine cell type labels
# ==========================================
# --- Parameter: split non-neuronal types by brain region? ---
split_by_brain <- FALSE  # Set to TRUE to analyze per brain region

# --- Check brain column exists ---
if (split_by_brain) {
  if (!"brain" %in% colnames(meta)) {
    stop("'brain' column not found in metadata. Available columns: ",
         paste(head(colnames(meta), 30), collapse = ", "))
  }
  print("Brain regions in metadata:")
  print(sort(table(meta$brain), decreasing = TRUE))
}

# --- Helper: check if neuron (ends with Glut/GABA, case-insensitive) ---
is_neuron <- function(label) {
  grepl("Glut$|GABA$|glut$|gaba$", label)
}

# --- Store original label ---
meta$class_id_label_original <- meta$class_id_label

# --- Redefine labels ---
meta$class_id_label <- sapply(seq_len(nrow(meta)), function(i) {
  lbl <- meta$class_id_label_original[i]
  
  # Neuron: keep unchanged (brain never appended for neurons)
  if (is_neuron(lbl)) {
    return(lbl)
  }
  
  sub <- meta$Final_subclass_name[i]
  
  # OPC-Oligo → OPC or Oligo
  if (grepl("OPC-Oligo", lbl, ignore.case = TRUE)) {
    if (!is.na(sub) && grepl("OPC", sub, ignore.case = TRUE) && !grepl("Oligo", sub, ignore.case = TRUE)) {
      new_lbl <- gsub("OPC-Oligo", "OPC", lbl, ignore.case = TRUE)
    } else {
      new_lbl <- gsub("OPC-Oligo", "Oligo", lbl, ignore.case = TRUE)
    }
  }
  # Immune → Micro or Other Immune
  else if (grepl("Immune", lbl, ignore.case = TRUE)) {
    if (!is.na(sub) && grepl("Microglia", sub, ignore.case = TRUE)) {
      new_lbl <- gsub("Immune", "Micro", lbl, ignore.case = TRUE)
    } else {
      new_lbl <- gsub("Immune", "Other Immune", lbl, ignore.case = TRUE)
    }
  } else {
    new_lbl <- lbl
  }
  
  # Optionally append brain region for non-neurons
  if (split_by_brain && !is_neuron(lbl)) {
    br <- meta$brain[i]
    if (!is.na(br) && br != "") {
      new_lbl <- paste(new_lbl, br, sep = "_")
    }
  }
  
  return(new_lbl)
})

print("Redefined cell type labels:")
print(sort(unique(meta$class_id_label)))
print(paste("Total cell types for analysis:", length(unique(meta$class_id_label))))

if (split_by_brain) {
  print("Cell counts per redefined type (non-neuron) × brain region:")
  nn_idx <- !is_neuron(meta$class_id_label_original)
  print(sort(table(meta$class_id_label[nn_idx])))
}

[markdown]
## 4. Differential Analysis
For each **cell type × month** combination:
1. Filter regulons ≥10% activation in AD **or** WT
2. Mann-Whitney U test on AUC scores
3. BH-FDR correction

[code]
# ==========================================
# 4.1 Main loop: cell_type × month → Mann-Whitney U
# ==========================================
cell_types <- unique(meta$class_id_label)
months     <- unique(meta$month)
results_list <- list()

for (ct in cell_types) {
  for (m in months) {
    
    # Subset cells for this cell_type + month
    cells_ct_m <- meta$cell_id[meta$class_id_label == ct & meta$month == m]
    if (length(cells_ct_m) < 20) next
    
    meta_sub <- meta[cells_ct_m, ]
    act_sub  <- activity[cells_ct_m, ]
    auc_sub  <- AUC_gene[cells_ct_m, ]
    
    cells_AD <- meta_sub$cell_id[meta_sub$genotype == "AD"]
    cells_WT <- meta_sub$cell_id[meta_sub$genotype == "WT"]
    if (length(cells_AD) < 3 | length(cells_WT) < 3) next
    
    # Activation ratio per group
    tf_AD <- act_sub[cells_AD, regulon_mapping$TF_name, drop = FALSE]
    tf_WT <- act_sub[cells_WT, regulon_mapping$TF_name, drop = FALSE]
    ratio_AD <- colMeans(tf_AD == 1, na.rm = TRUE)
    ratio_WT <- colMeans(tf_WT == 1, na.rm = TRUE)
    
    # Keep TFs active in ≥10% of either group
    keep_idx <- which((ratio_AD >= 0.1 | ratio_WT >= 0.1) %in% TRUE)
    if (length(keep_idx) == 0) next
    
    valid_tfs     <- regulon_mapping$TF_name[keep_idx]
    valid_auc     <- regulon_mapping$AUC_name[keep_idx]
    ratio_AD_kept <- ratio_AD[keep_idx]
    ratio_WT_kept <- ratio_WT[keep_idx]
    
    # Test each valid regulon
    res <- lapply(seq_along(valid_auc), function(i) {
      reg  <- valid_auc[i]
      tf   <- valid_tfs[i]
      
      val_AD <- as.numeric(auc_sub[cells_AD, reg])
      val_WT <- as.numeric(auc_sub[cells_WT, reg])
      wt     <- wilcox.test(val_AD, val_WT, exact = FALSE)
      
      mean_AD <- mean(val_AD, na.rm = TRUE)
      mean_WT <- mean(val_WT, na.rm = TRUE)
      log2FC  <- log2((mean_AD + 1e-5) / (mean_WT + 1e-5))
      
      data.frame(
        CellType        = ct,
        Month           = m,
        Regulon         = reg,
        TF              = tf,
        Active_Ratio_AD = ratio_AD_kept[i],
        Active_Ratio_WT = ratio_WT_kept[i],
        Mean_AD         = mean_AD,
        Mean_WT         = mean_WT,
        log2FC          = log2FC,
        p_val           = wt$p.value,
        stringsAsFactors = FALSE
      )
    })
    
    res_df <- do.call(rbind, res)
    res_df$p_val_adj <- p.adjust(res_df$p_val, method = "BH")
    results_list[[paste(ct, m, sep = "_")]] <- res_df
  }
}

# Aggregate & filter significant
final_diff_GRN <- do.call(rbind, results_list)
rownames(final_diff_GRN) <- NULL
significant_GRN <- final_diff_GRN %>% filter(p_val_adj < 0.05)

print(paste("Total tests:", nrow(final_diff_GRN),
            "| Significant (FDR<0.05):", nrow(significant_GRN)))
head(significant_GRN, 8)

[markdown]
## 5. Results: Coverage & Export

[code]
# ==========================================
# 5.1 How many cell types does each regulon span?
# ==========================================
regulon_counts <- significant_GRN %>%
  group_by(Regulon, TF) %>%
  summarise(Num_CellTypes = n_distinct(CellType), .groups = "drop")

ggplot(regulon_counts, aes(x = Num_CellTypes)) +
  geom_bar(fill = "#8DA0CB", color = "black", alpha = 0.85, width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = seq(1, max(regulon_counts$Num_CellTypes), by = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(
    x        = "Number of cell types with significant difference",
    y        = "Count of regulons",
    title    = "Regulon specificity vs ubiquity",
    subtitle = "How many cell types each differentially-active regulon spans"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

[code]
# ==========================================
# 5.2 Export
# ==========================================
# Broadly-dysregulated regulons (≥20 cell types)
print("Regulons spanning ≥20 cell types:")
print(regulon_counts %>% filter(Num_CellTypes >= 20))

# Save
# write.csv(significant_GRN,
#           '/home1/yzhang/Project_SN_Multiome/Analysis/Scenic+/Version_2_mouse/Differential_GRN_0515.csv',
#           row.names = FALSE)

[markdown]
## 6. Cross-Brain-Region Comparison
Only runs when `split_by_brain <- TRUE`. Parses the `CellType` column to extract base cell type and brain region, then compares which regulons are shared vs brain-region-specific.

[code]
# ==========================================
# 6.1 Cross-brain-region comparison
# ==========================================
if (split_by_brain) {
  
  # Parse CellType into base type and brain region
  # e.g. "31 OPC_FC" → BaseCellType = "31 OPC", BrainRegion = "FC"
  significant_GRN <- significant_GRN %>%
    mutate(
      BaseCellType = str_replace(CellType, "_[^_]+$", ""),
      BrainRegion  = str_extract(CellType, "[^_]+$")
    )
  
  # ========================================
  # 6.1.1 Regulon overlap across brain regions per base cell type
  # ========================================
  regulon_brain_span <- significant_GRN %>%
    group_by(BaseCellType, Month, Regulon, TF) %>%
    summarise(
      N_BrainRegions = n_distinct(BrainRegion),
      BrainRegions   = paste(sort(unique(BrainRegion)), collapse = ", "),
      .groups = "drop"
    )
  
  # Summary: shared (>=2 brain regions) vs region-specific regulons
  summary_shared <- regulon_brain_span %>%
    group_by(BaseCellType) %>%
    summarise(
      Total_Sig_Regulons    = n(),
      Shared_Across_Regions = sum(N_BrainRegions >= 2),
      Region_Specific       = sum(N_BrainRegions == 1),
      .groups = "drop"
    ) %>%
    arrange(desc(Shared_Across_Regions))
  
  print("=== 1. Shared vs brain-region-specific regulons per cell type ===")
  print(summary_shared)
  
  # ========================================
  # 6.1.2 Shared regulons (significant in ≥2 brain regions)
  # ========================================
  shared_regulons <- regulon_brain_span %>%
    filter(N_BrainRegions >= 2) %>%
    arrange(BaseCellType, Month, desc(N_BrainRegions))
  
  print("=== 2. Regulons significant in ≥2 brain regions (same cell type & month) ===")
  print(shared_regulons)
  
  # ========================================
  # 6.1.3 Sig regulon count matrix: cell type × brain region
  # ========================================
  region_counts <- significant_GRN %>%
    count(BaseCellType, BrainRegion, name = "N_Sig_Regulons") %>%
    pivot_wider(
      names_from  = BrainRegion,
      values_from = N_Sig_Regulons,
      values_fill = 0
    )
  
  print("=== 3. Significant regulon counts: cell type × brain region ===")
  print(as.data.frame(region_counts))
  
  # ========================================
  # 6.1.4 Pairwise Jaccard overlap between brain regions per cell type
  # ========================================
  cell_types_base <- unique(significant_GRN$BaseCellType)
  overlap_list <- list()
  
  for (bt in cell_types_base) {
    sig_bt <- significant_GRN %>% filter(BaseCellType == bt)
    regions <- unique(sig_bt$BrainRegion)
    if (length(regions) < 2) next
    
    for (i in 1:(length(regions) - 1)) {
      for (j in (i + 1):length(regions)) {
        r1 <- regions[i]
        r2 <- regions[j]
        set1 <- sig_bt$Regulon[sig_bt$BrainRegion == r1]
        set2 <- sig_bt$Regulon[sig_bt$BrainRegion == r2]
        inter <- length(intersect(set1, set2))
        un    <- length(union(set1, set2))
        jac   <- ifelse(un > 0, round(inter / un, 3), 0)
        
        overlap_list[[length(overlap_list) + 1]] <- data.frame(
          BaseCellType = bt,
          Region1 = r1,
          Region2 = r2,
          N_Shared = inter,
          Jaccard  = jac,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(overlap_list) > 0) {
    overlap_df <- do.call(rbind, overlap_list)
    print("=== 4. Pairwise Jaccard overlap of sig regulons between brain regions ===")
    print(overlap_df %>% arrange(desc(Jaccard)))
  }
  
} else {
  print("split_by_brain = FALSE — skipping cross-brain-region comparison.")
  print("Set split_by_brain <- TRUE in Section 3.5 and re-run to enable.")
}

# ==========================================
# 【Seurat_R_Pipeline】来源: 2026-05-16.md
# 合并日期: 2026-05-16
# ==========================================
library(seurat)
data=read.RDS(data.rds)