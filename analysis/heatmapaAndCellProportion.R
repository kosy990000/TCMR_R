suppressPackageStartupMessages({
  library(Seurat)               # NormalizeData, ScaleData
  library(dplyr)                # 데이터 처리
  library(tidyr)                # pivot_longer
  library(tibble)               # rownames_to_column
  library(ggplot2)              # 시각화
  library(scales)               # 퍼센트 스케일
  library(ComplexHeatmap)       # Heatmap
  library(circlize)             # colorRamp2
  library(ggsci)                # pal_d3
  library(jsonlite)             # JSON
  library(grid)                 # grid
})

config <- fromJSON("test_protein.json", simplifyVector = FALSE)

run_heatmap_and_cellProportion <- function(
    seurat_path,
    output_dir,
    cell_types,
    n_factor,
    group_var = "Region_group",
    group_order = c(
      "glom_CTRL_pre", "glom_CTRL_post",
      "glom_TCMR_pre", "glom_TCMR_post",
      "tubulointer_CTRL_pre", "tubulointer_CTRL_post",
      "tubulointer_TCMR_pre", "tubulointer_TCMR_post"
    )
) {
  message("==== Start: run_heatmap_and_cellProportion ====")
  
  figure_dir <- file.path(output_dir, "figure")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  weight_path <- file.path(output_dir, "factor-weight", paste0("weights_Factor", n_factor, ".csv"))
  
  # 1. Load Object & Setup Metadata
  seurat_obj <- readRDS(seurat_path)
  DefaultAssay(seurat_obj) <- "Protein"
  
  cell_types <- unlist(cell_types)
  group_order <- unlist(group_order)
  
  # 메타데이터 추출 및 Factor 설정 (정규화 전에도 수행 가능)
  meta <- seurat_obj@meta.data
  meta[[group_var]] <- factor(meta[[group_var]], levels = group_order)
  
  #----------------------------------------------#
  # 1. Cell proportion bar plot (No Normalization Needed)
  #----------------------------------------------#
  message("Calculating cell proportions...")
  
  celltype_cols <- cell_types
  n_cell_types <- length(celltype_cols)
  
  cell_meta <- meta %>%
    select(all_of(celltype_cols), all_of(group_var)) %>%
    pivot_longer(cols = all_of(celltype_cols),
                 names_to = "celltype",
                 values_to = "value")
  
  cell_prop <- cell_meta %>%
    group_by(!!sym(group_var), celltype) %>%
    summarise(total_value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(!!sym(group_var)) %>%
    mutate(freq = total_value / sum(total_value))
  
  cell_prop[[group_var]] <- factor(cell_prop[[group_var]], levels = group_order)
  
  my_palette <- pal_d3("category20")(n_cell_types)
  p_cell <- ggplot(cell_prop, aes(x = !!sym(group_var), y = freq, fill = celltype)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = my_palette) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Proportion (%)", fill = "Cell Type") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
  
  ggsave(file.path(figure_dir, "cell_prop.pdf"), p_cell, dpi = 300, width = 12, height = 8)
  message("Cell proportion plot saved.")
  
  #----------------------------------------------#
  # 2. Factor–cell proportion correlation (No Normalization Needed)
  #----------------------------------------------#
  message("Computing factor-celltype correlation...")
  
  factor_cols <- grep("^Factor", colnames(meta), value = TRUE)
  cor_mat <- cor(meta[, factor_cols], meta[, celltype_cols],
                 use = "pairwise.complete.obs", method = "pearson")
  
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("Factor") %>%
    pivot_longer(-Factor, names_to = "CellType", values_to = "Correlation")
  
  p2 <- ggplot(cor_df, aes(x = CellType, y = Factor, fill = Correlation)) +
    geom_tile(color = "grey90", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limit = c(-1, 1), name = "Pearson R") +
    labs(title = "Correlation between Factors and Cell proportions",
         x = "Cell Type", y = "Factor") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(figure_dir, "factor_celltype_correlation_heatmap.pdf"),
         p2, width = 8, height = 5, dpi = 300)
  message("Correlation heatmap saved.")
  
  #----------------------------------------------#
  # *** Data Normalization & Scaling Step ***
  # (Heatmap을 그리기 위해 여기서 수행)
  #----------------------------------------------#
  message("Performing Normalization (CLR) and Scaling for Heatmap...")
  
  if (IsGlobal(seurat_obj, assay = "Protein", slot = "data") || 
      length(GetAssayData(seurat_obj, assay = "Protein", layer = "data")) == 0) {
    
    message(" -> Normalization not found. Running NormalizeData (CLR)...")
    seurat_obj <- NormalizeData(
      seurat_obj,
      assay = "Protein",
      normalization.method = "CLR",
      margin = 2
    )
  } else {
    message(" -> Normalized data found. Skipping NormalizeData.")
  }
  
  seurat_obj <- ScaleData(
    seurat_obj,
    assay = "Protein"
  )
  
  #----------------------------------------------#
  # 3. Factor weight-based gene heatmap
  #----------------------------------------------#
  message("Drawing factor weight-based gene heatmap...")
  
  weight_list <- read.csv(weight_path)
  weight_sorted <- weight_list[order(weight_list$value, decreasing = TRUE), ]
  selected_genes <- c(head(weight_sorted$feature, 20), tail(weight_sorted$feature, 20))
  
  # 이제 ScaleData가 완료되었으므로 scale.data 레이어 접근 가능
  mat <- GetAssayData(
    seurat_obj,
    assay = "Protein",
    layer = "scale.data"
  )[selected_genes, ]
  
  # Heatmap 정렬을 위해 meta 정보를 다시 활용
  mat <- mat[, order(meta[[group_var]])]
  region_group <- meta[[group_var]][order(meta[[group_var]])]
  names(region_group) <- colnames(mat)
  
  region_colors <- structure(pal_d3("category20")(length(group_order)), names = group_order)
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ha <- HeatmapAnnotation(
    Group = region_group,
    col = list(Group = region_colors),
    annotation_name_side = "right"
  )
  
  p_heatmap <- Heatmap(
    mat, name = "Z-score", top_annotation = ha, col = col_fun,
    cluster_rows = TRUE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE
  )
  
  pdf(file.path(figure_dir, "heatmap_top_bottom20.pdf"), width = 10, height = 8)
  draw(p_heatmap)
  dev.off()
  message("Gene expression heatmap saved.")
  message("==== Completed: run_heatmap_and_cellProportion ====")
}


# Loop execution
for (params in config$factor_analysis$params_list) {
  output_dir <- params$output_dir
  
  message("========================================")
  message("Running factor analysis with parameters:")
  message("output_dir   : ", output_dir)
  message("seurat_path  : ", file.path(output_dir, "seurat_obj_raw.rds"))
  message("group_var    : ", as.character(params$group_var))
  message("n_factor     : ", params$factor_num)
  message("========================================")
  
  run_heatmap_and_cellProportion(
    seurat_path = file.path(output_dir, "seurat_obj_raw.rds"),
    output_dir = params$output_dir,
    group_var =  as.character(params$group_var),
    group_order = unlist(params$group_order),
    cell_types = unlist(params$cell_types),
    n_factor = params$factor_num 
  )
}
