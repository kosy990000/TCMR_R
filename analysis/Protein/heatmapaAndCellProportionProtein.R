suppressPackageStartupMessages({
  library(Seurat)               # NormalizeData, ScaleData
  library(DESeq2)               # DESeq2 VST normalization
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

config <- fromJSON("analysis/config/config_visualization_protein_new.json", simplifyVector = FALSE)

run_heatmap_and_cellProportion <- function(
    seurat_path,
    output_dir,
    cell_types,
    n_factor,
    weight_base_dir,
    normalization_method = "CLR",  # "CLR" or "DESeq"
    group_var = "Region_group",
    group_order = c(
      "glom_CTRL_pre", "glom_CTRL_post",
      "glom_TCMR_pre", "glom_TCMR_post",
      "tubulointer_CTRL_pre", "tubulointer_CTRL_post",
      "tubulointer_TCMR_pre", "tubulointer_TCMR_post"
    )
) {
  message("==== Start: run_heatmap_and_cellProportion ====")
  message(paste0("Normalization method: ", normalization_method))

  figure_dir <- file.path(output_dir, "figure")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  weight_path <- file.path(weight_base_dir, "factor-weight", paste0("weights_Factor", n_factor, ".csv"))

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
  # *** Normalization Step ***
  # (Heatmap을 그리기 위해 여기서 수행)
  #----------------------------------------------#
  if (normalization_method == "CLR") {
    message(">>> Performing CLR normalization for Protein data...")

    seurat_obj <- NormalizeData(
      seurat_obj,
      assay = "Protein",
      normalization.method = "CLR",
      margin = 2  # cell 에대해 정규화
    )
    message(" -> CLR normalization completed.")

  } else if (normalization_method == "DESeq") {
    message(">>> Performing DESeq2 VST normalization for Protein data...")

    # Protein count 데이터 추출
    expression_data <- GetAssayData(seurat_obj, assay = "Protein", layer = "counts")

    # --- DESeq2 variance stabilization ---
    message("Variance stabilizing transformation...")
    dds <- DESeqDataSetFromMatrix(
      countData = expression_data,
      colData = meta,
      design = ~1
    )
    rm(expression_data); gc()

    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vst_mat <- assay(vsd)
    rm(dds, vsd); gc()

    # VST 정규화된 데이터를 Seurat 객체의 data layer에 저장
    seurat_obj[["Protein"]]@layers$data <- vst_mat
    message(" -> VST normalization completed.")

  } else {
    stop("normalization_method must be either 'CLR' or 'DESeq'")
  }

  # Z-score scaling (Seurat ScaleData 사용)
  message(" -> Running ScaleData for Z-score transformation...")
  seurat_obj <- ScaleData(
    seurat_obj,
    assay = "Protein"
  )
  message(" -> Scaling completed.")

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


# Single execution with config
vis_config <- config$visualization

message("========================================")
message("Running heatmap and cell proportion visualization:")
message("seurat_path    : ", vis_config$seurat_path)
message("output_dir     : ", vis_config$output_dir)
message("weight_base_dir: ", vis_config$weight_base_dir)
message("group_var      : ", vis_config$group_var)
message("n_factor       : ", vis_config$n_factor)
message("========================================")

run_heatmap_and_cellProportion(
  seurat_path = vis_config$seurat_path,
  output_dir = vis_config$output_dir,
  weight_base_dir = vis_config$weight_base_dir,
  group_var = vis_config$group_var,
  group_order = unlist(vis_config$group_order),
  cell_types = unlist(vis_config$cell_types),
  n_factor = vis_config$n_factor,
  normalization_method = vis_config$normalization_method
)
