suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(jsonlite)
})

config <- fromJSON("analysis/config/config_no_con_new.json", simplifyVector = FALSE)

output_dir  <- file.path(config$paths$base_output_dir, "Deseq4000Factor4no_con_new")
seurat_path <- file.path(output_dir, "seurat_obj_raw.rds")
figure_dir  <- file.path(output_dir, "figure")
n_factors   <- 4
cell_types  <- unlist(config$fixed_params$cell_types)

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

seurat_obj <- readRDS(seurat_path)
meta_data  <- seurat_obj@meta.data

factor_cols   <- paste0("Factor", 1:n_factors)
factor_cols   <- factor_cols[factor_cols %in% colnames(meta_data)]
celltype_cols <- cell_types[cell_types %in% colnames(meta_data)]

if (length(factor_cols) == 0)  stop("Factor 컬럼이 메타데이터에 없습니다.")
if (length(celltype_cols) == 0) stop("유효한 cell type 컬럼이 없습니다.")

cor_mat <- cor(meta_data[, factor_cols], meta_data[, celltype_cols],
               use = "pairwise.complete.obs", method = "pearson")

cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "CellType", values_to = "Correlation")

p_heatmap <- ggplot(cor_df, aes(x = CellType, y = Factor, fill = Correlation)) +
  geom_tile(color = "grey90", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, limit = c(-1, 1), name = "Pearson R") +
  labs(title = "Correlation between Factors and Cell proportions (no_con, Factor4)",
       x = "Cell Type", y = "Factor") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

out_path <- file.path(figure_dir, "factor_celltype_correlation_heatmap_corr3.pdf")
ggsave(out_path, p_heatmap, width = 10, height = 5, dpi = 300)
message("Saved: ", out_path)
