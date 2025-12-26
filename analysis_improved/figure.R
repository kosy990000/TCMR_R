suppressPackageStartupMessages({
  library(MOFA2)
  library(Seurat)
  library(sva)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(matrixStats)
  library(RColorBrewer)
  library(purrr)
  library(jsonlite)
})

config <- fromJSON("config.json", simplifyVector = FALSE)

figure_MOFA <- function(
    mofa_path,
    seurat_path,
    output_dir = "MOFA_Output",
    n_factors = 4,
    x_var_list = c("Region_group"),
    cor_vars = c("T_cells")
) {
  message("==== Plotting MOFA Results ====")
  seurat_obj <- readRDS(seurat_path)
  figure_dir <- file.path(output_dir, "figure")


  meta_data <- seurat_obj@meta.data


  # 순서 변경
  meta_data$Region_group <- factor(meta_data$Region_group, levels = c(
    "glom_CTRL_pre", "glom_CTRL_post",
    "glom_TCMR_pre", "glom_TCMR_post",
    "tubulointer_CTRL_pre", "tubulointer_CTRL_post",
    "tubulointer_TCMR_pre", "tubulointer_TCMR_post"
  ))




  # --------------------------------------
  plot_all_factors <- function(meta_data, n_factors, fill_palette = "Set3", x_var, set_ncol = 2) {

    #long_df 로 factor score 변경
    long_df <- meta_data %>%
      pivot_longer(cols = starts_with("Factor"),
                   names_to = "Factor", values_to = "Score") %>%
      filter(Factor %in% paste0("Factor", 1:n_factors))
    n_groups <- length(unique(long_df[[x_var]]))
    palette_expanded <- colorRampPalette(brewer.pal(8, fill_palette))(n_groups)


    ggplot(long_df, aes(x = !!sym(x_var), y = Score, fill = !!sym(x_var))) +
      geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 0.9) +
      scale_fill_manual(values = palette_expanded) +
      facet_wrap(~Factor, ncol = set_ncol, scales = "free_y") +
      labs(
        x = x_var,
        y = "Factor Score"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),  # 글씨 크고 두껍게
        axis.title.x = element_text(size = 14, face = "bold"),                         # x축 제목도 강조
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
      )
  }


  for (x_var in x_var_list) {
    message(x_var)
    p <- plot_all_factors(
      meta_data = meta_data,
      n_factors = n_factors,
      x_var = x_var
    )
    ggsave(file.path(figure_dir, paste0("Factor_boxplot_", x_var, ".pdf")), p)
  }


  p3 <- plot_all_factors(meta_data, n_factor = 3, x_var = "Region_group", set_ncol = 1)
  out_path <- file.path(figure_dir, paste0("Factor_boxplot_agg3", ".pdf"))
  ggsave(out_path, p3, width = 6, height = 10)



}



for (params in config$factor_analysis$params_list) {
  output_dir <- params$output_dir
  figure_MOFA(
    seurat_path = file.path(output_dir, "seurat_obj_raw.rds"),
    mofa_path = file.path(output_dir, "mofa_model_raw.rds"),
    output_dir = params$output_dir,
    n_factors = params$n_factors,
    x_var_list = unlist(params$x_var_list),
    cor_vars = unlist(params$cor_vars)
  )
}

