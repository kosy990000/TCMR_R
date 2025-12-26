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


  # 상관분석
  # --------------------------------------
  cor_summary_table <- function(meta_df, n_factors, vars) {
    cor_results <- list()
    for (i in 1:n_factors) {
      factor_col <- paste0("Factor", i)
      if (!factor_col %in% names(meta_df)) next

      for (var in vars) {
        if (!var %in% names(meta_df)) next
        res <- cor.test(meta_df[[factor_col]], meta_df[[var]], method = "pearson")
        corr <- round(res$estimate, 3)
        p_val <-  signif(res$p.value, 3)
        # ------------------------------------- #
        cor_results[[length(cor_results) + 1]] <- data.frame(
          Factor = factor_col,
          Variable = var,
          Correlation = corr,
          P_value = p_val
        )
        p <- ggplot(meta_df, aes_string(x = var, y = factor_col, color = "Region_group")) +
          geom_point(alpha = 0.8, size = 1) +
          geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
          labs(
            title = paste0(factor_col, " vs ", var," proportion" ," (R=", corr, ", p=", p_val, ")"),
            x = paste0(var, " proportion"),
            y = paste0(factor_col, " score"),
            color = "Category"
          ) +
          theme_minimal(base_size = 13)


        ggsave(
          filename = paste0(figure_dir, "/", factor_col, "_", var, ".pdf"),
          plot = p,
          width = 8,
          height = 6,
          dpi = 300
        )
      }
    }
    cor_table <- do.call(rbind, cor_results)
    cor_table <- cor_table[order(cor_table$Factor, cor_table$Variable), ]
    rownames(cor_table) <- NULL
    return(cor_table)
  }

  cor_result_df <- cor_summary_table(meta_data, n_factors, cor_vars)
  write.csv(cor_result_df, file.path(output_dir, "correlation_summary.csv"), row.names = FALSE)
  message("plots and correlation summary saved.")

}

params_list <- list(
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "result/CLSProNorfactor4",
    top_n = 4000,
    n_factors = 4,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  ),
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "result/CLSProNorfactor6",
    top_n = 4000,
    n_factors = 6,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  ),
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "result/CLSProNorfactor8",
    top_n = 4000,
    n_factors = 8,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  )
)
# -------------------------------------------------------------------
# 2) 반복 실행
# -------------------------------------------------------------------
for (i in seq_along(params_list)) {

  cat("\n============================\n")
  cat("Running MOFA job", i, "\n")
  cat("============================\n\n")

  p <- params_list[[i]]
  output_dir <- p$output_dir
  figure_MOFA(
    seurat_path = file.path(output_dir, "seurat_obj_raw.rds"),
    mofa_path = file.path(output_dir, "mofa_model_raw.rds"),
    output_dir = p$output_dir,
    n_factors = p$n_factors,
    x_var_list = c("Region_group"),
    cor_vars = c("T_cells")
  )
}

##------------------------------------------------------------




