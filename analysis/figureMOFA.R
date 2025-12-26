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

# Load configuration
config <- fromJSON("analysis_improved/config/config.json", simplifyVector = FALSE)


##============================================================================##
## 함수: figure_MOFA
## 설명: MOFA 분석 결과를 시각화하는 함수
##       - Factor weight 저장 및 시각화
##       - Factor score boxplot
##       - Factor-변수 상관관계 분석 및 scatter plot
##       - Factor-세포타입 상관관계 heatmap
##============================================================================##
figure_MOFA <- function(
    mofa_path,          # MOFA 모델 객체 경로 (.rds 파일)
    seurat_path,        # Seurat 객체 경로 (.rds 파일)
    output_dir = "MOFA_Output",  # 결과 저장 디렉토리
    n_factors = 4,      # 분석할 Factor 개수
    x_var_list = c("Region_group"),  # Boxplot에서 x축으로 사용할 변수 리스트
    cor_vars = c("T_cells"),  # 상관관계 분석할 변수 리스트
    cell_types = NULL,  # 세포 타입 컬럼명 리스트 (Factor-celltype correlation용)
    view_name = "RNA"   # MOFA view 이름 ("RNA" 또는 "Protein")
) {
  message("==== Plotting MOFA Results ====")

  # -------------------------------------------------------------------------
  # 1. 데이터 로딩 및 디렉토리 설정
  # -------------------------------------------------------------------------
  # MOFA 모델 객체 로드
  mofa_obj <- readRDS(mofa_path)
  # Seurat 객체 로드 (Factor score가 메타데이터에 포함되어 있음)
  seurat_obj <- readRDS(seurat_path)

  # 출력 디렉토리 생성
  figure_dir <- file.path(output_dir, "figure")  # Figure 저장 폴더
  factor_weight_dir <- file.path(output_dir, "factor-weight")  # Factor weight 저장 폴더
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(factor_weight_dir, showWarnings = FALSE, recursive = TRUE)

  # -------------------------------------------------------------------------
  # 2. Factor Weights 저장 (CSV)
  # -------------------------------------------------------------------------
  # MOFA 모델에서 각 Feature(유전자)의 weight를 추출
  # weight가 클수록 해당 factor에서 중요한 유전자임
  message("Saving factor weights...")
  weights_df <- get_weights(mofa_obj, view = view_name, as.data.frame = TRUE)

  # 각 Factor별로 weight를 절대값 기준 내림차순 정렬하여 CSV로 저장
  for (factor_name in unique(weights_df$factor)) {
    sub_df <- filter(weights_df, factor == factor_name) %>%
      arrange(desc(abs(value)))  # 절대값 기준 내림차순 정렬 (중요한 유전자가 위로)
    write.csv(sub_df,
              file = file.path(factor_weight_dir, paste0("weights_", factor_name, ".csv")),
              row.names = FALSE)
  }

  # -------------------------------------------------------------------------
  # 3. Factor Weights 시각화 (Top 20 features)
  # -------------------------------------------------------------------------
  # 각 Factor의 상위 20개 중요 유전자를 barplot으로 시각화
  message("Plotting factor weights...")
  for (i in seq_len(n_factors)) {
    p <- plot_top_weights(mofa_obj, factors = i, nfeatures = 20)
    ggsave(file.path(figure_dir, paste0("factor_weights_Factor", i, ".pdf")), p, dpi = 300)
  }

  # -------------------------------------------------------------------------
  # 4. 메타데이터 준비
  # -------------------------------------------------------------------------
  # Seurat 객체에서 메타데이터 추출 (Factor score + 세포 정보 포함)
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

  # --------------------------------------
  # Factor-celltype correlation heatmap
  # 모든 Factor와 모든 cell type 간의 상관관계를 한눈에 보기 위한 heatmap
  # --------------------------------------
  if (!is.null(cell_types) && length(cell_types) > 0) {
    message("Creating Factor-celltype correlation heatmap...")

    # Factor 컬럼들 추출 (Factor1, Factor2, ...)
    factor_cols <- grep("^Factor", colnames(meta_data), value = TRUE)

    # cell_types 중에서 실제로 meta_data에 존재하는 컬럼만 선택
    celltype_cols <- cell_types[cell_types %in% colnames(meta_data)]

    if (length(celltype_cols) > 0) {
      # 상관계수 행렬 계산: Factor x Cell type
      cor_mat <- cor(meta_data[, factor_cols], meta_data[, celltype_cols],
                     use = "pairwise.complete.obs", method = "pearson")

      # 상관계수 행렬을 long format으로 변환 (ggplot용)
      cor_df <- as.data.frame(cor_mat) %>%
        rownames_to_column("Factor") %>%
        pivot_longer(-Factor, names_to = "CellType", values_to = "Correlation")

      # Heatmap 그리기
      p_heatmap <- ggplot(cor_df, aes(x = CellType, y = Factor, fill = Correlation)) +
        geom_tile(color = "grey90", linewidth = 0.5) +  # 타일 그리기
        geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3, fontface = "bold") +  # 상관계수 값 표시
        scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",  # 파란색-흰색-빨간색 gradient
                             midpoint = 0, limit = c(-1, 1), name = "Pearson R") +
        labs(title = "Correlation between Factors and Cell proportions",
             x = "Cell Type", y = "Factor") +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
              axis.text.y = element_text(face = "bold"),
              plot.title = element_text(hjust = 0.5, face = "bold"))

      # PDF로 저장
      ggsave(file.path(figure_dir, "factor_celltype_correlation_heatmap.pdf"),
             p_heatmap, width = 8, height = 5, dpi = 300)
      message("Factor-celltype correlation heatmap saved.")
    } else {
      message("Warning: No valid cell type columns found in metadata.")
    }
  }

}

# Build params_list from config
params_list <- lapply(config$run_list, function(run) {
  list(
    output_dir = file.path(config$paths$base_output_dir, run$name),
    n_factors = run$n_factors
  )
})

# Run figure generation for each configuration
for (i in seq_along(params_list)) {
  cat("\n============================\n")
  cat("Generating figures for job", i, "/", length(params_list), "\n")
  cat("Output:", params_list[[i]]$output_dir, "\n")
  cat("============================\n\n")

  p <- params_list[[i]]

  figure_MOFA(
    seurat_path = file.path(p$output_dir, "seurat_obj_raw.rds"),
    mofa_path = file.path(p$output_dir, "mofa_model_raw.rds"),
    output_dir = p$output_dir,
    n_factors = p$n_factors,
    x_var_list = unlist(config$fixed_params$x_var_list),
    cor_vars = unlist(config$fixed_params$cor_vars),
    cell_types = unlist(config$fixed_params$cell_types),
    view_name = config$fixed_params$view_name
  )
}

