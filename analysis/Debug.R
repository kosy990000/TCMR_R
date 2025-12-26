
suppressPackageStartupMessages({
  library(MOFA2)
  library(Seurat)
  library(sva)
  library(dplyr)
  library(tidyr)
  library(DESeq2)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(matrixStats)
  library(RColorBrewer)
  library(purrr)
})

run_MOFA_analysis <- function(
    seurat_path,
    output_dir = "MOFA_Output",
    n_factors = 4,
    subset_region = NULL,
    x_var_list = c("Region_group", "ROI_clusters"),
    cor_vars = c("nCount_RNA", "T_cells")
) {
  # --------------------------------------
  # load_package
  # --------------------------------------
  
  message("Seurat 객체 불러오는 중...")
  seurat_obj <- readRDS(seurat_path)
  
  # --------------------------------------
  #  선택적 subset
  # --------------------------------------
  if (!is.null(subset_region)) {
    seurat_obj <- subset(seurat_obj, subset = Region == subset_region)
    message(paste0("Subset 적용: Region == ", subset_region))
  }
  
  # --------------------------------------
  #  데이터 준비
  # --------------------------------------
  message("RNA counts 추출 및 메타데이터 구성 중...")
  expression_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  meta <- seurat_obj@meta.data[, c("Group", "Region_group", "Region", "ScanLabel", "Condition")]
  meta$sample <- colnames(expression_data)
  
  
  root_dir <- getwd()
  output_dir <- file.path(root_dir, output_dir)
  figure_dir <- file.path(output_dir, "figure")
  # 없으면 만들기
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(figure_dir, showWarnings = FALSE)
  
  
  dds <- DESeqDataSetFromMatrix(
    countData = expression_data,
    colData   = meta,
    design    = ~ 1
  )
  
  # 정규화 및 분산 안정화 변환
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  vst_mat <- assay(vsd)
  
  top4000 <- order(rowVars(vst_mat), decreasing = TRUE)[1:4000]
  vst_mat <- vst_mat[top4000, ]
  
  # --------------------------------------
  #  MOFA용 데이터 변환
  # --------------------------------------
  message("MOFA용 데이터 변환 중...")
  dt <- as.data.frame(vst_mat) %>%
    tibble::rownames_to_column("feature") %>%
    pivot_longer(cols = -feature, names_to = "sample", values_to = "value") %>%
    left_join(meta, by = "sample") %>%
    # # nogroup
    mutate(view = "RNA_counts", group = Condition) %>% # group을 condition으로 CTRL/TCMR
    select(sample, group, feature, view, value)
  
  rm(expression_data)
  
  # --------------------------------------
  #  MOFA 객체 생성 및 설정
  # --------------------------------------
  message(" MOFA 객체 준비 중...")
  mofa_obj <- create_mofa_from_df(dt)
  rm(dt)
  
  data_opts <- get_default_data_options(mofa_obj)
  model_opts <- get_default_model_options(mofa_obj)
  train_opts <- get_default_training_options(mofa_obj)
  
  # seed 및 factor 추출
  model_opts$num_factors <- n_factors
  train_opts$seed <- 2025
  
  
  mofa_obj <- prepare_mofa(mofa_obj,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)
  
  # --------------------------------------
  #  학습
  # --------------------------------------
  message(paste0("MOFA 모델 학습 시작 (Factor 개수: ", n_factors, ")"))
  outfile <- file.path(output_dir, paste0("MOFA_factor", n_factors, ".hdf5"))
  mofa_obj <- run_mofa(mofa_obj, use_basilisk = TRUE, outfile = outfile)
  message("MOFA 학습 완료 및 모델 저장:", outfile)
  
  # --------------------------------------
  message("결과 시각화 및 저장 중...")
  
  # factor 별 weigh 저장 상위 20개
  for (i in seq_len(n_factors)) {
    p <- plot_top_weights(mofa_obj, factors = i, nfeatures = 20)
    ggsave(file.path(figure_dir, paste0("factor_weights_Factor", i, ".jpg")),
           p, dpi = 300)
  }
  
  # Factor score 병합
  factors_df <- get_factors(mofa_obj, factors = "all", as.data.frame = TRUE)
  factors_wide <- factors_df %>%
    select(sample, factor, value) %>%
    pivot_wider(names_from = factor, values_from = value)
  
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("sample")
  
  meta_joined <- left_join(meta, factors_wide, by = "sample") %>%
    tibble::column_to_rownames("sample")
  
  seurat_obj@meta.data <- meta_joined
  
  
  # 순서 변경
  meta_joined$Region_group <- factor(meta_joined$Region_group, levels = c(
    "glom_CTRL_pre", "glom_CTRL_post",
    "glom_TCMR_pre", "glom_TCMR_post",
    "tubulointer_CTRL_pre", "tubulointer_CTRL_post",
    "tubulointer_TCMR_pre", "tubulointer_TCMR_post"
  ))
  
  # --------------------------------------
  # 8 Factor boxplot 함수
  # --------------------------------------
  plot_all_factors <- function(meta_joined, n_factors, fill_palette = "Set3", x_var, set_ncol = 2) {
    
    #long_df 로 factor score 변경
    long_df <- meta_joined %>%
      pivot_longer(cols = starts_with("Factor"),
                   names_to = "Factor", values_to = "Score") %>%
      filter(Factor %in% paste0("Factor", 1:n_factors))
    
    
    
    n_groups <- length(unique(long_df[[x_var]]))
    palette_expanded <- colorRampPalette(brewer.pal(8, fill_palette))(n_groups)
    
    
    ggplot(long_df, aes(x = !!sym(x_var), y = Score, fill = !!sym(x_var))) +
      geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.7) +   # 🔹 폭 줄이기
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
  
  
  # --------------------------------------
  #  Factor별 그룹 비교 플롯 저장
  # --------------------------------------
  for (x_var in x_var_list) {
    p <- plot_all_factors(meta_joined, n_factors = n_factors, x_var = x_var)
    out_path <- file.path(figure_dir, paste0("Factor_boxplot_", x_var, ".pdf"))
    ggsave(out_path, p)
    message("저장 완료 ", out_path)
  }
  
  
  p3 <- plot_all_factors(meta_joined, n_factor = 3, x_var = "Region_group", set_ncol = 1)
  out_path <- file.path(figure_dir, paste0("Factor_boxplot_agg3", ".pdf"))
  ggsave(out_path, p3, width = 6, height = 10)
  
  
  # --------------------------------------
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
  
  cor_result_df <- cor_summary_table(meta_joined, n_factors = n_factors, vars = cor_vars)
  write.csv(cor_result_df, file.path(output_dir, "correlation_summary.csv"), row.names = FALSE)
  message("상관분석 결과 저장 완료 → correlation_summary.csv")
  
  # --------------------------------------
  # Seurat 객체 저장
  # --------------------------------------
  saveRDS(seurat_obj, file = file.path(output_dir, "complete_ROI.rds"))
  message("Seurat 객체 저장 완료 → complete_ROI.rds")
  
  message(" 분석 완료: 결과 폴더 → ", output_dir)
}



# --------------------------------------
# 실행할 파라미터 목록 정의
# --------------------------------------

params_list <- list(
  
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "result/Deseq4000Factor4con",
    top_n = 4000,
    n_factors = 4,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  ),
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "result/Deseq4000Factor6con",
    top_n = 4000,
    n_factors = 6,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  ),
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "result/Deseq4000Factor8con",
    top_n = 4000,
    n_factors = 8,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  )
)


# --------------------------------------
# (모든 실행에서 동일하게 적용됨)
# --------------------------------------
x_var_list <- c("Region_group", "ROI_clusters", "ScanLabel" ) # sample = ScanLabel
cor_vars <- c("T_cells")

# --------------------------------------
# 연속 실행
# --------------------------------------
walk(params_list, function(p) {
  cat("\n=====================================\n")
  cat("실행 시작:", p$output_dir, "\n")
  cat("Factor:", p$n_factors, 
      " | Subset:", ifelse(is.null(p$subset_region), "None", p$subset_region), "\n")
  cat("=====================================\n")
  
  run_MOFA_analysis(
    seurat_path = p$seurat_path,
    output_dir = p$output_dir,
    n_factors = p$n_factors,
    subset_region = p$subset_region,
    x_var_list = x_var_list,
    cor_vars = cor_vars
  )
  
  cat("완료:", p$output_dir, "\n\n")
})
