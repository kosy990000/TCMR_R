library(Seurat)
library(DESeq2)
library(MOFA2)
library(dplyr)
library(tidyr)
library(tibble)
library(jsonlite)
library(Matrix)
library(ggplot2)




run_MOFA <- function(
    seurat_path,
    top_n = 4000,
    output_dir = "MOFA_Output",
    n_factors = 4,
    subset_region = NULL,
    MOFA_group = "Condition"
) {
  # message("==== RUN MOFA ====")

  # --- Load and subset ---
  seurat_obj <- readRDS(seurat_path)
  message(paste0("Subset 이전 셀 수: ", ncol(seurat_obj)))

  if (!is.null(subset_region)) {
    seurat_obj <- subset(seurat_obj, subset = Region == subset_region)
    message(paste0("Subset 적용: Region == ", subset_region))
  }

  message(paste0("Subset 이후 셀 수: ", ncol(seurat_obj)))

  # --- Expression matrix & metadata ---
  # protein
  expression_data <- GetAssayData(seurat_obj, assay = "Protein", layer = "counts")
  meta <- seurat_obj@meta.data
  meta$sample <- colnames(expression_data)
  
  figure_dir <- file.path(output_dir, "figure")
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(figure_dir, showWarnings = FALSE)
  #protein
  
  #------------------------------------------------------------------------#
  # --- DESeq2 variance stabilization ---
  message("Variance stabilizing transformation...")
  dds <- DESeqDataSetFromMatrix(
    countData = expression_data,
    colData = meta,
    design = ~1
  )
  rm(expression_data); gc()  #
  
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  vst_mat <- assay(vsd)
  rm(dds, vsd); gc()
  
  gc()
#------------------------------------------------------------------------#
  # --- Prepare MOFA input ---
  message("Creating MOFA input dataframe...")
  dt <- as.data.frame(vst_mat) %>%
    tibble::rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "value") %>%
    left_join(meta, by = "sample") %>%
    mutate(view = "Protein", group = .data[[MOFA_group]]) %>%
    select(sample, group, feature, view, value)
  rm(vst_mat); gc()  #

  # --- Create and train MOFA ---
  message("Preparing MOFA model...")
  mofa_obj <- create_mofa_from_df(dt)
  rm(dt); gc()  #

  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$seed <- 2025

  mofa_obj <- prepare_mofa(
    mofa_obj,
    data_options = get_default_data_options(mofa_obj),
    model_options = model_opts,
    training_options = train_opts
  )
  rm(model_opts, train_opts); gc()


  outfile <- file.path(output_dir, paste0("MOFA_factor", n_factors, ".hdf5"))
  message("Running MOFA...")
  mofa_obj <- run_mofa(mofa_obj, use_basilisk = TRUE, outfile = outfile)
  message("MOFA model trained and saved: ", outfile)

  #confirm mofa obj
  plot_data_overview(mofa_obj)

#------------- save Factor weight and Factor figure ----------------------------#

  figure_dir <- file.path(output_dir, "figure")
  factor_weight_dir <- file.path(output_dir, "factor-weight")
  dir.create(figure_dir, showWarnings = FALSE)
  dir.create(factor_weight_dir, showWarnings = FALSE)

  view_name <- "Protein"
  weights_df <- get_weights(mofa_obj, view = view_name, as.data.frame = TRUE)

#------------------ save factor weight lst and factor weight figure ------------#
  for (factor_name in unique(weights_df$factor)) {
    sub_df <- filter(weights_df, factor == factor_name)
    write.csv(sub_df,
              file = file.path(factor_weight_dir, paste0("weights_", factor_name, ".csv")),
              row.names = FALSE)
  }

  # --- Factor weights plot ---
  for (i in seq_len(n_factors)) {
    p <- plot_top_weights(mofa_obj, factors = i, nfeatures = 20)
    ggsave(file.path(figure_dir, paste0("factor_weights_Factor", i, ".pdf")), p, dpi = 300)
  }



#------------------------------------------------------------------------------#
  # --- Attach factor scores to Seurat ---
  message("Extracting and merging factor scores...")
  factors_df <- get_factors(mofa_obj, factors = "all", as.data.frame = TRUE)

  factors_wide <- factors_df %>%
    dplyr::select(sample, factor, value) %>%
    tidyr::pivot_wider(names_from = factor, values_from = value)
  rm(factors_df); gc()

  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("sample")

  meta_joined <- dplyr::left_join(meta, factors_wide, by = "sample") %>%
    tibble::column_to_rownames("sample")
  rm(factors_wide); gc()

  meta_joined <- meta_joined[match(colnames(seurat_obj), rownames(meta_joined)), ]
  seurat_obj@meta.data <- meta_joined
  rm(meta, meta_joined); gc()

  # --- Save Seurat + MOFA ---
  saveRDS(seurat_obj, file.path(output_dir, "seurat_obj_raw.rds"))
  saveRDS(mofa_obj, file.path(output_dir, "mofa_model_raw.rds"))
  rm(seurat_obj, mofa_obj); gc()

  message("Core analysis complete. Ready for visualization.")
}
params_list <- list(
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "data/results/DeseqProNorfactor4",
    top_n = 4000,
    n_factors = 4,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  ),
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "data/results/DeseqProNorfactor6",
    top_n = 4000,
    n_factors = 6,
    subset_region = "tubulointer",
    MOFA_group = "Condition"
  ),
  list(
    seurat_path = "raw-data/raw_ROI_data.rds",
    output_dir = "data/results/DeseqProNorfactor8",
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

  run_MOFA(
    seurat_path = p$seurat_path,
    output_dir = p$output_dir,
    top_n = p$top_n,
    n_factors = p$n_factors,
    subset_region = p$subset_region,
    MOFA_group = p$MOFA_group
  )
}



#for (params in config$factor_analysis$params_list) {
#  run_MOFA(
#    seurat_path = params$seurat_path,
#    output_dir = params$output_dir,
#    top_n = params$top_n,
#    n_factors = params$n_factors,
#    subset_region = params$subset_region,
#    MOFA_group = params$MOFA_group
#  )
#}
