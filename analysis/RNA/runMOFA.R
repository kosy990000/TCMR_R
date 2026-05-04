library(Seurat)
library(DESeq2)
library(MOFA2)
library(dplyr)
library(tidyr)
library(tibble)
library(jsonlite)
library(Matrix)
library(ggplot2)

# Load configuration
config <- fromJSON("analysis/config/config_no_con_new.json", simplifyVector = FALSE)


run_MOFA <- function(
    seurat_path,
    top_n = 4000,
    output_dir = "MOFA_Output",
    n_factors = 4,
    subset_region = NULL,
    MOFA_group = "Condition"
) {
  message("==== RUN MOFA ====")

  # --- Load and subset ---
  seurat_obj <- readRDS(seurat_path)
  message(paste0("Subset 이전 셀 수: ", ncol(seurat_obj)))

  if (!is.null(subset_region)) {
    seurat_obj <- subset(seurat_obj, subset = Region == subset_region)
    message(paste0("Subset 적용: Region == ", subset_region))
  }

  message(paste0("Subset 이후 셀 수: ", ncol(seurat_obj)))

  # --- Expression matrix & metadata ---
  expression_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  # protein
  #expression_data <- GetAssayData(seurat_obj, assay = "Protein", slot = "counts")
  meta <- seurat_obj@meta.data
  meta$sample <- colnames(expression_data)

  dir.create(output_dir, showWarnings = FALSE)

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

  # --- Feature selection ---
  message("Selecting top variable genes...")
  topN <- order(rowVars(vst_mat), decreasing = TRUE)[1:top_n]
  vst_mat <- vst_mat[topN, ]
  gc()
#------------------------------------------------------------------------#
  # --- Prepare MOFA input ---
  message("Creating MOFA input dataframe...")
  dt <- as.data.frame(vst_mat) %>%
    tibble::rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "value") %>%
    left_join(meta, by = "sample") %>%
    mutate(view = "RNA", group = .data[[MOFA_group]]) %>%
    select(sample, group, feature, view, value)
  rm(vst_mat); gc()  #

  # --- Create and train MOFA ---
  message("Preparing MOFA model...")
  mofa_obj <- create_mofa_from_df(dt)
  rm(dt); gc()  #

  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$seed <- config$fixed_params$seed

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

  message("MOFA model training complete. Model saved to: ", output_dir)
}

# Build params_list from config
params_list <- lapply(config$run_list, function(run) {
  list(
    seurat_path = config$paths$seurat_data,
    output_dir = file.path(config$paths$base_output_dir, run$name),
    top_n = config$fixed_params$top_n,
    n_factors = run$n_factors,
    subset_region = config$fixed_params$subset_region,
    MOFA_group = config$fixed_params$MOFA_group
  )
})

# Run MOFA for each configuration
for (i in seq_along(params_list)) {
  cat("\n============================\n")
  cat("Running MOFA job", i, "/", length(params_list), "\n")
  cat("Output:", params_list[[i]]$output_dir, "\n")
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

