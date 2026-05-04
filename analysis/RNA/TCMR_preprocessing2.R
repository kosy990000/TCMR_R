library(Seurat)
library(DESeq2)
library(MOFA2)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)

# RNA, Protein 모두 tubulointer 만 진행
#--------------------------------RNA-------------------------------------#
# --- DESeq2 variance stabilization ---
# parameter 4000 genes, Deseq기반 정규화
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
#----------------------------------MOFA------------------------------------#
# --- Prepare MOFA input ---
# 파라미터 MOFA_Group은 CTRL과 TCMR = 두 그룹, seed 2025
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

#--------------------------------- protein -----------------------------------#
# --- Normalization ---
# 모든 protein에 대해 진행, CLR기반인 경우 seurat_CLR 사용, 아니면 Deseq기반
if (normalization_method == "CLR") {
  message("Applying CLR normalization...")
  seurat_obj <- NormalizeData(
    seurat_obj,
    assay = "Protein",
    normalization.method = "CLR",
    margin = 2  # cell 에대해 정규화
  )
  expression_data <- GetAssayData(seurat_obj, assay = "Protein", layer = "data")
  meta$sample <- colnames(expression_data)
  vst_mat <- as.matrix(expression_data)   # CLR-normalized matrix
  rm(expression_data); gc()
  
} else if (normalization_method == "DESeq") {
  message("Applying DESeq2 variance stabilizing transformation...")
  expression_data <- GetAssayData(seurat_obj, assay = "Protein", layer = "counts")
  meta$sample <- colnames(expression_data)
  
  # DESeq2 variance stabilization
  dds <- DESeqDataSetFromMatrix(
    countData = expression_data,
    colData = meta,
    design = ~1
  )
  rm(expression_data); gc()
  
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  vst_mat <- assay(vsd)
  rm(dds, vsd); gc()
  
} else {
  stop("normalization_method must be either 'CLR' or 'DESeq'")
}
#----------------------------mofa내부 파라미터는 RNA와 같음 ---------------------#

