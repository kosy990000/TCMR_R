library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(Matrix) # 희소 행렬 변환용

# 1. Load H5ad 
root_dir <- getwd()
data_dir <- file.path(root_dir, "raw-data")
data_name <- "cell_proportion_ROI_no_human.h5ad"
data_path <- file.path(data_dir, data_name)

sce <- readH5AD(data_path)

# [최적화] 일반 행렬을 희소 행렬(Sparse Matrix)로 변환 (메모리 절약)
# 'X'가 없으면 'counts' 등을 확인해야 할 수도 있습니다.
if ("X" %in% assayNames(sce)) {
  counts_mat <- assay(sce, "X")
} else {
  stop("Error: 'X' assay가 없습니다. assayNames(sce)를 확인하세요.")
}

# Seurat은 sparse matrix를 좋아합니다.
if (!is(counts_mat, "sparseMatrix")) {
  counts_mat <- as(counts_mat, "CsparseMatrix")
}

meta <- as.data.frame(colData(sce))

# [검증 1] 세포 매칭 (작성하신 부분 그대로 유지 - 아주 좋습니다)
if (ncol(counts_mat) != nrow(meta)) {
  stop("Error: 세포 개수가 일치하지 않습니다!")
}

if (!identical(colnames(counts_mat), rownames(meta))) {
  warning("세포 이름 순서가 다릅니다. 메타데이터 순서에 맞춰 정렬합니다.")
  counts_mat <- counts_mat[, rownames(meta)]
} else {
  print("Pass: RNA 데이터와 메타데이터의 매핑이 정확합니다.")
}

# 2. Seurat Object 생성 (RNA)
seurat_obj <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta,
  min.cells = 0,
  min.features = 0
)

# 3. Protein 데이터 처리 (핵심 수정 구간)
if ("protein" %in% reducedDimNames(sce)) {
  
  print("Found 'protein' in reducedDims...")
  
  # (1) 데이터 추출
  raw_protein <- reducedDim(sce, "protein")
  
  # (2) 전치 (Cells x Proteins -> Proteins x Cells)
  protein_data <- t(raw_protein) 
  
  # (3) [안전장치] 단백질 이름 찾기 (우선순위 로직 적용)
  # Priority 1: 이미 컬럼 이름이 붙어있는 경우 (가장 이상적)
  if (!is.null(colnames(raw_protein))) {
    rownames(protein_data) <- colnames(raw_protein)
    print(" -> 단백질 이름을 reducedDim 컬럼명에서 가져왔습니다.")
    
    # Priority 2: metadata(uns)에 따로 저장된 경우 (작성하신 코드)
  } else if (!is.null(metadata(sce)$protein_features)) {
    protein_names <- metadata(sce)$protein_features
    if(length(protein_names) == nrow(protein_data)){
      rownames(protein_data) <- protein_names
      print(" -> 단백질 이름을 metadata에서 찾아 매핑했습니다.")
    } else {
      warning("경고: metadata의 단백질 이름 개수와 데이터 행 개수가 맞지 않습니다!")
    }
    
    # Priority 3: 이름이 아예 없는 경우 (임시 이름 부여 - 에러 방지)
  } else {
    warning("경고: 단백질 이름을 찾을 수 없습니다! 'Prot_1', 'Prot_2'...로 임의 지정합니다.")
    rownames(protein_data) <- paste0("Prot_", 1:nrow(protein_data))
  }
  
  # (4) Assay 추가
  seurat_obj[["Protein"]] <- CreateAssayObject(counts = protein_data)
  print("성공: 단백질 데이터가 'Protein' Assay로 추가되었습니다.")
  
} else {
  print("알림: 이 데이터에는 'protein' 정보가 reducedDim에 없습니다.")
}

# 4. 메타데이터 가공
seurat_obj$Group <- paste(seurat_obj$Condition, seurat_obj$Timepoint, sep = "_")
seurat_obj$Region_group <- paste(seurat_obj$Region, seurat_obj$Group, sep = "_")

# 5. 저장
saveRDS(seurat_obj, file.path(data_dir, "raw_ROI_data_no_human.rds"))

# 6. 최종 검증 (Seurat V5 호환성 고려)
print("=== 검증 시작 ===")

# RNA 확인 (GetAssayData 권장)
rna_check <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[1:5, 1:5]
print("RNA Counts (Top 5x5):")
print(rna_check)

# Protein 확인
if ("Protein" %in% names(seurat_obj@assays)) {
  prot_check <- GetAssayData(seurat_obj, assay = "Protein", layer = "counts")[1:5, 1:5]
  print("Protein Counts (Top 5x5):")
  print(prot_check)
}

