library(Seurat)
library(dplyr)
library(tidyr)
# 1. 파일 로드
# (이미 로드하셨다면 이 부분은 건너뛰세요)
obj_raw <- readRDS("raw-data/raw_ROI_data.rds")
obj_cluster <- readRDS("raw-data/ROI_cluster.rds")

# Check dimensions (Features x Cells)
dim(obj_check)

# Check available assays
Assays(obj_check)

# See the active assay
DefaultAssay(obj_check)
# ==============================================================================
# 함수: 데이터 상태 정밀 진단 (SCT 유무 + 정수 여부 확인)
# ==============================================================================
check_data_status <- function(obj, name) {
  cat("\n======================================================\n")
  cat(paste0("  파일 진단: ", name, "\n"))
  cat("======================================================\n")
  
  # 1. Assay 확인
  assays <- Assays(obj)
  cat("1. 보유 Assay:", paste(assays, collapse = ", "), "\n")
  
  # 2. SCT 여부
  if ("SCT" %in% assays) {
    cat("2. SCTransform 적용됨: [YES] (주의: DefaultAssay가 SCT인지 확인 필요)\n")
  } else {
    cat("2. SCTransform 적용됨: [NO] (Clean Data일 가능성 높음)\n")
  }
  
  # 3. RNA Counts 정밀 분석
  if ("RNA" %in% assays) {
    # Seurat 버전에 따른 데이터 추출
    counts <- tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
    
    # 샘플링 확인 (속도를 위해)
    sample_vals <- as.vector(counts[1:100, 1:10]) # 일부만 추출
    sample_vals <- sample_vals[sample_vals > 0]   # 0 제외
    sample_vals <- head(sample_vals, 50)          # 50개만 확인
    
    is_integer <- all(sample_vals %% 1 == 0)
    max_val <- max(sample_vals, na.rm = TRUE)
    
    cat("3. RNA Counts 데이터 타입:\n")
    if (is_integer) {
      cat("   ✅ [PASS] 정수(Integer)입니다. (Raw Data)\n")
    } else {
      cat("   ❌ [FAIL] 소수점(Float)입니다! (이미 정규화됨)\n")
    }
    cat(paste0("   (참고: 샘플 최댓값 = ", max_val, ")\n"))
  }
}

# ==============================================================================
# 실행: 두 파일 상태 확인
# ==============================================================================
check_data_status(obj_raw, "FILE 1 (Raw Data)")
check_data_status(obj_cluster, "FILE 2 (Cluster Data)")


# ==============================================================================
# [핵심] 두 파일의 Raw Count 값 완전 일치 여부 검증
# ==============================================================================
cat("\n======================================================\n")
cat("  [최종 검증] 두 파일의 데이터 값이 동일한가?\n")
cat("======================================================\n")

# 1. 공통된 세포와 유전자 찾기 (교집합)
common_cells <- intersect(colnames(obj_raw), colnames(obj_cluster))
common_genes <- intersect(rownames(obj_raw), rownames(obj_cluster))

cat(paste0("- 공통 세포 수: ", length(common_cells), "\n"))
cat(paste0("- 공통 유전자 수: ", length(common_genes), "\n"))

if (length(common_cells) == 0) {
  stop("🚨 두 파일 간에 공통된 세포 이름이 하나도 없습니다! 비교 불가능.")
}

# 2. 공통 영역의 데이터 추출
# (주의: 행/열 순서를 강제로 맞춰야 정확한 비교가 가능)
mat_raw <- GetAssayData(obj_raw, assay = "RNA", layer = "counts")[common_genes, common_cells]
mat_cluster <- GetAssayData(obj_cluster, assay = "RNA", layer = "counts")[common_genes, common_cells]

# 3. 값 비교 (identical 함수 사용)
if (identical(mat_raw, mat_cluster)) {
  cat("\n🎉 [결과: 완벽 일치] 🎉\n")
  cat("두 파일의 RNA counts 값은 토씨 하나 안 틀리고 똑같습니다.\n")
  cat("-> 결론: 'ROI_cluster.rds'를 써도 데이터 값 자체는 원본과 동일합니다.\n")
  cat("-> 다만, 세포 목록(Subset)이 다를 수 있으니 그 점만 주의하세요.\n")
} else {
  cat("\n🚨 [결과: 불일치] 🚨\n")
  cat("공통된 세포임에도 불구하고 값이 다릅니다!\n")
  
  # 얼마나 다른지 확인
  diff_vals <- sum(mat_raw != mat_cluster)
  cat(paste0("-> 다른 값의 개수: ", diff_vals, "\n"))
  cat("-> 원인: ROI_cluster 파일 저장 시, 정규화된 값으로 덮어씌워졌거나 가공되었을 수 있습니다.\n")
  cat("-> 조치: 무조건 'raw_ROI_data.rds'를 사용하세요.\n")
}

cat("=== 전체 유전자(행) 개수 비교 ===\n")
cat("Raw 파일 유전자 수:    ", nrow(obj_raw), "\n")
cat("Cluster 파일 유전자 수:", nrow(obj_cluster), "\n")

cat("\n=== nCount_RNA와 실제 합계 비교 ===\n")
# Cluster 파일의 첫 번째 세포 실제 합계 계산
real_sum <- sum(GetAssayData(obj_cluster, assay="RNA", layer="counts")[,1])
meta_sum <- obj_cluster@meta.data$nCount_RNA[1]

cat("Cluster 파일의 첫 세포 - 실제 Matrix 합계:", real_sum, "\n")
cat("Cluster 파일의 첫 세포 - 적혀있는 Metadata:", meta_sum, "\n")

if (real_sum != meta_sum) {
  cat("👉 결론: Metadata(nCount_RNA)가 현재 데이터와 싱크가 안 맞습니다.\n")
  cat("   (과거의 기록이거나, 유전자가 필터링된 후 업데이트되지 않은 것입니다.)\n")
} else {
  cat("👉 결론: 합계는 맞습니다. 그렇다면 Raw 파일과 유전자 목록 자체가 다른 것입니다.\n")
}



compare_EVERYTHING <- function(path_raw, path_cluster) {
  
  message(">>> 📂 데이터 로딩 중... (시간이 조금 걸릴 수 있습니다)")
  obj_raw <- readRDS(path_raw)
  obj_cluster <- readRDS(path_cluster)
  
  # 메타데이터 추출 (행 이름을 컬럼으로 변환하여 안전하게 비교)
  meta_raw <- obj_raw@meta.data %>% tibble::rownames_to_column("CellID")
  meta_cluster <- obj_cluster@meta.data %>% tibble::rownames_to_column("CellID")
  
  cat("\n======================================================\n")
  cat("🔍 [1단계] 세포(행) 구성 비교: 누가 빠졌나?\n")
  cat("======================================================\n")
  
  cells_raw <- meta_raw$CellID
  cells_cluster <- meta_cluster$CellID
  
  common_cells <- intersect(cells_raw, cells_cluster)
  only_raw <- setdiff(cells_raw, cells_cluster)     # Raw에만 있는 세포 (삭제된 애들)
  only_cluster <- setdiff(cells_cluster, cells_raw) # Cluster에만 있는 세포 (이상한 경우)
  
  cat(paste0("1. Raw 파일 총 세포 수    : ", length(cells_raw), "\n"))
  cat(paste0("2. Cluster 파일 총 세포 수: ", length(cells_cluster), "\n"))
  cat(paste0("3. 교집합(공통) 세포 수   : ", length(common_cells), "\n"))
  
  if (length(only_raw) > 0) {
    cat(paste0("\n🚨 [차이 발견] Raw 데이터에만 존재하는 세포가 ", length(only_raw), "개 있습니다.\n"))
    cat("   -> (해석) QC 과정에서 Low Quality 세포들이 제거된 것으로 보입니다.\n")
    cat("   -> 이 '제거된 세포들' 때문에 DESeq2 정규화 배경이 달라져 결과가 변한 것입니다.\n")
  } else {
    cat("\n✅ 세포 목록이 완벽하게 동일합니다. (QC 필터링 차이 없음)\n")
  }
  
  if (length(only_cluster) > 0) {
    cat(paste0("\n⚠️ [주의] Cluster 데이터에 갑자기 생긴 세포가 ", length(only_cluster), "개 있습니다.\n"))
    cat("   -> 세포 이름이 변경되었거나, 다른 데이터가 합쳐졌을 수 있습니다.\n")
  }
  
  cat("\n======================================================\n")
  cat("🔍 [2단계] 타겟(컬럼) 값 전수 비교\n")
  cat("   (공통된 ", length(common_cells), "개 세포에 대해서만 비교합니다)\n")
  cat("======================================================\n")
  
  # 공통 세포 기준으로 데이터 정렬 (순서 맞추기)
  df_raw <- meta_raw %>% filter(CellID %in% common_cells) %>% arrange(CellID)
  df_cluster <- meta_cluster %>% filter(CellID %in% common_cells) %>% arrange(CellID)
  
  # 공통 컬럼 찾기
  cols_raw <- colnames(df_raw)
  cols_cluster <- colnames(df_cluster)
  common_cols <- intersect(cols_raw, cols_cluster)
  
  # 컬럼 비교 루프
  diff_report <- list()
  
  for (col in common_cols) {
    if (col == "CellID") next
    
    # 값 추출 (문자열로 변환하여 비교)
    val_raw <- as.character(df_raw[[col]])
    val_cluster <- as.character(df_cluster[[col]])
    
    # NA 처리 (NA는 문자열 "NA"로 치환해서 비교)
    val_raw[is.na(val_raw)] <- "NA_VALUE"
    val_cluster[is.na(val_cluster)] <- "NA_VALUE"
    
    if (all(val_raw == val_cluster)) {
      # cat(paste0("✅ [일치] ", col, "\n")) # 너무 많으면 생략
    } else {
      # 다른 개수 세기
      diff_count <- sum(val_raw != val_cluster)
      cat(paste0("❌ [불일치] 컬럼명: '", col, "' -> ", diff_count, "개의 세포에서 값이 다름!\n"))
      
      # 예시 출력
      diff_idx <- which(val_raw != val_cluster)[1] # 첫 번째 다른 곳
      cat(paste0("    예시 (Cell: ", df_raw$CellID[diff_idx], ")\n"))
      cat(paste0("      - Raw 파일 값    : ", val_raw[diff_idx], "\n"))
      cat(paste0("      - Cluster 파일 값: ", val_cluster[diff_idx], "\n\n"))
    }
  }
  
  # 한쪽에만 있는 컬럼 확인
  raw_only_cols <- setdiff(cols_raw, cols_cluster)
  cluster_only_cols <- setdiff(cols_cluster, cols_raw)
  
  if(length(raw_only_cols) > 0) {
    cat("\n[참고] Raw 파일에만 있는 컬럼:\n")
    print(raw_only_cols)
  }
  if(length(cluster_only_cols) > 0) {
    cat("\n[참고] Cluster 파일에만 있는 컬럼 (새로 추가된 분석 결과 등):\n")
    print(cluster_only_cols)
  }
  
  message("\n>>> 비교 완료.")
}

# -----------------------------------------------------------
# 🚀 실행 (경로만 본인 것으로 수정해서 돌리세요)
# -----------------------------------------------------------
compare_EVERYTHING(
  path_raw = "raw-data/raw_ROI_data.rds", 
  path_cluster = "raw-data/ROI_cluster.rds"
)


