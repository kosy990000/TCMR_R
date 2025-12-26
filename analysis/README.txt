================================================================================
MOFA Analysis - Improved Version
================================================================================

프로젝트 구조
================================================================================

analysis_improved/
├── config/                          # 설정 파일 디렉토리
│   ├── config.json                  # MOFA 실행 설정 (RNA 분석용)
│   └── config_visualization.json    # Heatmap & Cell Proportion 시각화 설정
│
├── RNA/                             # RNA 분석 스크립트
│   ├── runMOFA.R                    # MOFA 모델 실행
│   ├── figureMOFA.R                 # MOFA 결과 시각화
│   └── heatmapaAndCellProportion.R  # Heatmap & Cell Proportion 생성
│
├── Protein/                         # Protein 분석 스크립트
│   ├── runMOFAProtein.R             # MOFA 모델 실행 (Protein)
│   ├── runMOFAProteinDeseq.R        # MOFA + DESeq 실행 (Protein)
│   ├── figureMOFAProtein.R          # MOFA 결과 시각화 (Protein)
│   └── heatmapaAndCellProportionProtein.R  # Heatmap & Cell Proportion (Protein)
│
├── figure.R                         # 기타 figure 생성
└── multiRunMofa.R                   # 여러 MOFA 모델 동시 실행


파일 설명
================================================================================

[RNA 폴더]
----------
1. runMOFA.R
   - MOFA 모델 학습 및 실행
   - config/config.json 사용
   - 여러 factor 수(4, 6, 8)로 반복 실행
   - 출력: seurat_obj_raw.rds, mofa_model_raw.rds

2. figureMOFA.R
   - MOFA 결과 시각화
   - config/config.json 사용
   - 생성 결과:
     * Factor weights (CSV & 시각화)
     * Factor score boxplot
     * Factor-변수 상관관계 scatter plot
     * Factor-celltype 상관관계 heatmap

3. heatmapaAndCellProportion.R
   - Cell proportion 및 Gene expression heatmap 생성
   - config/config_visualization.json 사용
   - 한 번만 실행 (반복 불필요)
   - 생성 결과:
     * Cell proportion bar plot
     * Gene expression heatmap (top/bottom 20 genes)


[Protein 폴더]
--------------
- RNA 폴더와 동일한 구조, Protein assay 사용


[Config 폴더]
-------------
1. config.json
   - MOFA 분석 전체 설정
   - 경로, 파라미터, 실행 리스트 포함
   - 여러 factor 수로 반복 실행 설정

2. config_visualization.json
   - Heatmap & Cell Proportion 시각화 설정
   - 단일 실행용
   - seurat_path: raw-data에서 로드
   - weight_base_dir: result 폴더에서 factor weight 로드


실행 방법
================================================================================

[1단계] MOFA 모델 실행
Rscript analysis_improved/RNA/runMOFA.R

[2단계] MOFA 결과 시각화
Rscript analysis_improved/RNA/figureMOFA.R

[3단계] Heatmap & Cell Proportion 생성
Rscript analysis_improved/RNA/heatmapaAndCellProportion.R


주요 변경 사항
================================================================================

1. 폴더 구조 재구성
   - RNA/Protein 분리
   - Config 파일 별도 관리

2. Cell Proportion & Heatmap 분리
   - figureMOFA.R: Factor-celltype correlation heatmap 포함
   - heatmapaAndCellProportion.R: Cell proportion bar plot & Gene expression heatmap만 포함

3. Config 기반 실행
   - 하드코딩 제거
   - JSON 설정 파일로 파라미터 관리

4. Visualization 최적화
   - 중복 제거 (한 번만 실행)
   - Raw data에서 seurat 로드
   - Result에서 weight 로드


출력 디렉토리 구조
================================================================================

improve_result/
├── Deseq4000Factor4con/
│   ├── seurat_obj_raw.rds
│   ├── mofa_model_raw.rds
│   ├── factor-weight/
│   │   └── weights_Factor*.csv
│   └── figure/
│       ├── factor_weights_Factor*.pdf
│       ├── Factor_boxplot_*.pdf
│       ├── Factor*_*.pdf (scatter plots)
│       └── factor_celltype_correlation_heatmap.pdf
│
├── Deseq4000Factor6con/
│   └── (동일 구조)
│
├── Deseq4000Factor8con/
│   └── (동일 구조)
│
└── visualization/
    └── figure/
        ├── cell_prop.pdf
        └── heatmap_top_bottom20.pdf


참고 사항
================================================================================

- RNA 분석이 기본 설정
- Protein 분석 시 Protein 폴더의 스크립트 사용
- Config 파일 수정으로 파라미터 변경 가능
- Factor 수는 config.json의 run_list에서 조정

================================================================================
