#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 필요한 라이브러리 임포트
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.colors import LinearSegmentedColormap

# 한글 폰트 설정 (필요한 경우)
plt.rcParams['font.family'] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False

# 설정: 기본 디렉토리
base_dir = '/home/swr0460/daeseong'
data_dir = os.path.join(base_dir, 'data')
figures_dir = os.path.join(base_dir, 'figures')
output_dir = os.path.join(base_dir, 'output')

# 디렉토리 생성
os.makedirs(figures_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Scanpy의 기본 figure 경로를 설정
sc.settings.figdir = figures_dir

# 재현성을 위한 랜덤 시드 설정
np.random.seed(42)

def run_downstream_analysis():
    """
    QC 이후의 다운스트림 분석을 수행합니다.
    - 정규화 및 스케일링
    - 고변동 유전자 선택
    - 차원 축소 (PCA, UMAP, t-SNE)
    - 클러스터링
    - 마커 유전자 분석
    """
    
    print("QC 완료된 데이터 로딩 중...")
    # QC 완료된 데이터 로드
    try:
        adata = sc.read_h5ad(os.path.join(output_dir, 'filtered_data.h5ad'))
        print(f"QC 완료된 데이터 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    except FileNotFoundError:
        print("QC 완료된 데이터 파일을 찾을 수 없습니다. 원본 데이터를 로드합니다.")
        # 원본 데이터 로드 및 기본 QC 적용
        adata = sc.read_10x_mtx(
            os.path.join(data_dir, 'filtered_feature_bc_matrix/'),
            var_names='gene_symbols',
            cache=True
        )
        
        # 기본 QC 필터링
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        # QC 필터링 적용
        adata = adata[(adata.obs.n_genes_by_counts >= 200) &
                      (adata.obs.n_genes_by_counts <= 2500) &
                      (adata.obs.pct_counts_mt <= 5), :].copy()
        
        print(f"원본 데이터 로드 및 기본 QC 적용 완료: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    
    # 1. 정규화 및 로그 변환
    print("데이터 정규화 및 로그 변환 중...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 2. 고변동 유전자 선택
    print("고변동 유전자 선택 중...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # 고변동 유전자 시각화
    sc.pl.highly_variable_genes(adata, save='highly_variable_genes.pdf', show=False)
    
    print(f"선택된 고변동 유전자 수: {np.sum(adata.var.highly_variable)}")
    
    # 고변동 유전자만 선택 (선택사항)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    
    # 3. 데이터 스케일링
    print("데이터 스케일링 중...")
    sc.pp.scale(adata_hvg, max_value=10)
    
    # 4. 주성분 분석 (PCA)
    print("PCA 수행 중...")
    sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
    
    # PCA 결과 시각화
    sc.pl.pca(adata_hvg, color='n_genes_by_counts', save='pca.pdf', show=False)
    sc.pl.pca_variance_ratio(adata_hvg, n_pcs=50, save='pca_variance_ratio.pdf', show=False)
    
    # 5. 이웃 그래프 계산
    print("이웃 그래프 계산 중...")
    sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=30)
    
    # 6. UMAP 차원 축소
    print("UMAP 차원 축소 수행 중...")
    sc.tl.umap(adata_hvg)
    
    # 7. t-SNE 차원 축소 (선택사항)
    print("t-SNE 차원 축소 수행 중...")
    sc.tl.tsne(adata_hvg, n_pcs=30)
    
    # 8. 클러스터링
    print("클러스터링 수행 중...")
    # Leiden 알고리즘을 사용한 클러스터링
    sc.tl.leiden(adata_hvg, resolution=0.5)
    
    # 다양한 resolution 값으로 클러스터링 시도 (선택사항)
    resolutions = [0.3, 0.5, 0.8, 1.0]
    for res in resolutions:
        sc.tl.leiden(adata_hvg, resolution=res, key_added=f'leiden_res{res}')
    
    # 9. 차원 축소 결과 시각화
    print("차원 축소 결과 시각화 중...")
    
    # UMAP 시각화 - 클러스터링 결과
    sc.pl.umap(adata_hvg, color=['leiden'], palette='tab20', save='umap_clusters.pdf', show=False)
    
    # UMAP 시각화 - QC 지표
    sc.pl.umap(adata_hvg, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
              save='umap_qc_metrics.pdf', show=False)
    
    # t-SNE 시각화 - 클러스터링 결과
    sc.pl.tsne(adata_hvg, color=['leiden'], palette='tab20', save='tsne_clusters.pdf', show=False)
    
    # 다양한 resolution 값으로 클러스터링 결과 비교
    sc.pl.umap(adata_hvg, color=[f'leiden_res{res}' for res in resolutions], 
              wspace=0.5, save='umap_resolutions.pdf', show=False)
    
    # 10. 마커 유전자 분석
    print("마커 유전자 분석 중...")
    sc.tl.rank_genes_groups(adata_hvg, 'leiden', method='wilcoxon')
    
    # 마커 유전자 결과 시각화
    sc.pl.rank_genes_groups(adata_hvg, n_genes=25, sharey=False, save='marker_genes.pdf', show=False)
    
    # 상위 마커 유전자 히트맵
    sc.pl.rank_genes_groups_heatmap(adata_hvg, n_genes=10, groupby='leiden', 
                                   save='marker_genes_heatmap.pdf', show=False)
    
    # 상위 마커 유전자 점도표
    sc.pl.rank_genes_groups_dotplot(adata_hvg, n_genes=5, groupby='leiden', 
                                   save='marker_genes_dotplot.pdf', show=False)
    
    # 상위 마커 유전자 바이올린 플롯
    sc.pl.rank_genes_groups_violin(adata_hvg, n_genes=5, 
                                  save='marker_genes_violin.pdf', show=False)
    
    # 11. 마커 유전자 테이블 저장
    marker_genes = pd.DataFrame()
    for i in range(len(adata_hvg.obs['leiden'].cat.categories)):
        markers_df = sc.get.rank_genes_groups_df(adata_hvg, group=str(i))
        markers_df['cluster'] = i
        marker_genes = pd.concat([marker_genes, markers_df])
    
    marker_genes.to_csv(os.path.join(output_dir, 'marker_genes.csv'), index=False)
    
    # 12. 클러스터별 세포 수 계산 및 시각화
    cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
    cluster_counts.name = '세포 수'
    cluster_counts.index.name = '클러스터'
    
    cluster_counts.to_csv(os.path.join(output_dir, 'cluster_cell_counts.csv'))
    
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(x=cluster_counts.index, y=cluster_counts.values)
    ax.set_title('클러스터별 세포 수')
    ax.set_xlabel('클러스터')
    ax.set_ylabel('세포 수')
    
    # 막대 위에 숫자 표시
    for i, count in enumerate(cluster_counts.values):
        ax.text(i, count + 5, str(count), ha='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'cluster_cell_counts.pdf'))
    plt.savefig(os.path.join(figures_dir, 'cluster_cell_counts.png'), dpi=300)
    plt.close()
    
    # 13. 클러스터 간 상관관계 분석 (선택사항)
    print("클러스터 간 상관관계 분석 중...")
    
    # 클러스터별 평균 발현량 계산
    sc.tl.dendrogram(adata_hvg, groupby='leiden')
    
    # 클러스터 간 상관관계 히트맵
    sc.pl.dendrogram(adata_hvg, groupby='leiden', save='cluster_dendrogram.pdf', show=False)
    
    # 14. 결과 저장
    print("분석 결과 저장 중...")
    adata_hvg.write(os.path.join(output_dir, 'downstream_results.h5ad'))
    
    # 원본 데이터에 클러스터 정보 추가 (선택사항)
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.obsm['X_tsne'] = adata_hvg.obsm['X_tsne']
    
    # 전체 데이터 저장
    adata.write(os.path.join(output_dir, 'complete_results.h5ad'))
    
    print(f"\n다운스트림 분석이 완료되었습니다!")
    print(f"- 시각화 결과: {figures_dir}")
    print(f"- 분석 결과 데이터: {os.path.join(output_dir, 'downstream_results.h5ad')}")
    print(f"- 마커 유전자 목록: {os.path.join(output_dir, 'marker_genes.csv')}")
    print(f"- 클러스터별 세포 수: {os.path.join(output_dir, 'cluster_cell_counts.csv')}")

if __name__ == "__main__":
    run_downstream_analysis() 