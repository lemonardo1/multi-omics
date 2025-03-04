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

def run_trajectory_analysis():
    """
    세포 궤적 분석을 수행합니다.
    - 의사시간(pseudotime) 분석
    - 세포 분화 경로 시각화
    - 시간에 따른 유전자 발현 변화 분석
    """
    
    print("주석 완료된 데이터 로딩 중...")
    # 주석 완료된 데이터 로드
    try:
        adata = sc.read_h5ad(os.path.join(output_dir, 'annotated_results.h5ad'))
        print(f"주석 완료된 데이터 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    except FileNotFoundError:
        try:
            adata = sc.read_h5ad(os.path.join(output_dir, 'downstream_results.h5ad'))
            print(f"클러스터링 완료된 데이터 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
        except FileNotFoundError:
            try:
                adata = sc.read_h5ad(os.path.join(output_dir, 'complete_results.h5ad'))
                print(f"전체 데이터 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
            except FileNotFoundError:
                print("분석 완료된 데이터 파일을 찾을 수 없습니다. 이전 분석 스크립트를 먼저 실행해주세요.")
                return
    
    # 1. 데이터 준비
    print("궤적 분석을 위한 데이터 준비 중...")
    
    # 이미 정규화, 스케일링, 차원 축소가 완료되었는지 확인
    if 'X_pca' not in adata.obsm:
        print("PCA 결과가 없습니다. 기본 전처리를 수행합니다.")
        # 기본 전처리 수행
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_hvg = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata_hvg, max_value=10)
        sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
        sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=30)
        sc.tl.umap(adata_hvg)
        adata = adata_hvg
    
    # 2. 확산 의사시간(Diffusion Pseudotime) 계산
    print("확산 의사시간(DPT) 계산 중...")
    
    # 시작점(root cell) 선택 - 여기서는 첫 번째 클러스터의 중심을 사용
    # 실제 분석에서는 생물학적 지식을 바탕으로 시작점을 선택해야 함
    if 'leiden' in adata.obs:
        # 첫 번째 클러스터의 세포들 중에서 선택
        root_cells = adata.obs[adata.obs['leiden'] == '0'].index
        if len(root_cells) > 0:
            root_cell = root_cells[0]
            # 루트 세포의 인덱스 찾기
            root_index = np.where(adata.obs.index == root_cell)[0][0]
            print(f"시작점으로 클러스터 0의 세포를 선택: {root_cell}")
        else:
            # 임의의 세포 선택
            root_cell = adata.obs.index[0]
            root_index = 0
            print(f"시작점으로 임의의 세포를 선택: {root_cell}")
    else:
        # 임의의 세포 선택
        root_cell = adata.obs.index[0]
        root_index = 0
        print(f"시작점으로 임의의 세포를 선택: {root_cell}")
    
    # 확산 맵 계산
    sc.tl.diffmap(adata)
    
    # 루트 세포 인덱스 설정
    adata.uns['iroot'] = root_index
    
    # 확산 의사시간 계산
    sc.tl.dpt(adata, n_dcs=10, n_branchings=0)
    
    # 3. 의사시간 시각화
    print("의사시간 시각화 중...")
    
    # UMAP에 의사시간 시각화
    sc.pl.umap(adata, color=['dpt_pseudotime'], cmap='viridis', 
              save='pseudotime.pdf', show=False)
    
    # 4. PAGA (Partition-based Graph Abstraction) 분석
    print("PAGA 분석 수행 중...")
    
    # PAGA 계산
    if 'leiden' in adata.obs:
        sc.tl.paga(adata, groups='leiden')
        
        # PAGA 그래프 시각화
        sc.pl.paga(adata, threshold=0.05, layout='fr', save='paga.pdf', show=False)
        
        # PAGA 초기화된 UMAP 시각화
        sc.tl.umap(adata, init_pos='paga')
        sc.pl.umap(adata, color=['leiden', 'dpt_pseudotime'], 
                  save='paga_initialized_umap.pdf', show=False)
    else:
        print("클러스터 정보(leiden)가 없어 PAGA 분석을 건너뜁니다.")
    
    # 5. 의사시간에 따른 유전자 발현 변화 분석
    print("의사시간에 따른 유전자 발현 변화 분석 중...")
    
    # 의사시간에 따라 발현이 변하는 상위 유전자 선택
    try:
        # log1p 'base' 키 문제 해결
        if 'log1p' in adata.uns and 'base' not in adata.uns['log1p']:
            adata.uns['log1p'] = {'base': None}
        
        sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
        
        # 상위 차별 발현 유전자 추출
        if 'rank_genes_groups' in adata.uns:
            top_genes = []
            for group in adata.uns['rank_genes_groups']['names'].dtype.names:
                top_genes.extend(adata.uns['rank_genes_groups']['names'][group][:5])
            
            # 중복 제거
            top_genes = list(set(top_genes))
            
            # 의사시간에 따른 상위 유전자 발현 시각화
            sc.pl.diffmap(adata, color=['dpt_pseudotime'] + top_genes[:5], 
                        save='top_genes_diffmap.pdf', show=False)
            
            # 의사시간에 따른 유전자 발현 히트맵
            adata.obs['dpt_pseudotime_bins'] = pd.cut(adata.obs['dpt_pseudotime'], bins=10)
            sc.pl.heatmap(adata, var_names=top_genes[:20], groupby='dpt_pseudotime_bins', 
                        save='pseudotime_heatmap.pdf', show=False)
    except Exception as e:
        print(f"유전자 발현 변화 분석 중 오류 발생: {e}")
        print("유전자 발현 변화 분석을 건너뜁니다.")
    
    # 6. 궤적 분기점 분석 (선택사항)
    print("궤적 분기점 분석 중...")
    
    # 분기점 계산
    try:
        sc.tl.dpt(adata, n_dcs=10, n_branchings=1)
        
        # 분기점 시각화
        if 'dpt_groups' in adata.obs:
            sc.pl.umap(adata, color=['dpt_groups'], 
                      save='dpt_groups.pdf', show=False)
            
            # 분기점별 차별 발현 유전자 분석
            sc.tl.rank_genes_groups(adata, 'dpt_groups', method='wilcoxon')
            sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='dpt_groups', 
                                          save='branch_markers_heatmap.pdf', show=False)
    except Exception as e:
        print(f"분기점 분석 중 오류 발생: {e}")
        print("분기점 분석을 건너뜁니다.")
    
    # 7. Velocity 분석 준비 (RNA velocity 데이터가 있는 경우)
    print("RNA velocity 분석 준비 중...")
    print("참고: 실제 RNA velocity 분석을 위해서는 spliced/unspliced 정보가 포함된 데이터가 필요합니다.")
    
    # 8. 결과 저장
    print("궤적 분석 결과 저장 중...")
    
    # Convert the interval column to strings (if it exists) to avoid serialization issues
    if 'dpt_pseudotime_bins' in adata.obs.columns:
        adata.obs['dpt_pseudotime_bins'] = adata.obs['dpt_pseudotime_bins'].astype(str)
        
    adata.write(os.path.join(output_dir, 'trajectory_results.h5ad'))
    
    # 9. 의사시간 정보 CSV 저장
    pseudotime_df = pd.DataFrame({
        'cell_id': adata.obs.index,
        'pseudotime': adata.obs['dpt_pseudotime']
    })
    
    if 'leiden' in adata.obs:
        pseudotime_df['cluster'] = adata.obs['leiden']
    
    if 'cell_type' in adata.obs:
        pseudotime_df['cell_type'] = adata.obs['cell_type']
    
    if 'dpt_groups' in adata.obs:
        pseudotime_df['branch'] = adata.obs['dpt_groups']
    
    pseudotime_df.to_csv(os.path.join(output_dir, 'pseudotime_data.csv'), index=False)
    
    # 10. 추가 시각화: 의사시간에 따른 세포 분포
    plt.figure(figsize=(10, 6))
    sns.histplot(data=adata.obs, x='dpt_pseudotime', bins=50)
    plt.title('의사시간에 따른 세포 분포')
    plt.xlabel('의사시간')
    plt.ylabel('세포 수')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'pseudotime_distribution.pdf'))
    plt.savefig(os.path.join(figures_dir, 'pseudotime_distribution.png'), dpi=300)
    plt.close()
    
    # 11. 클러스터별 의사시간 분포 (바이올린 플롯)
    if 'leiden' in adata.obs:
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=adata.obs, x='leiden', y='dpt_pseudotime')
        plt.title('클러스터별 의사시간 분포')
        plt.xlabel('클러스터')
        plt.ylabel('의사시간')
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'cluster_pseudotime_violin.pdf'))
        plt.savefig(os.path.join(figures_dir, 'cluster_pseudotime_violin.png'), dpi=300)
        plt.close()
    
    # 12. 세포 유형별 의사시간 분포 (바이올린 플롯)
    if 'cell_type' in adata.obs:
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=adata.obs, x='cell_type', y='dpt_pseudotime')
        plt.title('세포 유형별 의사시간 분포')
        plt.xlabel('세포 유형')
        plt.ylabel('의사시간')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'cell_type_pseudotime_violin.pdf'))
        plt.savefig(os.path.join(figures_dir, 'cell_type_pseudotime_violin.png'), dpi=300)
        plt.close()
    
    print(f"\n세포 궤적 분석이 완료되었습니다!")
    print(f"- 시각화 결과: {figures_dir}")
    print(f"- 궤적 분석 결과 데이터: {os.path.join(output_dir, 'trajectory_results.h5ad')}")
    print(f"- 의사시간 데이터: {os.path.join(output_dir, 'pseudotime_data.csv')}")

if __name__ == "__main__":
    run_trajectory_analysis() 


"""
# 시각화 결과 파일들
# 1. 의사시간 UMAP 시각화 
'pseudotime.pdf'

# 2. PAGA 그래프 시각화
'paga.pdf'

# 3. PAGA 초기화된 UMAP 시각화 
'paga_initialized_umap.pdf'

# 4. 상위 유전자 발현 시각화
'top_genes_diffmap.pdf'

# 5. 의사시간에 따른 유전자 발현 히트맵
'pseudotime_heatmap.pdf'

# 6. 궤적 분기점 시각화
'dpt_groups.pdf'

# 7. 분기점별 차별 발현 유전자 히트맵
'branch_markers_heatmap.pdf'

# 8. 의사시간에 따른 세포 분포 그래프
'pseudotime_distribution.pdf'
'pseudotime_distribution.png'

# 9. 클러스터별 의사시간 분포 바이올린 플롯
'cluster_pseudotime_violin.pdf'
'cluster_pseudotime_violin.png'

# 10. 세포 유형별 의사시간 분포 바이올린 플롯
'cell_type_pseudotime_violin.pdf'
'cell_type_pseudotime_violin.png'


데이터 결과 파일들

# 1. 궤적 분석 결과 AnnData 객체
'trajectory_results.h5ad'

# 2. 의사시간 정보가 담긴 CSV 파일
'pseudotime_data.csv'
"""