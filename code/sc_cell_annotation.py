#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 필요한 라이브러리 임포트
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
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

# 세포 유형별 마커 유전자 정의 (예시)
# 실제 분석에서는 해당 조직/세포 유형에 맞는 마커 유전자로 수정 필요
CELL_MARKERS = {
    'T 세포': ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD4', 'IL7R'],
    'B 세포': ['CD79A', 'CD79B', 'MS4A1', 'CD19'],
    'NK 세포': ['NCAM1', 'NKG7', 'KLRD1', 'KLRF1'],
    '단핵구': ['CD14', 'LYZ', 'CST3', 'FCGR3A', 'MS4A7'],
    '대식세포': ['MARCO', 'MSR1', 'CD68', 'ITGAM'],
    '수지상세포': ['FCER1A', 'CLEC10A', 'CLEC9A', 'LILRA4', 'CLEC4C'],
    '상피세포': ['EPCAM', 'KRT8', 'KRT18', 'KRT19'],
    '섬유아세포': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'],
    '내피세포': ['PECAM1', 'VWF', 'CDH5', 'CLDN5'],
    '적혈구': ['HBA1', 'HBA2', 'HBB'],
    '혈소판': ['PPBP', 'PF4']
}

def run_cell_annotation():
    """
    클러스터링된 세포에 대한 세포 유형 주석(annotation)을 수행합니다.
    - 마커 유전자 기반 세포 유형 주석
    - 세포 유형별 발현 패턴 시각화
    - 세포 유형 분포 분석
    """
    
    print("클러스터링 완료된 데이터 로딩 중...")
    # 클러스터링 완료된 데이터 로드
    try:
        adata = sc.read_h5ad(os.path.join(output_dir, 'downstream_results.h5ad'))
        print(f"클러스터링 완료된 데이터 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    except FileNotFoundError:
        try:
            adata = sc.read_h5ad(os.path.join(output_dir, 'complete_results.h5ad'))
            print(f"전체 데이터 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
        except FileNotFoundError:
            print("클러스터링 완료된 데이터 파일을 찾을 수 없습니다. sc_downstream.py를 먼저 실행해주세요.")
            return
    
    # 1. 마커 유전자 발현 확인
    print("마커 유전자 발현 확인 중...")
    
    # 마커 유전자 목록 평탄화 및 중복 제거
    all_markers = list(set([gene for genes in CELL_MARKERS.values() for gene in genes]))
    
    # 데이터에 존재하는 마커 유전자만 필터링
    existing_markers = [gene for gene in all_markers if gene in adata.var_names]
    missing_markers = [gene for gene in all_markers if gene not in adata.var_names]
    
    if missing_markers:
        print(f"주의: 다음 마커 유전자가 데이터에 없습니다: {', '.join(missing_markers)}")
    
    print(f"분석에 사용할 마커 유전자 수: {len(existing_markers)}")
    
    # 2. 클러스터별 마커 유전자 발현 시각화
    print("클러스터별 마커 유전자 발현 시각화 중...")
    
    # 도트플롯 - 클러스터별 마커 유전자 발현
    if existing_markers:
        sc.pl.dotplot(adata, var_names=existing_markers, groupby='leiden', 
                     dendrogram=True, save='marker_genes_by_cluster.pdf', show=False)
    
    # 3. 세포 유형 스코어 계산
    print("세포 유형 스코어 계산 중...")
    
    # 각 세포 유형별 스코어 계산
    cell_type_scores = {}
    for cell_type, markers in CELL_MARKERS.items():
        # 데이터에 존재하는 마커만 사용
        existing_cell_markers = [gene for gene in markers if gene in adata.var_names]
        if existing_cell_markers:
            # 해당 세포 유형의 마커 유전자 평균 발현량 계산
            sc.tl.score_genes(adata, gene_list=existing_cell_markers, score_name=f'{cell_type}_score')
            cell_type_scores[cell_type] = existing_cell_markers
    
    # 4. 세포 유형 스코어 시각화
    print("세포 유형 스코어 시각화 중...")
    
    # 세포 유형 스코어 컬럼 목록
    score_columns = [key for key in adata.obs.columns if key.endswith('_score')]
    
    if score_columns:
        # UMAP에 세포 유형 스코어 시각화
        sc.pl.umap(adata, color=score_columns, ncols=3, save='cell_type_scores.pdf', show=False)
    
    # 5. 클러스터별 세포 유형 할당
    print("클러스터별 세포 유형 할당 중...")
    
    # 클러스터별 평균 세포 유형 스코어 계산
    cluster_scores = pd.DataFrame(index=adata.obs['leiden'].cat.categories)
    
    for score_col in score_columns:
        cluster_means = adata.obs.groupby('leiden')[score_col].mean()
        cluster_scores[score_col] = cluster_means
    
    # 각 클러스터에 가장 높은 스코어를 가진 세포 유형 할당
    cluster_scores['assigned_cell_type'] = cluster_scores[score_columns].idxmax(axis=1)
    cluster_scores['assigned_cell_type'] = cluster_scores['assigned_cell_type'].str.replace('_score', '')
    
    # 결과 저장
    cluster_scores.to_csv(os.path.join(output_dir, 'cluster_cell_type_scores.csv'))
    
    # 6. 세포 유형 할당 결과 시각화
    print("세포 유형 할당 결과 시각화 중...")
    
    # 클러스터별 할당된 세포 유형 히트맵
    plt.figure(figsize=(12, 8))
    sns.heatmap(cluster_scores[score_columns], cmap='viridis', annot=True, fmt='.2f')
    plt.title('클러스터별 세포 유형 스코어')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'cluster_cell_type_heatmap.pdf'))
    plt.savefig(os.path.join(figures_dir, 'cluster_cell_type_heatmap.png'), dpi=300)
    plt.close()
    
    # 7. 세포에 세포 유형 정보 추가
    print("세포에 세포 유형 정보 추가 중...")
    
    # 클러스터 ID를 기반으로 세포 유형 할당
    adata.obs['cell_type'] = adata.obs['leiden'].astype(str).map(cluster_scores['assigned_cell_type'])
    
    # 8. 세포 유형별 UMAP 시각화
    print("세포 유형별 UMAP 시각화 중...")
    
    # 세포 유형별 색상 지정 - NumPy 배열 대신 딕셔너리로 변환
    cell_types = adata.obs['cell_type'].unique()
    n_cell_types = len(cell_types)
    colors = plt.cm.tab20(np.linspace(0, 1, n_cell_types))
    
    # 세포 유형과 색상을 매핑하는 딕셔너리 생성
    cell_type_color_dict = {
        cell_type: colors[i].tolist() for i, cell_type in enumerate(cell_types)
    }
    
    # UMAP에 세포 유형 시각화 - 색상 딕셔너리 전달
    sc.pl.umap(adata, color='cell_type', palette=cell_type_color_dict, 
              save='cell_types.pdf', show=False)
    
    # 9. 세포 유형별 마커 유전자 발현 시각화
    print("세포 유형별 마커 유전자 발현 시각화 중...")
    
    # 세포 유형별 마커 유전자 발현 도트플롯
    if existing_markers:
        sc.pl.dotplot(adata, var_names=existing_markers, groupby='cell_type', 
                     dendrogram=True, save='marker_genes_by_cell_type.pdf', show=False)
    
    # 10. 세포 유형별 분포 분석
    print("세포 유형별 분포 분석 중...")
    
    # 세포 유형별 세포 수 계산
    cell_type_counts = adata.obs['cell_type'].value_counts().sort_values(ascending=False)
    cell_type_counts.name = '세포 수'
    cell_type_counts.index.name = '세포 유형'
    
    # 결과 저장
    cell_type_counts.to_csv(os.path.join(output_dir, 'cell_type_counts.csv'))
    
    # 세포 유형별 분포 시각화
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(x=cell_type_counts.index, y=cell_type_counts.values)
    ax.set_title('세포 유형별 세포 수')
    ax.set_xlabel('세포 유형')
    ax.set_ylabel('세포 수')
    plt.xticks(rotation=45, ha='right')
    
    # 막대 위에 숫자 표시
    for i, count in enumerate(cell_type_counts.values):
        ax.text(i, count + 5, str(count), ha='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'cell_type_distribution.pdf'))
    plt.savefig(os.path.join(figures_dir, 'cell_type_distribution.png'), dpi=300)
    plt.close()
    
    # 11. 세포 유형별 비율 파이 차트
    plt.figure(figsize=(10, 10))
    plt.pie(cell_type_counts.values, labels=cell_type_counts.index, autopct='%1.1f%%', 
           shadow=True, startangle=90)
    plt.axis('equal')
    plt.title('세포 유형 비율')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'cell_type_pie_chart.pdf'))
    plt.savefig(os.path.join(figures_dir, 'cell_type_pie_chart.png'), dpi=300)
    plt.close()
    
    # 12. 세포 유형 정보를 포함한 결과 저장
    print("세포 유형 정보를 포함한 결과 저장 중...")
    adata.write(os.path.join(output_dir, 'annotated_results.h5ad'))
    
    # 13. 세포 유형 주석 정보 JSON 저장
    cell_type_info = {
        'cell_types': {ct: list(markers) for ct, markers in cell_type_scores.items()},
        'cluster_assignments': cluster_scores['assigned_cell_type'].to_dict(),
        'cell_type_counts': cell_type_counts.to_dict()
    }
    
    with open(os.path.join(output_dir, 'cell_type_annotation.json'), 'w', encoding='utf-8') as f:
        json.dump(cell_type_info, f, ensure_ascii=False, indent=2)
    
    print(f"\n세포 유형 주석 분석이 완료되었습니다!")
    print(f"- 시각화 결과: {figures_dir}")
    print(f"- 주석 결과 데이터: {os.path.join(output_dir, 'annotated_results.h5ad')}")
    print(f"- 클러스터별 세포 유형 스코어: {os.path.join(output_dir, 'cluster_cell_type_scores.csv')}")
    print(f"- 세포 유형별 세포 수: {os.path.join(output_dir, 'cell_type_counts.csv')}")
    print(f"- 세포 유형 주석 정보: {os.path.join(output_dir, 'cell_type_annotation.json')}")

if __name__ == "__main__":
    run_cell_annotation() 