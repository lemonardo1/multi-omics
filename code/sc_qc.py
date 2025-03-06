#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 필요한 라이브러리 임포트
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

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

def run_qc():
    """단일세포 RNA-seq 데이터에 대한 품질 관리(QC) 분석을 수행합니다."""
    
    print("데이터 로딩 중...")
    # 10x Genomics 데이터 읽기
    adata = sc.read_10x_mtx(
        os.path.join(data_dir, 'filtered_feature_bc_matrix/'),
        var_names='gene_symbols',
        cache=True
    )
    
    print(f"데이터 로딩 완료: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    
    # 기본적인 전처리 단계
    print("기본 필터링 적용 중...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # 품질 지표 계산
    print("QC 지표 계산 중...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 미토콘드리아 유전자 주석
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # QC 지표 요약 통계
    qc_summary = pd.DataFrame({
        '세포 수': [adata.n_obs],
        '유전자 수': [adata.n_vars],
        '세포당 평균 유전자 수': [adata.obs['n_genes_by_counts'].mean()],
        '세포당 중앙값 유전자 수': [adata.obs['n_genes_by_counts'].median()],
        '세포당 평균 UMI 수': [adata.obs['total_counts'].mean()],
        '세포당 중앙값 UMI 수': [adata.obs['total_counts'].median()],
        '평균 미토콘드리아 비율(%)': [adata.obs['pct_counts_mt'].mean()]
    }).T.rename(columns={0: '값'})
    
    print("\nQC 요약 통계:")
    print(qc_summary)
    
    # QC 요약 통계 저장
    qc_summary.to_csv(os.path.join(output_dir, 'qc_summary.csv'))
    
    # 품질 지표 시각화
    print("QC 시각화 생성 중...")
    
    # 1. 히스토그램 - 세포당 유전자 수와 미토콘드리아 비율
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    sns.histplot(data=adata.obs, x='n_genes_by_counts', bins=60, ax=axs[0])
    axs[0].set_title('세포당 유전자 수 분포')
    axs[0].set_xlabel('유전자 수')
    axs[0].set_ylabel('세포 수')
    
    sns.histplot(data=adata.obs, x='pct_counts_mt', bins=60, ax=axs[1])
    axs[1].set_title('미토콘드리아 비율 분포')
    axs[1].set_xlabel('미토콘드리아 비율 (%)')
    axs[1].set_ylabel('세포 수')
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'qc_histograms.pdf'))
    plt.savefig(os.path.join(figures_dir, 'qc_histograms.png'), dpi=300)
    plt.close()
    
    # 2. 산점도 - 유전자 수 vs UMI 수, 색상은 미토콘드리아 비율
    plt.figure(figsize=(8, 6))
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', 
                 save='gene_count_vs_mt.pdf', show=False)
    plt.savefig(os.path.join(figures_dir, 'gene_count_vs_mt.png'), dpi=300)
    plt.close()
    
    # 3. 바이올린 플롯 - 주요 QC 지표
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True, save='qc_violins.pdf', show=False)
    
    # 4. 추가 바이올린 플롯 - 더 자세한 시각화
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    
    sns.violinplot(data=adata.obs, y='n_genes_by_counts', ax=axs[0])
    axs[0].set_title('세포당 유전자 수')
    axs[0].set_ylabel('유전자 수')
    
    sns.violinplot(data=adata.obs, y='total_counts', ax=axs[1])
    axs[1].set_title('세포당 UMI 수')
    axs[1].set_ylabel('UMI 수')
    
    sns.violinplot(data=adata.obs, y='pct_counts_mt', ax=axs[2])
    axs[2].set_title('미토콘드리아 비율')
    axs[2].set_ylabel('미토콘드리아 비율 (%)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'qc_violins_detailed.pdf'))
    plt.savefig(os.path.join(figures_dir, 'qc_violins_detailed.png'), dpi=300)
    plt.close()
    
    # 5. 상관관계 히트맵
    qc_vars = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    corr = adata.obs[qc_vars].corr()
    
    plt.figure(figsize=(8, 7))
    sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1, center=0)
    plt.title('QC 지표 간 상관관계')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'qc_correlation.pdf'))
    plt.savefig(os.path.join(figures_dir, 'qc_correlation.png'), dpi=300)
    plt.close()
    
    # QC 기준에 따른 필터링 결과 시각화
    print("필터링 기준 시각화 중...")
    
    # 필터링 기준 설정
    min_genes = 200
    max_genes = 2500
    max_mt_pct = 5
    
    # 필터링 전후 세포 수 계산
    total_cells = adata.n_obs
    cells_pass_min_genes = np.sum(adata.obs['n_genes_by_counts'] >= min_genes)
    cells_pass_max_genes = np.sum(adata.obs['n_genes_by_counts'] <= max_genes)
    cells_pass_mt = np.sum(adata.obs['pct_counts_mt'] <= max_mt_pct)
    cells_pass_all = np.sum((adata.obs['n_genes_by_counts'] >= min_genes) & 
                           (adata.obs['n_genes_by_counts'] <= max_genes) & 
                           (adata.obs['pct_counts_mt'] <= max_mt_pct))
    
    # 필터링 결과 저장
    filter_results = pd.DataFrame({
        '기준': ['총 세포 수', 
                f'최소 {min_genes}개 이상 유전자 발현', 
                f'최대 {max_genes}개 이하 유전자 발현', 
                f'미토콘드리아 비율 {max_mt_pct}% 이하',
                '모든 기준 통과'],
        '세포 수': [total_cells, cells_pass_min_genes, cells_pass_max_genes, cells_pass_mt, cells_pass_all],
        '비율(%)': [100, 
                  cells_pass_min_genes/total_cells*100, 
                  cells_pass_max_genes/total_cells*100, 
                  cells_pass_mt/total_cells*100,
                  cells_pass_all/total_cells*100]
    })
    
    filter_results.to_csv(os.path.join(output_dir, 'filtering_results.csv'), index=False)
    
    print("\n필터링 결과:")
    print(filter_results)
    
    # 필터링 결과 시각화
    plt.figure(figsize=(10, 6))
    sns.barplot(data=filter_results, x='기준', y='비율(%)')
    plt.xticks(rotation=45, ha='right')
    plt.title('QC 필터링 기준별 통과 세포 비율')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'filtering_results.pdf'))
    plt.savefig(os.path.join(figures_dir, 'filtering_results.png'), dpi=300)
    plt.close()
    
    # 필터링 적용
    print("QC 필터링 적용 중...")
    adata_filtered = adata[(adata.obs.n_genes_by_counts >= min_genes) &
                          (adata.obs.n_genes_by_counts <= max_genes) &
                          (adata.obs.pct_counts_mt <= max_mt_pct), :].copy()
    
    print(f"필터링 후 남은 세포: {adata_filtered.n_obs} (원래 {adata.n_obs}에서 {adata_filtered.n_obs/adata.n_obs*100:.1f}%)")
    
    # 필터링 전후 비교 시각화
    fig, axs = plt.subplots(2, 3, figsize=(18, 12))
    
    # 필터링 전
    sns.violinplot(data=adata.obs, y='n_genes_by_counts', ax=axs[0, 0])
    axs[0, 0].set_title('필터링 전: 세포당 유전자 수')
    axs[0, 0].axhline(y=min_genes, color='r', linestyle='--')
    axs[0, 0].axhline(y=max_genes, color='r', linestyle='--')
    
    sns.violinplot(data=adata.obs, y='total_counts', ax=axs[0, 1])
    axs[0, 1].set_title('필터링 전: 세포당 UMI 수')
    
    sns.violinplot(data=adata.obs, y='pct_counts_mt', ax=axs[0, 2])
    axs[0, 2].set_title('필터링 전: 미토콘드리아 비율')
    axs[0, 2].axhline(y=max_mt_pct, color='r', linestyle='--')
    
    # 필터링 후
    sns.violinplot(data=adata_filtered.obs, y='n_genes_by_counts', ax=axs[1, 0])
    axs[1, 0].set_title('필터링 후: 세포당 유전자 수')
    
    sns.violinplot(data=adata_filtered.obs, y='total_counts', ax=axs[1, 1])
    axs[1, 1].set_title('필터링 후: 세포당 UMI 수')
    
    sns.violinplot(data=adata_filtered.obs, y='pct_counts_mt', ax=axs[1, 2])
    axs[1, 2].set_title('필터링 후: 미토콘드리아 비율')
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'before_after_filtering.pdf'))
    plt.savefig(os.path.join(figures_dir, 'before_after_filtering.png'), dpi=300)
    plt.close()
    
    # 필터링된 데이터 저장
    adata_filtered.write(os.path.join(output_dir, 'filtered_data.h5ad'))
    
    print(f"\nQC 분석이 완료되었습니다!")
    print(f"- QC 시각화: {figures_dir}")
    print(f"- 필터링된 데이터: {os.path.join(output_dir, 'filtered_data.h5ad')}")
    print(f"- QC 통계: {os.path.join(output_dir, 'qc_summary.csv')}")
    print(f"- 필터링 결과: {os.path.join(output_dir, 'filtering_results.csv')}")

if __name__ == "__main__":
    run_qc() 
