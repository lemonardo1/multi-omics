#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 필요한 라이브러리 임포트
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import anndata
from scipy import stats
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
import gseapy as gp
from adjustText import adjust_text
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# 설정: 기본 디렉토리
base_dir = '/home/swr0460/daeseong'
data_dir = os.path.join(base_dir, 'data')
figures_dir = os.path.join(base_dir, 'figures')
output_dir = os.path.join(base_dir, 'output')
advanced_dir = os.path.join(figures_dir, 'advanced')

# 디렉토리 생성
os.makedirs(advanced_dir, exist_ok=True)

# Scanpy의 기본 figure 경로를 설정
sc.settings.figdir = advanced_dir

# 재현성을 위한 랜덤 시드 설정
np.random.seed(42)

# 한글 폰트 설정 (필요한 경우)
plt.rcParams['font.family'] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False

def run_advanced_analysis():
    """
    세포 궤적 분석 결과에 대한 추가 분석을 수행합니다.
    - 의사시간 기반 유전자 발현 상관관계 분석
    - 전사 인자 활성화 분석
    - 세포 상태 전환 분석
    - 심화 시각화
    """
    
    print("\n===== 세포 궤적 분석 결과 심화 분석 시작 =====")
    
    # 1. 이전 분석 결과 로드
    print("이전 궤적 분석 결과 로드 중...")
    try:
        adata = sc.read_h5ad(os.path.join(output_dir, 'trajectory_results.h5ad'))
        print(f"궤적 분석 결과 로드 성공: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    except FileNotFoundError:
        print("궤적 분석 결과 파일이 없습니다. 먼저 sc_trajectory.py를 실행해주세요.")
        return
    
    # 2. 의사시간과 유전자 발현 간의 상관관계 분석
    print("의사시간과 유전자 발현 간의 상관관계 분석 중...")
    
    # 의사시간 벡터 준비
    pseudotime = adata.obs['dpt_pseudotime'].values
    
    # 상관관계 계산을 위한 데이터프레임 생성
    gene_names = adata.var_names
    expression_matrix = adata.X.toarray() if isinstance(adata.X, np.ndarray) is False else adata.X
    
    # 스피어만 상관계수 계산
    correlation_scores = []
    p_values = []
    
    for i in range(expression_matrix.shape[1]):
        gene_expr = expression_matrix[:, i]
        corr, p_val = spearmanr(pseudotime, gene_expr)
        correlation_scores.append(corr)
        p_values.append(p_val)
    
    # 결과를 데이터프레임으로 정리
    correlation_df = pd.DataFrame({
        'gene': gene_names,
        'correlation': correlation_scores,
        'p_value': p_values
    })
    
    # 다중검정 보정 (Benjamini-Hochberg)
    correlation_df['adjusted_p_value'] = stats.false_discovery_control(correlation_df['p_value'])
    
    # 통계적으로 유의미한 상관관계를 가진 유전자 필터링
    significant_genes = correlation_df[(correlation_df['adjusted_p_value'] < 0.05)]
    
    # 상관관계 방향에 따라 분류
    positive_corr_genes = significant_genes[significant_genes['correlation'] > 0].sort_values('correlation', ascending=False)
    negative_corr_genes = significant_genes[significant_genes['correlation'] < 0].sort_values('correlation', ascending=True)
    
    print(f"의사시간과 양의 상관관계를 가진 유전자: {len(positive_corr_genes)} 개")
    print(f"의사시간과 음의 상관관계를 가진 유전자: {len(negative_corr_genes)} 개")
    
    # 상위 상관관계 유전자 저장
    correlation_df.to_csv(os.path.join(output_dir, 'pseudotime_gene_correlation.csv'), index=False)
    positive_corr_genes.head(200).to_csv(os.path.join(output_dir, 'positive_corr_genes.csv'), index=False)
    negative_corr_genes.head(200).to_csv(os.path.join(output_dir, 'negative_corr_genes.csv'), index=False)
    
    # 상위 상관관계 유전자 시각화
    top_n_genes = 10
    
    # 상위 양의 상관관계 유전자
    if len(positive_corr_genes) > 0:
        top_pos_genes = positive_corr_genes.head(top_n_genes)['gene'].tolist()
        
        # 각 유전자에 대해 개별적으로 플롯 생성
        for i, gene in enumerate(top_pos_genes):
            sc.pl.scatter(adata, x='dpt_pseudotime', y=gene, 
                        color='dpt_pseudotime', color_map='viridis', 
                        save=f'top_positive_corr_gene_{i}_{gene}.pdf', show=False)
        
        # 모든 상위 유전자에 대한 시각화 (matplotlib 사용)
        fig, axes = plt.subplots(5, 2, figsize=(14, 20))
        axes = axes.flatten()
        
        for i, gene in enumerate(top_pos_genes):
            if i < len(axes):
                gene_idx = np.where(adata.var_names == gene)[0][0]
                x = adata.obs['dpt_pseudotime'].values
                y = adata.X.toarray()[:, gene_idx] if not isinstance(adata.X, np.ndarray) else adata.X[:, gene_idx]
                
                axes[i].scatter(x, y, c=x, cmap='viridis', s=10, alpha=0.5)
                axes[i].set_title(gene)
                axes[i].set_xlabel('의사시간')
                axes[i].set_ylabel('발현량')
        
        plt.tight_layout()
        plt.savefig(os.path.join(advanced_dir, 'top_positive_corr_genes_combined.pdf'))
        plt.savefig(os.path.join(advanced_dir, 'top_positive_corr_genes_combined.png'), dpi=300)
        plt.close()
    
    # 상위 음의 상관관계 유전자
    if len(negative_corr_genes) > 0:
        top_neg_genes = negative_corr_genes.head(top_n_genes)['gene'].tolist()
        
        # 각 유전자에 대해 개별적으로 플롯 생성
        for i, gene in enumerate(top_neg_genes):
            sc.pl.scatter(adata, x='dpt_pseudotime', y=gene, 
                        color='dpt_pseudotime', color_map='viridis', 
                        save=f'top_negative_corr_gene_{i}_{gene}.pdf', show=False)
        
        # 모든 상위 유전자에 대한 시각화 (matplotlib 사용)
        fig, axes = plt.subplots(5, 2, figsize=(14, 20))
        axes = axes.flatten()
        
        for i, gene in enumerate(top_neg_genes):
            if i < len(axes):
                gene_idx = np.where(adata.var_names == gene)[0][0]
                x = adata.obs['dpt_pseudotime'].values
                y = adata.X.toarray()[:, gene_idx] if not isinstance(adata.X, np.ndarray) else adata.X[:, gene_idx]
                
                axes[i].scatter(x, y, c=x, cmap='viridis', s=10, alpha=0.5)
                axes[i].set_title(gene)
                axes[i].set_xlabel('의사시간')
                axes[i].set_ylabel('발현량')
        
        plt.tight_layout()
        plt.savefig(os.path.join(advanced_dir, 'top_negative_corr_genes_combined.pdf'))
        plt.savefig(os.path.join(advanced_dir, 'top_negative_corr_genes_combined.png'), dpi=300)
        plt.close()
    
    # 3. 상위 상관관계 유전자들의 GO 분석 (gseapy 사용)
    try:
        print("상위 상관관계 유전자들의 GO 분석 중...")
        
        # 양의 상관관계 유전자 GO 분석
        if len(positive_corr_genes) >= 20:
            pos_gene_list = positive_corr_genes.head(100)['gene'].tolist()
            
            # GO 분석 수행
            pos_enr = gp.enrichr(gene_list=pos_gene_list,
                              organism='Human',
                              gene_sets=['GO_Biological_Process_2021'],
                              outdir=os.path.join(advanced_dir, 'gseapy_pos'))
            
            # 결과 시각화 (최신 버전 gseapy 지원)
            try:
                # 먼저 직접 plot 메서드를 시도
                pos_enr.plot(top_term=10, figsize=(12, 6), title='양의 상관관계 유전자 GO 분석')
                plt.tight_layout()
                plt.savefig(os.path.join(advanced_dir, 'positive_corr_GO.pdf'))
                plt.savefig(os.path.join(advanced_dir, 'positive_corr_GO.png'), dpi=300)
                plt.close()
            except AttributeError:
                # AttributeError 발생 시 결과 데이터프레임을 직접 시각화
                print("gseapy 버전이 다른 것 같습니다. 대체 시각화 방법 사용...")
                
                # 결과가 데이터프레임인지 확인
                if isinstance(pos_enr, pd.DataFrame):
                    result_df = pos_enr
                    # 상위 10개 결과 선택
                    if 'Adjusted P-value' in result_df.columns:
                        top_results = result_df.sort_values('Adjusted P-value').head(10)
                    elif 'P-value' in result_df.columns:
                        top_results = result_df.sort_values('P-value').head(10)
                    else:
                        top_results = result_df.head(10)
                    
                    # 시각화할 열 선택
                    term_col = 'Term' if 'Term' in top_results.columns else ('GO_Biological_Process_2021' if 'GO_Biological_Process_2021' in top_results.columns else top_results.columns[0])
                    pval_col = 'Adjusted P-value' if 'Adjusted P-value' in top_results.columns else ('P-value' if 'P-value' in top_results.columns else top_results.columns[1])
                    
                    # 시각화
                    plt.figure(figsize=(12, 6))
                    plt.barh(top_results[term_col], -np.log10(top_results[pval_col]))
                    plt.xlabel('-log10(P-value)')
                    plt.ylabel('GO Term')
                    plt.title('양의 상관관계 유전자 GO 분석')
                    plt.tight_layout()
                    plt.savefig(os.path.join(advanced_dir, 'positive_corr_GO.pdf'))
                    plt.savefig(os.path.join(advanced_dir, 'positive_corr_GO.png'), dpi=300)
                    plt.close()
                else:
                    print("GO 분석 결과가 예상한 형식이 아닙니다.")
                    
                    # 결과를 CSV로 저장
                    if hasattr(pos_enr, 'results'):
                        pos_enr.results.to_csv(os.path.join(output_dir, 'positive_corr_GO_results.csv'))
                        print(f"GO 분석 결과를 {os.path.join(output_dir, 'positive_corr_GO_results.csv')}에 저장했습니다.")
        
        # 음의 상관관계 유전자 GO 분석
        if len(negative_corr_genes) >= 20:
            neg_gene_list = negative_corr_genes.head(100)['gene'].tolist()
            
            # GO 분석 수행
            neg_enr = gp.enrichr(gene_list=neg_gene_list,
                              organism='Human',
                              gene_sets=['GO_Biological_Process_2021'],
                              outdir=os.path.join(advanced_dir, 'gseapy_neg'))
            
            # 결과 시각화 (최신 버전 gseapy 지원)
            try:
                # 먼저 직접 plot 메서드를 시도
                neg_enr.plot(top_term=10, figsize=(12, 6), title='음의 상관관계 유전자 GO 분석')
                plt.tight_layout()
                plt.savefig(os.path.join(advanced_dir, 'negative_corr_GO.pdf'))
                plt.savefig(os.path.join(advanced_dir, 'negative_corr_GO.png'), dpi=300)
                plt.close()
            except AttributeError:
                # AttributeError 발생 시 결과 데이터프레임을 직접 시각화
                print("gseapy 버전이 다른 것 같습니다. 대체 시각화 방법 사용...")
                
                # 결과가 데이터프레임인지 확인
                if isinstance(neg_enr, pd.DataFrame):
                    result_df = neg_enr
                    # 상위 10개 결과 선택
                    if 'Adjusted P-value' in result_df.columns:
                        top_results = result_df.sort_values('Adjusted P-value').head(10)
                    elif 'P-value' in result_df.columns:
                        top_results = result_df.sort_values('P-value').head(10)
                    else:
                        top_results = result_df.head(10)
                    
                    # 시각화할 열 선택
                    term_col = 'Term' if 'Term' in top_results.columns else ('GO_Biological_Process_2021' if 'GO_Biological_Process_2021' in top_results.columns else top_results.columns[0])
                    pval_col = 'Adjusted P-value' if 'Adjusted P-value' in top_results.columns else ('P-value' if 'P-value' in top_results.columns else top_results.columns[1])
                    
                    # 시각화
                    plt.figure(figsize=(12, 6))
                    plt.barh(top_results[term_col], -np.log10(top_results[pval_col]))
                    plt.xlabel('-log10(P-value)')
                    plt.ylabel('GO Term')
                    plt.title('음의 상관관계 유전자 GO 분석')
                    plt.tight_layout()
                    plt.savefig(os.path.join(advanced_dir, 'negative_corr_GO.pdf'))
                    plt.savefig(os.path.join(advanced_dir, 'negative_corr_GO.png'), dpi=300)
                    plt.close()
                else:
                    print("GO 분석 결과가 예상한 형식이 아닙니다.")
                    
                    # 결과를 CSV로 저장
                    if hasattr(neg_enr, 'results'):
                        neg_enr.results.to_csv(os.path.join(output_dir, 'negative_corr_GO_results.csv'))
                        print(f"GO 분석 결과를 {os.path.join(output_dir, 'negative_corr_GO_results.csv')}에 저장했습니다.")
    
    except Exception as e:
        print(f"GO 분석 중 오류 발생: {e}")
        print("GO 분석을 건너뜁니다.")
    
    # 4. 의사시간에 따른 세포 상태 전환 분석 (변화율 계산)
    print("의사시간에 따른 세포 상태 전환 분석 중...")
    
    # 의사시간 구간 생성 (20개 구간)
    n_bins = 20
    adata.obs['pt_bins'] = pd.cut(adata.obs['dpt_pseudotime'], bins=n_bins)
    
    # 구간별 평균 발현량 계산 - 수정된 방식
    bin_means = pd.DataFrame(index=adata.obs['pt_bins'].cat.categories, columns=gene_names)
    
    for bin_cat in adata.obs['pt_bins'].cat.categories:
        # 해당 구간에 속하는 세포 인덱스 가져오기
        cell_indices = np.where(adata.obs['pt_bins'] == bin_cat)[0]
        
        # 해당 세포들의 평균 발현량 계산
        if len(cell_indices) > 0:
            bin_expr_mean = np.mean(expression_matrix[cell_indices, :], axis=0)
            bin_means.loc[bin_cat] = bin_expr_mean
    
    # 구간별 발현 변화율 계산
    expression_change_rate = pd.DataFrame()
    for i in range(1, bin_means.shape[0]):
        # 이전 구간과 현재 구간의 발현량 차이
        change = bin_means.iloc[i] - bin_means.iloc[i-1]
        # 변화율 계산
        change_rate = change / (bin_means.iloc[i-1] + 1e-10)  # 0으로 나누기 방지
        expression_change_rate[f'bin_{i-1}_to_{i}'] = change_rate
    
    # 각 구간에서의 변화율 합계 계산
    total_change_per_bin = expression_change_rate.abs().sum()
    
    # 변화율이 가장 큰 구간 식별
    transition_points = total_change_per_bin.sort_values(ascending=False).head(3)
    print("주요 세포 상태 전환 지점:")
    for bin_name, change_value in transition_points.items():
        bin_idx = int(bin_name.split('_')[1])
        start_pt = bin_means.index[bin_idx].left
        end_pt = bin_means.index[bin_idx].right
        print(f"  - 의사시간 {start_pt:.3f}~{end_pt:.3f} 구간 (변화율: {change_value:.2f})")
    
    # 세포 상태 전환 시각화
    plt.figure(figsize=(12, 6))
    plt.bar(total_change_per_bin.index, total_change_per_bin.values)
    plt.xticks(rotation=45)
    plt.title('의사시간 구간별 유전자 발현 변화율')
    plt.xlabel('의사시간 구간')
    plt.ylabel('유전자 발현 변화 총량')
    plt.tight_layout()
    plt.savefig(os.path.join(advanced_dir, 'expression_change_rate.pdf'))
    plt.savefig(os.path.join(advanced_dir, 'expression_change_rate.png'), dpi=300)
    plt.close()
    
    # 5. 전환 지점에서 가장 크게 변화하는 유전자 식별
    for transition_bin in transition_points.index:
        bin_idx = int(transition_bin.split('_')[1])
        
        # 해당 구간에서 변화량이 큰 유전자 선택
        genes_change = expression_change_rate[transition_bin].abs().sort_values(ascending=False)
        top_changing_genes = genes_change.head(20)
        
        # 결과 저장
        pd.DataFrame({
            'gene': top_changing_genes.index,
            'change_rate': top_changing_genes.values,
            'change_direction': np.where(expression_change_rate.loc[top_changing_genes.index, transition_bin] > 0, 'up', 'down')
        }).to_csv(os.path.join(output_dir, f'transition_genes_{transition_bin}.csv'), index=False)
        
        # 전환 지점 유전자 시각화
        plt.figure(figsize=(10, 8))
        
        # 각 유전자의 의사시간에 따른 발현 곡선
        for gene in top_changing_genes.head(10).index:
            # 20개 구간의 중간점 계산
            bin_midpoints = [bin_idx.mid for bin_idx in bin_means.index]
            gene_expression = bin_means[gene].values
            
            # 곡선 그리기
            plt.plot(bin_midpoints, gene_expression, label=gene, marker='o', linewidth=2)
        
        # 전환 지점 표시
        transition_midpoint = (bin_means.index[bin_idx].left + bin_means.index[bin_idx].right) / 2
        plt.axvline(x=transition_midpoint, color='red', linestyle='--', alpha=0.5,
                    label=f'전환 지점 ({bin_means.index[bin_idx].left:.2f}~{bin_means.index[bin_idx].right:.2f})')
        
        plt.title(f'전환 지점 {transition_bin}에서 변화하는 상위 유전자')
        plt.xlabel('의사시간')
        plt.ylabel('유전자 발현 (평균)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(advanced_dir, f'transition_genes_{transition_bin}.pdf'))
        plt.savefig(os.path.join(advanced_dir, f'transition_genes_{transition_bin}.png'), dpi=300)
        plt.close()
    
    # 6. 분기점 별 세포 경로 특성화 (분기점이 있는 경우)
    if 'dpt_groups' in adata.obs:
        print("분기점 별 세포 경로 특성화 중...")
        
        # 분기점 별 세포 수 확인
        branch_counts = adata.obs['dpt_groups'].value_counts()
        print(f"분기점 별 세포 수: {branch_counts.to_dict()}")
        
        # 분기점 별 세포 유형 분포 시각화 (cell_type이 있는 경우)
        if 'cell_type' in adata.obs:
            branch_cell_types = pd.crosstab(adata.obs['dpt_groups'], adata.obs['cell_type'], normalize='index')
            
            plt.figure(figsize=(14, 8))
            branch_cell_types.plot(kind='bar', stacked=True, colormap='tab20')
            plt.title('분기점 별 세포 유형 분포')
            plt.xlabel('분기점 그룹')
            plt.ylabel('비율')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(advanced_dir, 'branch_cell_type_distribution.pdf'))
            plt.savefig(os.path.join(advanced_dir, 'branch_cell_type_distribution.png'), dpi=300)
            plt.close()
        
        # 분기점 시각화 (3D UMAP)
        try:
            if 'X_umap' in adata.obsm:
                # 3D UMAP 계산 (기존 2D UMAP 기반)
                sc.tl.umap(adata, n_components=3, init_pos='paga')
                
                # 3D 시각화
                fig = plt.figure(figsize=(12, 10))
                ax = fig.add_subplot(111, projection='3d')
                
                # 분기점 별로 다른 색상 적용
                for branch in adata.obs['dpt_groups'].unique():
                    mask = adata.obs['dpt_groups'] == branch
                    ax.scatter(
                        adata.obsm['X_umap'][mask, 0],
                        adata.obsm['X_umap'][mask, 1],
                        adata.obsm['X_umap'][mask, 2],
                        label=f'Branch {branch}',
                        alpha=0.7,
                        s=30
                    )
                
                ax.set_title('3D UMAP with Branching Points')
                ax.set_xlabel('UMAP1')
                ax.set_ylabel('UMAP2')
                ax.set_zlabel('UMAP3')
                ax.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(advanced_dir, 'umap3d_branches.pdf'))
                plt.savefig(os.path.join(advanced_dir, 'umap3d_branches.png'), dpi=300)
                plt.close()
        except Exception as e:
            print(f"3D UMAP 시각화 중 오류 발생: {e}")
    
    # 7. 특정 유전자 세트(마커/조절자)의 발현 패턴 분석 (사용자 정의 마커 사용)
    print("특정 유전자 세트의 발현 패턴 분석 중...")
    
    # 예시 마커 유전자 (실제 분석에서는 세포 유형에 맞게 조정 필요)
    # 줄기세포 관련 마커
    stem_markers = ['POU5F1', 'NANOG', 'SOX2', 'KLF4', 'MYC']
    
    # 분화 관련 마커 (다양한 세포 유형)
    diff_markers = ['PAX6', 'FOXG1', 'GATA6', 'GATA4', 'SOX17', 'PDX1', 'NKX2-1']
    
    # 전사 인자
    tf_markers = ['RUNX1', 'RUNX2', 'CEBPA', 'CEBPB', 'FOXA2', 'FOXO1', 'STAT3']
    
    # 발현 패턴 시각화
    all_markers = stem_markers + diff_markers + tf_markers
    
    # 데이터셋에 존재하는 마커만 선택
    existing_markers = [gene for gene in all_markers if gene in adata.var_names]
    
    if existing_markers:
        print(f"데이터셋에서 {len(existing_markers)}/{len(all_markers)} 마커 유전자를 찾았습니다.")
        
        # 의사시간에 따른 마커 발현 시각화
        fig = plt.figure(figsize=(15, 10))
        gs = GridSpec(3, 1, height_ratios=[1, 1, 1], figure=fig)
        
        # 줄기세포 마커
        stem_existing = [gene for gene in stem_markers if gene in adata.var_names]
        if stem_existing:
            ax1 = fig.add_subplot(gs[0])
            for gene in stem_existing:
                gene_idx = np.where(adata.var_names == gene)[0][0]
                x = adata.obs['dpt_pseudotime'].values
                y = expression_matrix[:, gene_idx]
                
                # 발현 값에 대한 스무딩 (이동 평균)
                df = pd.DataFrame({'x': x, 'y': y})
                df = df.sort_values('x')
                df['y_smooth'] = df['y'].rolling(window=50, min_periods=10).mean()
                
                ax1.scatter(df['x'], df['y'], alpha=0.1, s=5)
                ax1.plot(df['x'], df['y_smooth'], linewidth=2, label=gene)
            
            ax1.set_title('줄기세포 마커 발현 패턴')
            ax1.set_xlabel('의사시간')
            ax1.set_ylabel('발현량')
            ax1.legend()
        
        # 분화 마커
        diff_existing = [gene for gene in diff_markers if gene in adata.var_names]
        if diff_existing:
            ax2 = fig.add_subplot(gs[1])
            for gene in diff_existing:
                gene_idx = np.where(adata.var_names == gene)[0][0]
                x = adata.obs['dpt_pseudotime'].values
                y = expression_matrix[:, gene_idx]
                
                # 발현 값에 대한 스무딩 (이동 평균)
                df = pd.DataFrame({'x': x, 'y': y})
                df = df.sort_values('x')
                df['y_smooth'] = df['y'].rolling(window=50, min_periods=10).mean()
                
                ax2.scatter(df['x'], df['y'], alpha=0.1, s=5)
                ax2.plot(df['x'], df['y_smooth'], linewidth=2, label=gene)
            
            ax2.set_title('분화 마커 발현 패턴')
            ax2.set_xlabel('의사시간')
            ax2.set_ylabel('발현량')
            ax2.legend()
        
        # 전사 인자
        tf_existing = [gene for gene in tf_markers if gene in adata.var_names]
        if tf_existing:
            ax3 = fig.add_subplot(gs[2])
            for gene in tf_existing:
                gene_idx = np.where(adata.var_names == gene)[0][0]
                x = adata.obs['dpt_pseudotime'].values
                y = expression_matrix[:, gene_idx]
                
                # 발현 값에 대한 스무딩 (이동 평균)
                df = pd.DataFrame({'x': x, 'y': y})
                df = df.sort_values('x')
                df['y_smooth'] = df['y'].rolling(window=50, min_periods=10).mean()
                
                ax3.scatter(df['x'], df['y'], alpha=0.1, s=5)
                ax3.plot(df['x'], df['y_smooth'], linewidth=2, label=gene)
            
            ax3.set_title('전사 인자 발현 패턴')
            ax3.set_xlabel('의사시간')
            ax3.set_ylabel('발현량')
            ax3.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(advanced_dir, 'marker_genes_expression.pdf'))
        plt.savefig(os.path.join(advanced_dir, 'marker_genes_expression.png'), dpi=300)
        plt.close()
    else:
        print("데이터셋에서 지정된 마커 유전자를 찾을 수 없습니다.")
    
    # 8. 의사시간에 따른 주성분 분석 (PCA)
    print("의사시간에 따른 주성분 분석 중...")
    
    # 의사시간 구간 (10개)
    pt_bins = 10
    adata.obs['pt_bins_pca'] = pd.cut(adata.obs['dpt_pseudotime'], bins=pt_bins, labels=False)
    
    # 구간별 데이터 추출 및 PCA 수행
    pca = PCA(n_components=2)
    
    # 결과 저장용 데이터프레임
    pca_df = pd.DataFrame(columns=['PC1', 'PC2', 'Pseudotime_Bin'])
    
    for bin_idx in range(pt_bins):
        bin_cells = adata.obs['pt_bins_pca'] == bin_idx
        if sum(bin_cells) > 10:  # 최소 10개 이상의 세포가 있는 경우만 계산
            bin_pca = pca.fit_transform(expression_matrix[bin_cells])
            
            # 결과 저장
            bin_pca_df = pd.DataFrame({
                'PC1': bin_pca[:, 0],
                'PC2': bin_pca[:, 1],
                'Pseudotime_Bin': bin_idx
            })
            pca_df = pd.concat([pca_df, bin_pca_df], ignore_index=True)
    
    # 의사시간 구간별 PCA 시각화
    plt.figure(figsize=(10, 8))
    
    # 구간별로 다른 색상 적용
    cmap = plt.cm.viridis
    norm = plt.Normalize(pca_df['Pseudotime_Bin'].min(), pca_df['Pseudotime_Bin'].max())
    
    for bin_idx in sorted(pca_df['Pseudotime_Bin'].unique()):
        bin_data = pca_df[pca_df['Pseudotime_Bin'] == bin_idx]
        
        plt.scatter(
            bin_data['PC1'], 
            bin_data['PC2'],
            c=cmap(norm(bin_idx)),
            label=f'Bin {bin_idx}',
            alpha=0.7,
            s=50
        )
    
    plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), label='의사시간 구간')
    plt.title('의사시간 구간별 유전자 발현 PCA')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.tight_layout()
    plt.savefig(os.path.join(advanced_dir, 'pseudotime_pca.pdf'))
    plt.savefig(os.path.join(advanced_dir, 'pseudotime_pca.png'), dpi=300)
    plt.close()
    
    # 9. 궤적 시각화 개선 (PAGA + 의사시간)
    if 'paga' in adata.uns and 'connectivities' in adata.uns['paga']:
        print("개선된 PAGA 시각화 생성 중...")
        
        # PAGA 연결성 가져오기
        paga_connectivities = adata.uns['paga']['connectivities'].toarray()
        
        # 강한 연결만 유지
        threshold = 0.05
        paga_strong = paga_connectivities.copy()
        paga_strong[paga_strong < threshold] = 0
        
        # 클러스터별 의사시간 평균 계산
        if 'leiden' in adata.obs:
            cluster_pt = adata.obs.groupby('leiden')['dpt_pseudotime'].mean().sort_values()
            
            # PAGA 그래프와 의사시간을 결합한 시각화
            plt.figure(figsize=(12, 10))
            
            # 노드 위치 정의 (원형 배치)
            n_clusters = len(cluster_pt)
            angles = np.linspace(0, 2*np.pi, n_clusters, endpoint=False)
            
            radius = 5
            pos = {cluster: (radius * np.cos(angle), radius * np.sin(angle)) 
                  for cluster, angle in zip(cluster_pt.index, angles)}
            
            # 노드 크기 계산 (클러스터 크기에 비례)
            cluster_sizes = adata.obs['leiden'].value_counts()
            node_size_scale = 2000  # 스케일 조정
            node_sizes = {cluster: np.sqrt(size) * node_size_scale / np.sqrt(cluster_sizes.max()) 
                         for cluster, size in cluster_sizes.items()}
            
            # 노드 그리기
            for cluster in cluster_pt.index:
                plt.scatter(
                    pos[cluster][0], 
                    pos[cluster][1],
                    s=node_sizes[cluster],
                    c=[plt.cm.viridis(cluster_pt[cluster] / cluster_pt.max())],
                    edgecolors='black',
                    linewidth=1,
                    zorder=3
                )
                
                # 클러스터 라벨 추가
                plt.text(
                    pos[cluster][0], 
                    pos[cluster][1],
                    cluster,
                    ha='center', 
                    va='center',
                    fontsize=12,
                    fontweight='bold',
                    color='black',
                    zorder=4
                )
            
            # 엣지 그리기
            for i, source in enumerate(cluster_pt.index):
                for j, target in enumerate(cluster_pt.index):
                    if i != j and paga_strong[int(source), int(target)] > 0:
                        # 연결 강도에 비례하는 선 굵기
                        linewidth = paga_strong[int(source), int(target)] * 10
                        
                        # 의사시간에 따른 방향성 화살표
                        if cluster_pt[source] <= cluster_pt[target]:
                            # 순방향 (낮은 의사시간에서 높은 의사시간으로)
                            plt.annotate(
                                "",
                                xy=pos[target],
                                xytext=pos[source],
                                arrowprops=dict(
                                    arrowstyle="-|>",
                                    lw=linewidth,
                                    color='gray',
                                    alpha=0.7,
                                    connectionstyle="arc3,rad=0.1"
                                ),
                                zorder=2
                            )
            
            # 의사시간 컬러바 추가
            sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(0, cluster_pt.max()))
            sm.set_array([])
            cbar = plt.colorbar(sm, label='의사시간')
            
            plt.title('PAGA 그래프 + 의사시간 (화살표는 분화 방향)')
            plt.axis('equal')
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(os.path.join(advanced_dir, 'enhanced_paga_pseudotime.pdf'))
            plt.savefig(os.path.join(advanced_dir, 'enhanced_paga_pseudotime.png'), dpi=300)
            plt.close()
    
    print(f"\n======== 심화 분석 완료 ========")
    print(f"심화 분석 결과가 다음 위치에 저장되었습니다: {advanced_dir}")

if __name__ == "__main__":
    run_advanced_analysis()


"""
<think>사용자가 요청한 것은 세포 궤적 분석 결과에 대한 추가 분석 코드입니다. 기존 코드를 참고하여 추가 분석을 위한 코드를 생성하겠습니다. 추가 분석으로는 다음과 같은 내용을 포함할 수 있습니다:

1. 의사시간과 상관관계가 높은 유전자 식별 및 GO 분석
2. 특정 분화 단계에서의 전사 인자 활성화 분석
3. 의사시간에 따른 세포 상태 전환 지점 식별
4. 궤적 분석 결과 시각화 개선 및 추가 시각화
5. 알려진 마커 유전자와의 비교

이러한 분석을 수행하는 코드를 작성하겠습니다. 기존 코드에 이어 실행할 수 있는 별도의 스크립트로 작성하겠습니다.</think># 세포 궤적 분석을 위한 추가 분석 코드

아래 코드는 이전에 수행한 궤적 분석 결과를 바탕으로 추가 분석을 수행합니다. 새 파일(예: `sc_trajectory_advanced.py`)로 저장하여 사용하세요.


## 설명

이 코드는 이전에 실행한 `sc_trajectory.py`의 결과를 가져와 다음과 같은 심화 분석을 수행합니다:

1. **의사시간-유전자 발현 상관관계 분석**:
   - 각 유전자와 의사시간 간의 스피어만 상관계수 계산
   - 양의 상관관계(의사시간에 따라 증가) 및 음의 상관관계(의사시간에 따라 감소) 유전자 식별
   - 가장 유의미한 상관관계를 가진 유전자들의 발현 패턴 시각화

2. **GO 분석**:
   - 상관관계가 높은 유전자들의 생물학적 기능 분석
   - 의사시간에 따라 변화하는 생물학적 과정 이해 

3. **세포 상태 전환점 분석**:
   - 의사시간에 따른 유전자 발현 변화율 계산
   - 급격한 발현 변화가 일어나는 전환점 식별
   - 전환점에서 변화하는 주요 유전자 식별

4. **분기점별 세포 특성화**:
   - 서로 다른 분화 경로 간의 차이점 분석
   - 분기점별 세포 유형 분포 비교
   - 3D UMAP을 통한 분기점 시각화

5. **마커 유전자 분석**:
   - 알려진 마커 유전자의 발현 패턴 분석
   - 줄기세포 마커, 분화 마커, 전사인자 등의 발현 변화 추적

6. **개선된 시각화**:
   - 의사시간 구간별 주성분 분석(PCA)
   - 강화된 PAGA 시각화 (의사시간과 분화 방향 표시)

이 코드를 실행하면 보다 심층적인 세포 궤적 분석 결과를 얻을 수 있으며, 다양한 시각화 자료와 함께 분화 과정의 상세한 이해가 가능합니다.
"""