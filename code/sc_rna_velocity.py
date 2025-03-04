#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 사용법:
# python sc_rna_velocity.py --input code/results.h5ad

"""
sc_rna_velocity.py: scRNA-seq 데이터를 이용한 RNA velocity 분석 스크립트

이 스크립트는 scVelo와 scanpy를 사용하여 단일세포 RNA-seq 데이터로부터 
RNA velocity를 계산하고 시각화합니다. 세포 분화의 방향성과 동역학을 추론하는 데 
사용됩니다.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import scvelo as scv
import cellrank as cr
from scipy import io
from scipy.sparse import csr_matrix
import argparse
from datetime import datetime

# 경고 메시지 무시
import warnings
warnings.filterwarnings('ignore')

# 시각화 설정
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120, facecolor='white', frameon=True)
scv.settings.set_figure_params(dpi=120, facecolor='white', frameon=True)
scv.settings.verbosity = 3
scv.settings.presenter_view = True

def setup_directories():
    """결과 저장을 위한 디렉토리 설정"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join("output", f"velocity_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)
    
    figures_dir = os.path.join("figures", f"velocity_{timestamp}")
    os.makedirs(figures_dir, exist_ok=True)
    
    return output_dir, figures_dir

def parse_arguments():
    """명령줄 인수 파싱"""
    parser = argparse.ArgumentParser(description='RNA velocity 분석 스크립트')
    
    parser.add_argument('--input', type=str, required=True, 
                        help='입력 AnnData 객체 (.h5ad 파일) 경로')
    parser.add_argument('--loom', type=str, default=None,
                        help='velocyto에서 생성된 .loom 파일 경로 (선택 사항)')
    parser.add_argument('--spliced_layer', type=str, default=None,
                        help='Spliced counts가 저장된 레이어명 (기본: X)')
    parser.add_argument('--unspliced_layer', type=str, default=None,
                        help='Unspliced counts가 저장된 레이어명 (기본: None)')
    parser.add_argument('--mode', type=str, default='stochastic',
                        choices=['deterministic', 'stochastic', 'dynamical'],
                        help='RNA velocity 모드 (deterministic, stochastic, dynamical)')
    parser.add_argument('--embed', type=str, default='umap',
                        help='사용할 임베딩 (umap, tsne, pca 등)')
    parser.add_argument('--cluster_key', type=str, default='leiden',
                        help='클러스터 정보가 있는 obs 컬럼명')
    parser.add_argument('--n_jobs', type=int, default=8,
                        help='사용할 CPU 코어 수')
    
    return parser.parse_args()

def load_data(args):
    """데이터 로드 및 전처리"""
    print(f"[INFO] {args.input} 파일 로드 중...")
    adata = sc.read_h5ad(args.input)
    
    # .loom 파일이 제공된 경우
    if args.loom:
        print(f"[INFO] {args.loom} loom 파일 로드 중...")
        ldata = scv.read(args.loom, cache=True)
        
        # 공통 바코드 및 유전자 찾기
        common_barcodes = np.intersect1d(adata.obs_names, ldata.obs_names)
        common_genes = np.intersect1d(adata.var_names, ldata.var_names)
        
        # 데이터 서브셋
        adata = adata[common_barcodes, common_genes].copy()
        ldata = ldata[common_barcodes, common_genes].copy()
        
        # spliced/unspliced 정보 추가
        adata.layers["spliced"] = ldata.layers["spliced"]
        adata.layers["unspliced"] = ldata.layers["unspliced"]
        
    # 레이어가 명시적으로 지정된 경우
    elif args.spliced_layer and args.unspliced_layer:
        print(f"[INFO] 지정된 레이어 사용: spliced={args.spliced_layer}, unspliced={args.unspliced_layer}")
        adata.layers["spliced"] = adata.layers[args.spliced_layer]
        adata.layers["unspliced"] = adata.layers[args.unspliced_layer]
    
    # 사전 정제된 counts가 없는 경우, 자동 추론
    else:
        print("[INFO] spliced/unspliced counts 추론 중...")
        scv.pp.filter_and_normalize(adata)
        scv.pp.moments(adata)
        try:
            scv.tl.recover_dynamics(adata)
        except Exception as e:
            print(f"[WARNING] 동적 모델 복구 중 오류: {e}")
    
    return adata

def preprocess_velocity(adata):
    """RNA velocity 분석을 위한 전처리"""
    print("[INFO] RNA velocity 전처리 중...")
    
    # 품질 관리 기본값 (필요시 조정)
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)
    
    # 모멘트 계산 (가중치 있는 PCA 기반 속도)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    
    return adata

def compute_velocity(adata, mode='stochastic'):
    """RNA velocity 계산"""
    print(f"[INFO] {mode} 모드로 RNA velocity 계산 중...")
    
    if mode == 'deterministic':
        scv.tl.velocity(adata, mode='deterministic')
    elif mode == 'stochastic':
        scv.tl.velocity(adata, mode='stochastic')
    elif mode == 'dynamical':
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')
    
    # Velocity 그래프 및 임베딩 속도 계산
    scv.tl.velocity_graph(adata)
    
    return adata

def visualize_velocity(adata, embed, cluster_key, output_dir, figures_dir):
    """RNA velocity 시각화"""
    print("[INFO] RNA velocity 시각화 중...")
    
    # 저장 경로 설정
    velocity_grid_path = os.path.join(figures_dir, f"velocity_grid_{embed}.png")
    velocity_stream_path = os.path.join(figures_dir, f"velocity_stream_{embed}.png")
    velocity_embedding_path = os.path.join(figures_dir, f"velocity_embedding_{embed}.png")
    
    # 속도 그리드 플롯
    scv.pl.velocity_embedding_grid(adata, basis=embed, color=cluster_key, 
                                  save=velocity_grid_path, title='RNA Velocity Grid')
    
    # 속도 스트림 플롯
    scv.pl.velocity_embedding_stream(adata, basis=embed, color=cluster_key, 
                                    save=velocity_stream_path, title='RNA Velocity Stream')
    
    # 속도 화살표 플롯
    scv.pl.velocity_embedding(adata, basis=embed, color=cluster_key, 
                             save=velocity_embedding_path, title='RNA Velocity Arrows')
    
    # 속도 그래프 신뢰도 분석
    scv.tl.velocity_confidence(adata)
    scv.pl.velocity_embedding_confidence(adata, basis=embed, 
                                        save=os.path.join(figures_dir, "velocity_confidence.png"))
    
    # 유전자 수준의 phase portrait
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:12]
    scv.pl.scatter(adata, basis=top_genes[:6], ncols=3, color=cluster_key, 
                  save=os.path.join(figures_dir, "velocity_phase_portrait1.png"))
    
    if len(top_genes) > 6:
        scv.pl.scatter(adata, basis=top_genes[6:12], ncols=3, color=cluster_key, 
                      save=os.path.join(figures_dir, "velocity_phase_portrait2.png"))

def perform_advanced_analysis(adata, embed, cluster_key, output_dir, figures_dir):
    """세포 궤적 및 상위 기여 유전자 등 고급 분석"""
    print("[INFO] 고급 velocity 분석 수행 중...")
    
    # PAGA 궤적 분석
    scv.tl.paga(adata, groups=cluster_key)
    scv.pl.paga(adata, basis=embed, size=50, alpha=0.1, min_edge_width=2, node_size_scale=1.5,
               save=os.path.join(figures_dir, "velocity_paga.png"))
    
    # 속도에 가장 크게 기여하는 유전자 식별
    scv.tl.rank_velocity_genes(adata, groupby=cluster_key, min_corr=0.3)
    
    # 상위 속도 유전자 시각화
    df = scv.get_df(adata, 'rank_velocity_genes/names')
    cluster_keys = list(df.keys())
    
    # 클러스터 그룹별로 상위 속도 유전자 저장
    for i, cluster in enumerate(cluster_keys):
        if i >= 6:  # 최대 6개 클러스터만 시각화
            break
        genes = df[cluster][:5]  # 각 클러스터의 상위 5개 유전자
        
        try:
            scv.pl.scatter(adata, genes, color=cluster_key, legend_loc='none', size=80,
                          save=os.path.join(figures_dir, f"top_velocity_genes_cluster{cluster}.png"))
        except Exception as e:
            print(f"[WARNING] 클러스터 {cluster}의 유전자 시각화 오류: {e}")
    
    # 속도 슈도시간 계산 (세포 발달 순서)
    try:
        scv.tl.velocity_pseudotime(adata)
        scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', 
                      save=os.path.join(figures_dir, "velocity_pseudotime.png"))
    except Exception as e:
        print(f"[WARNING] 슈도시간 계산 오류: {e}")
    
    # 잠재적 출발 세포 식별
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,
                  save=os.path.join(figures_dir, "latent_time.png"))
    
    # 속도 기반 마커 유전자 식별 및 차별 발현 분석
    top_velocity_genes = scv.tl.rank_velocity_genes(adata, groupby=cluster_key, n_genes=10)
    
    # 결과 저장
    adata.write(os.path.join(output_dir, "adata_with_velocity.h5ad"))
    pd.DataFrame(df).to_csv(os.path.join(output_dir, "velocity_genes.csv"))
    
    return adata

def run_cellrank_analysis(adata, cluster_key, output_dir, figures_dir):
    """CellRank를 사용한 세포 운명 예측"""
    try:
        print("[INFO] CellRank 분석 수행 중...")
        
        # CellRank 객체 초기화 (RNA velocity 기반)
        cr_knn = cr.tl.transition_matrix.KernelTransition(adata, weight_connectivities=0.2)
        cr_knn.compute_transition_matrix()
        
        # 마르코프 체인 초기화
        g_knn = cr.tl.estimators.GPCCA(cr_knn)
        
        # 고유값 분석
        g_knn.compute_eigendecomposition()
        cr.pl.eigendecomposition(g_knn, save=os.path.join(figures_dir, "cellrank_eigendecomposition.png"))
        
        # 매크로스테이트 식별
        g_knn.compute_macrostates(n_states=5)  # 상태 수는 데이터에 따라 조정 필요
        cr.pl.macrostates(g_knn, basis='umap', save=os.path.join(figures_dir, "cellrank_macrostates.png"))
        
        # 흡수 확률 및 운명 예측
        g_knn.compute_absorption_probabilities()
        cr.pl.absorption_probabilities(g_knn, basis='umap', save=os.path.join(figures_dir, "cellrank_absorption.png"))
        
        # 상위 운명 유전자 식별
        g_knn.compute_lineage_drivers()
        cr.pl.lineage_drivers(g_knn, save=os.path.join(figures_dir, "cellrank_lineage_drivers.png"))
        
        # 결과 저장
        with open(os.path.join(output_dir, "cellrank_results.txt"), "w") as f:
            f.write("CellRank Analysis Results\n")
            f.write("-------------------------\n")
            f.write(f"Number of terminal states: {len(g_knn.terminal_states_names)}\n")
            f.write(f"Terminal states: {g_knn.terminal_states_names}\n")
            
        return g_knn
    
    except Exception as e:
        print(f"[WARNING] CellRank 분석 오류: {e}")
        return None

def main():
    """메인 실행 함수"""
    # 인수 파싱
    args = parse_arguments()
    
    # 출력 디렉토리 설정
    output_dir, figures_dir = setup_directories()
    
    # 데이터 로드
    adata = load_data(args)
    print(f"[INFO] 데이터 로드 완료: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    
    # 전처리
    adata = preprocess_velocity(adata)
    
    # RNA velocity 계산
    adata = compute_velocity(adata, mode=args.mode)
    
    # 시각화
    visualize_velocity(adata, args.embed, args.cluster_key, output_dir, figures_dir)
    
    # 고급 분석
    adata = perform_advanced_analysis(adata, args.embed, args.cluster_key, output_dir, figures_dir)
    
    # CellRank 분석 (선택적)
    try:
        g_knn = run_cellrank_analysis(adata, args.cluster_key, output_dir, figures_dir)
    except:
        print("[WARNING] CellRank 분석을 건너뜁니다. (라이브러리 문제 또는 호환성 문제)")
    
    print(f"[INFO] 분석 완료!")
    print(f"[INFO] 결과는 {output_dir}에 저장되었습니다.")
    print(f"[INFO] 시각화는 {figures_dir}에 저장되었습니다.")

if __name__ == "__main__":
    main() 