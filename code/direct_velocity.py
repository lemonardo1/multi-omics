#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
사용법:
python /home/swr0460/daeseong/code/direct_velocity.py --input /home/swr0460/daeseong/code/results.h5ad

direct_velocity.py: scvelo를 사용하여 10x Genomics 데이터에서 직접 RNA velocity 계산

이 스크립트는 10x Genomics 데이터를 읽고 전처리하여 RNA velocity 분석을 수행합니다.
velocyto 전처리 단계 없이 scvelo를 사용하여 RNA velocity 계산을 시도합니다.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import scvelo as scv
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
    parser = argparse.ArgumentParser(description='직접 RNA velocity 계산 스크립트')
    
    parser.add_argument('--input', type=str, required=True, 
                        help='입력 AnnData 객체 (.h5ad 파일) 또는 10x 데이터 디렉토리 경로')
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
    """데이터 로드"""
    print(f"[INFO] {args.input} 데이터 로드 중...")
    
    if args.input.endswith('.h5ad'):
        # .h5ad 파일 로드
        adata = sc.read_h5ad(args.input)
    else:
        # 10x Genomics 데이터 로드
        try:
            adata = sc.read_10x_mtx(args.input)
            print(f"[INFO] 10x 데이터 로드 완료: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
            
            # 기본 전처리
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            # PCA 및 Neighbors
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            sc.pp.pca(adata, n_comps=30)
            sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
            
            # UMAP 및 클러스터링 
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.8)
            
            print(f"[INFO] 기본 전처리 및 클러스터링 완료")
        except Exception as e:
            print(f"[ERROR] 10x 데이터 로드 실패: {e}")
            sys.exit(1)
    
    return adata

def process_velocity(adata, mode='stochastic', n_jobs=8):
    """scvelo를 사용한 RNA velocity 처리"""
    print(f"[INFO] RNA velocity 전처리 및 계산 중...")
    
    # spliced/unspliced 레이어 확인
    has_velocity_data = 'spliced' in adata.layers and 'unspliced' in adata.layers
    
    if not has_velocity_data:
        print("[INFO] spliced/unspliced 레이어가 없습니다. 데이터에서 추정을 시도합니다.")
        try:
            # 방법 1: 기본 레이어에서 동역학 복구 시도
            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            print(f"[INFO] 모멘트 계산 완료")
            
            # 동역학 모델 복구 시도
            scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
            print(f"[INFO] 동역학 모델 복구 완료")
            
            # RNA velocity 계산 (dynamical 모드)
            scv.tl.velocity(adata, mode='dynamical')
            print(f"[INFO] dynamical 모드로 RNA velocity 계산 완료")
            
        except Exception as e:
            print(f"[WARNING] 동역학 모델 복구 실패: {e}")
            print("[INFO] spliced/unspliced 데이터 생성을 시도합니다...")
            
            # 방법 2: 기존 데이터로 간단한 velocity 추정
            # PCA 기반 속도 계산 (La Manno et al. 방법)
            try:
                # neighbors가 없는 경우 계산
                if 'neighbors' not in adata.uns:
                    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
                
                # 속도 그래프 직접 계산
                scv.tl.velocity_pseudotime(adata)
                scv.tl.velocity_clusters(adata)
                scv.tl.rank_velocity_genes(adata, n_genes=10)
                
                print(f"[INFO] 속도 기반 유사시간 및 클러스터 계산 완료")
                return adata
                
            except Exception as e2:
                print(f"[ERROR] RNA velocity 계산 실패: {e2}")
                print("[INFO] 원시 데이터(10x)로 다시 시도하거나 velocyto 전처리를 수행하세요.")
                return adata
    else:
        # 기존 코드: spliced/unspliced 레이어가 있는 경우
        # RNA velocity 전처리
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        print(f"[INFO] 필터링 및 정규화 완료")
        
        # 로그 변환 및 모멘트 계산
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        print(f"[INFO] 모멘트 계산 완료")
    
    try:
        # 동역학 모델 복구 (dynamical 모드에만 필요)
        if mode == 'dynamical' and not adata.uns.get('recover_dynamics', False):
            scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
            print(f"[INFO] 동역학 모델 복구 완료")
        
        # RNA velocity 계산
        scv.tl.velocity(adata, mode=mode)
        print(f"[INFO] {mode} 모드로 RNA velocity 계산 완료")
        
        # Velocity 그래프 및 임베딩 속도 계산
        scv.tl.velocity_graph(adata)
        print(f"[INFO] Velocity 그래프 계산 완료")
        
        return adata
    
    except Exception as e:
        print(f"[ERROR] RNA velocity 계산 중 오류: {e}")
        
        # 가능한 경우 이전 단계까지의 결과 반환
        return adata

def visualize_velocity(adata, embed, cluster_key, output_dir, figures_dir):
    """RNA velocity 시각화"""
    print("[INFO] RNA velocity 시각화 중...")
    
    # 결과 경로 설정
    velocity_embedding_path = os.path.join(figures_dir, f"velocity_embedding_{embed}.png")
    velocity_stream_path = os.path.join(figures_dir, f"velocity_stream_{embed}.png")
    
    # velocity 계산 여부 확인
    has_velocity = 'velocity' in adata.uns
    has_velocity_graph = 'velocity_graph' in adata.uns
    
    try:
        # 기본 UMAP/임베딩 플롯
        sc.pl.embedding(adata, basis=embed, color=cluster_key, 
                       title=f'세포 클러스터 ({cluster_key})',
                       save=os.path.join(figures_dir, f"clusters_{embed}.png"))
        
        if has_velocity and has_velocity_graph:
            # 속도 임베딩 플롯
            scv.pl.velocity_embedding(adata, basis=embed, color=cluster_key, 
                                    title='RNA Velocity (arrows)',
                                    save=velocity_embedding_path)
            
            # 속도 스트림 플롯
            scv.pl.velocity_embedding_stream(adata, basis=embed, color=cluster_key, 
                                          title='RNA Velocity Stream',
                                          save=velocity_stream_path)
            
            # 유전자 레벨 속도 (상위 6개 유전자)
            if 'velocity_genes' in adata.var.columns:
                top_genes = adata.var['velocity_genes'].sort_values(ascending=False).index[:6]
                scv.pl.velocity(adata, top_genes, ncols=2, 
                              save=os.path.join(figures_dir, "top_velocity_genes.png"))
                
            print(f"[INFO] Velocity 시각화 완료.")
        else:
            # velocity 없는 경우 기본 시각화
            if 'velocity_pseudotime' in adata.obs:
                sc.pl.embedding(adata, basis=embed, color='velocity_pseudotime', 
                               color_map='viridis', 
                               title='Velocity Pseudotime',
                               save=os.path.join(figures_dir, "velocity_pseudotime.png"))
            
            sc.pl.embedding(adata, basis=embed, color=[cluster_key], 
                          title=f'세포 군집화 ({cluster_key})',
                          save=os.path.join(figures_dir, "clusters.png"))
            
            print(f"[WARNING] velocity 그래프가 생성되지 않아 기본 시각화만 수행합니다.")
        
        print(f"[INFO] 시각화 완료. 결과는 {figures_dir}에 저장되었습니다.")
    
    except Exception as e:
        print(f"[ERROR] 시각화 중 오류: {e}")
        print(f"[INFO] 시각화를 건너뜁니다.")

def main():
    """메인 실행 함수"""
    # 인수 파싱
    args = parse_arguments()
    
    # 출력 디렉토리 설정
    output_dir, figures_dir = setup_directories()
    
    # 데이터 로드
    adata = load_data(args)
    print(f"[INFO] 데이터 로드 완료: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
    
    # RNA velocity 계산
    adata = process_velocity(adata, mode=args.mode, n_jobs=args.n_jobs)
    
    # 시각화
    visualize_velocity(adata, args.embed, args.cluster_key, output_dir, figures_dir)
    
    # 결과 저장
    result_path = os.path.join(output_dir, "velocity_results.h5ad")
    adata.write(result_path)
    print(f"[INFO] 분석 결과가 {result_path}에 저장되었습니다.")

if __name__ == "__main__":
    main() 