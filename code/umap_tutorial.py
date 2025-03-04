import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import load_digits, load_iris
from umap import UMAP
from sklearn.preprocessing import StandardScaler

# 재현성을 위한 랜덤 시드 설정
np.random.seed(42)

def load_and_prepare_data(dataset='digits'):
    """
    데이터 로드 및 전처리 함수
    """
    if dataset == 'digits':
        # MNIST 손글씨 데이터 로드 (간단한 버전)
        digits = load_digits()
        X = digits.data
        y = digits.target
        data_name = "Digits"
    else:
        # iris 데이터 로드 (더 간단한 데이터셋)
        iris = load_iris()
        X = iris.data
        y = iris.target
        data_name = "Iris"
    
    # 데이터 스케일링
    X_scaled = StandardScaler().fit_transform(X)
    
    return X_scaled, y, data_name

def explore_umap_parameters(X, y, data_name):
    """
    다양한 UMAP 파라미터 탐색
    """
    # 다양한 파라미터 조합
    n_neighbors_list = [5, 15, 50]
    min_dist_list = [0.1, 0.5, 0.8]
    
    fig, axes = plt.subplots(len(n_neighbors_list), len(min_dist_list), 
                            figsize=(15, 15))
    
    for i, n_neighbors in enumerate(n_neighbors_list):
        for j, min_dist in enumerate(min_dist_list):
            # UMAP 모델 생성 및 학습
            umap_model = UMAP(n_neighbors=n_neighbors, 
                            min_dist=min_dist,
                            random_state=42)
            
            # 데이터 변환
            embedding = umap_model.fit_transform(X)
            
            # 결과 시각화
            axes[i,j].scatter(embedding[:, 0], embedding[:, 1], 
                            c=y, cmap='Spectral', s=5)
            axes[i,j].set_title(f'n_neighbors={n_neighbors}, min_dist={min_dist}')
            axes[i,j].set_xticks([])
            axes[i,j].set_yticks([])
    
    plt.suptitle(f'UMAP Parameter Exploration - {data_name} Dataset', 
                 fontsize=16, y=0.95)
    plt.tight_layout()
    plt.show()

def compare_with_different_metrics(X, y, data_name):
    """
    다양한 거리 메트릭을 사용한 UMAP 비교
    """
    metrics = ['euclidean', 'manhattan', 'cosine']
    fig, axes = plt.subplots(1, len(metrics), figsize=(15, 5))
    
    for i, metric in enumerate(metrics):
        umap_model = UMAP(metric=metric, random_state=42)
        embedding = umap_model.fit_transform(X)
        
        axes[i].scatter(embedding[:, 0], embedding[:, 1], 
                       c=y, cmap='Spectral', s=5)
        axes[i].set_title(f'Metric: {metric}')
        axes[i].set_xticks([])
        axes[i].set_yticks([])
    
    plt.suptitle(f'UMAP with Different Metrics - {data_name} Dataset', 
                 fontsize=16, y=1.05)
    plt.tight_layout()
    plt.show()

def main():
    # 두 가지 데이터셋에 대해 실험
    for dataset in ['digits', 'iris']:
        print(f"\nAnalyzing {dataset} dataset...")
        
        # 데이터 로드 및 전처리
        X, y, data_name = load_and_prepare_data(dataset)
        
        # 파라미터 탐색
        print("Exploring different UMAP parameters...")
        explore_umap_parameters(X, y, data_name)
        
        # 거리 메트릭 비교
        print("Comparing different distance metrics...")
        compare_with_different_metrics(X, y, data_name)

if __name__ == "__main__":
    main() 