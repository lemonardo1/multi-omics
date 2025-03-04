# 필요한 라이브러리 임포트
import scanpy as sc  # 단일세포 RNA-seq 분석을 위한 주요 라이브러리
import pandas as pd  # 데이터 처리를 위한 라이브러리
import numpy as np   # 수치 계산을 위한 라이브러리
import matplotlib.pyplot as plt  # 시각화를 위한 라이브러리
import seaborn as sns  # 통계적 데이터 시각화를 위한 라이브러리
import os

# 설정: 기본 디렉토리
base_dir = '/home/swr0460/daeseong'
data_dir = os.path.join(base_dir, 'data')
code_dir = os.path.join(base_dir, 'code')
figures_dir = os.path.join(base_dir, 'figures')
output_dir = os.path.join(base_dir, 'output')

# Scanpy의 기본 figure 경로를 설정 (UMAP 저장 파일 등)
sc.settings.figdir = figures_dir

# 재현성을 위한 랜덤 시드 설정
np.random.seed(42)

# 10x Genomics 데이터 읽기
# 10x Genomics는 단일세포 RNA 시퀀싱에서 가장 널리 사용되는 플랫폼 중 하나입니다.
adata = sc.read_10x_mtx(
    os.path.join(data_dir, 'filtered_feature_bc_matrix/'),  # .mtx 파일이 있는 디렉토리
    var_names='gene_symbols',            # 변수명으로 유전자 심볼 사용
    cache=True                           # 데이터 캐싱 활성화
)

# 기본적인 전처리 단계
# 최소 200개의 유전자가 발현된 세포만 선택
sc.pp.filter_cells(adata, min_genes=200)
# 최소 3개의 세포에서 발현된 유전자만 선택
sc.pp.filter_genes(adata, min_cells=3)

# 품질 지표 계산
# MT-: 미토콘드리아 유전자를 나타내며, 세포 스트레스나 사멸의 지표로 사용됨
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 미토콘드리아 유전자 주석
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 품질 지표 시각화
# 세포당 유전자 수와 미토콘드리아 비율의 분포를 확인
fig, axs = plt.subplots(1, 2, figsize=(15, 5))
sns.histplot(data=adata.obs, x='n_genes_by_counts', bins=60, ax=axs[0])
axs[0].set_title('세포당 유전자 수 분포')
sns.histplot(data=adata.obs, x='pct_counts_mt', bins=60, ax=axs[1])
axs[1].set_title('미토콘드리아 비율 분포')

os.makedirs(figures_dir, exist_ok=True)
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'new_qc_plots.pdf'))
plt.close()

# 품질 기준에 따른 세포 필터링
# 최소 200개 유전자를 가진 세포 선택
sc.pp.filter_cells(adata, min_genes=200)
# 너무 많은 유전자(2500개 이상)를 가진 세포는 제외 (더블릿 가능성)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# 미토콘드리아 비율이 5% 이상인 세포 제외 (사멸 세포 가능성)
adata = adata[adata.obs.pct_counts_mt < 5, :]

# 데이터 정규화
# 각 세포의 총 카운트를 10000으로 정규화
sc.pp.normalize_total(adata, target_sum=1e4)
# 로그 변환 (log(x + 1)) - 데이터의 분포를 더 정규분포에 가깝게 만듦
sc.pp.log1p(adata)

# 고변동 유전자 선택
# 세포 타입을 구분하는데 가장 유용한 유전자들을 선택
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 데이터 스케일링
# 각 유전자의 평균을 0, 분산을 1로 조정 (이상치의 영향을 줄이기 위해 최대값 10으로 제한)
sc.pp.scale(adata, max_value=10)

# PCA (주성분 분석) 수행
# 차원 축소를 통해 데이터의 주요 변동성을 포착
sc.tl.pca(adata, svd_solver='arpack')

# 이웃 그래프 계산
# 세포들 간의 유사성을 기반으로 그래프 구성
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP 차원 축소
# 고차원 데이터를 2차원으로 시각화하기 위한 비선형 차원 축소
sc.tl.umap(adata)

# Leiden 클러스터링
# 세포들을 비슷한 발현 패턴을 가진 그룹으로 분류
sc.tl.leiden(adata)

# UMAP 시각화
# 클러스터링 결과를 2차원 평면에 표시
sc.pl.umap(adata, color=['leiden'], save='new_clusters.pdf')

# 결과 저장
# 분석된 데이터를 h5ad 형식으로 저장
os.makedirs(code_dir, exist_ok=True)
adata.write(os.path.join(code_dir, 'new_results.h5ad'))

print("분석이 완료되었습니다! 결과가 '{}'에 저장되었습니다.".format(os.path.join(code_dir, 'new_results.h5ad'))) 
