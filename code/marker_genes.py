import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

# 저장된 AnnData 결과 파일 불러오기
adata = sc.read('code/results.h5ad')

# 클러스터 정보 확인 (예: 'leiden' 열 사용)
print('발견된 클러스터:', adata.obs['leiden'].unique())

# 마커 유전자 식별: 클러스터별로 통계적 유의성이 높은 유전자 찾기, wilcoxon 검정을 사용
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# 마커 유전자 결과 시각화: 각 클러스터 별 상위 25개 유전자 플롯 저장
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes.pdf')

def save_marker_genes(adata, output_prefix='marker_genes'):
    # output 디렉토리 생성 (없으면 생성)
    os.makedirs('output', exist_ok=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    for group in groups:
        df = pd.DataFrame({
            'gene': result['names'][group],
            'logfoldchanges': result['logfoldchanges'][group],
            'pvals': result['pvals'][group],
            'pvals_adj': result['pvals_adj'][group]
        })
        # CSV 파일을 output 디렉토리에 저장
        df.to_csv(f"output/{output_prefix}_{group}.csv", index=False)

save_marker_genes(adata)

print('마커 유전자 식별 및 저장이 완료되었습니다!') 