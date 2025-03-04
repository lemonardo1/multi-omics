import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# 저장된 AnnData 결과 파일 불러오기
adata = sc.read('code/results.h5ad')

# 'leiden' 클러스터를 기준으로 차별 발현 분석 수행 (t-test 사용)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# 분석 결과를 시각화하여 플롯 저장
# sharey=False: 모든 플롯에 대해 동일한 y축 범위 사용   
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_diff_expression.pdf')

# 차별 발현 분석 결과를 CSV 파일로 저장하는 함수
def save_diff_expression(adata, output_prefix='diff_expression'):
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

save_diff_expression(adata)

print('차별 발현 분석 수행 및 결과 저장이 완료되었습니다!') 