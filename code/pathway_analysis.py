"""
Pathway Analysis, Gene Ontology, 그리고 세포 상태 및 분화 경로 분석 코드

- 기존 분석 결과 (new_results.h5ad)를 기반으로 클러스터별 marker gene을 산출합니다.
- gseapy의 enrichr를 이용해 각 클러스터에 대한 GO & KEGG pathway enrichment 분석을 수행합니다.
- Diffusion pseudotime (DPT) 분석을 수행하여 세포 분화 경로를 시각화합니다.
- 분석 결과를 멋지게 시각화 (dotplot, UMAP, pseudotime plot) 합니다.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# gseapy: enrichr을 이용한 pathway/GO 분석 (pip install gseapy)
try:
    import gseapy as gp
except ImportError:
    raise ImportError("gseapy 모듈이 필요합니다. 'pip install gseapy'로 설치하세요.")

# 경로 설정: 기본 경로 /home/swr0460/daeseong/
base_dir = '/home/swr0460/daeseong'
data_dir = os.path.join(base_dir, 'data')
code_dir = os.path.join(base_dir, 'code')
figures_dir = os.path.join(base_dir, 'figures')
output_dir = os.path.join(base_dir, 'output')

# 결과 저장 디렉토리 확인 및 생성
os.makedirs(figures_dir, exist_ok=True)
os.makedirs(code_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Scanpy의 기본 figure 경로 설정 (UMAP 등 저장 위치)
sc.settings.figdir = figures_dir

# 분석 결과 파일 불러오기
adata_path = os.path.join(code_dir, 'new_results.h5ad')
adata = sc.read_h5ad(adata_path)
print("분석 결과를 불러왔습니다:", adata_path)

# 클러스터별 marker gene 산출 (클러스터 정보는 'leiden'으로 저장되어 있다고 가정)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', n_genes=100)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
plt.savefig(os.path.join(figures_dir, "rank_genes_groups.pdf"))
plt.close()
print("클러스터별 marker gene 분석 및 시각화를 저장했습니다.")

# 각 클러스터별 상위 50개 marker gene 리스트 준비 후 Enrichr 분석 수행
# 분석에는 GO Biological Process와 KEGG 2019 Human 라이브러리 사용
cluster_ids = adata.obs['leiden'].cat.categories

enrich_results = {}
for cluster in cluster_ids:
    # marker gene DataFrame 추출
    genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).loc[:, cluster].tolist()
    top_genes = list(pd.unique(genes))[:50]  # 중복 제거 후 상위 50개 선택
    
    print(f"클러스터 {cluster}의 상위 marker genes (총 {len(top_genes)}):", top_genes)
    
    # Enrichr 분석: GO Biological Process
    enr_go = gp.enrichr(gene_list=top_genes,
                        description=f"cluster_{cluster}_GO",
                        gene_sets='GO_Biological_Process_2018',
                        outdir=None,  # 저장하지 않고 바로 사용
                        cutoff=0.05)
    
    # KEGG pathway 분석
    enr_kegg = gp.enrichr(gene_list=top_genes,
                          description=f"cluster_{cluster}_KEGG",
                          gene_sets='KEGG_2019_Human',
                          outdir=None,
                          cutoff=0.05)
    
    enrich_results[cluster] = {'GO': enr_go.results, 'KEGG': enr_kegg.results}

    # dotplot 시각화 (상위 10개 항목)
    for key, df in zip(['GO', 'KEGG'], [enr_go.results, enr_kegg.results]):
        if not df.empty:
            df_plot = df.sort_values("Adjusted P-value").head(10)
            plt.figure(figsize=(8, 6))
            sns.scatterplot(data=df_plot, x='Combined Score', y=df_plot.index, size='Overlap',
                            sizes=(50, 300), hue='Adjusted P-value', palette="viridis", legend='brief')
            plt.title(f"Cluster {cluster} {key} Enrichment")
            plt.xlabel("Combined Score")
            plt.ylabel("Enriched Terms")
            plt.tight_layout()
            plot_filename = os.path.join(figures_dir, f"cluster_{cluster}_{key}_enrichment.pdf")
            plt.savefig(plot_filename)
            plt.close()
            print(f"클러스터 {cluster} {key} enrichment 시각화를 저장했습니다:", plot_filename)

# 세포 상태 및 분화 경로 분석: Diffusion pseudotime (DPT) & PAGA
# 이미 이웃 그래프와 UMAP이 계산되었으므로 해당 분석 사용
sc.tl.dpt(adata)
sc.pl.umap(adata, color=['dpt_pseudotime'], save='_dpt.pdf', show=False)
plt.close()

# PAGA 그래프 생성 및 시각화
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden'], save='_paga.pdf')
plt.close()

print("Diffusion pseudotime 및 PAGA 분석 결과를 시각화했습니다.")

# 분석 결과 저장 (Enrichment 결과는 CSV 파일로 저장)
for cluster in enrich_results:
    for lib in enrich_results[cluster]:
        df = enrich_results[cluster][lib]
        if not df.empty:
            output_file = os.path.join(output_dir, f"cluster_{cluster}_{lib}_enrichment.csv")
            df.to_csv(output_file, index=False)
            print(f"클러스터 {cluster} {lib} enrichment 결과 저장:", output_file)

print("전체 pathway 및 gene ontology, 세포 상태/분화 경로 분석을 완료하였습니다!")



"""
### 코드 설명

1. **경로 설정**  
   기본 디렉토리 \(`/home/swr0460/daeseong/`\)를 기준으로 데이터, 코드, figures, output 경로를 설정합니다.

2. **기존 분석 결과 불러오기**  
   `new_results.h5ad` 파일을 읽어와서 이후 분석에 사용합니다.

3. **클러스터별 marker gene 분석**  
   Scanpy의 `sc.tl.rank_genes_groups()`를 이용해 각 클러스터의 marker gene을 산출하고, 결과를 UMAP 및 PDF로 저장합니다.

4. **Enrichr를 이용한 Pathway/GO 분석**  
   각 클러스터의 상위 marker gene (중복 제거 후 상위 50개)을 대상으로, GO Biological Process와 KEGG 경로 분석을 수행합니다.  
   각 결과별 상위 10개 항목을 dotplot으로 시각화하여 figures 폴더에 저장합니다.

5. **세포 상태 및 분화 경로 분석**  
   Diffusion pseudotime (DPT) 분석과 PAGA를 사용하여 세포 분화 경로를 시각화합니다.

6. **분석 결과 저장**  
   Enrichment 분석 결과는 CSV 파일로 output 디렉토리에 저장됩니다.

이 코드를 실행하면 \(`/home/swr0460/daeseong/`\) 내의 각 디렉토리에 분석 결과와 시각화 파일들이 저장됩니다.

"""