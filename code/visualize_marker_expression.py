import scanpy as sc
import matplotlib.pyplot as plt

# 결과 파일 로드 (경로를 환경에 맞게 수정하세요)
adata = sc.read_h5ad("/home/swr0460/daeseong/code/results.h5ad")

# 주요 마커 유전자 목록 (필요에 따라 수정하세요)
marker_genes = ["CD3D", "NKG7", "MS4A1"]

# UMAP에 각 마커 유전자의 발현 패턴 시각화
sc.pl.umap(adata, color=marker_genes, save="_marker_expression.pdf", show=True) 