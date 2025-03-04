# RNA Velocity 분석 가이드

이 문서는 10x Genomics 단일세포 RNA-seq 데이터에서 RNA velocity 분석을 수행하는 방법을 설명합니다.

## 개요

RNA velocity는 세포의 미래 상태를 예측하여 발현 동역학과 전사체 상태 변화를 분석하는 방법입니다. 이 분석은 다음 단계로 구성됩니다:

1. **Velocyto 전처리**: BAM 파일에서 spliced/unspliced 카운트 추출
2. **Velocity 계산**: 추출된 카운트에서 RNA velocity 계산
3. **시각화 및 분석**: RNA velocity 결과 해석 및 시각화

## 1. Velocyto 전처리

### 필요 파일

- **BAM 파일**: 10x Genomics CellRanger의 `possorted_genome_bam.bam` 파일
- **GTF 파일**: 참조 유전체 어노테이션 파일 (예: `genes.gtf`)
- **마스크 파일**(선택사항): 반복 영역을 마스킹하기 위한 GTF 파일

### velocyto 설치

```bash
# conda 환경 사용 시
conda activate scRNA
pip install velocyto

# 필요한 의존성 설치
conda install -c conda-forge -y numpy scipy cython numba matplotlib scikit-learn h5py click
```

### velocyto 실행 (CLI)

```bash
# BAM 파일에서 직접 실행
velocyto run --samtools-memory=32G --samtools-threads=8 \
  -m /path/to/repeat_msk.gtf \
  /path/to/possorted_genome_bam.bam \
  /path/to/reference.gtf

# 10x Genomics 출력 디렉토리에서 실행
velocyto run10x -m /path/to/repeat_msk.gtf \
  /path/to/10x_genomics_output \
  /path/to/reference.gtf
```

### velocyto 스크립트 사용

제공된 스크립트를 사용하여 더 편리하게 velocyto 전처리를 수행할 수 있습니다:

```bash
# velocyto 전처리 스크립트 실행
python /home/swr0460/daeseong/code/velocyto_preprocessing.py \
  --bam /path/to/possorted_genome_bam.bam \
  --gtf /path/to/reference.gtf \
  --output_dir /path/to/output \
  --mask_file /path/to/repeat_msk.gtf \
  --convert_to_h5ad
```

## 2. loom 파일을 h5ad로 변환

velocyto는 결과를 `.loom` 형식으로 저장합니다. 이를 `.h5ad` 형식으로 변환할 수 있습니다:

```bash
# loom 파일을 h5ad로 변환
python /home/swr0460/daeseong/code/velocyto_loom_to_h5ad.py \
  --loom /path/to/velocyto_output.loom \
  --output /path/to/output.h5ad
```

## 3. RNA Velocity 분석

변환된 `.h5ad` 파일을 사용하여 RNA velocity 분석을 수행할 수 있습니다:

```bash
# RNA velocity 분석 실행
python /home/swr0460/daeseong/code/direct_velocity.py \
  --input /path/to/velocyto_output_velocity.h5ad \
  --mode dynamical
```

## Python에서 직접 분석

Python에서 직접 분석할 경우의 일반적인 워크플로우:

```python
import scvelo as scv
import scanpy as sc

# loom 파일 로드
adata = scv.read_loom('/path/to/velocyto_output.loom')

# 또는 변환된 h5ad 파일 로드
# adata = scv.read_h5ad('/path/to/velocyto_output_velocity.h5ad')

# 필터링 및 정규화
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# 모멘트 계산
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# RNA velocity 계산 (stochastic, deterministic, dynamical 중 선택)
scv.tl.recover_dynamics(adata)  # dynamical 모드에만 필요
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

# 시각화
scv.pl.velocity_embedding(adata, basis='umap', color='cell_type')
scv.pl.velocity_embedding_stream(adata, basis='umap', color='cell_type')
```

## 문제해결

### 1. 'spliced'/'unspliced' 카운트가 없는 경우

이 문제는 velocyto 전처리가 필요하다는 의미입니다. 위 가이드에 따라 velocyto를 실행하세요.

### 2. Velocyto 실행 오류

Velocyto 실행 중 메모리 부족 오류가 발생하면 `--samtools-memory` 옵션을 조정하세요.

### 3. BAM 파일 문제

BAM 파일이 CellRanger에서 생성한 `possorted_genome_bam.bam` 형식인지 확인하세요.

## 참고자료

- [Velocyto 공식 문서](http://velocyto.org/)
- [scVelo 공식 문서](https://scvelo.readthedocs.io/)
- [RNA Velocity 원논문](https://www.nature.com/articles/s41586-018-0414-6)
- [scVelo 논문](https://www.nature.com/articles/s41587-020-0591-3) 