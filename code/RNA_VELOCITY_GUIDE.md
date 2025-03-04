# RNA velocity 분석 가이드

이 가이드는 `sc_rna_velocity.py` 스크립트를 사용하여 단일세포 RNA-seq 데이터로부터 RNA velocity 분석을 수행하는 방법을 설명합니다.

## RNA velocity란?

RNA velocity는 단일세포 수준에서 세포의 미래 상태를 예측할 수 있는 강력한 도구입니다. 이 분석은 spliced(성숙) mRNA와 unspliced(전구체) mRNA의 비율을 측정하여 유전자 발현의 변화율을 추정합니다. 이를 통해 세포 발달 방향과 궤적을 예측할 수 있습니다.

## 필수 요구사항

이 스크립트를 실행하기 위해서는 다음 패키지가 필요합니다:
- Python 3.8 이상
- scanpy
- anndata
- scvelo
- cellrank
- loompy (선택 사항, .loom 파일 분석 시 필요)

모든 의존성 패키지는 다음 명령으로 설치할 수 있습니다:
```bash
pip install -r requirements.txt
```

## 데이터 준비

RNA velocity 분석을 위해서는 다음 두 가지 방법 중 하나로 데이터를 준비할 수 있습니다:

1. **Velocyto로 생성된 .loom 파일**: 이는 spliced/unspliced 정보가 이미 포함된 파일입니다.
   ```bash
   velocyto run10x -m repeat_msk.gtf my_10x_directory/
   ```

2. **AnnData 객체**: 이미 존재하는 AnnData 객체에 spliced/unspliced 계수가 별도의 레이어로 저장된 경우

## 스크립트 사용법

기본적인 실행 방법은 다음과 같습니다:

```bash
python sc_rna_velocity.py --input results.h5ad
```

### 주요 옵션

- `--input`: 입력 AnnData 객체 경로 (필수)
- `--loom`: Velocyto로 생성된 .loom 파일 경로 (선택 사항)
- `--spliced_layer`: Spliced counts가 저장된 레이어명
- `--unspliced_layer`: Unspliced counts가 저장된 레이어명
- `--mode`: RNA velocity 계산 모드 (deterministic, stochastic, dynamical)
- `--embed`: 시각화에 사용할 임베딩 (umap, tsne, pca)
- `--cluster_key`: 클러스터 정보가 있는 obs 컬럼명
- `--n_jobs`: 사용할 CPU 코어 수

### 예제 명령어

1. **기본 분석 (자동 spliced/unspliced 추론)**:
   ```bash
   python sc_rna_velocity.py --input results.h5ad --mode stochastic --embed umap --cluster_key leiden
   ```

2. **Loom 파일 사용**:
   ```bash
   python sc_rna_velocity.py --input results.h5ad --loom velocyto_results.loom --mode dynamical
   ```

3. **특정 레이어 지정**:
   ```bash
   python sc_rna_velocity.py --input results.h5ad --spliced_layer spliced --unspliced_layer unspliced
   ```

## 결과 해석

스크립트는 다음과 같은 주요 결과를 생성합니다:

1. **Velocity 임베딩 시각화**:
   - 그리드 플롯: 각 영역의 속도 방향을 보여줌
   - 스트림 플롯: 세포 궤적의 연속적인 흐름을 시각화
   - 화살표 플롯: 개별 세포의 속도 방향을 표시

2. **고급 분석 결과**:
   - PAGA 궤적 분석: 클러스터 간 연결성과 전환 확률
   - 상위 속도 유전자: 각 클러스터의 주요 속도 유전자
   - 속도 슈도시간: 세포 발달 순서
   - 잠재 시간 분석: 세포 발달 경로의 시작점과 종료점 예측

3. **CellRank 분석**:
   - 세포 운명 예측: 최종 세포 상태와 전환 확률
   - 상위 운명 유전자: 세포 운명 결정에 중요한 유전자

## 팁과 주의사항

- RNA velocity 분석은 계산 집약적일 수 있으므로, 대규모 데이터셋의 경우 서브샘플링을 고려하세요.
- 최상의 결과를 위해 고품질의 클러스터링과 UMAP/tSNE 임베딩이 필요합니다.
- Dynamical 모드는 가장 정확하지만 가장 많은 계산 자원을 필요로 합니다.
- 결과 해석 시 생물학적 맥락과 실험 디자인을 항상 고려하세요.

## 참고 문헌

1. La Manno, G. et al. (2018). RNA velocity of single cells. Nature, 560, 494-498.
2. Bergen, V. et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. Nature Biotechnology, 38, 1408-1414.
3. Lange, M. et al. (2022). CellRank for directed single-cell fate mapping. Nature Methods, 19, 159-170. 