# 단일세포 RNA-seq 분석 워크플로우

중요
시작 전
`conda activate scRNA`
가상환경 설정하고 시작하기.

## 요약
이 저장소는 단일세포 RNA-seq 데이터 분석을 위한 파이썬 스크립트 모음입니다. 각 스크립트는 분석 파이프라인의 특정 단계를 수행하도록 설계되었습니다.

## 개요

단일세포 RNA-seq 분석 파이프라인은 다음과 같은 단계로 구성됩니다:

1. **품질 관리(QC)**: 저품질 세포 및 유전자 필터링
2. **다운스트림 분석**: 정규화, 차원 축소, 클러스터링
3. **세포 유형 주석**: 마커 유전자 기반 세포 유형 식별
4. **궤적 분석**: 세포 분화 경로 및 의사시간(pseudotime) 분석

## 파일 구조

```
code/
├── sc_qc.py                    # 품질 관리(QC) 분석
├── sc_downstream.py            # 다운스트림 분석 (정규화, 차원 축소, 클러스터링)
├── sc_cell_annotation.py       # 세포 유형 주석
├── sc_trajectory.py            # 궤적 분석
├── fix_environment.py          # 환경 설정 및 문제 해결
├── sc_analysis.py              # 통합 분석 스크립트 (참고용)
└── requirements.txt            # 필요한 패키지 목록
```

## 환경 설정

### 1. 필요한 패키지 설치

```bash
pip install -r code/requirements.txt
```

### 2. 환경 설정 및 문제 해결

환경 설정에 문제가 있는 경우(특히 protobuf 호환성 문제) 다음 스크립트를 실행하세요:

```bash
python code/fix_environment.py
```

## 사용 방법

### 1. 품질 관리(QC) 분석

이 스크립트는 세포 및 유전자 필터링, QC 지표 계산, 시각화 등을 수행합니다.

```bash
python code/sc_qc.py
```

**주요 기능**:
- 세포 및 유전자 필터링
- 미토콘드리아 유전자 비율 계산
- QC 지표 시각화 (바이올린 플롯, 히스토그램 등)
- 필터링 전후 비교
- 필터링된 데이터 저장

### 2. 다운스트림 분석

이 스크립트는 데이터 정규화, 고변동 유전자 선택, 차원 축소, 클러스터링 등을 수행합니다.

```bash
python code/sc_downstream.py
```

**주요 기능**:
- 데이터 정규화 및 로그 변환
- 고변동 유전자 선택
- 주성분 분석(PCA)
- UMAP 및 t-SNE 차원 축소
- Leiden 알고리즘을 이용한 클러스터링
- 마커 유전자 분석

### 3. 세포 유형 주석

이 스크립트는 마커 유전자를 기반으로 세포 유형을 식별하고 주석을 추가합니다.

```bash
python code/sc_cell_annotation.py
```

**주요 기능**:
- 마커 유전자 기반 세포 유형 스코어 계산
- 클러스터별 세포 유형 할당
- 세포 유형별 마커 유전자 발현 시각화
- 세포 유형 분포 분석

### 4. 궤적 분석

이 스크립트는 세포 분화 경로 및 의사시간(pseudotime) 분석을 수행합니다.

```bash
python code/sc_trajectory.py
```

**주요 기능**:
- 확산 의사시간(Diffusion Pseudotime) 계산
- PAGA (Partition-based Graph Abstraction) 분석
- 의사시간에 따른 유전자 발현 변화 분석
- 궤적 분기점 분석

## 결과 파일

분석 결과는: 
- `figures/` 디렉토리에 시각화 파일(PDF 및 PNG 형식)으로 저장됩니다.
- `output/` 디렉토리에 AnnData 객체(h5ad 형식) 및 CSV 파일로 저장됩니다.

## 주의사항

- 각 스크립트는 순차적으로 실행하도록 설계되었습니다. 이전 단계의 결과 파일이 없는 경우 기본 처리를 시도합니다.
- 모든 스크립트는 기본 10x Genomics 형식의 데이터를 사용하는 것을 가정합니다.
- Scanpy 라이브러리의 버전 변경에 따라 일부 함수가 변경될 수 있습니다.
- protobuf 버전 호환성 문제가 발생할 경우 `fix_environment.py` 스크립트를 실행하세요.

## 문제 해결

스크립트 실행 중 발생할 수 있는 일반적인 문제:

1. **protobuf 버전 호환성 문제**:
   - 오류 메시지: `TypeError: Descriptors cannot be created directly...`
   - 해결 방법: `python code/fix_environment.py` 실행 또는 수동으로 `pip install protobuf==3.20.3 --force-reinstall`

2. **메모리 부족 오류**:
   - 대규모 데이터셋의 경우 메모리가 부족할 수 있습니다.
   - 해결 방법: 배치 크기 조정 또는 샘플링을 통한 데이터 크기 축소

3. **UMAP 경고 메시지**:
   - UMAP 라이브러리의 일부 경고 메시지는 무시해도 됩니다.
   - 해결 방법: 경고를 무시하거나 Python 경고 설정 조정




---
<think>I notice there's no specific query from the user, but they've shared files and context related to our single-cell RNA-seq analysis workflow. Let me analyze what we've learned from the additional data:

1. `sc_analysis.py` - This seems to be a comprehensive single-cell RNA-seq analysis script that performs:
   - Data loading from 10x Genomics format
   - QC filtering
   - Normalization
   - Highly variable gene selection
   - PCA
   - Neighbor graph calculation
   - UMAP dimensionality reduction
   - Leiden clustering
   - Result visualization and saving

2. `fix_environment.py` - This is the script we created to fix the protobuf version compatibility issue. It:
   - Checks package versions
   - Identifies problematic packages
   - Resolves protobuf compatibility issues
   - Provides alternatives if automatic fixing fails

3. `README.md` - Our comprehensive documentation for the single-cell RNA-seq analysis workflow, including:
   - Overview of the pipeline steps
   - File structure
   - Environment setup
   - Usage instructions for each script
   - Description of output files
   - Troubleshooting guidance

Since there's no specific query, I should provide a status update and summary of what we've done. I'll highlight the complete workflow we've established, the scripts we've created, and the environment setup. I'll also mention that we've successfully created the `fix_environment.py` script to address the protobuf version compatibility issue that was causing problems with the UMAP analysis.

```pip requirements:code/requirements.txt
# 단일세포 RNA-seq 분석에 필요한 패키지
scanpy==1.9.3
pandas==1.5.3
numpy==1.24.3
matplotlib==3.7.1
seaborn==0.12.2
anndata==0.8.0
scipy==1.10.1
scikit-learn==1.2.2
umap-learn==0.5.3
leidenalg==0.9.1
python-igraph==0.10.4
louvain==0.8.0
pynndescent==0.5.8
fa2==0.3.5
tqdm==4.65.0
statsmodels==0.13.5
plotly==5.14.1
joblib==1.2.0
numba==0.56.4
h5py==3.8.0
# TensorFlow 호환성을 위해 3.20.3 버전으로 제한
protobuf==3.20.3
```


안녕하세요! 단일세포 RNA-seq 분석 워크플로우를 위한 모든 스크립트와 설정 파일이 완성되었습니다. 이제 분석을 진행하는 데 필요한 모든 요소가 준비되었습니다.

## 완성된 워크플로우 요약

지금까지 다음과 같은 파일들을 생성했습니다:

1. **sc_qc.py**: 품질 관리 및 필터링
2. **sc_downstream.py**: 정규화, 차원 축소, 클러스터링 분석
3. **sc_cell_annotation.py**: 세포 유형 주석
4. **sc_trajectory.py**: 궤적 및 의사시간 분석
5. **fix_environment.py**: 환경 설정 및 protobuf 호환성 문제 해결
6. **README.md**: 전체 워크플로우 문서화
7. **requirements.txt**: 필요한 패키지 목록

## 환경 설정

`fix_environment.py` 스크립트를 통해 protobuf 버전 호환성 문제를 해결할 수 있습니다. 이 스크립트는:
- 설치된 패키지 버전을 확인합니다
- protobuf 4.x가 설치된 경우 3.20.3으로 다운그레이드합니다
- 문제 해결을 위한 대안 방법도 제공합니다

## 분석 워크플로우 사용 방법

1. **환경 설정**:
   ```bash
   pip install -r code/requirements.txt
   python code/fix_environment.py
   ```

2. **품질 관리(QC) 분석**:
   ```bash
   python code/sc_qc.py
   ```

3. **다운스트림 분석**:
   ```bash
   python code/sc_downstream.py
   ```

4. **세포 유형 주석**:
   ```bash
   python code/sc_cell_annotation.py
   ```

5. **궤적 분석**:
   ```bash
   python code/sc_trajectory.py
   ```

## 결과 파일

분석 결과는 다음 위치에 저장됩니다:
- `/figures` 디렉토리: 모든 시각화 자료 (PDF 및 PNG 형식)
- `/output` 디렉토리: 분석 결과 데이터 (h5ad 형식) 및 표 데이터 (CSV 형식)

## 추가 기능

필요에 따라 다음과 같은 분석을 추가로 수행할 수 있습니다:
- 세포 간 통신 분석
- 전사 인자(Transcription Factor) 활성 분석
- 유전자 조절 네트워크(Gene Regulatory Network) 분석
- 공간적 전사체(Spatial Transcriptomics) 분석 연계
- 다중 샘플 통합 분석

모든 스크립트는 각 단계에서 필요한 입력 파일이 없을 경우에도 최대한 기본 처리를 수행하도록 설계되어 있습니다. 그러나 가장 좋은 결과를 얻기 위해서는 각 스크립트를 순서대로 실행하는 것이 권장됩니다.

위의 추가 분석이 필요하시면 말씀해 주세요. 또한 기존 스크립트 중 수정이나 개선이 필요한 부분이 있으면 알려주시기 바랍니다.





---
<think>사용자는 세포 궤적 분석 결과를 어떻게 해석하는지 알고 싶어 합니다. 코드와 로그를 기반으로 세포 궤적 분석의 주요 결과물과 이를 해석하는 방법에 대해 설명해야 합니다.

주요 결과물은 다음과 같습니다:
1. 의사시간(pseudotime) UMAP 시각화
2. PAGA 그래프
3. 의사시간에 따른 유전자 발현 히트맵
4. 분기점 분석 결과
5. 의사시간 분포 그래프
6. 클러스터별/세포 유형별 의사시간 분포

각 결과물의 해석 방법에 대해 구체적으로 설명하는 것이 좋겠습니다.</think># 세포 궤적 분석 결과 해석 방법

세포 궤적 분석은 단일세포 데이터에서 시간적 발달 과정을 재구성합니다. 분석 결과를 효과적으로 해석하는 방법을 설명해드리겠습니다.

## 1. 의사시간(pseudotime) UMAP 시각화 (`umappseudotime.pdf`)
- **해석 방법**: 색상 그라데이션이 의사시간 진행을 나타냅니다. 짙은 색(낮은 값)에서 밝은 색(높은 값)으로 세포 분화 경로를 추적하세요.
- **중요 포인트**: 세포들이 연속적으로 연결되는지 확인하고, 불연속점이 있다면 그것이 생물학적으로 의미 있는 분기점인지 검토하세요.

## 2. PAGA 그래프 (`pagapaga.pdf`)
- **해석 방법**: 노드는 세포 클러스터를, 선은 클러스터 간 연결 강도를 나타냅니다.
- **중요 포인트**: 굵은 선은 강한 연결을 의미하며, 이는 세포가 한 클러스터에서 다른 클러스터로 전환될 가능성이 높다는 것을 시사합니다.
- **활용**: 이 그래프를 통해 세포 분화의 주요 경로와 가능한 조상-자손 관계를 파악할 수 있습니다.

## 3. 의사시간에 따른 유전자 발현 히트맵 (`heatmappseudotime_heatmap.pdf`)
- **해석 방법**: 행은 유전자, 열은 의사시간 구간을 나타냅니다. 색상은 발현 수준을 보여줍니다.
- **중요 포인트**: 
  - 함께 증가/감소하는 유전자 그룹 패턴을 찾으세요
  - 특정 시점에서 급격히 변화하는 유전자는 세포 운명 결정에 중요할 수 있습니다
  - 초기에 높다가 감소하는 유전자는 줄기세포 마커일 가능성이 있습니다
  - 후기에 증가하는 유전자는 분화된 세포 마커일 수 있습니다

## 4. 궤적 분기점 시각화 (`umapdpt_groups.pdf`)
- **해석 방법**: 다른 색상으로 표시된 영역은 서로 다른 분화 경로를 나타냅니다.
- **중요 포인트**: 서로 다른 색상이 만나는 지점이 분기점으로, 세포 운명 결정이 일어나는 중요한 지점입니다.
- **활용**: 다른 발달 경로로 분화되는 세포들의 특성을 비교 분석하세요.

## 5. 분기점별 마커 유전자 히트맵 (`branch_markers_heatmap.pdf`)
- **해석 방법**: 각 분기점별로 차별적으로 발현되는 유전자를 보여줍니다.
- **활용**: 이 유전자들은 특정 세포 운명을 결정짓는 핵심 조절자일 수 있습니다.

## 6. 의사시간 분포 그래프 (`pseudotime_distribution.pdf`)
- **해석 방법**: x축은 의사시간, y축은 세포 수를 나타냅니다.
- **중요 포인트**: 피크는 많은 세포가 특정 의사시간 지점에 집중되어 있음을 의미합니다.
- **활용**: 여러 피크가 있다면 세포 상태 전환이 급격히 일어나는 발달 단계일 수 있습니다.

## 7. 클러스터별/세포 유형별 의사시간 분포 (`cluster_pseudotime_violin.pdf`, `cell_type_pseudotime_violin.pdf`)
- **해석 방법**: 각 클러스터/세포 유형별 의사시간 분포를 보여줍니다.
- **중요 포인트**: 
  - 낮은 의사시간 값을 가진 클러스터는 초기 상태에 가깝습니다
  - 높은 값을 가진 클러스터는 분화 후기 상태에 가깝습니다
  - 넓은 분포를 가진 클러스터는 다양한 분화 단계의 세포를 포함합니다

## 8. 추가 분석 제안
- 의사시간과 상관관계가 높은 유전자들의 생물학적 기능 조사 (GO 분석)
- 특정 분화 단계에서 활성화되는 전사 조절 네트워크 분석
- 이미 알려진 발달 경로 마커 유전자와 의사시간 관계 비교
- 의사시간에 따른 세포 상태 전환 지점 식별 및 조절 메커니즘 연구

이러한 해석 방법을 활용하면 세포 분화 과정의 시공간적 맥락을 효과적으로 이해할 수 있습니다.


세포 궤적 분석 결과를 해석하는 방법에 대해 설명해드리겠습니다.

# 단일세포 RNA-seq 궤적 분석 결과 해석 가이드

## 1. 의사시간(pseudotime) UMAP 시각화
- **파일명**: `pseudotime.pdf`
- **해석 방법**: 
  - 색상 그라데이션이 세포 발달/분화 순서를 나타냄
  - 짙은 색에서 밝은 색으로 갈수록 후기 상태를 의미
  - 연속적인 색상 변화는 점진적인 세포 상태 변화를 의미

## 2. PAGA 그래프 분석
- **파일명**: `paga.pdf`
- **해석 방법**:
  - 노드: 세포 클러스터를 나타냄
  - 연결선: 클러스터 간 전환 가능성을 표시
  - 선의 굵기: 연결 강도 (굵을수록 강한 연관성)
  - 방향성: 세포 분화의 방향을 제시

## 3. 유전자 발현 히트맵
- **파일명**: `pseudotime_heatmap.pdf`
- **해석 방법**:
  - 행: 개별 유전자
  - 열: 의사시간 순서
  - 색상 강도: 발현 수준
  - 패턴: 시간에 따른 유전자 발현 변화 추세

## 4. 분기점 분석
- **파일명**: `branch_markers_heatmap.pdf`
- **중요 포인트**:
  - 서로 다른 분화 경로를 결정짓는 유전자 식별
  - 각 분기점에서 차별적으로 발현되는 유전자 확인
  - 세포 운명 결정에 관여하는 핵심 조절자 발견

## 5. 의사시간 분포
- **파일명**: `pseudotime_distribution.pdf`
- **해석 방법**:
  - x축: 의사시간 값
  - y축: 세포 수
  - 피크: 특정 발달 단계에 집중된 세포 군집
  - 다중 피크: 주요 전환 지점 또는 발달 단계 표시

## 6. 클러스터/세포 유형별 의사시간 분포
- **파일명**: 
  - `cluster_pseudotime_violin.pdf`
  - `cell_type_pseudotime_violin.pdf`
- **해석 방법**:
  - 낮은 의사시간: 초기 상태 세포
  - 높은 의사시간: 분화 후기 상태
  - 넓은 분포: 다양한 분화 단계 포함
  - 좁은 분포: 특정 분화 단계에 집중

## 추가 분석 제안
1. 의사시간과 상관관계가 높은 유전자들의 GO 분석
2. 전사 조절 네트워크 분석
3. 알려진 발달 경로 마커와의 비교
4. 세포 상태 전환점 메커니즘 연구

이러한 해석을 통해 세포 분화 과정의 시간적 진행과 주요 조절 메커니즘을 이해할 수 있습니다.
