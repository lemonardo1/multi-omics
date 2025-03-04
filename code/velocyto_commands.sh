#!/bin/bash
# velocyto 명령어 예제 모음

# 1. velocyto 설치
# conda 환경을 사용할 경우
conda_install_velocyto="
# scRNA 환경에 velocyto 설치
conda activate scRNA
pip install velocyto

# 필요한 추가 의존성 설치
conda install -c conda-forge -y numpy scipy cython numba matplotlib scikit-learn h5py click
"

# pip로 직접 설치할 경우 
pip_install_velocyto="
pip install velocyto
"

# 2. 10x genomics 데이터에 대한 velocyto 실행 예시
run_velocyto_10x="
# 기본 실행
velocyto run10x -m repeat_msk.gtf sample_directory reference.gtf

# 또는 일반적인 run 명령어 사용
velocyto run -b filtered_barcodes.tsv -o output_dir -m repeat_msk.gtf possorted_genome_bam.bam reference.gtf

# 샘플 이름 지정 예시
velocyto run -b filtered_barcodes.tsv -o output_dir --samtools-memory=32G --samtools-threads=8 -m repeat_msk.gtf possorted_genome_bam.bam reference.gtf
"

# 3. 예시: 표준 10x Genomics 데이터에 대한 velocyto 실행
example_command="
# 샘플 디렉토리 지정
SAMPLE_DIR=/path/to/10x_genomics_output

# GTF 파일 경로
REFERENCE_GTF=/path/to/reference/genes.gtf

# 반복 영역 마스크 파일 (선택사항)
MASK_GTF=/path/to/mask/repeat_msk.gtf

# velocyto 실행
velocyto run10x -m $MASK_GTF $SAMPLE_DIR $REFERENCE_GTF

# 또는, BAM 파일 직접 지정
BAM_FILE=/path/to/possorted_genome_bam.bam
velocyto run --samtools-memory=32G --samtools-threads=8 -m $MASK_GTF $BAM_FILE $REFERENCE_GTF
"

# 4. 작성한 Python 스크립트 실행 예시
run_python_scripts="
# velocyto 전처리 실행
python /home/swr0460/daeseong/code/velocyto_preprocessing.py --bam /path/to/possorted_genome_bam.bam --gtf /path/to/reference.gtf --output_dir /path/to/output --mask_file /path/to/repeat_msk.gtf --convert_to_h5ad

# loom 파일을 h5ad로 변환
python /home/swr0460/daeseong/code/velocyto_loom_to_h5ad.py --loom /path/to/velocyto_output.loom --output /path/to/output.h5ad

# 변환된 h5ad 파일로 RNA velocity 분석
python /home/swr0460/daeseong/code/direct_velocity.py --input /path/to/velocyto_output_velocity.h5ad --mode dynamical
"

# 출력
echo "velocyto 설치 명령어 (conda):"
echo "$conda_install_velocyto"
echo 
echo "velocyto 설치 명령어 (pip):"
echo "$pip_install_velocyto"
echo
echo "velocyto 실행 예시:"
echo "$run_velocyto_10x"
echo
echo "실제 샘플에 대한 velocyto 실행 예시:"
echo "$example_command"
echo
echo "작성한 Python 스크립트 실행 예시:"
echo "$run_python_scripts"
