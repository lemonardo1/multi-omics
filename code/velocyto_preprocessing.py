#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
velocyto_preprocessing.py: 10x Genomics 데이터에서 velocyto를 사용한 RNA velocity 전처리 

사용법:
python /home/swr0460/daeseong/code/velocyto_preprocessing.py --bam /path/to/possorted_genome_bam.bam --gtf /path/to/genes.gtf --output_dir /path/to/output --mask_file /path/to/mask.gtf (선택사항)

이 스크립트는 velocyto를 사용하여 10x Genomics 데이터에서 spliced/unspliced 카운트를 생성합니다.
결과로 생성된 loom 파일은 scvelo를 사용한 RNA velocity 분석에 사용할 수 있습니다.
"""

import os
import sys
import argparse
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvelo as scv
from datetime import datetime

# 경고 메시지 무시
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    """명령줄 인수 파싱"""
    parser = argparse.ArgumentParser(description='velocyto 전처리 스크립트')
    
    parser.add_argument('--bam', type=str, required=True, 
                        help='10x Genomics possorted_genome_bam.bam 파일 경로')
    parser.add_argument('--gtf', type=str, required=True,
                        help='참조 유전체 GTF 파일 경로')
    parser.add_argument('--mask_file', type=str, default=None,
                        help='마스킹할 영역의 GTF 파일 경로 (선택사항)')
    parser.add_argument('--output_dir', type=str, default=None,
                        help='결과 파일을 저장할 디렉토리 (기본값: BAM 파일과 같은 디렉토리)')
    parser.add_argument('--samtools_memory', type=str, default='32G',
                        help='samtools에 할당할 메모리 (기본값: 32G)')
    parser.add_argument('--samtools_threads', type=int, default=8,
                        help='samtools에 할당할 스레드 수 (기본값: 8)')
    parser.add_argument('--convert_to_h5ad', action='store_true',
                        help='loom 파일을 h5ad로 변환할지 여부')
    
    return parser.parse_args()

def setup_directories(args):
    """결과 저장을 위한 디렉토리 설정"""
    # 출력 디렉토리 설정
    if args.output_dir is None:
        args.output_dir = os.path.dirname(args.bam)
    
    # 출력 디렉토리 생성
    os.makedirs(args.output_dir, exist_ok=True)
    
    return args.output_dir

def run_velocyto(args, output_dir):
    """velocyto CLI 실행"""
    print(f"[INFO] velocyto 전처리 시작...")
    
    # BAM 파일 기반으로 샘플 이름 추출
    sample_name = os.path.basename(args.bam).replace(".bam", "")
    
    # velocyto 명령어 구성
    cmd = ['velocyto', 'run', 
           f'--samtools-memory={args.samtools_memory}',
           f'--samtools-threads={args.samtools_threads}',
           '-v']  # 자세한 출력
    
    # 마스크 파일이 지정된 경우 추가
    if args.mask_file:
        cmd.extend(['-m', args.mask_file])
    
    # BAM 파일, GTF 파일 추가
    cmd.extend([args.bam, args.gtf])
    
    print(f"[INFO] 실행 명령어: {' '.join(cmd)}")
    
    try:
        # velocyto 실행
        subprocess.run(cmd, check=True)
        
        # 기본 출력 위치: BAM 파일과 같은 디렉토리에 velocyto 하위 디렉토리
        velocyto_output_dir = os.path.join(os.path.dirname(args.bam), 'velocyto')
        if os.path.exists(velocyto_output_dir):
            loom_files = [f for f in os.listdir(velocyto_output_dir) if f.endswith('.loom')]
            if loom_files:
                loom_file = os.path.join(velocyto_output_dir, loom_files[0])
                print(f"[INFO] velocyto 처리 완료. 결과 파일: {loom_file}")
                return loom_file
        
        # 기본 위치에 없는 경우 다른 위치 확인
        alt_velocyto_dir = os.path.join(args.output_dir, 'velocyto')
        if os.path.exists(alt_velocyto_dir):
            loom_files = [f for f in os.listdir(alt_velocyto_dir) if f.endswith('.loom')]
            if loom_files:
                loom_file = os.path.join(alt_velocyto_dir, loom_files[0])
                print(f"[INFO] velocyto 처리 완료. 결과 파일: {loom_file}")
                return loom_file
        
        # 현재 디렉토리에서 loom 파일 찾기
        loom_files = [f for f in os.listdir('.') if f.endswith('.loom')]
        if loom_files:
            loom_file = os.path.abspath(loom_files[0])
            print(f"[INFO] velocyto 처리 완료. 결과 파일: {loom_file}")
            return loom_file
            
        print(f"[WARNING] loom 파일을 찾을 수 없습니다.")
        return None
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] velocyto 실행 실패: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] 예상치 못한 오류: {e}")
        sys.exit(1)

def convert_loom_to_h5ad(loom_file, output_dir):
    """loom 파일을 h5ad 형식으로 변환"""
    if loom_file is None:
        print(f"[WARNING] 변환할 loom 파일이 없습니다.")
        return None
        
    print(f"[INFO] loom 파일을 h5ad로 변환 중...")
    
    try:
        # loom 파일 로드
        adata = scv.read_loom(loom_file)
        
        # 파일 이름 설정
        h5ad_file = os.path.join(output_dir, 
                               os.path.basename(loom_file).replace('.loom', '_velocity.h5ad'))
        
        # 파일 저장
        adata.write(h5ad_file)
        
        print(f"[INFO] 변환 완료. h5ad 파일: {h5ad_file}")
        return h5ad_file
        
    except Exception as e:
        print(f"[ERROR] loom에서 h5ad로 변환 실패: {e}")
        return None

def check_velocyto_installation():
    """velocyto 설치 여부 확인"""
    try:
        subprocess.run(['velocyto', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except FileNotFoundError:
        return False

def main():
    """메인 실행 함수"""
    print("=" * 80)
    print("velocyto 전처리 스크립트")
    print("=" * 80)
    
    # velocyto 설치 확인
    if not check_velocyto_installation():
        print("[ERROR] velocyto가 설치되어 있지 않습니다.")
        print("다음 명령어로 설치해 주세요:")
        print("pip install velocyto")
        sys.exit(1)
    
    # 인수 파싱
    args = parse_arguments()
    
    # 파일 경로 확인
    if not os.path.exists(args.bam):
        print(f"[ERROR] BAM 파일을 찾을 수 없습니다: {args.bam}")
        sys.exit(1)
    
    if not os.path.exists(args.gtf):
        print(f"[ERROR] GTF 파일을 찾을 수 없습니다: {args.gtf}")
        sys.exit(1)
    
    if args.mask_file and not os.path.exists(args.mask_file):
        print(f"[ERROR] 마스크 파일을 찾을 수 없습니다: {args.mask_file}")
        sys.exit(1)
    
    # 출력 디렉토리 설정
    output_dir = setup_directories(args)
    print(f"[INFO] 출력 디렉토리: {output_dir}")
    
    # velocyto 실행
    loom_file = run_velocyto(args, output_dir)
    
    # h5ad로 변환 (옵션)
    if args.convert_to_h5ad and loom_file:
        h5ad_file = convert_loom_to_h5ad(loom_file, output_dir)
        
        if h5ad_file:
            print(f"\n[INFO] 변환된 h5ad 파일을 사용하여 RNA velocity 분석을 수행할 수 있습니다.")
            print(f"[INFO] 예시 명령어:")
            print(f"python /home/swr0460/daeseong/code/direct_velocity.py --input {h5ad_file}")
    
    print("\n[INFO] velocyto 전처리 완료.")
    
    # 다음 단계 안내
    if loom_file:
        print("\n== 다음 단계 ==")
        print("1. 생성된 loom 파일을 scvelo에서 사용:")
        print("   ```python")
        print("   import scvelo as scv")
        print(f"   adata = scv.read_loom('{loom_file}')")
        print("   scv.pp.filter_and_normalize(adata)")
        print("   scv.pp.moments(adata)")
        print("   scv.tl.velocity(adata)")
        print("   scv.tl.velocity_graph(adata)")
        print("   scv.pl.velocity_embedding(adata, basis='umap')")
        print("   ```")
        print("2. 또는 direct_velocity.py 스크립트 사용:")
        if args.convert_to_h5ad and h5ad_file:
            print(f"   python /home/swr0460/daeseong/code/direct_velocity.py --input {h5ad_file}")
        else:
            print(f"   python /home/swr0460/daeseong/code/direct_velocity.py --input {loom_file}")

if __name__ == "__main__":
    main() 