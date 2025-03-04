#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
velocyto_loom_to_h5ad.py: velocyto loom 파일을 h5ad 형식으로 변환하는 유틸리티

사용법:
python /home/swr0460/daeseong/code/velocyto_loom_to_h5ad.py --loom /path/to/velocyto.loom --output /path/to/output.h5ad

이 스크립트는 velocyto로 생성된 loom 파일을 scvelo 분석에 적합한 h5ad 형식으로 변환합니다.
"""

import os
import sys
import argparse
import scvelo as scv
from datetime import datetime

# 경고 메시지 무시
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    """명령줄 인수 파싱"""
    parser = argparse.ArgumentParser(description='velocyto loom 파일을 h5ad로 변환')
    
    parser.add_argument('--loom', type=str, required=True, 
                        help='velocyto로 생성된 loom 파일 경로')
    parser.add_argument('--output', type=str, default=None,
                        help='출력 h5ad 파일 경로 (기본값: loom 파일과 같은 디렉토리에 _velocity.h5ad)')
    
    return parser.parse_args()

def convert_loom_to_h5ad(loom_file, output_file=None):
    """loom 파일을 h5ad 형식으로 변환"""
    print(f"[INFO] '{loom_file}' 파일을 로드 중...")
    
    try:
        # loom 파일 로드
        adata = scv.read_loom(loom_file)
        print(f"[INFO] 데이터 로드 완료: {adata.shape[0]} 세포, {adata.shape[1]} 유전자")
        
        # 출력 파일 경로 설정
        if output_file is None:
            output_file = os.path.splitext(loom_file)[0] + "_velocity.h5ad"
        
        # 출력 디렉토리 확인 및 생성
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # 파일 저장
        print(f"[INFO] h5ad 파일로 저장 중: {output_file}")
        adata.write(output_file)
        
        print(f"[INFO] 변환 완료. h5ad 파일: {output_file}")
        return output_file
        
    except Exception as e:
        print(f"[ERROR] loom에서 h5ad로 변환 실패: {e}")
        return None

def main():
    """메인 실행 함수"""
    print("=" * 80)
    print("velocyto loom 파일을 h5ad로 변환")
    print("=" * 80)
    
    # 인수 파싱
    args = parse_arguments()
    
    # 파일 경로 확인
    if not os.path.exists(args.loom):
        print(f"[ERROR] loom 파일을 찾을 수 없습니다: {args.loom}")
        sys.exit(1)
    
    # loom 파일 변환
    h5ad_file = convert_loom_to_h5ad(args.loom, args.output)
    
    if h5ad_file:
        print("\n== 다음 단계 ==")
        print(f"1. 변환된 h5ad 파일을 direct_velocity.py 스크립트로 분석할 수 있습니다:")
        print(f"   python /home/swr0460/daeseong/code/direct_velocity.py --input {h5ad_file}")
        print("\n2. 또는 Python에서 직접 사용:")
        print("   ```python")
        print("   import scvelo as scv")
        print(f"   adata = scv.read_h5ad('{h5ad_file}')")
        print("   scv.pp.filter_and_normalize(adata)")
        print("   scv.pp.moments(adata)")
        print("   scv.tl.velocity(adata)")
        print("   scv.tl.velocity_graph(adata)")
        print("   scv.pl.velocity_embedding(adata, basis='umap')")
        print("   ```")

if __name__ == "__main__":
    main() 