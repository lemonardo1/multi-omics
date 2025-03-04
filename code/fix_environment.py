#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
단일세포 RNA-seq 분석 환경 설정 및 문제 해결 스크립트
- protobuf 버전 호환성 문제 해결
- 필요한 패키지 설치
- 환경 설정 체크
"""

import subprocess
import sys
import os
import pkg_resources

def check_packages():
    """
    필요한 패키지 버전을 확인하고 문제가 있는 경우 해결 방법을 제안합니다.
    """
    print("패키지 버전 확인 중...")
    
    # 설치된 패키지 목록 가져오기
    installed_packages = {pkg.key: pkg.version for pkg in pkg_resources.working_set}
    
    # 필요한 패키지 및 권장 버전
    required_packages = {
        'scanpy': '>=1.9.0',
        'pandas': '>=1.3.0',
        'numpy': '>=1.20.0',
        'matplotlib': '>=3.4.0',
        'seaborn': '>=0.11.0',
        'protobuf': '<=3.20.3',  # protobuf 버전 제한
        'umap-learn': '>=0.5.0',
    }
    
    # 버전 문제가 있는 패키지 목록
    problematic_packages = []
    
    for package, version_req in required_packages.items():
        if package in installed_packages:
            print(f"설치됨: {package} {installed_packages[package]}")
            
            # protobuf 버전 체크 (3.20.x 이하 권장)
            if package == 'protobuf' and installed_packages[package].startswith('4.'):
                problematic_packages.append(package)
                print(f"  ⚠️ 경고: protobuf 버전 {installed_packages[package]}는 TensorFlow와 호환성 문제가 있을 수 있습니다.")
                print("     권장 버전: 3.20.3 이하")
        else:
            print(f"설치되지 않음: {package}")
            problematic_packages.append(package)
    
    return problematic_packages

def fix_protobuf_issue():
    """
    protobuf 호환성 문제 해결
    """
    print("\nprotobuf 호환성 문제 해결 중...")
    
    try:
        # 현재 버전 확인
        protobuf_version = pkg_resources.get_distribution("protobuf").version
        print(f"현재 protobuf 버전: {protobuf_version}")
        
        if protobuf_version.startswith('4.'):
            print("호환성을 위해 protobuf 버전을 다운그레이드합니다 (3.20.3)...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", "protobuf==3.20.3", "--force-reinstall"])
            print("✅ protobuf 다운그레이드 완료")
        else:
            print("현재 protobuf 버전은 호환됩니다. 작업이 필요하지 않습니다.")
    
    except pkg_resources.DistributionNotFound:
        print("protobuf가 설치되어 있지 않습니다. 설치 중...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "protobuf==3.20.3"])
        print("✅ protobuf 설치 완료")
    
    except Exception as e:
        print(f"⚠️ 오류 발생: {e}")
        print("\n대체 해결 방법:")
        print("1. 환경 변수 설정: PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python")
        print("   (이 방법은 순수 Python 구현을 사용하므로 성능이 느립니다)")
        print("2. 수동으로 다음 명령어 실행:")
        print("   pip install protobuf==3.20.3 --force-reinstall")
        return False
    
    return True

def setup_environment():
    """
    전체 환경 설정
    """
    print("==== 단일세포 RNA-seq 분석 환경 설정 ====")
    
    # 패키지 확인
    problematic_packages = check_packages()
    
    # protobuf 문제 해결
    if 'protobuf' in problematic_packages or len(problematic_packages) > 0:
        success = fix_protobuf_issue()
        if not success:
            print("\n⚠️ 환경 설정 중 문제가 발생했습니다. 수동으로 해결해주세요.")
            return
    
    print("\n✅ 환경 설정이 완료되었습니다!")
    print("이제 분석 스크립트를 실행할 수 있습니다:")
    print("1. QC 분석: python code/sc_qc.py")
    print("2. 다운스트림 분석: python code/sc_downstream.py")
    print("3. 세포 유형 주석: python code/sc_cell_annotation.py")
    print("4. 궤적 분석: python code/sc_trajectory.py")

if __name__ == "__main__":
    setup_environment() 