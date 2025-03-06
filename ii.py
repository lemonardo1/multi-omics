import os
from pathlib import Path
import humanize

def find_large_files(directory, limit=20):
    files = []
    
    for path in Path(directory).rglob('*'):
        if path.is_file():
            size = os.path.getsize(path)
            files.append((size, path))
    
    # 크기순으로 정렬
    files.sort(reverse=True)
    
    print(f"\n가장 큰 {limit}개 파일:")
    print("-" * 80)
    print(f"{'크기':>10} {'경로':<70}")
    print("-" * 80)
    
    for size, path in files[:limit]:
        # 사람이 읽기 쉬운 형태로 파일 크기 변환
        readable_size = humanize.naturalsize(size)
        print(f"{readable_size:>10} {str(path):<70}")

# 실행
find_large_files('daeseong')
