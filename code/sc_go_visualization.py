#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 필요한 라이브러리 임포트
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import seaborn as sns
import sys
import matplotlib as mpl
import matplotlib.colors as mcolors
from adjustText import adjust_text

# GO 분석을 위한 goatools 임포트 (필요한 경우 설치)
try:
    from goatools import obo_parser
    from goatools.godag_plot import plot_gos, plot_results
    from goatools.gosubdag.gosubdag import GoSubDag
    from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot
except ImportError:
    print("goatools 라이브러리가 설치되지 않았습니다.")
    print("pip install goatools를 실행하여 설치할 수 있습니다.")

# 설정: 기본 디렉토리
base_dir = '/home/swr0460/daeseong'
data_dir = os.path.join(base_dir, 'data')
figures_dir = os.path.join(base_dir, 'figures')
output_dir = os.path.join(base_dir, 'output')
go_viz_dir = os.path.join(figures_dir, 'go_visualization')

# 디렉토리 생성
os.makedirs(go_viz_dir, exist_ok=True)

# 한글 폰트 설정 (필요한 경우)
plt.rcParams['font.family'] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False

def run_go_visualization():
    """
    GO 분석 결과에 대한 고급 시각화를 수행합니다.
    - GO 계층 구조 시각화
    - GO 네트워크 시각화
    - GO 용어 관계 그래프
    """
    
    print("\n===== GO 분석 결과 시각화 시작 =====")
    
    # GO 분석 결과 파일 확인
    pos_result_file = os.path.join(output_dir, 'positive_corr_GO_results.csv')
    neg_result_file = os.path.join(output_dir, 'negative_corr_GO_results.csv')
    
    if not (os.path.exists(pos_result_file) or os.path.exists(neg_result_file)):
        print("GO 분석 결과 파일이 없습니다. 먼저 sc_trajectory_advanced.py를 실행해주세요.")
        # 대신 기본 GO 용어 시각화 예제 실행
        print("기본 GO 용어 시각화 예제를 실행합니다.")
        create_go_network_example()
        return
    
    # 결과 파일 로드
    if os.path.exists(pos_result_file):
        print(f"양의 상관관계 GO 분석 결과 로드 중: {pos_result_file}")
        pos_go_results = pd.read_csv(pos_result_file)
        visualize_go_results(pos_go_results, "positive_correlation", go_viz_dir)
    
    if os.path.exists(neg_result_file):
        print(f"음의 상관관계 GO 분석 결과 로드 중: {neg_result_file}")
        neg_go_results = pd.read_csv(neg_result_file)
        visualize_go_results(neg_go_results, "negative_correlation", go_viz_dir)
    
    # GO 용어 계층 구조 시각화 (goatools 사용)
    try:
        print("GO 용어 계층 구조 시각화 시도 중...")
        # 수정된 부분: 반환 값 확인
        success = create_go_hierarchy_visualization()
        if not success:
            # 계층 구조 시각화가 실패한 경우에만 네트워크 예제 실행
            create_go_network_example()
    except Exception as e:
        print(f"GO 계층 구조 시각화 중 오류 발생: {e}")
        print("대체 GO 네트워크 시각화를 시도합니다.")
        create_go_network_example()

def create_go_network_example():
    """
    샘플 GO 네트워크 시각화를 생성합니다.
    (실제 GO 데이터가 없는 경우 예시로 사용)
    """
    print("GO 네트워크 시각화 예제 생성 중...")
    
    # 샘플 GO 용어 및 관계 정의
    go_terms = {
        "GO:0006351": "transcription, DNA-templated",
        "GO:0045944": "positive regulation of transcription by RNA polymerase II",
        "GO:0000122": "negative regulation of transcription by RNA polymerase II",
        "GO:0042438": "melanin biosynthetic process",
        "GO:0030318": "melanocyte differentiation",
        "GO:0008544": "epidermis development",
        "GO:0006928": "movement of cell or subcellular component",
        "GO:0048066": "developmental pigmentation"
    }
    
    # 네트워크 구조 정의 (상위-하위 관계)
    edge_connections = [
        ("GO:0006351", "GO:0045944"),
        ("GO:0006351", "GO:0000122"),
        ("GO:0030318", "GO:0042438"),
        ("GO:0008544", "GO:0030318"),
        ("GO:0006928", "GO:0008544"),
        ("GO:0048066", "GO:0042438")
    ]
    
    # 가상의 p-value 할당
    p_values = {
        "GO:0006351": 1e-5,
        "GO:0045944": 1e-6,
        "GO:0000122": 1e-4,
        "GO:0042438": 5e-5,
        "GO:0030318": 2e-5,
        "GO:0008544": 3e-4,
        "GO:0006928": 2e-4,
        "GO:0048066": 1e-3
    }
    
    # 그래프 생성
    G = nx.DiGraph()
    
    # 노드 추가
    for term_id, term_name in go_terms.items():
        G.add_node(term_id, 
                  name=term_name,
                  p_value=p_values[term_id],
                  significance=-np.log10(p_values[term_id]))
    
    # 엣지 추가
    for parent, child in edge_connections:
        G.add_edge(parent, child)
    
    # 시각화 - 수정된 부분
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # 노드 위치 계산 (계층적 레이아웃)
    pos = nx.spring_layout(G, k=0.5, iterations=50)
    
    # 노드 크기와 색상 설정 (p-value에 따라)
    node_sizes = [300 * G.nodes[node]['significance'] for node in G.nodes()]
    node_colors = [G.nodes[node]['significance'] for node in G.nodes()]
    
    # 노드 그리기
    nodes = nx.draw_networkx_nodes(G, pos, 
                                 node_size=node_sizes,
                                 node_color=node_colors,
                                 cmap=plt.cm.viridis,
                                 alpha=0.8,
                                 ax=ax)
    
    # 엣지 그리기 - 수정된 부분: 변수 이름을 변경하여 원래의 edges를 덮어쓰지 않도록 함
    edge_collection = nx.draw_networkx_edges(G, pos, 
                                  edge_color='gray',
                                  arrows=True,
                                  arrowsize=15,
                                  width=1.5,
                                  connectionstyle='arc3,rad=0.1',
                                  alpha=0.6,
                                  ax=ax)
    
    # 노드 라벨 그리기
    labels = {node: f"{go_terms[node]}\n({p_values[node]:.1e})" for node in G.nodes()}
    text_objects = nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_family='sans-serif', ax=ax)
    
    # 라벨 조정 (겹침 방지)
    adjust_text([text_objects[node] for node in G.nodes()], arrowprops=dict(arrowstyle='-', color='black'))
    
    # 컬러바 추가 - 수정된 부분
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(node_colors), vmax=max(node_colors)))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='-log10(p-value)')
    
    # 그래프 제목 및 레이아웃 설정
    ax.set_title('GO 용어 네트워크 시각화 예제')
    ax.axis('off')
    plt.tight_layout()
    
    # 저장
    plt.savefig(os.path.join(go_viz_dir, 'go_network_example.pdf'))
    plt.savefig(os.path.join(go_viz_dir, 'go_network_example.png'), dpi=300)
    plt.close()
    
    print(f"GO 네트워크 시각화 예제가 {go_viz_dir}에 저장되었습니다.")
    
    # 두 번째 시각화: 포커스 그래프와 컨텍스트 그래프
    create_focus_context_graph(go_terms, edge_connections, p_values)

def create_focus_context_graph(go_terms, edges, p_values):
    """
    이미지에서 보이는 것과 유사한 포커스 그래프와 컨텍스트 그래프를 생성합니다.
    """
    print("GO 포커스/컨텍스트 그래프 생성 중...")
    
    # 가상의 계층 레벨 설정
    go_levels = {
        "GO:0006351": 3,
        "GO:0045944": 4, 
        "GO:0000122": 4,
        "GO:0042438": 5,
        "GO:0030318": 4,
        "GO:0008544": 3,
        "GO:0006928": 2,
        "GO:0048066": 3
    }
    
    # 쿼리 용어(파란색) vs 조상 용어(보라색) 구분
    query_terms = ["GO:0045944", "GO:0000122", "GO:0042438", "GO:0030318"]
    
    # 그림 설정
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 10), gridspec_kw={'width_ratios': [1, 1]})
    
    # 1. 포커스 그래프 (왼쪽)
    G_focus = nx.DiGraph()
    
    # 노드 추가
    for term_id, term_name in go_terms.items():
        is_query = term_id in query_terms
        G_focus.add_node(term_id, 
                       name=term_name,
                       level=go_levels[term_id],
                       is_query=is_query)
    
    # 엣지 추가
    for parent, child in edges:
        G_focus.add_edge(parent, child)
    
    # 계층적 레이아웃을 구현하기 위한 노드 위치 설정
    pos_focus = {}
    x_spacing = 0.3
    level_counts = {}
    
    # 각 레벨별 노드 수 계산
    for node, data in G_focus.nodes(data=True):
        level = data['level']
        if level not in level_counts:
            level_counts[level] = 0
        level_counts[level] += 1
    
    # 노드 위치 설정
    for node, data in G_focus.nodes(data=True):
        level = data['level']
        count = level_counts[level]
        idx = sum(1 for n, d in G_focus.nodes(data=True) if d['level'] == level and list(G_focus.nodes).index(n) < list(G_focus.nodes).index(node))
        y = -level
        x = (idx - count/2) * x_spacing
        pos_focus[node] = (x, y)
    
    # 노드 색상 설정 (쿼리 vs 조상)
    node_colors = ['blue' if G_focus.nodes[node]['is_query'] else 'purple' for node in G_focus.nodes]
    
    # 그리기
    nx.draw_networkx_nodes(G_focus, pos_focus, 
                         ax=ax1,
                         node_size=300,
                         node_color=node_colors,
                         alpha=0.9)
    
    nx.draw_networkx_edges(G_focus, pos_focus, 
                         ax=ax1,
                         edge_color='gray',
                         arrows=True,
                         arrowsize=10,
                         width=1.0,
                         alpha=0.4)
    
    # 축 설정
    ax1.set_title('Focus Graph', fontsize=16)
    ax1.axis('off')
    
    # 레이어 눈금 표시
    unique_levels = sorted(set(data['level'] for _, data in G_focus.nodes(data=True)))
    y_ticks = [-level for level in unique_levels]
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels([str(level) for level in unique_levels])
    for y in y_ticks:
        ax1.axhline(y=y, color='lightgray', linestyle='--', alpha=0.5, zorder=-1)
    
    # 범례
    legend_elements = [
        Patch(facecolor='blue', label='query terms'),
        Patch(facecolor='purple', label='ancestor terms')
    ]
    ax1.legend(handles=legend_elements, loc='upper left')
    
    # 2. 컨텍스트 그래프 (오른쪽)
    # 각 레벨별 노드 수와 조상/자손 노드 수 표시
    levels = sorted(unique_levels)
    
    # 가상의 컨텍스트 데이터 생성 - 수정된 부분: 배열 길이 맞추기
    # levels 배열 길이에 맞춰 데이터 생성
    total_nodes = [50, 120, 340, 580]
    ancestral_nodes = [3, 8, 15, 25]
    
    # 배열 길이 확인 및 조정
    if len(levels) > len(total_nodes):
        # 부족한 데이터는 마지막 값으로 채우기
        total_nodes.extend([total_nodes[-1]] * (len(levels) - len(total_nodes)))
        ancestral_nodes.extend([ancestral_nodes[-1]] * (len(levels) - len(ancestral_nodes)))
    elif len(levels) < len(total_nodes):
        # 초과하는 데이터는 잘라내기
        total_nodes = total_nodes[:len(levels)]
        ancestral_nodes = ancestral_nodes[:len(levels)]
    
    context_data = {
        'level': levels,
        'total_nodes': total_nodes,
        'ancestral_nodes': ancestral_nodes
    }
    
    # 막대 그래프 그리기
    bar_width = 0.8
    bars = ax2.barh(context_data['level'], context_data['total_nodes'], 
                  height=bar_width, color='gray', label='other terms')
    
    # 조상/자손 노드 표시
    ax2.barh(context_data['level'], context_data['ancestral_nodes'], 
           height=bar_width, color='purple', label='ancestor / descendant terms')
    
    # 축 설정
    ax2.set_title('Context Graph', fontsize=16)
    ax2.set_xlabel('number of nodes per level')
    ax2.set_yticks(levels)
    ax2.set_yticklabels([str(level) for level in levels])
    
    # 레이블 추가 (노드 수와 조상 노드 수)
    for i, level in enumerate(context_data['level']):
        total = context_data['total_nodes'][i]
        ancestral = context_data['ancestral_nodes'][i]
        ax2.text(total + 50, level, f"{total} [{ancestral}]", va='center')
    
    # 범례
    ax2.legend(loc='upper right')
    
    # GO 용어 설명 추가
    ax1.text(0.05, -len(levels)-1, '\n'.join([
        f"i. {go_terms['GO:0045944']}",
        f"ii. {go_terms['GO:0000122']}",
        f"iii. {go_terms['GO:0006351']}",
        f"iv. {go_terms['GO:0042438']}",
        f"v. {go_terms['GO:0030318']}",
        f"vi. {go_terms['GO:0008544']}",
        f"vii. {go_terms['GO:0006928']}",
        f"viii. {go_terms['GO:0048066']}"
    ]), fontsize=9, va='top')
    
    # 노드에 로마자 라벨 추가
    roman_labels = {
        "GO:0045944": "i",
        "GO:0000122": "ii",
        "GO:0006351": "iii",
        "GO:0042438": "iv",
        "GO:0030318": "v",
        "GO:0008544": "vi",
        "GO:0006928": "vii",
        "GO:0048066": "viii"
    }
    
    for node, pos in pos_focus.items():
        if node in roman_labels:
            ax1.text(pos[0], pos[1], roman_labels[node], fontsize=10, 
                   ha='center', va='center', color='white')
    
    plt.tight_layout()
    plt.savefig(os.path.join(go_viz_dir, 'go_focus_context_graph.pdf'))
    plt.savefig(os.path.join(go_viz_dir, 'go_focus_context_graph.png'), dpi=300)
    plt.close()
    
    print(f"GO 포커스/컨텍스트 그래프가 {go_viz_dir}에 저장되었습니다.")

def visualize_go_results(go_results, result_type, output_dir):
    """
    GO 분석 결과를 시각화합니다.
    (실제 분석 결과를 사용)
    """
    # 결과 확인
    if go_results.empty:
        print(f"{result_type} GO 분석 결과가 비어 있습니다.")
        return
    
    print(f"{result_type} GO 분석 결과 시각화 중...")
    
    # 열 이름 추출
    p_value_col = next((col for col in go_results.columns if 'p-value' in col.lower() or 'pvalue' in col.lower()), None)
    term_col = next((col for col in go_results.columns if 'term' in col.lower() or 'name' in col.lower() or 'des' in col.lower()), go_results.columns[0])
    
    if p_value_col is None:
        print("P-value 열을 찾을 수 없습니다.")
        return
    
    # 상위 15개 GO 용어 선택
    top_go = go_results.sort_values(by=p_value_col).head(15)
    
    # GO 용어 이름 정리 (너무 긴 경우 자르기)
    top_go['short_term'] = top_go[term_col].apply(lambda x: x[:60] + '...' if len(str(x)) > 60 else x)
    
    # 시각화
    plt.figure(figsize=(12, 8))
    
    # -log10(p-value)로 변환
    log_p_values = -np.log10(top_go[p_value_col])
    
    # 막대 그래프
    bars = plt.barh(top_go['short_term'], log_p_values, color=sns.color_palette("viridis", len(top_go)))
    
    # 축 및 제목 설정
    plt.xlabel('-log10(p-value)')
    plt.ylabel('GO Term')
    plt.title(f'Top GO Terms ({result_type})')
    
    # p-value 표시
    for i, bar in enumerate(bars):
        plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, 
                f"{top_go.iloc[i][p_value_col]:.2e}", 
                va='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'go_{result_type}_top_terms.pdf'))
    plt.savefig(os.path.join(output_dir, f'go_{result_type}_top_terms.png'), dpi=300)
    plt.close()
    
    # GO 방사형 네트워크 시각화 생성
    create_go_radial_network(top_go, term_col, p_value_col, result_type, output_dir)

def create_go_radial_network(go_data, term_col, p_value_col, result_type, output_dir):
    """
    GO 용어의 방사형 네트워크 시각화를 생성합니다.
    """
    print(f"{result_type} GO 방사형 네트워크 생성 중...")
    
    G = nx.Graph()
    
    # 중심 노드 추가
    G.add_node('center', name='Center', p_value=1.0)
    
    # GO 용어 노드 추가
    for idx, row in go_data.iterrows():
        term_id = f"term_{idx}"
        G.add_node(term_id, 
                  name=str(row[term_col]),  # 문자열로 변환
                  p_value=row[p_value_col],
                  significance=-np.log10(row[p_value_col]))
        
        # 중심과 연결
        G.add_edge('center', term_id)
    
    # 시각화
    fig, ax = plt.subplots(figsize=(14, 14))  # 명시적으로 figure와 axes 생성
    
    # 방사형 레이아웃
    pos = nx.spring_layout(G, k=0.5)
    
    # 중심 노드는 중앙에 위치
    pos['center'] = np.array([0, 0])
    
    # 노드 크기와 색상 설정
    # 중심 노드
    nx.draw_networkx_nodes(G, pos, 
                         nodelist=['center'],
                         node_size=1000,
                         node_color='lightgray',
                         edgecolors='black',
                         ax=ax)  # ax 지정
    
    # GO 용어 노드
    term_nodes = [node for node in G.nodes() if node != 'center']
    node_sizes = [300 * G.nodes[node].get('significance', 1) for node in term_nodes]
    node_colors = [G.nodes[node].get('significance', 1) for node in term_nodes]
    
    nodes = nx.draw_networkx_nodes(G, pos, 
                                 nodelist=term_nodes,
                                 node_size=node_sizes,
                                 node_color=node_colors,
                                 cmap=plt.cm.viridis,
                                 alpha=0.8,
                                 ax=ax)  # ax 지정
    
    # 엣지 그리기
    nx.draw_networkx_edges(G, pos, 
                         width=1.5,
                         alpha=0.6,
                         edge_color='lightgray',
                         ax=ax)  # ax 지정
    
    # 노드 라벨 그리기
    # 중심 노드
    nx.draw_networkx_labels(G, pos, 
                          labels={'center': 'GO Terms'},
                          font_size=12,
                          font_weight='bold',
                          ax=ax)  # ax 지정
    
    # GO 용어 노드 (짧은 이름으로)
    term_labels = {node: str(G.nodes[node]['name'])[:20] + '...' if len(str(G.nodes[node]['name'])) > 20 else str(G.nodes[node]['name']) 
                 for node in term_nodes}
    nx.draw_networkx_labels(G, pos, 
                          labels=term_labels,
                          font_size=8,
                          ax=ax)  # ax 지정
    
    # 컬러바 추가 - 수정된 부분
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(node_colors), vmax=max(node_colors)))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='-log10(p-value)')  # ax 인자 추가
    
    # 그래프 제목 및 레이아웃 설정
    ax.set_title(f'GO 용어 방사형 네트워크 ({result_type})')  # plt.title 대신 ax.set_title 사용
    ax.axis('off')  # plt.axis 대신 ax.axis 사용
    plt.tight_layout()
    
    # 저장
    plt.savefig(os.path.join(output_dir, f'go_{result_type}_radial_network.pdf'))
    plt.savefig(os.path.join(output_dir, f'go_{result_type}_radial_network.png'), dpi=300)
    plt.close()

def create_go_hierarchy_visualization():
    """
    GO 계층 구조의 시각화를 생성합니다. (goatools 사용)
    이 함수는 goatools 라이브러리가 설치되어 있어야 합니다.
    """
    print("GO 계층 구조 시각화 생성 중...")
    
    try:
        # Graphviz 설치 확인
        try:
            import subprocess
            result = subprocess.run(['dot', '-V'], capture_output=True, text=True, check=False)
            if result.returncode != 0:
                print("Graphviz 'dot' 명령어를 찾을 수 없습니다. Graphviz를 설치해주세요.")
                print("대체 GO 네트워크 시각화를 실행합니다.")
                return False
        except FileNotFoundError:
            print("Graphviz 'dot' 명령어를 찾을 수 없습니다. Graphviz를 설치해주세요.")
            print("대체 GO 네트워크 시각화를 실행합니다.")
            return False
        
        # 1. GO OBO 파일 다운로드
        obo_file = os.path.join(data_dir, 'go-basic.obo')
        if not os.path.exists(obo_file):
            print(f"GO OBO 파일을 다운로드합니다: {obo_file}")
            import urllib.request
            urllib.request.urlretrieve('http://purl.obolibrary.org/obo/go/go-basic.obo', obo_file)
        
        # 2. GO DAG 파싱 - 수정된 부분: relationship=True 추가
        print("GO DAG 파싱 중...")
        go_dag = obo_parser.GODag(obo_file, optional_attrs=['relationship'])
        
        # 3. 관심있는 GO 용어 설정 (세포 분화 관련)
        go_ids = [
            'GO:0045944',  # positive regulation of transcription by RNA polymerase II
            'GO:0000122',  # negative regulation of transcription by RNA polymerase II
            'GO:0030318',  # melanocyte differentiation
            'GO:0042438',  # melanin biosynthetic process
            'GO:0008544',  # epidermis development
        ]
        
        # 4. GO SubDAG 생성
        print("GO SubDAG 생성 중...")
        go_subdag = GoSubDag(go_ids, go_dag, relationships=True, prt=None)
        
        # 5. Plot 객체 생성
        go_plot = GoSubDagPlot(go_subdag)
        
        # 6. 계층 구조 시각화
        print("GO 계층 구조 그리기...")
        output_file = os.path.join(go_viz_dir, 'GO_hierarchy')
        try:
            go_plot.plt_dag(output_file)
            print(f"GO 계층 구조가 {go_viz_dir}에 저장되었습니다.")
            return True
        except Exception as e:
            print(f"GO DAG 시각화 생성 오류: {str(e)}")
            print("대체 GO 네트워크 시각화를 실행합니다.")
            return False
            
    except Exception as e:
        print(f"GO 계층 구조 시각화 중 오류 발생: {str(e)}")
        print("대체 GO 네트워크 시각화를 실행합니다.")
        return False

if __name__ == "__main__":
    run_go_visualization()
