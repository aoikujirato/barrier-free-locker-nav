import streamlit as st
import geopandas as gpd
import pandas as pd
import networkx as nx
import folium
from streamlit_folium import st_folium
from shapely.geometry import Point, LineString
from streamlit_js_eval import get_geolocation

# 1. ページ設定とスタイル
st.set_page_config(page_title="ついでにロッカー", layout="wide")

st.markdown("""
    <style>
    .main-header {
        background-color: #2E7D32;
        padding: 20px;
        border-radius: 10px;
        text-align: center;
        margin-bottom: 10px;
    }
    .header-top { color: white; font-size: 1.2rem; font-weight: bold; margin-bottom: 5px; }
    .header-title { color: white; font-size: 1.8rem; font-weight: bold; margin: 0; }
    .sub-description {
        color: #333333; font-size: 1.3rem; font-weight: bold;
        text-align: center; margin-top: 15px; margin-bottom: 25px;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown('<div class="main-header"><div class="header-top">ついでにロッカー</div><div class="header-title">新宿西口駅：バリアフリー・ロッカー経由ナビ</div></div>', unsafe_allow_html=True)
st.markdown('<div class="sub-description">新宿西口駅：コインロッカー経由＆バリアフリールートの経路検索</div>', unsafe_allow_html=True)

message_placeholder = st.empty()
if 'calc_result' not in st.session_state:
    st.session_state.calc_result = None

if st.session_state.calc_result:
    message_placeholder.success("最適ルートを表示します (ルート案内はマップ下に表示)")
else:
    message_placeholder.info("左のサイドバーを開き、検索設定から出発地と目的地を選択して「経路を検索」を押してください")

# 2. ロジック関数とデータ読み込み
def get_rank_weight(rank_str, mode_text):
    try:
        if not isinstance(rank_str, str) or len(rank_str) < 3: return 100.0
        step_rank = rank_str[2].upper()
        if "バリアフリー" in mode_text:
            weights = {'S': 1.0, 'A': 1.0, 'B': 5.0, 'C': 50.0, 'Z': 9999.0, 'X': 9999.0}
        else:
            weights = {'S': 1.0, 'A': 1.0, 'B': 1.2, 'C': 2.0, 'Z': 5.0, 'X': 5.0}
        return weights.get(step_rank, 100.0)
    except: return 100.0

@st.cache_resource
def load_raw_data():
    links = gpd.read_file("link_merge.geojson")
    nodes = gpd.read_file("node_merge.geojson")
    lockers = pd.read_csv("lockers.csv")
    
    try:
        locations_df = pd.read_csv("location.csv")
        locations_df['display_label'] = locations_df.apply(
            lambda r: f"{r['name']} ({r['lat']:.5f}, {r['lon']:.5f})", axis=1
        )
    except:
        locations_df = pd.DataFrame()

    nodes['node_id'] = nodes['node_id'].astype(str)
    links['start_id'] = links['start_id'].astype(str)
    links['end_id'] = links['end_id'].astype(str)
    return links, nodes, lockers, locations_df

links_df, nodes_gdf, lockers_df, locations_df = load_raw_data()

def get_nearest_node(lat, lon, floor):
    same_floor = nodes_gdf[nodes_gdf['floor'] == floor]
    target_nodes = same_floor if not same_floor.empty else nodes_gdf
    dist = target_nodes.geometry.distance(Point(lon, lat))
    return str(target_nodes.iloc[dist.argmin()]['node_id'])

def build_graph(mode_text):
    G = nx.Graph()
    for _, node in nodes_gdf.iterrows():
        G.add_node(str(node['node_id']), floor=node['floor'], lat=node['lat'], lon=node['lon'])
    
    for _, row in links_df.iterrows():
        r_type = str(row.get('route_type', '1'))
        rank = row.get('rank', 'SSS')
        step = rank[2].upper() if len(rank) >= 3 else 'S'
        
        if "バリアフリー" in mode_text:
            if r_type == '6' or step in ['Z', 'X']: continue 
        
        dist = float(row.get('distance', 0))
        if dist == 0 and r_type == '4': dist = 5.0 
        w = get_rank_weight(rank, mode_text)
        G.add_edge(str(row['start_id']), str(row['end_id']), weight=dist * w, actual_dist=dist, geometry=row.geometry, route_type=r_type, rank=rank)
    
    nodes_gdf['pos_key'] = nodes_gdf.apply(lambda r: (round(r['lat'], 4), round(r['lon'], 4)), axis=1)
    for _, group in nodes_gdf.groupby('pos_key'):
        if len(group) > 1:
            node_list = group.to_dict('records')
            for i in range(len(node_list)):
                for j in range(i + 1, len(node_list)):
                    n1, n2 = node_list[i], node_list[j]
                    if n1['floor'] != n2['floor']:
                        G.add_edge(str(n1['node_id']), str(n2['node_id']), weight=20.0, actual_dist=5.0, route_type='4', rank='SSS', geometry=LineString([(n1['lon'], n1['lat']), (n2['lon'], n2['lat'])]))
    return G

#3. UI（サイドバー）
if 'gps_coords' not in st.session_state:
    st.session_state.gps_coords = None

with st.sidebar:
    st.header("🔍 検索設定")
    search_mode = st.radio("移動モードを選択", ("バリアフリーモード（段差をなるべく避ける）", "標準モード（多少の階段はOK）"))
    st.markdown("---")
    
    if not locations_df.empty:
        loc_options = ["", "現在地 (GPS)"] + locations_df['display_label'].tolist()
        start_label = st.selectbox("出発地を選択", loc_options, index=0)
        end_label = st.selectbox("目的地を選択", loc_options, index=0)

        # 1. 現在地が選ばれた時にGPSを起動
        if start_label == "現在地 (GPS)" or end_label == "現在地 (GPS)":
            # keyを使わずに呼び出し
            loc = get_geolocation()

            if loc:
                # 緯度・経度の取得
                lat = loc['coords']['latitude']
                lon = loc['coords']['longitude']
                
                # 前回の取得内容と異なる場合（または初回）のみ保存して再読み込み
                if st.session_state.gps_coords is None or \
                   abs(st.session_state.gps_coords['lat'] - lat) > 0.0001:
                    st.session_state.gps_coords = {"lat": lat, "lon": lon}
                    st.rerun() 
            else:
                st.info("位置情報の許可を待っています。ポップアップが出ない場合は、ブラウザのアドレスバー左にある「鍵マーク」から位置情報を許可し、ページを再読み込みしてください。")

        st.markdown("---")
        # （ロッカーサイズ選択部分はそのまま）
        st.markdown('<span style="font-size: 14px;">ロッカーサイズを選択</span>', unsafe_allow_html=True)
        selected_sizes = []
        if st.checkbox("小サイズ", value=True): selected_sizes.append('small')
        if st.checkbox("中サイズ"): selected_sizes.append('medium')
        if st.checkbox("大サイズ"): selected_sizes.append('large')
        st.markdown("---")

# --- 地点確定ロジック ---
        # 最初に「ボタンはまだ押されていない/表示しない」状態にしておく（重要！）
        search_btn = False
        start_node = None
        end_node = None
        is_gps_waiting = False 

        if start_label == "" or end_label == "":
            st.warning("出発地と目的地を選択してください。")
        else:
            # 1. 出発地の特定
            if start_label == "現在地 (GPS)":
                if st.session_state.gps_coords:
                    gps = st.session_state.gps_coords
                    start_node = get_nearest_node(gps['lat'], gps['lon'], 0)
                else:
                    is_gps_waiting = True
            else:
                start_info = locations_df[locations_df['display_label'] == start_label].iloc[0]
                start_node = get_nearest_node(start_info['lat'], start_info['lon'], start_info['floor'])

            # 2. 目的地の特定
            if end_label == "現在地 (GPS)":
                if st.session_state.gps_coords:
                    gps = st.session_state.gps_coords
                    end_node = get_nearest_node(gps['lat'], gps['lon'], 0)
                else:
                    is_gps_waiting = True
            else:
                end_info = locations_df[locations_df['display_label'] == end_label].iloc[0]
                end_node = get_nearest_node(end_info['lat'], end_info['lon'], end_info['floor'])

            # 3. ボタンの表示判定
            if is_gps_waiting:
                # GPS待ちのときはボタンを出さずにメッセージを出す
                st.info("現在地を測位しています。少々お待ちください...")
            elif start_node and end_node:
                # 両方のノードが確定した時だけボタンを生成し、結果を search_btn に入れる
                search_btn = st.button("この条件で経路を検索", use_container_width=True)


# 4. 計算と表示 
if 'calc_result' not in st.session_state:
    st.session_state.calc_result = None

if search_btn:
    with st.spinner("条件に合うロッカーとバリアフリー経路を計算中..."):
        G = build_graph(search_mode)
        best_cost = float('inf')

        if selected_sizes:
            filtered_lockers = lockers_df[lockers_df[selected_sizes].gt(0).any(axis=1)]
        else:
            filtered_lockers = lockers_df 

        if filtered_lockers.empty:
            st.error("選択されたサイズの空きがあるロッカーが見つかりませんでした。")
        else:
            final_path, final_locker, final_l_node = None, None, None
            for _, lkr in filtered_lockers.iterrows():
                try:
                    l_pt = Point(lkr['lon'], lkr['lat'])
                    lkr_floor = lkr.get('floor', 0)
                    same_floor_nodes = nodes_gdf[nodes_gdf['floor'] == lkr_floor]
                    target_nodes = same_floor_nodes if not same_floor_nodes.empty else nodes_gdf
                    distances = target_nodes.geometry.distance(l_pt)
                    l_node_id = str(target_nodes.iloc[distances.argmin()]['node_id'])
                    
                    cost1 = nx.shortest_path_length(G, start_node, l_node_id, weight='weight')
                    path1 = nx.shortest_path(G, start_node, l_node_id, weight='weight')
                    cost2 = nx.shortest_path_length(G, l_node_id, end_node, weight='weight')
                    path2 = nx.shortest_path(G, l_node_id, end_node, weight='weight')
                    
                    if (cost1 + cost2) < best_cost:
                        best_cost = cost1 + cost2
                        final_path = path1[:-1] + path2
                        final_locker = lkr
                        final_l_node = l_node_id
                except: continue

            if final_path:
                st.session_state.calc_result = {
                    "path": final_path, "locker": final_locker, 
                    "l_node": final_l_node, "mode": search_mode, "graph": G
                }
                st.success("最適ルートを表示中 (ルート案内はマップ下に表示)")
            else:
                st.error("経路が見つかりませんでした。")

# C. 描画処理
if st.session_state.calc_result:
    res = st.session_state.calc_result
    path, locker, l_node, current_G = res["path"], res["locker"], res["l_node"], res["graph"]

    # 1. 地図の初期化
    m = folium.Map(zoom_start=18)
    all_coordinates = []

    # 階層に応じた色定義
    def get_floor_color(f_val):
        try:
            f = int(round(float(f_val)))
            colors = {
                1: "#FF4B4B",   # 地上階（赤）
                0: "#3D91FF",   # 地上 0F (青)
                -1: "#4CAF50",  # B1F (緑)
                -2: "#FF9800",  # B2F (オレンジ)
                -3: "#9C27B0",  # B3F (紫)
                -4: "#212121"   # B4F以下 (黒)
            }
            return colors.get(f, "#808080")
        except: return "#808080"

    # --- 2. 地点の表示 ---

    start_node_info = nodes_gdf[nodes_gdf['node_id'] == path[0]].iloc[0]
    if start_label == "現在地 (GPS)":
        # 現在地（青い丸）
        folium.CircleMarker(
            location=[start_node_info.lat, start_node_info.lon],
            radius=7, color='white', weight=2, fill=True, fill_color='#1A73E8', fill_opacity=1,
            tooltip="現在地"
        ).add_to(m)
    else:
        # 通常の出発地ピン（緑）
        folium.Marker( [start_node_info.lat, start_node_info.lon], popup="出発地", icon=folium.Icon(color='green', icon='play', prefix='fa')).add_to(m)
    all_coordinates.append([start_node_info.lat, start_node_info.lon])
    
    # 目的地（赤ピン）
    end_node_info = nodes_gdf[nodes_gdf['node_id'] == path[-1]].iloc[0]
    folium.Marker( [end_node_info.lat, end_node_info.lon], popup="目的地", icon=folium.Icon(color='red', icon='flag', prefix='fa')).add_to(m)
    all_coordinates.append([end_node_info.lat, end_node_info.lon])

    # ロッカー（階層別の色でピンを表示）
    lkr_f_raw = float(locker.get('floor', 0))
    lkr_f_int = int(round(lkr_f_raw))
    # folium.Iconが対応している基本色にマッピング
    icon_color_map = {0: 'blue', -1: 'green', -2: 'orange', -3: 'purple', -4: 'black'}
    l_color = icon_color_map.get(lkr_f_int, 'red' if lkr_f_int >= 1 else 'gray')

    folium.Marker(
    [locker.lat, locker.lon], popup=f"ロッカー: {locker.get('name')} ({lkr_f_raw}F)", icon=folium.Icon(color=l_color, icon='briefcase', prefix='fa')).add_to(m)
    all_coordinates.append([locker.lat, locker.lon])

    # --- 3. 経路の描画（階層別カラー） ---
    for u, v in zip(path[:-1], path[1:]):
        edge = current_G.get_edge_data(u, v)
        if edge and 'geometry' in edge:
            f_raw = current_G.nodes[u].get('floor', 0)
            f_color = get_floor_color(f_raw)
            coords = [(lat, lon) for lon, lat in edge['geometry'].coords]
            folium.PolyLine(
                coords, color=f_color, weight=6, opacity=0.8,
                tooltip=f"階数: {f_raw}F"
            ).add_to(m)
            all_coordinates.extend(coords)

    # --- 4. 階層凡例の追加（HTML） ---
    legend_html = '''
         <div style="position: fixed; 
                     bottom: 50px; left: 50px; width: 140px; height: 160px; 
                     border:2px solid grey; z-index:9999; font-size:12px;
                     background-color:white; opacity: 0.9; padding: 10px; border-radius: 5px;">
         <b>階層凡例</b><br>
         <span style="color:#FF4B4B">■</span> 地上階<br>
         <span style="color:#3D91FF">■</span> 地上 0F<br>
         <span style="color:#4CAF50">■</span> 地下 B1F<br>
         <span style="color:#FF9800">■</span> 地下 B2F<br>
         <span style="color:#9C27B0">■</span> 地下 B3F<br>
         <span style="color:#212121">■</span> 地下 B4F~<br>
         </div>
         '''
    m.get_root().html.add_child(folium.Element(legend_html))

    # 全地点が収まるように調整
    if all_coordinates:
        m.fit_bounds(all_coordinates)
    
    # 地図のレンダリング
    st_folium(m, width=900, height=450, returned_objects=[], key="navigation_map")

    # ルート案内
    st.markdown("---")
    st.subheader("ルート案内")
    ROUTE_TYPE_NAMES = {'2':"動く歩道",'4':"エレベーター",'5':"エスカレーター",'6':"⚠️ 階段",'7':"スロープ"}
    guidance, total_m = [], 0
    for u, v in zip(path[:-1], path[1:]):
        if u == l_node:
            guidance.append({"msg": f"📍 ロッカー「{locker.get('name','未設定')}」に到着", "dist": 0, "is_locker": True, "small": locker.get('small',0), "medium": locker.get('medium',0), "large": locker.get('large',0)})
        edge = current_G.get_edge_data(u, v)
        if edge:
            d = edge.get('actual_dist', 0)
            total_m += d
            rt = edge.get('route_type', '1')
            msg = f"{ROUTE_TYPE_NAMES[rt]}を利用" if rt in ROUTE_TYPE_NAMES else "道なりに進む"
            if guidance and not guidance[-1].get('is_locker') and guidance[-1]['msg'] == msg:
                guidance[-1]['dist'] += d
            else:
                guidance.append({"msg": msg, "dist": d})

    st.write(f"**総移動距離: 約 {total_m:.1f}m**")
    for i, g in enumerate(guidance):
        dist_info = f" **({g['dist']:.1f}m)**" if g['dist'] > 0 else ""
        if g.get('is_locker'):
            st.info(f"**Step {i+1}: {g['msg']}**\n\n（ロッカーの個数：大 {g['large']}個, 中 {g['medium']}個, 小 {g['small']}個）")
        else:

            st.write(f"**Step {i+1}:** {g['msg']}{dist_info}")

# --- アプリの最下部（フッター） ---
st.markdown("---") # 区切り線

# 注意事項
st.caption("注意事項")
st.caption("""
- アプリケーションに利用している公共交通データは、公共交通オープンデータチャレンジにおいて提供された「歩行空間ネットワークデータ(都営地下鉄大江戸線 新宿西口駅)」(国土交通省)を加工して作成しております。(2025年12月末時点のデータ)
- また、コインロッカーや駅地点データは公開情報を参照して製作者が自ら作成したものを利用しております。
- 本アプリの案内は、実際の現地の状況や混雑、工事等により異なる場合があります。必ず現地の案内看板や係員の指示に従ってください。
- 表示内容について、公共交通事業者への直接の問い合わせは行わないでください。
- 歩きスマホは大変危険です。地図を確認する際は、周囲の安全を確認した上で立ち止まってご利用ください。
""")

# お問い合わせ（小さな文字で表示）
st.markdown(
    """
    <div style="text-align: center; color: #888888; font-size: 0.8rem; margin-top: 30px;">
        © 2026 ついでにロッカー<br>
        お問い合わせ: <a href="mailto:aoikujirato [at] gmail.com" style="color: #888888;">aoikujirato [at] gmail.com ※[at] を@に置き換えてください</a><br>
    </div>
    """,
    unsafe_allow_html=True
)














