# visualize clusters on a map
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from shapely import set_precision
import maptograph as mtg
import numpy as np
import nametocode as ntc
from itertools import product

def graph_laplacian_sparse(graph):
    indptr = [0]
    indices = []
    data = []

    for node in graph:
        deg = 0.0
        def add_edge(neighbour, weight):
            nonlocal deg
            indices.append(neighbour)
            data.append(-weight)
            deg += weight
        
        for neighbour in graph.neighbors(node):
            add_edge(neighbour, graph.get_edge_data(node, neighbour)["weight"])

        indices.append(node)
        data.append(deg)

        indptr.append(len(indices))

    return (np.asarray(indptr, dtype=np.int32),
            np.asarray(indices, dtype=np.int32),
            np.asarray(data, dtype=np.float64))

def sir_rhs(state, out, beta=0.8, gamma=0.2):
    S, I, R = state[0], state[1], state[2]
    dS, dI, dR = out[0], out[1], out[2]
    inf = beta * S * I
    dS[...] = -inf
    dI[...] =  inf - gamma * I
    dR[...] =  gamma * I

def main():
    dt = 0.02
    steps = 1000
    snapshots = 50
    beta, gamma = 0.8, 0.2

    GRID = 0.2   
    #-------------HELP FOR DISPLAYING KOREAN FONTS------------#
    plt.rcParams["font.family"] = 'NanumGothic'
    plt.rcParams['axes.unicode_minus'] = False

    for name in ["Noto Sans CJK KR", "Noto Sans KR", "NanumGothic", 
                 "Malgun Gothic", "Apple SD Gothic Neo"]:
        if any(ft.name == name for ft in fm.fontManager.ttflist):
            plt.rcParams["font.family"] = name
            break
    #---------------------------------------------------------#

    shapefile_path = './shp/202101/SEOUL_SIG.shp'
    seoul_map = gpd.read_file(shapefile_path).to_crs(5179)

    anchors = seoul_map.representative_point()
    seoul_map["x"], seoul_map["y"] = anchors.x, anchors.y

    seoul_map["geometry"] = seoul_map.geometry.buffer(0).apply(
            lambda geom: set_precision(geom, GRID)
            )

    name_code_df = pd.read_csv("code_lookup.csv")
    name_dict = ntc.code_dict(code_df = name_code_df)

    codes = name_code_df["sgg"].to_list()
    edge_weights = {pair : 2*i 
                     for i, pair in enumerate(product(codes, codes))}

    G = mtg.maptograph(seoul_map, mode = "neighbours", pairs=edge_weights)

    seoul = mtg.graphDisplay(G, seoul_map)
    seoul.interactive_graph()
    plt.show()

