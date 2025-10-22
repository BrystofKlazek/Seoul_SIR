# visualize clusters on a map
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from shapely import set_precision
import maptograph as mtg
import numpy as np
import nametocode as ntc
import argparse
from itertools import product


parser = argparse.ArgumentParser(description="My parser")
parser.add_argument('--map', action=argparse.BooleanOptionalAction)
parser.add_argument('graph_style', type = str, nargs = "?", default="simple")
args = parser.parse_args()


# Hangul font setting
plt.rcParams["font.family"] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False

for name in ["Noto Sans CJK KR", "Noto Sans KR", "NanumGothic", "Malgun Gothic", "Apple SD Gothic Neo"]:
    if any(ft.name == name for ft in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = name
        break



# 0) Load + project
shapefile_path = './shp/202101/SEOUL_SIG.shp'
seoul_map = gpd.read_file(shapefile_path).to_crs(5179)

# 1) One inside point for plotting
anchors = seoul_map.representative_point()
seoul_map["x"], seoul_map["y"] = anchors.x, anchors.y

GRID = 0.2   # 20 cm grid; adjust 0.1–0.5 if needed
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

