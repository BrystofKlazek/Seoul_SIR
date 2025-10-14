# visualize clusters on a map
from re import I
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from shapely import set_precision
import maptograph as mtg
import argparse

parser = argparse.ArgumentParser(description="My parser")
parser.add_argument('--draw_map', action=argparse.BooleanOptionalAction)
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
seoul_map["geometry"] = seoul_map.geometry.buffer(0).apply(lambda geom: set_precision(geom, GRID))

G = mtg.maptograph(seoul_map)

seoul = mtg.graphDisplay(G, seoul_map)
seoul.interactive_graph(args.draw_map)
plt.show()

