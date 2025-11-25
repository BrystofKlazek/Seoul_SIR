import sys
from pathlib import Path

import json
import time
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from shapely import set_precision
import maptograph as mtg
import numpy as np
import nametocode as ntc
import networkx as nx
import matrix_builders as mb

# add C solver
this_dir = Path(__file__).resolve().parent
sys.path.append(str(this_dir.parent / "differential_solver_C"))

import SIR_diffusion_wrap as sir  

NORM = 5000000
TIME_STEP = 0.25
NUM_DAYS = 31    
NUM_SNAPSHOTS = 150    
BETA, GAMMA = 0.03, 0.01
SEED_IDX = 11020

def build_hourly_weight_table(weights_by_hour_raw, norm=NORM):
    #The raw json is in the structure of
    #weights_by_hour_raw[hour_str][u_str][v_str] = weight
    #I turn that to 
    #hourly[(u, v)] = np.array with each hour, already divided by norm
    hourly = {} 

    for h_str, by_u in weights_by_hour_raw.items():
        h = int(h_str)
        for u_str, neighs in by_u.items():
            u = int(u_str)
            for v_str, w in neighs.items():
                v = int(v_str)
                key = (u, v)
                arr = hourly.get(key)
                if arr is None:
                    arr = np.zeros(24, dtype=np.float64)
                    hourly[key] = arr
                arr[h] = float(w) / norm

    return hourly


def make_edge_weight_fn(hourly_weights, dt_snapshot):
    # Here I build a small helper function to be returned- 
    #edge_weight_fn(time, u, v)  -> weight from u to v at that frame,
    #built using the precomputed 24-hour arrays.
    # I prebuild these functions for speed.

    def edge_weight_fn(frame, u, v):
        t = frame * dt_snapshot
        h = t % 24.0
        h0 = int(h)
        h1 = (h0 + 1) % 24
        frac = h - h0

        arr = hourly_weights.get((int(u), int(v)))
        if arr is None:
            return 0.0

        w0 = arr[h0]
        w1 = arr[h1]
        return (1.0 - frac) * w0 + frac * w1

    return edge_weight_fn


def sir_rhs(state, out, beta=0.4, gamma=0.1):
    S, I, R = state[0], state[1], state[2]
    dS, dI, dR = out[0], out[1], out[2]
    inf = beta * S * I
    dS[...] = -inf
    dI[...] = inf - gamma * I
    dR[...] = gamma * I


def build_hourly_tensor_for_C(hourly_weights_dict, G, n_hours=24):
    #Here, I construct an 24 dim array of n by n (n is size of graph)
    # arrays to pass into C, where it is weights between nodes by hour.
    node_index_attr = nx.get_node_attributes(G, "index")
    nodes = sorted(G.nodes(), key=lambda u: node_index_attr[u])

    n = len(nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}

    # allocating the tensor
    hourly_arr = np.zeros((n_hours, n, n), dtype=np.float64)

    # Fill with rates from the dict - the keys are (u,v) codes
    for (u, v), arr in hourly_weights_dict.items():
        if u in node_to_idx and v in node_to_idx:
            i = node_to_idx[u]
            j = node_to_idx[v]
            hourly_arr[:, i, j] = arr

    return hourly_arr, nodes, node_to_idx


def main():
    # time step etc
    dt = TIME_STEP
    steps = int(24*NUM_DAYS/dt)
    snapshots = NUM_SNAPSHOTS
    beta, gamma = BETA, GAMMA

    GRID = 0.2

    # korean fonts (pick first one that exists)
    plt.rcParams["font.family"] = "NanumGothic"
    plt.rcParams["axes.unicode_minus"] = False
    for name in ["Noto Sans CJK KR", "Noto Sans KR", "NanumGothic",
                 "Malgun Gothic", "Apple SD Gothic Neo"]:
        if any(ft.name == name for ft in fm.fontManager.ttflist):
            plt.rcParams["font.family"] = name
            break

    # map + centroids
    shapefile_path = "./shp/202101/SEOUL_SIG.shp"
    seoul_map = gpd.read_file(shapefile_path).to_crs(5179)

    anchors = seoul_map.representative_point()
    seoul_map["x"], seoul_map["y"] = anchors.x, anchors.y

    GRID = 0.2
    seoul_map["geometry"] = seoul_map.geometry.buffer(0).apply(
        lambda geom: set_precision(geom, GRID)
    )

    # code table and OD data
    name_code_df = pd.read_csv("code_lookup.csv")
    name_dict = ntc.code_dict(code_df=name_code_df)

    codes = name_code_df["sgg"].to_list()
    with open("weights_03.json", encoding="utf-8") as f:
        weights_by_hour_raw = json.load(f)

    # compress per-hour json into arrays - the pipeline to the array we want
    hourly_weights = build_hourly_weight_table(weights_by_hour_raw, norm=NORM)

    # for the graph itself just use hour 0
    # We iterate over all arrs gotten with the key (u, v)
    edge_weights0 = {}
    for (u, v), arr in hourly_weights.items():
        w0 = arr[0]
        if w0 > 0.0:
            edge_weights0[(u, v)] = w0

    # build graph with static weights = hour 0 
    # It is prebuilt so that we know the size of it for the hourly tensors,
    # and for the animation. 
    G = mtg.maptograph(seoul_map, mode="from_file", pairs=edge_weights0)

    # ---- NEW: build dense hourly tensor for C solver ----
    hourly_tensor, nodes_order, node_to_idx = build_hourly_tensor_for_C(
        hourly_weights, G, n_hours=24
    )

    n = len(nodes_order)

    # initial S, I, R in the SAME node ordering as hourly_tensor
    S0 = np.ones(n, np.float64)
    I0 = np.zeros(n, np.float64)
    R0 = np.zeros(n, np.float64)

    # seed - district with code SEED_IDX (e.g. 11010)
    seed_code = SEED_IDX
    seed_idx = node_to_idx[seed_code]

    I0[seed_idx] = 0.1
    S0[seed_idx] -= I0[seed_idx]

    x0_graph = np.stack([S0, I0, R0], axis=0)

    # run Euler on the graph with hourly dense Laplacian + interpolation
    print("Solving SIR on Seoul graph (dense hourly, Euler, reaction-diffusion SIR)...")
    t0 = time.time()
    out_graph = sir.euler_solve_graph_rd(
        hourly_tensor,
        x0_graph,
        0.0,          # t0_hours, I start at time 0
        dt,
        steps,
        snapshots,
        rhs_fn=sir_rhs,
        params={"beta": beta, "gamma": gamma}
    )
    print(f"Graph solve done in {time.time() - t0:.2f}s, "
          f"snaps = {out_graph.shape[0]}")

    out = out_graph
    print("first snapshot:", out[0, :, 0])
    print("middle snapshot:", out[out.shape[0] // 2, :, 0])
    print("last snapshot:", out[-1, :, 0])
    print("I min/max over all:", out[:, 1, :].min(), out[:, 1, :].max())

    values_ts = out_graph
    var_names = ["S", "I", "R"]

    # time spacing between snapshots (0 to total_time) for the animation
    T = values_ts.shape[0]
    total_time = dt * steps
    dt_snapshot = total_time / max(T - 1, 1)

    # function for time-dependent edge weights - passed to animation
    edge_weight_fn = make_edge_weight_fn(hourly_weights, dt_snapshot)

    seoul = mtg.graphDisplay(G, seoul_map)
    seoul.interactive_graph()

    seoul.start_animation(
        values_ts=values_ts,
        dt_snapshot=dt_snapshot,
        var_names=var_names,
        var="I",
        edge_weight_fn=edge_weight_fn,
    )
    plt.show()

if __name__ == "__main__":
    main()

