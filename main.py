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
from itertools import product
import networkx as nx

this_dir = Path(__file__).resolve().parent
sys.path.append(str(this_dir.parent / "differential_solver_C"))

import SIR_diffusion_wrap as sir  


NORM = 100000

def interpolated_weight(weights_by_hour, u, v, when):
    """
    Linearly interpolate weight u->v at time `when` (in hours),
    using the hourly JSON data. Works whether u, v are str or int.
    Missing entries are treated as 0.
    """
    h_float = float(when) % 24.0
    h0 = int(h_float)
    h1 = (h0 + 1) % 24
    frac = h_float - h0

    h0_str = str(h0)
    h1_str = str(h1)

    u_key = str(u)
    v_key = str(v)

    # hour 0 weight
    hour0 = weights_by_hour.get(h0_str, {})
    u_dict0 = hour0.get(u_key, {})
    w0 = float(u_dict0.get(v_key, 0.0))

    # hour 1 weight
    hour1 = weights_by_hour.get(h1_str, {})
    u_dict1 = hour1.get(u_key, {})
    w1 = float(u_dict1.get(v_key, 0.0))

    w = (1.0 - frac) * w0 + frac * w1
    return w/NORM

def build_laplacian_timeseries_for_C(G, weights_by_hour,
                                     n_steps, dt, t0=0.0):
    """
    Build per-step CSR Laplacians for the graph G, using time-dependent
    weights from weights_by_hour and linear interpolation in time.

    Returns
    -------
    indptr_list  : list of length n_steps, each is np.int32, shape (n+1,)
    indices_list : list of length n_steps, each is np.int32, shape (nnz,)
    data_list    : list of length n_steps, each is np.float64, shape (nnz,)
    nodes        : list of node IDs, order used in the CSR
    """
    node_index = nx.get_node_attributes(G, "index")
    nodes = sorted(G.nodes(), key=lambda n: node_index[n])
    n = len(nodes)

    # --- neighbour pattern: use both successors and predecessors ---
    neighbours_idx = []
    for node in nodes:
        outs = set(G.successors(node))
        ins  = set(G.predecessors(node))
        neighs = sorted(node_index[nb] for nb in outs | ins)
        neighbours_idx.append(neighs)

    # if there are really no edges, nnz will be exactly n (diagonal only)
    indptr = [0]
    indices = []
    for i, neighs in enumerate(neighbours_idx):
        indices.extend(neighs)  # off-diagonals
        indices.append(i)       # diagonal
        indptr.append(len(indices))

    indptr_arr = np.asarray(indptr, dtype=np.int32)
    indices_arr = np.asarray(indices, dtype=np.int32)
    nnz = indices_arr.size

    print(f"Graph has {n} nodes, CSR nnz per step = {nnz}")

    # pattern is the same for every step → reuse same arrays
    indptr_list  = [indptr_arr]  * n_steps
    indices_list = [indices_arr] * n_steps

    data_list = []

    time = float(t0)
    for step in range(n_steps):
        data_step = np.empty(nnz, dtype=float)
        pos = 0
        for i, node in enumerate(nodes):
            deg = 0.0
            for j in neighbours_idx[i]:
                nb = nodes[j]

                # symmetric weight between node and nb, from OD file
                w = 0.0
                if G.has_edge(node, nb):
                    w += interpolated_weight(weights_by_hour, node, nb, time)
                if G.has_edge(nb, node):
                    w += interpolated_weight(weights_by_hour, nb, node, time)

                data_step[pos] = -w
                pos += 1
                deg += w

            data_step[pos] = deg
            pos += 1

        data_list.append(data_step.astype(np.float64, copy=False))
        time = (time + dt) % 24.0

    return indptr_list, indices_list, data_list, nodes

def sir_rhs(state, out, beta=0.4, gamma=0.1):
    S, I, R = state[0], state[1], state[2]
    dS, dI, dR = out[0], out[1], out[2]
    inf = beta * S * I
    dS[...] = -inf
    dI[...] =  inf - gamma * I
    dR[...] =  gamma * I


def main():
    dt = 0.01
    steps = 3000
    snapshots = 500
    beta, gamma = 0.8, 0.2

    DEFAULT_W = 0.02
    GRID = 0.2

    plt.rcParams["font.family"] = 'NanumGothic'
    plt.rcParams['axes.unicode_minus'] = False
    for name in ["Noto Sans CJK KR", "Noto Sans KR", "NanumGothic",
                 "Malgun Gothic", "Apple SD Gothic Neo"]:
        if any(ft.name == name for ft in fm.fontManager.ttflist):
            plt.rcParams["font.family"] = name
            break

    shapefile_path = './shp/202101/SEOUL_SIG.shp'
    seoul_map = gpd.read_file(shapefile_path).to_crs(5179)

    anchors = seoul_map.representative_point()
    seoul_map["x"], seoul_map["y"] = anchors.x, anchors.y

    seoul_map["geometry"] = seoul_map.geometry.buffer(0).apply(
        lambda geom: set_precision(geom, GRID)
    )

    name_code_df = pd.read_csv("code_lookup.csv")
    name_dict = ntc.code_dict(code_df=name_code_df)

    codes = name_code_df["sgg"].to_list()
    with open("weights_03.json", encoding="utf-8") as f:
        weights_by_hour = json.load(f)

    edge_weights0 = {}
    for u, neighs in weights_by_hour["0"].items():
        u_code = int(u)                 # convert to int
        for v, w in neighs.items():
            v_code = int(v)             # convert to int
            edge_weights0[(u_code, v_code)] = w / NORM

    G = mtg.maptograph(seoul_map, mode="from_file", pairs=edge_weights0)

    # build per-step CSR Laplacians
    indptr_list, indices_list, data_list, nodes = build_laplacian_timeseries_for_C(
        G, weights_by_hour, steps, dt
    )
    n = len(nodes)
    print(f"Graph has {n} nodes, CSR nnz at step 0 = {data_list[0].size}")

    # initial SIR state
    S0 = np.ones(n, np.float64)
    I0 = np.zeros(n, np.float64)
    R0 = np.zeros(n, np.float64)

    # seed: highest-degree district
    seed_idx = max(range(n), key=lambda k: G.degree[nodes[k]])
    I0[seed_idx] = 0.1
    S0[seed_idx] -= I0[seed_idx]

    x0_graph = np.stack([S0, I0, R0], axis=0)

    print("Solving SIR on Seoul graph (CSR, Euler, reaction–diffusion SIR)…")
    t0 = time.time()
    out_graph = sir.euler_solve_graph_csr_rd(
        indptr_list, indices_list, data_list, x0_graph,
        dt, steps, snapshots,
        rhs_fn=sir_rhs, params={"beta": beta, "gamma": gamma}
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

    # build an array (T, E) of weights for each frame, in the same edge order:
    edges = list(G.edges())
    E = len(edges)
    T = out_graph.shape[0]  # number of snapshots / frames

    total_time = dt * steps
    dt_snapshot = total_time / max(T - 1, 1)   # full 0..total_time range
    def edge_weight_fn(frame, u, v):
        # physical time of this frame, same as your timeseries plots
        t = frame * dt_snapshot
        return interpolated_weight(weights_by_hour, u, v, t)

    seoul = mtg.graphDisplay(G, seoul_map)
    seoul.interactive_graph()

    seoul.start_animation(
        values_ts=values_ts,
        dt_snapshot=dt_snapshot,
        var_names=var_names,
        var="I",
        edge_weight_fn=edge_weight_fn,   # <-- this is the only new thing
    )
    plt.show()

if __name__ == "__main__":
    main()

