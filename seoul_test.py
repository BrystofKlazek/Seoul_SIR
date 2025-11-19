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

import SIR_diffusion_wrap as sir  



def graph_laplacian_sparse(graph, weight_key="weight"):
    node_index = nx.get_node_attributes(graph, 'index')
    nodes = sorted(graph.nodes(), key=lambda n: node_index[n])  # indices 0..n-1
    n = len(nodes)
    
    # sanity check
    idx_vals = set(node_index.values())
    if idx_vals != set(range(n)):
        raise ValueError(
            "Node 'index' attributes must be exactly 0..n-1 and contiguous "
            f"(got {sorted(idx_vals)} for n={n})"
        )

    indptr = [0]
    indices = []
    data = []

    for node in nodes:
        i = node_index[node]
        deg = 0.0
        for neighbour in graph.neighbors(node):
            j = node_index[neighbour]
            w = graph.get_edge_data(node, neighbour)[weight_key]
            indices.append(j)
            data.append(-w)
            deg += w

        indices.append(i)
        data.append(deg)
        indptr.append(len(indices))

    return (np.asarray(indptr, dtype=np.int32),
            np.asarray(indices, dtype=np.int32),
            np.asarray(data, dtype=np.float64),
            nodes)

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
    name_dict = ntc.code_dict(code_df=name_code_df)

    codes = name_code_df["sgg"].to_list()
    edge_weights = {pair: 0.1 
                    for i, pair in enumerate(product(codes, codes))}

    G = mtg.maptograph(seoul_map, mode="neighbours", pairs=edge_weights)

    print("Building CSR Laplacian for Seoul graph…")
    indptr, indices, data, nodes = graph_laplacian_sparse(G)
    n = len(nodes)
    print(f"Graph has {n} nodes, CSR nnz = {len(data)}")

    #------ CONSTRUCTING SIR ON THE GRAPH ------#
    # S0: everyone susceptible
    S0 = np.ones(n, np.float64)
    # I0: single initially infected node 
    #(For test pick the highest-degree district)
    I0 = np.zeros(n, np.float64)
    seed_idx = max(range(n), key=lambda k: G.degree[nodes[k]])
    I0[seed_idx] = 0.1
    # R0: nobody recovered yet
    R0 = np.zeros(n, np.float64)

    x0_graph = np.stack([S0, I0, R0], axis=0) 

    #----- CALLING THE SOLVER ON THE GRAPH -----#
    print("Solving SIR on Seoul graph (CSR, Euler, reaction–diffusion SIR)…")
    t0 = time.time()
    out_graph = sir.euler_solve_graph_csr_rd(
        indptr, indices, data, x0_graph,
        dt, steps, snapshots,
        rhs_fn=sir_rhs, params={"beta": beta, "gamma": gamma}
    )
    print(f"Graph solve done in {time.time() - t0:.2f}s, "
          f"snaps = {out_graph.shape[0]}")   

    out = out_graph

    print("first snapshot:", out[0, :, 0])
    print("middle snapshot:", out[out.shape[0]//2, :, 0])
    print("last snapshot:", out[-1, :, 0])
    print("I min/max over all:", out[:, 1, :].min(), out[:, 1, :].max())

    values_ts = out_graph                   
    var_names = ["S", "I", "R"]

    T = values_ts.shape[0]
    total_time = dt * steps
    dt_snapshot = total_time / max(T - 1, 1)   # full 0..total_time range
    seoul = mtg.graphDisplay(G, seoul_map)
    seoul.interactive_graph()   # your existing click + flow highlighting

    seoul.attach_timeseries(
        values_ts,
        dt_snapshot=dt_snapshot,
        var_names=var_names,
    )

    seoul.start_animation(var="I")

    plt.show()
if __name__ == "__main__":
    main()

