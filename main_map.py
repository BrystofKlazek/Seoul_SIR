import json
import time
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from shapely import set_precision
import maptograph_map as mtg
import numpy as np
import nametocode as ntc
import networkx as nx
from matplotlib.lines import Line2D
import SIR_diffusion_wrap as sir

# Simulation parameters
NORM            = 9600000 * 4  # normalization factor for flow weights, gotten roughly
                               #from current soul population multiplied by four for the current scope
TIME_STEP       = 0.025        # euler step in hours
NUM_DAYS        = 40
NUM_SNAPSHOTS   = 100          # number of frames saved from the simulation
BETA, DELTA, GAMMA = 0.08, 0.04, 0.02
SEED_IDX        = 11050        # district code where infection starts
DISPLAY_MODE    = "animation"  # "animation" or "slider"
FLOW_TOP_K      = 12           # number of neighbours shown in the flow window


# Builds a dict of {(u,v): array of hourly weights} from the raw JSON structure.
# The raw JSON is nested as {hour: {u: {v: weight}}}, this flattens it into
# per-edge arrays of length n_hours. Weights are normalized by NORM.
def build_hourly_weight_table(weights_by_hour_raw, norm=NORM):
    if not weights_by_hour_raw:
        return {}

    hour_indices = [int(h_str) for h_str in weights_by_hour_raw.keys()]
    n_hours = max(hour_indices) + 1

    hourly: dict[tuple[int, int], np.ndarray] = {}

    for h_str, by_u in weights_by_hour_raw.items():
        h = int(h_str)
        for u_str, neighs in by_u.items():
            u = int(u_str)
            for v_str, w in neighs.items():
                v = int(v_str)
                key = (u, v)

                arr = hourly.get(key)
                if arr is None:
                    arr = np.zeros(n_hours, dtype=np.float64)
                    hourly[key] = arr

                arr[h] += float(w) / norm

    return hourly


# Returns a function that interpolates the flow weight between u and v at time t (hours).
# The weights are periodic with period_hours (168 for weekly data) for now, can be later reworked.
# Linear interpolation is used between the two surrounding hourly values.
def make_edge_weight_fn(hourly_weights):
    if not hourly_weights:
        def edge_weight_fn(_, __, ___):
            return 0.0
        return edge_weight_fn

    sample_key   = next(iter(hourly_weights))
    period_hours = len(hourly_weights[sample_key])

    def edge_weight_fn(t, u, v):
        # Wrap time into one period, then interpolate between surrounding hours
        h    = t % float(period_hours)
        h0   = int(h)
        h1   = (h0 + 1) % period_hours
        frac = h - h0

        arr = hourly_weights.get((int(u), int(v)))
        if arr is None:
            return 0.0

        return (1.0 - frac) * arr[h0] + frac * arr[h1]

    return edge_weight_fn


# Right-hand side of the SEIR equations for one node - this can be easily rewritten
# to fit other model types, if the structure is kept.
# state and out are arrays of shape (K, N) - number of vals and number of nodes - passed in to the C solver.
def seir_rhs(state, out, beta=BETA, gamma=GAMMA, delta=DELTA):
    S, E, I, R     = state[0], state[1], state[2], state[3]
    dS, dE, dI, dR = out[0],   out[1],   out[2],   out[3]
    inf     = beta * S * I
    dS[...] = -inf
    dE[...] =  inf - delta * E
    dI[...] =  delta * E - gamma * I
    dR[...] =  gamma * I


# Right-hand side of the simpler SIR equations (no exposed compartment).
#Unused for now
def sir_rhs(state, out, beta=BETA, gamma=GAMMA):
    S, I, R        = state[0], state[1], state[2]
    dS, dI, dR     = out[0],   out[1],   out[2]
    inf     = beta * S * I
    dS[...] = -inf
    dI[...] =  inf - gamma * I
    dR[...] =  gamma * I


# Converts the per-edge hourly weight dict into a (n_hours, n, n) tensor
# that the C solver expects. Nodes are sorted by their graph index attribute
# to ensure consistent ordering between the graph and the simulation array.
def build_hourly_tensor_for_C(hourly_weights_dict, G, n_hours=168):
    node_index_attr = nx.get_node_attributes(G, "index")
    nodes           = sorted(G.nodes(), key=lambda u: node_index_attr[u])

    n           = len(nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}

    hourly_arr = np.zeros((n_hours, n, n), dtype=np.float64)

    for (u, v), arr in hourly_weights_dict.items():
        if u in node_to_idx and v in node_to_idx:
            hourly_arr[:, node_to_idx[u], node_to_idx[v]] = arr

    return hourly_arr, nodes, node_to_idx


# Plots the hourly in/out flow weights for the seed node and saves to a PNG.
#This is for inspection if all works well.
# Shows the top TOP_K neighbours by peak flow, one color per neighbour.
def plot_seed_weights(hourly_weights, seed_code, output_path="seed_weights.png"):
    if not hourly_weights:
        print("No hourly_weights, nothing to plot.")
        return

    any_arr = next(iter(hourly_weights.values()))
    n_hours = len(any_arr)
    t       = np.arange(n_hours, dtype=float)

    # Separate flows into outgoing from seed and incoming to seed
    outgoing, incoming = {}, {}
    for (u, v), arr in hourly_weights.items():
        arr = np.asarray(arr, dtype=float)
        if u == seed_code and v != seed_code:
            outgoing[v] = arr
        if v == seed_code and u != seed_code:
            incoming[u] = arr

    if not outgoing and not incoming:
        print(f"No flows involving seed node {seed_code}, nothing to plot.")
        return

    TOP_K        = 30
    all_labels   = sorted(set(outgoing.keys()) | set(incoming.keys()))
    # Score each neighbour by its peak flow in either direction
    label_scores = {
        lab: max(outgoing.get(lab, np.array([0.0])).max(),
                 incoming.get(lab, np.array([0.0])).max())
        for lab in all_labels
    }
    shown_labels = sorted(all_labels, key=lambda x: label_scores[x], reverse=True)[:TOP_K]

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    label_color = {lab: color_cycle[i % len(color_cycle)] for i, lab in enumerate(shown_labels)}

    fig, (ax_out, ax_in) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    for lab in shown_labels:
        if lab in outgoing:
            ax_out.plot(t, outgoing[lab], color=label_color[lab], linewidth=1)
    ax_out.set_title(f"Outgoing weights from {seed_code}")
    ax_out.set_ylabel("weight")

    for lab in shown_labels:
        if lab in incoming:
            ax_in.plot(t, incoming[lab], color=label_color[lab], linewidth=1)
    ax_in.set_title(f"Incoming weights into {seed_code}")
    ax_in.set_ylabel("weight")

    # Label x axis with day names if the data covers full days, else it just puts the hour index
    if n_hours % 24 == 0:
        days      = n_hours // 24
        xticks    = np.arange(0, n_hours + 1, 24)
        day_names = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
        ax_in.set_xticks(xticks)
        ax_in.set_xticklabels([day_names[d % 7] for d in range(days + 1)])
        ax_in.set_xlabel("hour of week (day at 0:00)")
    else:
        ax_in.set_xlabel("hour index")

    legend_handles = [
        Line2D([0], [0], color=label_color[lab], lw=1, label=str(lab))
        for lab in shown_labels
    ]
    fig.legend(handles=legend_handles, labels=[str(l) for l in shown_labels],
               loc="center right", borderaxespad=0.0, frameon=False, fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(right=0.90)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved seed weight plot to {output_path}")


def main():
    dt             = TIME_STEP
    steps          = int(24 * NUM_DAYS / dt)
    beta, gamma, delta = BETA, GAMMA, DELTA

    # Korean font setup - tries several common fonts, falls back to default
    plt.rcParams["font.family"]        = "NanumGothic"
    plt.rcParams["axes.unicode_minus"] = False
    for name in ["Noto Sans CJK KR", "Noto Sans KR", "NanumGothic",
                 "Malgun Gothic", "Apple SD Gothic Neo"]:
        if any(ft.name == name for ft in fm.fontManager.ttflist):
            plt.rcParams["font.family"] = name
            break

    # Load and prepare the Seoul district shapefile
    seoul_map = gpd.read_file("./shp/202101/SEOUL_SIG.shp").to_crs(5179)
    anchors   = seoul_map.representative_point()
    seoul_map["x"], seoul_map["y"] = anchors.x, anchors.y
    seoul_map["geometry"] = seoul_map.geometry.buffer(0).apply(
        lambda geom: set_precision(geom, 0.2)
    )

    # Load name-to-code lookup and mobility weights
    name_dict = ntc.code_dict(code_df=pd.read_csv("code_lookup.csv"))
    with open("weights_03.json", encoding="utf-8") as f:
        weights_by_hour_raw = json.load(f)
    hourly_weights = build_hourly_weight_table(weights_by_hour_raw, norm=NORM)

    # Build graph - only edges that have nonzero flow at hour 0 are added.
    # The actual weight values stored on edges are unused after construction,
    # edge_weight_fn handles the real time-varying weights during display.
    edge_weights0 = {(u, v): arr[0] for (u, v), arr in hourly_weights.items() if arr[0] > 0.0}
    G = mtg.graph_construction(gdf=seoul_map, mode="from_file", pairs=edge_weights0)

    # Build tensor for C solver - shape (n_hours, n_nodes, n_nodes)
    hourly_tensor, nodes_order, node_to_idx = build_hourly_tensor_for_C(
        hourly_weights, G, n_hours=168
    )

    # Initial conditions - everyone susceptible except seed node
    n  = len(nodes_order)
    S0 = np.ones(n,  np.float64)
    E0 = np.zeros(n, np.float64)
    I0 = np.zeros(n, np.float64)
    R0 = np.zeros(n, np.float64)

    seed_idx      = node_to_idx[SEED_IDX]
    I0[seed_idx]  = 0.01
    S0[seed_idx] -= 0.01

    x0_graph = np.stack([S0, E0, I0, R0], axis=0)

    # Run simulation
    print("Solving equation on Seoul graph...")
    t0        = time.time()
    out_graph = sir.euler_solve_graph_rd(
        hourly_tensor, x0_graph, 0.0, dt, steps, NUM_SNAPSHOTS,
        rhs_fn=seir_rhs, params={"beta": beta, "gamma": gamma, "delta": delta}
    )
    print(f"Graph solve done in {time.time() - t0:.2f}s, snaps = {out_graph.shape[0]}")

    values_ts = out_graph
    var_names = ["S", "E", "I", "R"]

    # dt_snapshot is hours per simulation snapshot, used to build the time axis
    total_time  = dt * steps
    dt_snapshot = total_time / max(NUM_SNAPSHOTS - 1, 1)

    # edge_weight_fn(t, u, v) returns interpolated flow weight at time t in hours
    edge_weight_fn = make_edge_weight_fn(hourly_weights)

    plot_seed_weights(hourly_weights, SEED_IDX, "seed_weights.png")

    # Build and display interactive map
    seoul = mtg.graphDisplay(G, seoul_map)
    seoul.interactive_graph()

    if DISPLAY_MODE == "animation":
        seoul.start_animation(
            values_ts=values_ts, dt_snapshot=dt_snapshot,
            var_names=var_names, var="I", edge_weight_fn=edge_weight_fn,
        )
    elif DISPLAY_MODE == "slider":
        seoul.show_frame_slider(
            values_ts=values_ts, dt_snapshot=dt_snapshot,
            var_names=var_names, var="I", edge_weight_fn=edge_weight_fn,
            slider_label="frame",
        )
    else:
        raise ValueError(f"Unknown DISPLAY_MODE {DISPLAY_MODE!r}")

    seoul.attach_flow_timeseries(top_k=FLOW_TOP_K)
    plt.show()


if __name__ == "__main__":
    main()
