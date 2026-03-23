"""
predict_transport.py
====================
Exact linearised prediction of population distribution u_i(t) on a directed
graph under the transport equation

    du/dt = -L(t) u

given only the time-varying edge weight data (no simulation required).

Theory
------
The exact solution for a single node, treating neighbour values as fixed at
their initial condition u_j ≈ 1 (linearisation), is:

    u_i(t) = 1 + ∫_0^t  e^{-(W_i(t) - W_i(s))}  δ_i(s)  ds

where:
    W_i(t) = ∫_0^t w_i^out(τ) dτ     cumulative outflow
    δ_i(t) = w_i^in(t) - w_i^out(t)  instantaneous flow imbalance

This is computed efficiently as:

    u_i(t) = 1 + e^{-W_i(t)}  ∫_0^t  e^{W_i(s)}  δ_i(s)  ds
           = 1 + e^{-W_i(t)}  · cumsum[ e^{W_i} · δ_i · dt ]

which requires only a single pass over the time axis per node.

Usage
-----
    python predict_transport.py

Outputs
-------
    prediction_overlay.png   — predicted vs simulated u_i(t) for top nodes
    prediction_error.png     — absolute error over time
    prediction_scatter.png   — predicted vs simulated peak heights (scatter)
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ── configuration ─────────────────────────────────────────────────────────────
WEIGHTS_PATH = Path("weights_03.json")
NORM         = 9_600_000 * 5    # normalisation divisor from main_fixed.py
DT           = 0.5              # time step in hours
N_PERIODS    = 6                # number of weekly periods to simulate/predict
TOP_N_NODES  = 4                # nodes to show in overlay plot
N_COLS       = 2
# ─────────────────────────────────────────────────────────────────────────────


def load_weights(path, norm):
    """
    Load weights_03.json and return:
        hourly  : dict {(u,v): np.ndarray shape (n_hours,)}  — normalised flows
        n_hours : int
    """
    with open(path, encoding="utf-8") as f:
        raw = json.load(f)

    n_hours = max(int(h) for h in raw) + 1
    hourly  = {}

    for h_str, by_u in raw.items():
        h = int(h_str)
        for u_str, neighs in by_u.items():
            u = int(u_str)
            for v_str, w in neighs.items():
                v   = int(v_str)
                key = (u, v)
                if key not in hourly:
                    hourly[key] = np.zeros(n_hours)
                hourly[key][h] += float(w) / norm

    return hourly, n_hours


def build_node_flows(hourly, n_hours):
    """
    Aggregate edge flows into per-node in/out arrays.

    Returns
    -------
    nodes   : sorted list of node ids
    idx     : dict {node_id: integer index}
    w_in    : (n_nodes, n_hours)  — total hourly inflow rate
    w_out   : (n_nodes, n_hours)  — total hourly outflow rate
    """
    nodes = sorted({n for pair in hourly for n in pair})
    idx   = {n: i for i, n in enumerate(nodes)}
    n     = len(nodes)

    w_in  = np.zeros((n, n_hours))
    w_out = np.zeros((n, n_hours))

    for (u, v), arr in hourly.items():
        w_out[idx[u]] += arr
        w_in[idx[v]]  += arr

    return nodes, idx, w_in, w_out


def build_weight_tensor(hourly, nodes, idx, n_hours):
    """W[h, i, j] = flow from node j to node i at hour h."""
    n = len(nodes)
    W = np.zeros((n_hours, n, n))
    for (u, v), arr in hourly.items():
        W[:, idx[v], idx[u]] = arr   # destination i = v, source j = u
    return W


def simulate(W, w_out, n_hours, dt, n_steps):
    """
    Euler integration of  du/dt = -L(t) u.

    (Lu)_i = w_out_i(t) u_i - Σ_j W[t, i, j] u_j

    Returns trajectory of shape (n_steps+1, n_nodes).
    """
    n    = W.shape[1]
    u    = np.ones(n)
    traj = np.empty((n_steps + 1, n))
    traj[0] = u

    for step in range(n_steps):
        th = (step * dt) % n_hours
        h0 = int(th);  h1 = (h0 + 1) % n_hours;  a = th - h0

        W_t    = (1 - a) * W[h0]    + a * W[h1]           # (n, n)
        wout_t = (1 - a) * w_out[:, h0] + a * w_out[:, h1]  # (n,)

        u = u + dt * (W_t @ u - wout_t * u)
        traj[step + 1] = u

    return traj


def predict(w_in, w_out, n_hours, dt, n_steps):
    """
    Exact linearised prediction:

        u_i(t) ≈ 1 + e^{-W_i(t)}  cumsum[ e^{W_i(s)} δ_i(s) ds ]

    All operations are vectorised over nodes.

    Returns
    -------
    u_pred  : (n_steps+1, n_nodes)
    """
    n_nodes = w_in.shape[0]

    # Build time-interpolated arrays at each step: shape (n_steps+1, n_nodes)
    step_arr = np.arange(n_steps + 1)
    th_arr   = (step_arr * dt) % n_hours
    h0_arr   = th_arr.astype(int)
    h1_arr   = (h0_arr + 1) % n_hours
    a_arr    = th_arr - h0_arr

    # delta_t[step, node] = w_in_i(t) - w_out_i(t)
    delta_t = ((1 - a_arr[:, None]) * (w_in  - w_out)[:, h0_arr].T +
                    a_arr[:, None]  * (w_in  - w_out)[:, h1_arr].T)

    # w_out_t[step, node]
    wout_t  = ((1 - a_arr[:, None]) * w_out[:, h0_arr].T +
                    a_arr[:, None]  * w_out[:, h1_arr].T)

    # W_i(t) = cumulative outflow (trapezoid via cumsum * dt)
    W_cum = np.cumsum(wout_t * dt, axis=0)           # (n_steps+1, n_nodes)

    # Integral: e^{-W(t)} * cumsum[ e^{W(s)} * delta(s) * ds ]
    eW_delta_dt      = np.exp(W_cum) * delta_t * dt  # integrand
    cumulative_int   = np.cumsum(eW_delta_dt, axis=0)
    u_pred           = 1.0 + np.exp(-W_cum) * cumulative_int

    return u_pred


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_overlay(t_arr, traj, u_pred, nodes, top_indices, out_path):
    ncols = N_COLS
    nrows = (len(top_indices) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 3.5 * nrows), sharex=True)

    for ax, idx in zip(axes.flat, top_indices):
        ax.plot(t_arr, traj[:, idx],   color="steelblue", lw=1.4,
                label=r"$N_i$")
        ax.plot(t_arr, u_pred[:, idx], color="tomato",    lw=1.0, ls="--",
                label=r"$\widehat{N}_i$")
        ax.set_title(f"node {nodes[idx]}", fontsize=12)
        ax.set_ylabel("population", fontsize=12)
        ax.tick_params(labelsize=10)

    for ax in list(axes.flat)[len(top_indices):]:
        ax.set_visible(False)

    axes.flat[0].legend(fontsize=12)
    fig.supxlabel("time (hours)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_error(t_arr, traj, u_pred, nodes, top_indices, out_path):
    ncols = N_COLS
    nrows = (len(top_indices) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 3 * nrows), sharex=True)

    for ax, idx in zip(axes.flat, top_indices):
        err = np.abs(traj[:, idx] - u_pred[:, idx])
        ax.plot(t_arr, err, color="tomato", lw=0.9)
        ax.set_title(f"node {nodes[idx]}", fontsize=12)
        ax.set_ylabel("|error|", fontsize=11)
        ax.tick_params(labelsize=10)

    for ax in list(axes.flat)[len(top_indices):]:
        ax.set_visible(False)

    fig.supxlabel("time (hours)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_scatter(traj, u_pred, nodes, out_path):
    """Scatter of predicted vs simulated values at every timestep, every node."""
    fig, ax = plt.subplots(figsize=(6, 6))

    # Subsample for clarity (every 4th step)
    sim_vals  = traj[::4].ravel()
    pred_vals = u_pred[::4].ravel()

    ax.scatter(sim_vals, pred_vals, s=3, alpha=0.15, color="steelblue",
               rasterized=True)

    lo = min(sim_vals.min(), pred_vals.min()) - 0.001
    hi = max(sim_vals.max(), pred_vals.max()) + 0.001
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, label=r"ideal $\widehat{N}_i(t) = N_i(t)$")

    r = np.corrcoef(sim_vals, pred_vals)[0, 1]
    rmse = np.sqrt(np.mean((sim_vals - pred_vals) ** 2))
    ax.set_xlabel(r"$\widehat{N}_i(t)$", fontsize=12)
    ax.set_ylabel(r"$N_i(t)$", fontsize=12)
    ax.legend(fontsize=11)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")
    print(f"RMSE = {rmse}, r = {r}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading weights ...")
    hourly, n_hours = load_weights(WEIGHTS_PATH, NORM)
    print(f"  {n_hours} hours,  {len(hourly)} directed edges")

    nodes, idx, w_in, w_out = build_node_flows(hourly, n_hours)
    n = len(nodes)
    print(f"  {n} nodes")

    n_steps = int(N_PERIODS * n_hours / DT)
    t_arr   = np.arange(n_steps + 1) * DT
    print(f"\nSimulating {n_steps} steps  (dt={DT}h,  "
          f"{N_PERIODS} periods = {N_PERIODS * n_hours}h) ...")

    W_tensor = build_weight_tensor(hourly, nodes, idx, n_hours)
    traj     = simulate(W_tensor, w_out, n_hours, DT, n_steps)
    print(f"  trajectory shape: {traj.shape}")

    print("Computing exact linearised prediction ...")
    u_pred = predict(w_in, w_out, n_hours, DT, n_steps)
    print(f"  prediction shape: {u_pred.shape}")

    # ── statistics ────────────────────────────────────────────────────────────
    rmse     = np.sqrt(np.mean((traj - u_pred) ** 2))
    rel_rmse = rmse / traj.std()
    r        = np.corrcoef(traj.ravel(), u_pred.ravel())[0, 1]

    print(f"\nGlobal RMSE              : {rmse:.6f}")
    print(f"Relative RMSE (÷ std)    : {rel_rmse:.4f}")
    print(f"Pearson r                : {r:.4f}")

    print("\nPer-node RMSE:")
    for i, nd in enumerate(nodes):
        node_rmse = np.sqrt(np.mean((traj[:, i] - u_pred[:, i]) ** 2))
        print(f"  node {nd}: {node_rmse:.6f}")

    # ── plots ─────────────────────────────────────────────────────────────────
    variability = traj.std(axis=0)
    top_idx     = np.argsort(variability)[::-1][:TOP_N_NODES].tolist()
    print(f"\nTop {TOP_N_NODES} nodes by variability: {[nodes[i] for i in top_idx]}")

    out_dir = Path(".")
    plot_overlay(t_arr, traj, u_pred, nodes, top_idx,
                 out_dir / "prediction_overlay.png")
    plot_error  (t_arr, traj, u_pred, nodes, top_idx,
                 out_dir / "prediction_error.png")
    plot_scatter(traj, u_pred, nodes,
                 out_dir / "prediction_scatter.png")


if __name__ == "__main__":
    main()
