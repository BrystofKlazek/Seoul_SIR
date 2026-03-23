"""
make_table_errors.py
====================
Generates table_errors.tex — per-district RMSE by weekly block
and OLS error trend (slope, r, p-value).

Requires:  weights_03.json  in the same directory.
Output:    table_errors.tex
LaTeX preamble needs:  \\usepackage{booktabs}
"""

import json
import numpy as np
from scipy import stats
from pathlib import Path

WEIGHTS_PATH = Path("weights_03.json")
NORM         = 9_600_000 * 5
DT           = 0.5
N_PERIODS    = 6
n_hours      = 168
OUTPUT_PATH  = Path("table_errors.tex")

# ── load weights ──────────────────────────────────────────────────────────────
with open(WEIGHTS_PATH, encoding="utf-8") as f:
    raw = json.load(f)

hourly = {}
for h_str, by_u in raw.items():
    h = int(h_str)
    for u_str, neighs in by_u.items():
        u = int(u_str)
        for v_str, w in neighs.items():
            v   = int(v_str)
            key = (u, v)
            if key not in hourly:
                hourly[key] = np.zeros(n_hours)
            hourly[key][h] += float(w) / NORM

nodes = sorted({nn for pair in hourly for nn in pair})
idx   = {nn: i for i, nn in enumerate(nodes)}
n     = len(nodes)
w_in  = np.zeros((n, n_hours))
w_out = np.zeros((n, n_hours))
for (u, v), arr in hourly.items():
    w_out[idx[u]] += arr
    w_in[idx[v]]  += arr

W_tensor = np.zeros((n_hours, n, n))
for (u, v), arr in hourly.items():
    W_tensor[:, idx[v], idx[u]] = arr

# ── simulate ──────────────────────────────────────────────────────────────────
n_steps = int(N_PERIODS * n_hours / DT)
t_arr   = np.arange(n_steps + 1) * DT

u_sim = np.ones(n)
traj  = np.empty((n_steps + 1, n))
traj[0] = u_sim.copy()
for step in range(n_steps):
    th     = (step * DT) % n_hours
    h0     = int(th); h1 = (h0 + 1) % n_hours; a = th - h0
    W_t    = (1 - a) * W_tensor[h0] + a * W_tensor[h1]
    wout_t = (1 - a) * w_out[:, h0] + a * w_out[:, h1]
    u_sim  = u_sim + DT * (W_t @ u_sim - wout_t * u_sim)
    traj[step + 1] = u_sim.copy()

# ── predict ───────────────────────────────────────────────────────────────────
step_arr   = np.arange(n_steps + 1)
th_arr     = (step_arr * DT) % n_hours
h0_arr     = th_arr.astype(int)
h1_arr     = (h0_arr + 1) % n_hours
a_arr      = th_arr - h0_arr
delta      = w_in - w_out
delta_t    = ((1 - a_arr[:, None]) * delta[:, h0_arr].T +
                  a_arr[:, None]  * delta[:, h1_arr].T)
wout_t_all = ((1 - a_arr[:, None]) * w_out[:, h0_arr].T +
                  a_arr[:, None]  * w_out[:, h1_arr].T)
W_cum      = np.cumsum(wout_t_all * DT, axis=0)
eW_delta   = np.exp(W_cum) * delta_t * DT
u_pred     = 1.0 + np.exp(-W_cum) * np.cumsum(eW_delta, axis=0)

# ── error statistics ──────────────────────────────────────────────────────────
err     = np.abs(traj - u_pred)
block   = int(n_hours / DT)
n_weeks = N_PERIODS

week_rmse = np.zeros((n, n_weeks))
for b in range(n_weeks):
    s = b * block
    e = (b + 1) * block
    week_rmse[:, b] = np.sqrt(np.mean(err[s:e] ** 2, axis=0))

global_rmse = np.sqrt(np.mean(err ** 2, axis=0))

# ── OLS trend on weekly block means (6 points, honest degrees of freedom) ────
week_indices = np.arange(1, n_weeks + 1, dtype=float)  # 1..6
slopes = np.zeros(n)
pvals  = np.zeros(n)
rvals  = np.zeros(n)
for i in range(n):
    sl, _, rv, pv, _ = stats.linregress(week_indices, week_rmse[i])
    slopes[i] = sl
    pvals[i]  = pv
    rvals[i]  = rv

# ── build table ───────────────────────────────────────────────────────────────
lines = []
lines.append(r'\begin{table}[ht]')
lines.append(r'\centering')
lines.append(r'\small')
lines.append(r'\begin{tabular}{@{}lrrrrrrrrrrr@{}}')
lines.append(r'\toprule')
lines.append(
    r'District'
    r' & \multicolumn{7}{c}{RMSE\,($\times 10^{-4}$)}'
    r' & Slope ($\times10^{-4}$\,week$^{-1}$)'
    r' & $r$'
    r' & $p$ \\'
)
lines.append(r'\midrule')
lines.append(r' & W1 & W2 & W3 & W4 & W5 & W6 & Overall & & & & \\')
lines.append(r'\midrule')

for i, nd in enumerate(nodes):
    sig = r'$^*$' if pvals[i] < 0.05 else ''
    row = f'{nd}'
    for b in range(n_weeks):
        row += f' & {week_rmse[i, b]*1e4:.1f}'
    row += f' & {global_rmse[i]*1e4:.1f}'
    row += f' & {slopes[i]*1e4:.2f}'
    row += f' & {rvals[i]:.3f}'
    row += f' & {pvals[i]:.3f}{" " + sig if sig else ""}'
    row += r' \\'
    lines.append(row)

lines.append(r'\midrule')
mean_row = r'\textit{Mean}'
for b in range(n_weeks):
    mean_row += f' & {week_rmse[:, b].mean()*1e4:.2f}'
mean_row += f' & {global_rmse.mean()*1e4:.2f}'
mean_row += r' & & & \\'
lines.append(mean_row)

lines.append(r'\bottomrule')
lines.append(r'\end{tabular}')
lines.append(
    r'\caption{Per-district RMSE of the exact linearised prediction '
    r'against the full simulation broken down by weekly block '
    r'(columns W1--W6, $\times 10^{-4}$), '
    r'overall six-week RMSE, '
    r'and OLS linear trend fitted to the six weekly block means. '
    r'Slope in units of $\times 10^{-4}$ RMSE per week ($n=6$); '
    r'asterisks mark significance at $p < 0.05$.}'
)
lines.append(r'\label{tab:error_summary}')
lines.append(r'\end{table}')

with open(OUTPUT_PATH, 'w') as f:
    f.write('\n'.join(lines))

print(f"Written: {OUTPUT_PATH}")
