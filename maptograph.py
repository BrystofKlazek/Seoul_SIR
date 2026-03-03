import geopandas as gpd

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
import matplotlib.patheffects as pe

import networkx as nx

from shapely import set_precision
from shapely.geometry import Point

import nametocode as ntc

import numpy as np

FRAME_SKIP = 1

def maptograph(gdf, MIN_BORDER=25, label="SIG_KOR_NM",
              pairs=None, mode="all", weight_fn=None, name_dict=None):

    if name_dict is None:
        name_dict = ntc.code_dict()

    gdf = gdf.copy()
    gdf["geometry"] = gdf.geometry.buffer(0).apply(lambda geom: set_precision(geom, 0.2))
    sindex = gdf.sindex

    G = nx.DiGraph()

    # add nodes once
    for i, row in gdf.iterrows():
        code = name_dict.lookup(row[label])
        G.add_node(code, name=row[label], pos=(row["x"], row["y"]), code=code, index=i)

    def get_w(u, v):
        if weight_fn is not None:
            return float(weight_fn(u, v))
        if pairs is not None:
            w = pairs.get((u, v), 0.0)
            return 0.0 if w is None else float(w)
        return 1.0

    if mode == "all":
        nodes = list(G.nodes())
        for a, u in enumerate(nodes):
            for v in nodes[a+1:]:
                w_uv = get_w(u, v)
                w_vu = get_w(v, u)
                if w_uv: G.add_edge(u, v, weight=w_uv)
                if w_vu: G.add_edge(v, u, weight=w_vu)

    elif mode == "neighbours":
        for i, row_i in gdf.iterrows():
            u = name_dict.lookup(row_i[label])
            for j in sindex.intersection(row_i.geometry.bounds):
                if j <= i: 
                    continue
                row_j = gdf.iloc[j]
                L = float(row_j.geometry.boundary.intersection(row_i.geometry.boundary).length)
                if L >= MIN_BORDER:
                    v = name_dict.lookup(row_j[label])
                    w_uv = get_w(u, v); w_vu = get_w(v, u)
                    if w_uv: G.add_edge(u, v, weight=w_uv)
                    if w_vu: G.add_edge(v, u, weight=w_vu)

    elif mode == "from_file":
        if pairs:
            for (u, v), w in pairs.items():
                if u in G and v in G and w:
                    G.add_edge(u, v, weight=float(w))
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    return G


def intercept_check(x_coord, y_coord, polygons):
        pt = Point(x_coord, y_coord)
        cand = polygons.sindex.query(pt)
        deep_check = polygons.iloc[cand]
        
        # The following is a bit of magic with pandas masks, because
        # .covers method from geopandas returns a bool dataframe indexed
        # the same way as the geopandas one and states wether the condition
        # is true or not.
        
        mask = deep_check.covers(pt)
        containing_polygon = mask[mask]
        if containing_polygon.empty == False:
            return containing_polygon.index[0]
        else:
            return None

#I SHOULD WRITE COMMENTS FOR THIS - TODO
class graphDisplay:
    def __init__(self, graph, gdt=None, max_flow=None, map_name_col="SIG_KOR_NM",
                 gif_name=None):
        self.map = gdt
        self.graph = graph
        self.fig = plt.figure()
        gs = GridSpec(nrows=1, ncols=2, figure=self.fig, width_ratios=[4.8, 1.6])
        
        self.ax = self.fig.add_subplot(gs[0, 0])      
        self.ax_flow = self.fig.add_subplot(gs[0, 1])  
        self.ax_flow.axis("off")
       
        self._title_text = self.ax.text(
            0.5, 1.01, "",
            transform=self.ax.transAxes,
            ha="center", va="bottom",
            fontsize=10,
            animated=True
        )

        self.fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.08, wspace=0.02)
        self.cid = None
        self.position = nx.get_node_attributes(self.graph, 'pos')
        self.node_idxs = nx.get_node_attributes(self.graph, 'index')
        self._node_list = np.array(list(self.graph.nodes()), dtype=object)
        xy = np.array([self.position[n] for n in self._node_list], dtype=float)
        self._node_x = xy[:, 0]
        self._node_y = xy[:, 1]

        self.map_name_col = map_name_col

        self.max_flow = None
        if max_flow is None or max_flow <= 0.0:
            self.max_flow = 1.0
        else:
            self.max_flow = max_flow

        self.save_name = gif_name
        self._sel_marker = None
        self._legend_neighs = []   

        self.node_patches = {}

        self.values_ts = None
        self.var_names = None

        self.selected_node = None
        self.ts_fig = None
        self.ax_ts = None
        self.lines_vars = []
        self.time_marker = None

        self.dt_snapshot = None
        self.t_vec = None
        self.cbar = None

        self._anim = None

        self.edge_weight_fn = None

        self._flow_text_code = []
        self._flow_text_in = []
        self._flow_text_out = []
        self._flow_max_neigh = 30
        self._init_flow_panel()

        # Colormaps for in/out magnitudes
        self.in_weight_norm = mcolors.Normalize(vmin=0.0, vmax=self.max_flow, clip=True)
        self.out_weight_norm = mcolors.Normalize(vmin=0.0, vmax=self.max_flow, clip=True)
        self.in_weight_cmap = cm.get_cmap("viridis")
        self.out_weight_cmap = cm.get_cmap("viridis")

        self.flow_update_enabled = False   # OFF by default (faster)
        self._current_frame = 0

        self.fig.canvas.draw_idle()

    def _on_key(self, event):
        if event.key == "a":
            self.flow_update_enabled = not self.flow_update_enabled
            state = "ON" if self.flow_update_enabled else "OFF"
            if self.selected_node is not None:
                self.ax_flow.set_title(f"Flows: {self.selected_node}  (dynamic {state})",
                                   fontsize=9, loc="left", pad=6)
            self.fig.canvas.draw_idle()

        if event.key == "u":
            if self.selected_node is not None and self.edge_weight_fn is not None:
                self._update_legend_for_frame(self._current_frame)
                self.fig.canvas.draw_idle()

    def _init_flow_panel(self):
        self.ax_flow.axis("off")

        self.ax_flow.text(
            0.00, 1.00, "Neighbour",
            transform=self.ax_flow.transAxes,
            ha="left", va="top",
            fontsize=8, fontweight="bold"
        )
        self.ax_flow.text(
            0.52, 1.00, "In",
            transform=self.ax_flow.transAxes,
            ha="right", va="top",
            fontsize=8, fontweight="bold"
        )
        self.ax_flow.text(
            0.98, 1.00, "Out",
            transform=self.ax_flow.transAxes,
            ha="right", va="top",
            fontsize=8, fontweight="bold"
        )

        self._flow_text_code = []
        self._flow_text_in = []
        self._flow_text_out = []
        self.ax_flow.set_facecolor((1, 1, 1, 0.85))  
    
        MAX_NEIGH = self._flow_max_neigh
        dy = 1.0 / (MAX_NEIGH + 1)

        for row in range(MAX_NEIGH):
            y = 1.0 - (row + 1) * dy

            row_bbox = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)

            txt_code = self.ax_flow.text(
                0.00, y, "",
                transform=self.ax_flow.transAxes,
                ha="left", va="center",
                fontsize=7,
                family="monospace",
                visible=False,
            )
            txt_in = self.ax_flow.text(
                0.52, y, "",
                transform=self.ax_flow.transAxes,
                ha="right", va="center",
                fontsize=7,
                family="monospace",
                visible=False,
            )
            txt_out = self.ax_flow.text(
                0.98, y, "",
                transform=self.ax_flow.transAxes,
                ha="right", va="center",
                fontsize=7,
                family="monospace",
                visible=False,
            )

            self._flow_text_code.append(txt_code)
            self._flow_text_in.append(txt_in)
            self._flow_text_out.append(txt_out)

    def _format_flow_label(self, v, w_in, w_out):
        code_str = str(v)
        return f"{code_str}   ← {w_in:8.3g}   → {w_out:8.3g}"

    # SELECTION + FLOW PANEL     
    def _update_selection(self, node):
        self._sel_marker.center = self.graph.nodes[node]["pos"]

        edges = []
        for nb in self.graph.neighbors(node):
            w_out = self.graph.get_edge_data(node, nb)["weight"]
            w_in  = self.graph.get_edge_data(nb, node)["weight"]
            edges.append((node, nb, w_in, w_out))

        edges.sort(key=lambda e: e[2] + e[3], reverse=True)
        self._legend_neighs = [v for (_, v, _, _) in edges]

        if not edges:
            for i in range(self._flow_max_neigh):
                self._flow_text_code[i].set_visible(False)
                self._flow_text_in[i].set_visible(False)
                self._flow_text_out[i].set_visible(False)

            self.ax_flow.set_title(f"Flows: {node}", fontsize=9, loc="left", pad=6)
            self.selected_node = node
            self.fig.canvas.draw_idle()
            return

        in_vals  = [w_in  for (_, _, w_in,  _) in edges]
        out_vals = [w_out for (_, _, _, w_out) in edges]
        max_in  = max(in_vals)  if max(in_vals)  > 0 else 1.0
        max_out = max(out_vals) if max(out_vals) > 0 else 1.0
        self.in_weight_norm.vmin = 0.0
        self.in_weight_norm.vmax = max_in
        self.out_weight_norm.vmin = 0.0
        self.out_weight_norm.vmax = max_out

        MAX_NEIGH = self._flow_max_neigh
        n_show = min(len(edges), MAX_NEIGH)

        for i in range(n_show):
            _, v, w_in, w_out = edges[i]

            self._flow_text_code[i].set_text(f"{v}")
            self._flow_text_in[i].set_text(f"in {w_in:.3g}")
            self._flow_text_out[i].set_text(f"out {w_out:.3g}")

            self._flow_text_in[i].set_color(self.in_weight_cmap(self.in_weight_norm(w_in)))
            self._flow_text_out[i].set_color(self.out_weight_cmap(self.out_weight_norm(w_out)))

            self._flow_text_code[i].set_visible(True)
            self._flow_text_in[i].set_visible(True)
            self._flow_text_out[i].set_visible(True)

        for i in range(n_show, MAX_NEIGH):
            self._flow_text_code[i].set_visible(False)
            self._flow_text_in[i].set_visible(False)
            self._flow_text_out[i].set_visible(False)

        self.ax_flow.set_title(f"Flows from {node}", fontsize=8, loc="left")

        self.selected_node = node
        if self.values_ts is not None and self.ax_ts is not None:
            self._update_timeseries_lines()

        if self._anim is not None:
            try:
                self._anim._init_draw()  # refresh blit background after changing visible texts
            except Exception:
                pass

        self.fig.canvas.draw_idle()

    def _on_click_graph(self, event):
        if event.inaxes is not self.ax or event.xdata is None or event.ydata is None:
            return
        dx = self._node_x - event.xdata
        dy = self._node_y - event.ydata
        idx = int(np.argmin(dx*dx + dy*dy))
        node = self._node_list[idx]
        self._update_selection(node)    

    def _flow_artists(self):
        return self._flow_text_code + self._flow_text_in + self._flow_text_out

    def interactive_graph(self):
        self.draw_graph()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click_graph)
        self.fig.canvas.mpl_connect('key_press_event', self._on_key)  

    # STATIC GRAPH DRAWING 
    def draw_graph(self):
        R = 1000
        label_dy = 1.2 * R   # shift labels above the node

        self._sel_marker = Circle((float('nan'), float('nan')),
                                  radius=R,
                                  fill=False,
                                  edgecolor='crimson',
                                  linewidth=2.5,
                                  zorder=80)
        self._sel_marker.set_animated(True)
        self.ax.add_patch(self._sel_marker)

        for node in self.graph:
            circ = Circle(self.position[node],
                          radius=R,
                          color="white",
                          zorder=50)
            circ.set_animated(True)
            self.ax.add_patch(circ)
            self.node_patches[node] = circ

            x, y = self.position[node]
            self.ax.text(
                x, y + label_dy, f"{node}",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")],
                ha='center', va='bottom', zorder=100, fontsize=7
            )

        self.ax.autoscale_view()
        self.ax.set_aspect('equal')
        self.ax.set_autoscale_on(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)

    # TIME–SERIES 
    def attach_timeseries(self, values_ts, dt_snapshot,
                          var_names=None):
        vals = np.asarray(values_ts)

        T, K, N = vals.shape

        if var_names is None:
            var_names = [f"var{k}" for k in range(K)]

        self.values_ts = vals       
        self.var_names = list(var_names)

        if dt_snapshot is None or dt_snapshot <= 0:
            self.dt_snapshot = None
            self.t_vec = np.arange(T, dtype=float)
        else:
            self.dt_snapshot = float(dt_snapshot)
            self.t_vec = np.arange(T, dtype=float) * self.dt_snapshot

        if self.ax_ts is None or self.ts_fig is None:
            self.ts_fig, self.ax_ts = plt.subplots()

        self.ax_ts.clear()
        self.ax_ts.set_xlabel("time")
        self.ax_ts.set_ylabel("value")

        if self.selected_node is None:
            self.selected_node = list(self.graph.nodes())[0]

        idx = self.node_idxs[self.selected_node]
        node_vals = vals[:, :, idx]

        self.lines_vars = []
        for k in range(K):
            line, = self.ax_ts.plot(self.t_vec, node_vals[:, k],
                                    label=self.var_names[k])
            self.lines_vars.append(line)

        self.time_marker = self.ax_ts.axvline(self.t_vec[0],
                                              color="k",
                                              linestyle="--",
                                              linewidth=1)

        self.ax_ts.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            frameon=True
        )
        self.ts_fig.subplots_adjust(right=0.78) 

        self.ax_ts.set_title(f"values at node {self.selected_node}")
        self.ax_ts.set_xlim(self.t_vec[0], self.t_vec[-1])
        self.ax_ts.relim()
        self.ax_ts.autoscale_view()
        self.ts_fig.canvas.draw_idle()

    def _update_timeseries_lines(self):
        if self.values_ts is None or not self.lines_vars:
            return
        if self.selected_node not in self.node_idxs:
            return

        idx = self.node_idxs[self.selected_node]
        vals = self.values_ts[:, :, idx]

        for k, line in enumerate(self.lines_vars):
            line.set_ydata(vals[:, k])

        self.ax_ts.set_title(f"values at node {self.selected_node}")
        if self.t_vec is not None:
            self.ax_ts.set_xlim(self.t_vec[0], self.t_vec[-1])
        self.ax_ts.relim()
        self.ax_ts.autoscale_view()
        self.ts_fig.canvas.draw_idle()

    def _update_legend_for_frame(self, frame):
        if (
            self.edge_weight_fn is None
            or self.selected_node is None
            or not self._legend_neighs
            or not self._flow_text_in
        ):
            return

        u = self.selected_node

        rows = []
        for v in self._legend_neighs:
            w_out = self.edge_weight_fn(frame, u, v)
            w_in  = self.edge_weight_fn(frame, v, u)

            w_out = 0.0 if w_out is None else float(w_out)
            w_in  = 0.0 if w_in  is None else float(w_in)

            rows.append((w_in + w_out, v, w_in, w_out))

        if not rows:
            return

        rows.sort(key=lambda x: x[0], reverse=True)

        self._legend_neighs = [v for _, v, _, _ in rows]

        MAX_NEIGH = self._flow_max_neigh
        n_show = min(len(rows), MAX_NEIGH)

        for i in range(n_show):
            _, v, w_in, w_out = rows[i]

            self._flow_text_code[i].set_text(f"{v}")
            self._flow_text_in[i].set_text(f"in {w_in:.3g}")
            self._flow_text_out[i].set_text(f"out {w_out:.3g}")

            self._flow_text_in[i].set_color(self.in_weight_cmap(self.in_weight_norm(w_in)))
            self._flow_text_out[i].set_color(self.out_weight_cmap(self.out_weight_norm(w_out)))

            self._flow_text_code[i].set_visible(True)
            self._flow_text_in[i].set_visible(True)
            self._flow_text_out[i].set_visible(True)

        for i in range(n_show, MAX_NEIGH):
            self._flow_text_code[i].set_visible(False)
            self._flow_text_in[i].set_visible(False)
            self._flow_text_out[i].set_visible(False)    

    def start_animation(self, var_names, values_ts, dt_snapshot, var=0,
                        interval=60, cmap_name="plasma", fps=10,
                        edge_weight_fn=None):

        vals = np.asarray(values_ts)
        if vals.ndim == 2:
            vals = vals[:, np.newaxis, :]
        if vals.ndim != 3:
            raise ValueError("values_ts must have shape (T, K, N) or (T, N)")

        T, K, N = vals.shape
        self.attach_timeseries(values_ts, dt_snapshot, var_names)

        self.edge_weight_fn = edge_weight_fn

        if isinstance(var, str):
            if self.var_names is None or var not in self.var_names:
                raise ValueError(f"Unknown variable name {var!r}")
            k = self.var_names.index(var)
        else:
            k = int(var)
            if not (0 <= k < K):
                raise ValueError(f"var index {k} out of range [0, {K})")

        field_ts = vals[:, k, :]

        vmin = float(field_ts.min())
        vmax = float(field_ts.max())
        if vmax <= vmin:
            vmax = vmin + 1e-12

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        cmap = cm.get_cmap(cmap_name)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        self.cbar = self.fig.colorbar(
            sm,
            ax=self.ax,
            orientation="horizontal",
            fraction=0.05,
            pad=0.10,
        )

        if self.var_names is not None:
            self.cbar.set_label(self.var_names[k])

        frame0 = field_ts[0]
        for node in self.graph:
            j = self.node_idxs[node]
            col = cmap(norm(frame0[j]))
            self.node_patches[node].set_facecolor(col)

        def update(frame):
            vals_frame = field_ts[frame]
            for node in self.graph:
                j = self.node_idxs[node]
                self.node_patches[node].set_facecolor(cmap(norm(vals_frame[j])))

            self._current_frame = frame

            if self.flow_update_enabled and self.edge_weight_fn is not None and self.selected_node is not None:
                self._update_legend_for_frame(frame)

            if self.time_marker is not None and self.t_vec is not None:
                t = self.t_vec[frame]
                self.time_marker.set_xdata([t, t])
                if self.ax_ts is not None and (frame % FRAME_SKIP == 0):
                    self.ax_ts.figure.canvas.draw_idle()

            if self.t_vec is not None:
                t = self.t_vec[frame]
                self._title_text.set_text(f"{self.var_names[k]}(t = {t:.3f}), frame {frame}")            
            else:
                self._title_text.set_text(f"{self.var_names[k]}, frame {frame}")

            artists = list(self.node_patches.values()) + [self._sel_marker, self._title_text]
            if self.flow_update_enabled:
                artists += self._flow_artists()
            return artists

        self._anim = FuncAnimation(self.fig, update,
                                   frames=T,
                                   interval=interval,
                                   blit = True)

        if self.save_name is not None:
            from matplotlib.animation import PillowWriter
            writer = PillowWriter(fps=fps)
            self._anim.save(self.save_name, writer=writer)

        return self._anim



