import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
import matplotlib.patheffects as pe

from matplotlib.collections import LineCollection
import random

import networkx as nx

from shapely import set_precision
from shapely.geometry import Point

import nametocode as ntc

import numpy as np


#Funciton turning GEOPANDAS map with labels into a graph and turning the label names
#into a code
def maptograph(map, MIN_BORDER = 25, label="SIG_KOR_NM", pairs={}, mode="all",
               name_dict = ntc.code_dict()): 
    
    #GRID to snap to (so detect edges works)
    GRID = 0.2   
    map["geometry"] = map.geometry.buffer(0).apply(
            lambda geom: set_precision(geom, GRID)
            )

    #Search spatial index so that searching is faster 
    #(builds a tree for searching in the tree)
    sindex = map.sindex
    G = nx.DiGraph()
   
    #If we want to build graph from neighbours
    if mode=="neighbours":
        for i, row_i in map.iterrows():
            i_code = name_dict.lookup(row_i[label])
            #Adding nodes to the graph. Each node is indexed by the SIG code, 
            #then has also the
            #District name, position and a numerical index i. 
            #The code is there two times, just
            #for redundancy.
            G.add_node(i_code, name=row_i[label], pos=(row_i["x"], row_i["y"]),
                       code = i_code, index = i)
            
            #Searching in the sindex for intersections. 
            #Maybe the query method would be better
            #but this is enough. Sindex creates nested bounding boxes 
            #around the given geometry
            #the interestion is not 100% correct, but we then double 
            #check the results to find
            #true neighbours
            for j in sindex.intersection(row_i.geometry.bounds):
                if j <= i:
                    continue
                row_j = map.iloc[j]
                inter = row_j.geometry.boundary.intersection(row_i.geometry.boundary)
                L = float(inter.length)
                if L >= MIN_BORDER:
                    j_code = name_dict.lookup(row_j[label])
                    pair_ij = (i_code, j_code)
                    pair_ji = (j_code, i_code)
                    G.add_edge(i_code, j_code, weight=pairs[pair_ij])
                    G.add_edge(j_code, i_code, weight=pairs[pair_ji])
   
   #The same logic like the graph from neighbours but without 
   #needing to check intersections
    if mode=="all":
        for i, row_i in map.iterrows():
            i_code = name_dict.lookup(row_i[label])
            G.add_node(i_code, name=row_i[label], pos=(row_i["x"], row_i["y"]),
                       code = i_code, index = i)
            for j, row_j in map.iterrows():
                if j <= i:
                    continue
                j_code = name_dict.lookup(row_j[label])
                pair_ij = (i_code, j_code)
                pair_ji = (j_code, i_code)
                G.add_edge(i_code, j_code, weight=pairs[pair_ij])
                G.add_edge(j_code, i_code, weight=pairs[pair_ji])

    if mode=="from_file":
        # nodes: same as neighbours mode, but simpler
        for i, row_i in map.iterrows():
            i_code = name_dict.lookup(row_i[label])
            G.add_node(
                i_code,
                name=row_i[label],
                pos=(row_i["x"], row_i["y"]),
                code=i_code,
                index=i,
            )

        # edges: exactly the pairs you pass in
        if pairs is not None:
            for (u, v), w in pairs.items():
                # only add if both endpoints exist in the map
                if u in G and v in G:
                    G.add_edge(u, v, weight=w)
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
            idx = containing_polygon.index[0]
        else:
            idx = 0
        return idx

#I SHOULD WRITE COMMENTS FOR THIS - TODO
class graphDisplay:
    def __init__(self, graph, gdt=None, map_name_col="SIG_KOR_NM",
                 gif_name=None):
        self.map = gdt
        self.graph = graph
        self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(right=0.80)
        self.cid = None
        self.position = nx.get_node_attributes(self.graph, 'pos')
        self.code = nx.get_node_attributes(self.graph, 'code')
        self.node_idxs = nx.get_node_attributes(self.graph, 'index')
        self.node_dict = dict(zip(list(self.node_idxs), list(self.code)))
        self.map_name_col = map_name_col

        self.in_norm = mcolors.Normalize(vmin=None, vmax=None, clip=True)
        self.out_norm = mcolors.Normalize(vmin=None, vmax=None, clip=True)
        self.save_name = gif_name

        self._blit_enabled = False
        self._sel_marker = None
        self._line_cache = {}
        self._label_cache = {}
        self._legend = None
        self._legend_neighs = []   
        self._legend_texts = []   

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

    def _blit_draw(self):
        if not self._blit_enabled:
            self.fig.canvas.draw_idle()

    # Dummy line used only as a legend handle (never drawn on axes)
    def _draw_or_get_line(self, u, v, in_out: bool):
        key = (u, v, in_out)
        ln = self._line_cache.get(key)
        if ln is None:
            ln, = self.ax.plot([], [], alpha=0.0)  # invisible
            self._line_cache[key] = ln
        return ln

    def _set_or_get_label(self, u, v):
        key = (u, v)
        lb = self._label_cache.get(key)
        if lb is None:
            x1, y1 = self.graph.nodes[u]["pos"]
            x2, y2 = self.graph.nodes[v]["pos"]
            mid_x, mid_y = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
            ang = np.atan2(y2 - y1, x2 - x1)
            if ang > np.pi / 2:
                ang -= np.pi
            if ang < -np.pi / 2:
                ang += np.pi
            n_vec = (-np.sin(ang), np.cos(ang))
            mult = 250
            lb = self.ax.text(mid_x + mult * n_vec[0], mid_y + mult * n_vec[1],
                              "",
                              ha="center", va="bottom",
                              rotation=0, rotation_mode="anchor", fontsize=8,
                              color="black", zorder=2)
            self._label_cache[key] = lb
        return lb

    def _expand_norm(self, norm: mcolors.Normalize, w: float):
        if norm.vmin is None or norm.vmax is None:
            norm.vmin = w
            norm.vmax = w + 1e-12
        else:
            changed = False
            if w < norm.vmin:
                norm.vmin = w
                changed = True
            if w > norm.vmax:
                norm.vmax = w
                changed = True
            if changed and norm.vmax == norm.vmin:
                norm.vmax = norm.vmin + 1e-12

    # SELECTION + LEGEND (called on mouse click, NOT every frame)
    def _update_selection(self, node):
        # move the red selection ring
        self._sel_marker.center = self.graph.nodes[node]["pos"]

        # collect in/out weights to all neighbours (using current graph attrs)
        edges = []
        for nb in self.graph.neighbors(node):
            w_out = self.graph.get_edge_data(node, nb)["weight"]
            w_in  = self.graph.get_edge_data(nb, node)["weight"]
            edges.append((node, nb, w_in, w_out))

        #Sort neighbours by total flow (largest first)
        # THIS IS COMPLETLY OPTIONAL SORTING, BUT MAKES IT A BIT EASIER
        # TO ORIENT
        edges.sort(key=lambda e: e[2] + e[3], reverse=True)

        # remember neighbour order for the dynamic updates later
        # (The legend displays neighbours and the weight values)
        self._legend_neighs = [v for (_, v, _, _) in edges]

        # build handles + initial labels
        handles = []
        labels = []
        for u, v, w_in, w_out in edges:
            ln = self._draw_or_get_line(u, v, True)  # dummy handle
            label = f"{self.code[v]}: in {w_in:.3g}, out {w_out:.3g}"
            handles.append(ln)
            labels.append(label)

        if self._legend is not None:
            self._legend.remove()

        self._legend = self.ax.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            frameon=False,
            fontsize=7,   
        )

        # store text objects so we can update them fast later
        # (IN THE TIME EVOLUTION, WE JUST CHANGE THE TEXT AND DONT
        # REDRAW THE WHOLE LEGEND ITSELF)
        self._legend_texts = self._legend.get_texts()

        self.selected_node = node

        # update the SIR time–series window shown only when node changes
        if self.values_ts is not None and self.ax_ts is not None:
            self._update_timeseries_lines()

        self.fig.canvas.draw_idle()

    # CLICK HANDLER
    def _on_click_graph(self, event):
        if event.inaxes is not self.ax or event.xdata is None or event.ydata is None:
            return

        # find closest node
        nodes = list(self.graph)
        px, py = event.xdata, event.ydata
        d2 = []
        for node in nodes:
            x, y = self.position[node]
            dx = px - x
            dy = py - y
            d2.append(dx * dx + dy * dy)

        idx = nodes[int(np.argmin(d2))]
        self._update_selection(idx)

    def interactive_graph(self):
        self.draw_graph()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click_graph)

    # STATIC GRAPH DRAWING 
    def draw_graph(self):
        self._sel_marker = Circle((float('nan'), float('nan')),
                                  radius=1000,
                                  fill=False,
                                  edgecolor='crimson',
                                  linewidth=2.5,
                                  zorder=80)
        self.ax.add_patch(self._sel_marker)

        for node in self.graph:
            circ = Circle(self.position[node],
                          radius=1000,
                          color="white",
                          zorder=50)
            self.ax.add_patch(circ)
            self.node_patches[node] = circ

            self.ax.text(self.position[node][0], self.position[node][1],
                         f"{self.code[node]}",
                         path_effects=[
                             pe.withStroke(linewidth=2, foreground="white")],
                         ha='center', va='center', zorder=100, fontsize=9)

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

        self.ax_ts.legend(loc="upper right")
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

    # LEGEND UPDATE PER FRAME (uses edge_weight_fn that was built in main)
    def _update_legend_for_frame(self, frame):
        if (
            self.edge_weight_fn is None
            or self.selected_node is None
            or self._legend is None
            or not self._legend_texts
        ):
            return

        u = self.selected_node
        for v, text in zip(self._legend_neighs, self._legend_texts):
            w_out = self.edge_weight_fn(frame, u, v)
            w_in  = self.edge_weight_fn(frame, v, u)
            text.set_text(f"{self.code[v]}: in {w_in:.3g}, out {w_out:.3g}")

    # ANIMATION
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

            # update legend text using time-dependent weights
            if self.edge_weight_fn is not None and self.selected_node is not None:
                self._update_legend_for_frame(frame)

            # move time marker only
            if self.time_marker is not None and self.t_vec is not None:
                t = self.t_vec[frame]
                self.time_marker.set_xdata([t, t])
                if self.ax_ts is not None:
                    self.ax_ts.figure.canvas.draw_idle()

            if self.t_vec is not None:
                t = self.t_vec[frame]
                self.ax.set_title(f"{self.var_names[k]}(t = {t:.3f}), frame {frame}")
            else:
                self.ax.set_title(f"{self.var_names[k]}(frame {frame})")

            return list(self.node_patches.values()) + (
                [self.time_marker] if self.time_marker is not None else []
            )

        self._anim = FuncAnimation(self.fig, update,
                                   frames=T,
                                   interval=interval,
                                   blit=False)

        if self.save_name is not None:
            from matplotlib.animation import PillowWriter
            writer = PillowWriter(fps=fps)
            self._anim.save(self.save_name, writer=writer)

        return self._anim



