from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
import matplotlib.patheffects as pe
from matplotlib.widgets import RadioButtons

import geopandas as gpd
import networkx as nx
from shapely import set_precision
from shapely.geometry import Point
import nametocode as ntc
import numpy as np

FLOW_TOP_K = 12

# Graph construction from a map file/data frame
def graph_construction(gdf=None, df=None, MIN_BORDER=25, label="SIG_KOR_NM",
               pairs=None, mode="all", weight_fn=None, name_dict=None):
    if name_dict is None:
        name_dict = ntc.code_dict()

    if gdf is not None:
        working_df = gdf.copy()
        working_df["geometry"] = working_df.geometry.buffer(0).apply(
            lambda g: set_precision(g, 0.2)
        )
        sindex = working_df.sindex
    
    elif df is not None:
        working_df = df.copy()

    else:
        raise ValueError(f"No data to convert to graph supplied!")

    G = nx.DiGraph()
    has_label = label in working_df.columns

    for i, row in working_df.iterrows():
        if has_label:
            code = name_dict.lookup(row[label])
            name = row[label]
        else:
            code = i
            name = str(i)
        G.add_node(code, name=name, pos=(row["x"], row["y"]), code=code, index=i)    

    def get_w(u, v):
        if weight_fn is not None:
            return float(weight_fn(u, v))
        if pairs is not None:
            w = pairs.get((u, v), 0.0)
            return 0.0 if w is None else float(w)
        return 1.0

    #Following are listed all possible graph building styles, but they are not
    #useful anymore
    if mode == "all":
        nodes = list(G.nodes())
        for a, u in enumerate(nodes):
            for v in nodes[a + 1:]:
                w_uv, w_vu = get_w(u, v), get_w(v, u)
                if w_uv: G.add_edge(u, v, weight=w_uv)
                if w_vu: G.add_edge(v, u, weight=w_vu)

    elif mode == "neighbours":
        if gdf is None: 
            raise TypeError(f"GeoDataFrame expected for nieghbour building")
        for i, row_i in working_df.iterrows():
            u = name_dict.lookup(row_i[label])
            for j in sindex.intersection(row_i.geometry.bounds):
                if j <= i:
                    continue
                row_j = gdf.iloc[j]
                L = float(row_j.geometry.boundary
                          .intersection(row_i.geometry.boundary).length)
                if L >= MIN_BORDER:
                    v = name_dict.lookup(row_j[label])
                    w_uv, w_vu = get_w(u, v), get_w(v, u)
                    if w_uv: G.add_edge(u, v, weight=w_uv)
                    if w_vu: G.add_edge(v, u, weight=w_vu)

    #This is the one that is used, the previous modes were here just for testing
    #And I leave them be if I needed to do some of that testing again.
    elif mode == "from_file":
        if pairs:
            for (u, v), w in pairs.items():
                if u in G and v in G and w:
                    G.add_edge(u, v, weight=float(w))
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    return G


# Geometry helper for clicking - returns GeoDataFrame row-index of the polygon
#which has inside it (x,y) or returns None if the point is in no polygon
def intercept_check(x, y, polygons):
    pt = Point(x, y)
    candidates = polygons.sindex.query(pt)
    hits = polygons.iloc[candidates]
    mask = hits.covers(pt)
    containing = mask[mask]
    return containing.index[0] if not containing.empty else None


# Interactive display class
# Three separate figure windows are created on demand:
#    1. Main map figure     – animation of evolution or slider of the selected variable.
#    2. Timeseries figure   – per-node time series of all the simulated variables.
#    3. Flow figure         – displays the underlying flow rates between the selected node and
#                             its top-K neighbours, based on the top K specified at the top of
#                             this file.
#
# Typical call sequence is as follows:
#    city_graph = graphDisplay(G, city_map)
#    city_graph.interactive_graph()        
#    city_graph.start_animation(...)        
#    city_graph.attach_flow_timeseries(...)  
#    plt.show()

class graphDisplay:
    
    # Construction
    def __init__(self, graph, gdt=None, map_name_col="SIG_KOR_NM", gif_name=None):
        self.graph        = graph
        self.map          = gdt.copy() if gdt is not None else None
        self.map_name_col = map_name_col
        self.save_name    = gif_name

        # node look-ups with consistent indexing 
        self.node_idxs = nx.get_node_attributes(graph, "index")

        # main figure - the map of evolution with slider or animation
        self.fig      = plt.figure()
        gs            = GridSpec(2, 1, figure=self.fig,
                                 height_ratios=[1.0, 10.0])
        self.ax       = self.fig.add_subplot(gs[1, 0])
        self.ax_title = self.fig.add_subplot(gs[0, 0])
        self.ax_title.axis("off")
        self.fig.subplots_adjust(left=0.02, right=0.98,
                                 top=0.95, bottom=0.08)

        # The text is a separate subplot because previously I used blitting 
        # and putting the text as a separate entity fixed many issues. Blitting
        #is no longer used but the separate text remains.
        self._title_text = self.ax_title.text(
            0.5, 0.5, "",
            transform=self.ax_title.transAxes,
            ha="center", va="top", fontsize=10,
            clip_on=True, zorder=200,
            bbox=dict(facecolor="white", edgecolor="none",
                      alpha=0.9, pad=2.0),
        )

        # map polygon collection (filled in by draw_graph) - usually used from the GPD
        self._map_collection       = None
        self._prev_selected_idx    = None
        self._edge_colors          = None
        self._edge_widths          = None

        # colorbar storage
        self.cbar             = None
        self._field_ts        = None   # (T, N) array for the displayed var - T is snapshots in time and N is number of nodes
        self._field_var_name  = None
        self._field_norm      = None
        self._field_cmap      = None

        # Values from simulation to display
        self.values_ts  = None   # (T, K, N) array,  so all values, K is number of calculated vals. _field_ts selects from this
        self.var_names  = None
        self.dt_snapshot = None
        self.t_vec      = None

        # variable timeseries figure vars (per-node variables display window)
        self.ts_fig      = None
        self.ax_ts       = None
        self.lines_vars  = []
        self.time_marker = None

        # flow Timeseries figure vars (display of weights from the dataset in time)
        self.flow_ts_fig          = None
        self.ax_flow_ts           = None
        #Radio button windows - small individual axes 
        self.ax_flow_radio_neigh  = None
        self.ax_flow_radio_dir    = None
        # The buttons themselves
        self.radio_flow_neigh     = None
        self.radio_flow_dir       = None
        self.flow_time_marker     = None
        self._flow_series_cache   = {}
        self._selected_flow_neighbor = None
        self._selected_flow_dir   = "both"
        self._flow_ts_top_k       = FLOW_TOP_K

        # interaction state saved       
        self.selected_node   = None
        self.edge_weight_fn  = None
        self._current_frame  = 0
        #These are here so they do not get garpage collected
        self._anim           = None
        self._slider         = None
        self._slider_ax      = None

    # Map drawing method. Renders the map into self.ax. A GeoDataFrame map is needed.
    def draw_graph(self):
        if self.map is None:
            raise ValueError("graphDisplay requires a GDF map.")

        self.ax.clear()
        self.map = self.map.copy()

        #Previously I did not build the code column in main and did it here, so I leave it here 
        # incase I wanted to move back to the previous style
        if "code" not in self.map.columns:
            nd = ntc.code_dict()
            self.map["code"] = self.map[self.map_name_col].map(nd.lookup)

        # Plotting of the base map 
        self.map.plot(ax=self.ax, color="white",
                      edgecolor="black", linewidth=0.8, zorder=10)
        self.ax.set_aspect("equal")

        # Gets the last object added to the axis
        self._map_collection = self.ax.collections[-1]
        n = len(self.map)
        #Set all polygon faces white and edges black
        self._edge_colors = np.tile([[0., 0., 0., 1.]], (n, 1))
        self._edge_widths = np.full(n, 0.8)
        self._map_collection.set_edgecolor(self._edge_colors)
        self._map_collection.set_linewidth(self._edge_widths)

        # Setting of code text at representative points
        for _, row in self.map.iterrows():
            x = row["x"] if "x" in row.index else row.geometry.representative_point().x
            y = row["y"] if "y" in row.index else row.geometry.representative_point().y
            self.ax.text(
                x, y, f"{int(row['code'])}",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")],
                ha="center", va="center", zorder=40, fontsize=6,
            )

        # Setting of display artributes
        self.ax.autoscale_view()
        self.ax.set_aspect("equal")
        self.ax.set_autoscale_on(False)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)

    #Drawing of graph with handling of mouse-clicking, draws the under graphs and calls the event handler
    def interactive_graph(self):
        self.draw_graph()
        self.fig.canvas.mpl_connect("button_press_event", self._on_click_graph)

    # Click event handler
    def _on_click_graph(self, event):
        #Checking if click is in correct location
        if event.inaxes is not self.ax or event.xdata is None:
            return

        #Finding clicked polygon 
        poly_idx = intercept_check(event.xdata, event.ydata,
                                   self.map[["geometry"]])
        if poly_idx is None:
            return
        
        node = int(self.map.iloc[poly_idx]["code"])
        self._update_selection(node)

    #Highlit selected node and refresh panels dependent on it
    def _update_selection(self, node):
        # if there was any previous selection, turn it back to black and white 
        if self._prev_selected_idx is not None:
            self._edge_colors[self._prev_selected_idx] = [0., 0., 0., 1.]
            self._edge_widths[self._prev_selected_idx] = 0.8

        #Change color of current selected node 
        idx = self.node_idxs.get(node) 
        if idx is not None:
            self._edge_colors[idx] = [1., 0., 0., 1.]
            self._edge_widths[idx] = 2.4

        #Now update all needed variables
        self._prev_selected_idx = idx
        self._map_collection.set_edgecolor(self._edge_colors)
        self._map_collection.set_linewidth(self._edge_widths)

        self.selected_node = node

        # Updating of timeseries, vars and flow
        if self.values_ts is not None and self.ts_fig is not None:
            self._update_timeseries_lines()

        if self.flow_ts_fig is not None:
            self._rebuild_flow_series_cache()
            self._build_flow_controls()
            self._draw_selected_flow_timeseries()

        # If animation is happening, throw away cashed data so all works, else redraw
        if self._anim is not None:
            self._anim._init_draw()
        else:
            self.fig.canvas.draw_idle()

    # Simulation data showing a timeseries window, stores (T, K, N) simulation output and opens the timeseries window
    def attach_timeseries(self, values_ts, dt_snapshot, var_names=None):
        vals = np.asarray(values_ts) #If not already np array 
        T, K, N = vals.shape

        # Storage of values 
        self.values_ts = vals
        self.var_names = var_names if var_names is not None \
            else [f"var{k}" for k in range(K)]
        #Turn into list, ensures valid data
        self.var_names = list(self.var_names)

        #If valid timestep was supplied, builds axis using this timestep, else just numbers it
        if dt_snapshot and dt_snapshot > 0:
            self.dt_snapshot = float(dt_snapshot)
            self.t_vec = np.arange(T, dtype=float) * self.dt_snapshot
        else:
            self.dt_snapshot = None
            self.t_vec = np.arange(T, dtype=float)

       #If timeseries of vars nonexistant, create it 
        if self.ts_fig is None:
            self.ts_fig, self.ax_ts = plt.subplots()

        # Set of default node 
        if self.selected_node is None:
            self.selected_node = next(iter(self.node_idxs))

        # Index of selected node for further working
        idx = self.node_idxs[self.selected_node]
        node_vals = vals[:, :, idx]   # (T, K)

        #If I add an option to rerun later, I add clearing
        self.ax_ts.clear()

        #Initialization of the window and adding of values
        self.lines_vars = []
        for k in range(K):
            line, = self.ax_ts.plot(self.t_vec, node_vals[:, k],
                                    label=self.var_names[k])
            self.lines_vars.append(line)

        #Current time marker
        self.time_marker = self.ax_ts.axvline(
            self.t_vec[0], color="k", linestyle="--", linewidth=1
        )

        #Labels of the window
        self.ax_ts.set_xlabel("time")
        self.ax_ts.set_ylabel("value")
        self.ax_ts.set_title(f"values at node {self.selected_node}")
        self.ax_ts.set_xlim(self.t_vec[0], self.t_vec[-1])
        self.ax_ts.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0),
                          borderaxespad=0.0, frameon=True)
        self.ts_fig.subplots_adjust(right=0.78)
        self.ts_fig.canvas.draw_idle()

    #Redraw timeseries for the currently selected node via updating artists.
    # the time marker is moveed by the _render_frame()
    def _update_timeseries_lines(self):

        #Update of timeseries for the selected node
        idx  = self.node_idxs[self.selected_node]
        vals = self.values_ts[:, :, idx]
        for k, line in enumerate(self.lines_vars):
            line.set_ydata(vals[:, k])

        self.ax_ts.set_title(f"values at node {self.selected_node}")
        self.ax_ts.relim()
        self.ax_ts.autoscale_view()
        self.ts_fig.canvas.draw_idle()

    # Field rendering shared by animation and slider - call attach_timeseries, set up colormap/colorbar.
    def _prepare_field_view(self, var_names, values_ts, dt_snapshot,
                            var=0, cmap_name="plasma", edge_weight_fn=None):
        vals = np.asarray(values_ts)
        #Expands (T, N) to (T, K, N) with K = 1
        if vals.ndim == 2:
            vals = vals[:, np.newaxis, :]
        if vals.ndim != 3:
            raise ValueError("values_ts must be (T, K, N) or (T, N)")

        T, K, N = vals.shape
        self.attach_timeseries(vals, dt_snapshot, var_names)
        self.edge_weight_fn = edge_weight_fn

        # resolve variable index - turn name into index or ensure number as index
        if isinstance(var, str):
            if var not in self.var_names:
                raise ValueError(f"Unknown variable {var!r}")
            k = self.var_names.index(var)
        else:
            k = int(var)
            if not 0 <= k < K:
                raise ValueError(f"var index {k} out of range [0, {K})")

        #Values for display in the main map - so the selected variable
        field_ts = vals[:, k, :]
        vmin, vmax = float(field_ts.min()), float(field_ts.max())
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        cmap = cm.get_cmap(cmap_name)

        # rebuild colorbar, if this somehow gets called second time sometimes later
        if self.cbar is not None:
            self.cbar.remove()
        # Individual colorbar, because there is no object that the colormap would be called at
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        self.cbar = self.fig.colorbar(sm, ax=self.ax,
                                      orientation="horizontal",
                                      fraction=0.05, pad=0.10)
        self.cbar.set_label(self.var_names[k])

        self._field_ts       = field_ts
        self._field_var_name = self.var_names[k]
        self._field_norm     = norm
        self._field_cmap     = cmap
        
    #Update map colors and time markers for given frame from animation or slider 
    def _render_frame(self, frame, redraw_ts=False):
        frame      = int(frame)
        vals_frame = self._field_ts[frame]

        facecolors = []
        for _, row in self.map.iterrows():
            node = row.get("code")
            if node is None or node not in self.node_idxs:
                facecolors.append((1., 1., 1., 1.))
            else:
                j = self.node_idxs[node]
                facecolors.append(
                    self._field_cmap(self._field_norm(vals_frame[j]))
                )
        self._map_collection.set_facecolor(facecolors)

        self._current_frame = frame
        t = self.t_vec[frame] if self.t_vec is not None else frame

        if self.time_marker is not None:
            self.time_marker.set_xdata([t, t])
            if redraw_ts and self.ax_ts is not None:
                self.ax_ts.figure.canvas.draw_idle()

        if self.flow_time_marker is not None:
            self.flow_time_marker.set_xdata([t, t])
            if redraw_ts and self.flow_ts_fig is not None:
                self.flow_ts_fig.canvas.draw_idle()

        self._title_text.set_text(
            f"{self._field_var_name}  t = {t:.3f} h  (frame {frame})"
        )

    # Logic of animation
    def start_animation(self, var_names, values_ts, dt_snapshot, var=0,
                        interval=60, cmap_name="plasma", fps=10,
                        edge_weight_fn=None):
        self._prepare_field_view(var_names, values_ts, dt_snapshot,
                                     var, cmap_name, edge_weight_fn)
        self._render_frame(0)
        T = values_ts.shape[0]

        def update(frame):
            self._render_frame(frame, redraw_ts=True)

        self._anim = FuncAnimation(self.fig, update, frames=T,
                                   interval=interval, blit=False)

        if self.save_name is not None:
            from matplotlib.animation import PillowWriter
            self._anim.save(self.save_name, writer=PillowWriter(fps=fps))

        return self._anim

    # Logic of slider
    def show_frame_slider(self, var_names, values_ts, dt_snapshot, var=0,
                          cmap_name="plasma", edge_weight_fn=None,
                          slider_label="frame"):
        self._prepare_field_view(var_names, values_ts, dt_snapshot,
                                     var, cmap_name, edge_weight_fn)

        T = values_ts.shape[0]

        if self._slider_ax is None:
            self.fig.subplots_adjust(bottom=0.14)
            self._slider_ax = self.fig.add_axes([0.18, 0.03, 0.55, 0.03])
        else:
            self._slider_ax.clear()

        self._slider = Slider(self._slider_ax, slider_label,
                              valmin=0, valmax=T - 1,
                              valinit=0, valstep=1)
        self._render_frame(0)

        def on_change(val):
            self._render_frame(int(val), redraw_ts=True)
            self.fig.canvas.draw_idle()

        self._slider.on_changed(on_change)
        self.fig.canvas.draw_idle()
        return self._slider

    # Flow timeseries window attachment - Open (or if needed refresh) the flow timeseries figure.
    def attach_flow_timeseries(self, top_k=FLOW_TOP_K):
        self._flow_ts_top_k = top_k

        if self.flow_ts_fig is None:
            self.flow_ts_fig         = plt.figure()
            self.ax_flow_ts          = self.flow_ts_fig.add_axes(
                [0.08, 0.12, 0.70, 0.78])
            self.ax_flow_radio_neigh = self.flow_ts_fig.add_axes(
                [0.82, 0.46, 0.15, 0.38])
            self.ax_flow_radio_dir   = self.flow_ts_fig.add_axes(
                [0.82, 0.22, 0.15, 0.16])

        if self.selected_node is None:
            self.selected_node = list(self.graph.nodes())[0]

        self._rebuild_flow_series_cache()
        self._build_flow_controls()
        self._draw_selected_flow_timeseries()



    # Precompute hourly in/out flow series for top-K neighbours of selected node
    def _rebuild_flow_series_cache(self):
        self._flow_series_cache = {}

        # Nothing to compute if prerequisites are missing
        if (self.selected_node is None
                or self.edge_weight_fn is None
                or self.t_vec is None):
            return

        u = self.selected_node

        # Dense hourly time axis, can be finer than simulation snapshots for accurate curves
        # Currently it works at the resolution of hours, because finer does not currently make sense
        # with the weights linearly interpolated between hours
        total_time = float(self.t_vec[-1])
        flow_t     = np.arange(max(1, int(np.ceil(total_time)) + 1), dtype=float)
        neighs = []
        for v in self.graph.nodes():
            if v == u:
                continue
            # Only consider nodes that have any edge with u
            if not (self.graph.has_edge(u, v) or self.graph.has_edge(v, u)):
                continue

            out_s = np.array([self.edge_weight_fn(t, u, v) for t in flow_t], dtype=float)
            in_s  = np.array([self.edge_weight_fn(t, v, u) for t in flow_t], dtype=float)
            # Score by peak flow. This is used to select top-K
            score = max(float(out_s.max()), float(in_s.max()))
            if score > 0:
                neighs.append((score, v, flow_t, in_s, out_s))

        # Keep only top-K neighbours by peak flow for the flow window
        neighs.sort(key=lambda x: x[0], reverse=True)
        for _, v, ft, in_s, out_s in neighs[:self._flow_ts_top_k]:
            self._flow_series_cache[v] = {"t": ft, "in": in_s, "out": out_s}

        # If previously selected neighbour is no longer in top-K reset to first in the list
        if self._flow_series_cache:
            if self._selected_flow_neighbor not in self._flow_series_cache:
                self._selected_flow_neighbor = next(iter(self._flow_series_cache))
        
        else:
            self._selected_flow_neighbor = None

    #Build (or rebuild) the radio buttons for the flow window
    def _build_flow_controls(self):
        if self.ax_flow_radio_neigh is None:
            return

        self.ax_flow_radio_neigh.clear()
        self.ax_flow_radio_dir.clear()

        neigh_labels = [str(v) for v in self._flow_series_cache]
        if not neigh_labels:
            self.ax_flow_radio_neigh.text(
                0.1, 0.5, "No flows",
                transform=self.ax_flow_radio_neigh.transAxes)
            self.ax_flow_radio_neigh.set_axis_off()
            self.ax_flow_radio_dir.set_axis_off()
            self.radio_flow_neigh = self.radio_flow_dir = None
            return

        active_neigh = neigh_labels.index(str(self._selected_flow_neighbor))
        self.ax_flow_radio_neigh.set_title("Neighbour", fontsize=9)
        self.radio_flow_neigh = RadioButtons(
            self.ax_flow_radio_neigh, neigh_labels, active=active_neigh)

        dir_labels   = ["both", "out", "in"]
        active_dir   = dir_labels.index(self._selected_flow_dir)
        self.ax_flow_radio_dir.set_title("Direction", fontsize=9)
        self.radio_flow_dir = RadioButtons(
            self.ax_flow_radio_dir, dir_labels, active=active_dir)

        def on_neigh(label):
            try:
                self._selected_flow_neighbor = int(label)
            except ValueError:
                self._selected_flow_neighbor = label
            self._draw_selected_flow_timeseries()

        def on_dir(label):
            self._selected_flow_dir = label
            self._draw_selected_flow_timeseries()

        self.radio_flow_neigh.on_clicked(on_neigh)
        self.radio_flow_dir.on_clicked(on_dir)

    #Plot in/out series for the currently selected neighbour.
    def _draw_selected_flow_timeseries(self):
        if self.ax_flow_ts is None:
            return

        self.ax_flow_ts.clear()
        u, v = self.selected_node, self._selected_flow_neighbor

        if u is None or v is None or v not in self._flow_series_cache:
            self.ax_flow_ts.set_title("No selected flow")
            self.flow_ts_fig.canvas.draw_idle()
            return

        cache    = self._flow_series_cache[v]
        flow_t = cache["t"]
        in_s     = cache["in"]
        out_s    = cache["out"]

        if self._selected_flow_dir in ("both", "out"):
            self.ax_flow_ts.plot(flow_t, out_s, label=f"{u} → {v}")
        if self._selected_flow_dir in ("both", "in"):
            self.ax_flow_ts.plot(flow_t, in_s, linestyle="--",
                                 label=f"{v} → {u}")

        self.ax_flow_ts.set_title(f"Flow  {u} ↔ {v}")
        self.ax_flow_ts.set_xlabel("time (h)")
        self.ax_flow_ts.set_ylabel("weight")
        self.ax_flow_ts.legend(loc="upper left")

        # vertical marker at current animation frame
        if self.t_vec is not None:
            t = self.t_vec[self._current_frame]
            self.flow_time_marker = self.ax_flow_ts.axvline(
                t, color="k", linestyle="--", linewidth=1)

        self.flow_ts_fig.canvas.draw_idle()
