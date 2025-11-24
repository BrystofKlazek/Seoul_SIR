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

    #Search spatial index so that searching is faster (builds a tree for searching in the tree)
    sindex = map.sindex
    G = nx.DiGraph()
   
    #If we want to build graph from neighbours
    if mode=="neighbours":
        for i, row_i in map.iterrows():
            i_code = name_dict.lookup(row_i[label])
            #Adding nodes to the graph. Each node is indexed by the SIG code, then has also the
            #District name, position and a numerical index i. The code is there two times, just
            #for redundancy.
            G.add_node(i_code, name=row_i[label], pos=(row_i["x"], row_i["y"]),
                       code = i_code, index = i)
            
            #Searching in the sindex for intersections. Maybe the query method would be better
            #but this is enough. Sindex creates nested bounding boxes around the given geometry
            #the interestion is not 100% correct, but we then double check the results to find
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
   
   #The same logic like the graph from neighbours but without needing to check intersections
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


#A class to display the graph
class graphDisplay:
    # self.map - map to display (if not None)
    # self.graph - graph to display
    # self.graph and self.ax - stuff for plotting using matplotlib
    # self.cid - what to do incase of an event with matplotlib 
    # self.position - a dict mapping node to its position
    # self.node_idxs - a dict mapping a node to its numerical idx 
    # self.node_dict - a reverse dict mapping an node to its index
    # self.in_cmap - colormap wor out edges, self.out_smap - the same for out
    # self.save_name - where to save gif of evolution (if wanted)
    # self._sel_marker - selection marker (a red ring)
    # self._line_cache - a helper line cache for matplotlib for the interaction
    # with clicks
    # self._label_cache - a helper label cache in the same vein
    # self._legend - a helper for a label 

    # self.node_patches - a helper for nodes (mainly for color of evolution)     

    # self.values_ts - values into the timeseries graph       
    # self.var_names - names of vars calculated on the graph         

    # self.selected_node - What node is currently selected? Helper mainly for the timeseries.   
    # self.ts_fig - matplotlib stuff for plotting the timeseries
    # self.ax_ts - in the same vein         
    # self.lines_vars - plotted values of the vars     
    # self.time_marker - displayed vertical lines of current time in timeseries     

    # self.dt_snapshot - helper to calculate current time     
    # self.t_vec - sim           
    # self.cbar - cbar displaying what values each node has          

    # self._anim - helper for if animation is to be or not        

    def __init__(self, graph, gdt = None, map_name_col = "SIG_KOR_NM",
                 gif_name = None):
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
        
        self.in_cmap  = mcolors.LinearSegmentedColormap.from_list(
            "inmap", [
                "#005f73",  
                "#0a9396",  
                "#3b82f6", 
                "#312e81",  
            ])

        self.out_cmap = mcolors.LinearSegmentedColormap.from_list(
            "outmap", [
                "#ffb703", 
                "#fb8500",
                "#e63946",
                "#7f1d1d",
            ])        

        self.in_norm  = mcolors.Normalize(vmin=None, vmax=None, clip=True)
        self.out_norm = mcolors.Normalize(vmin=None, vmax=None, clip=True)
        self.save_name = gif_name

        #For further speedup down the line it seems promising to add blitting if
        #it is needed to make this part of the code run faster. Depends
        #on further implementation. This is left as an opportunity. But I hope
        #it wont be needed 
        self._blit_enabled = False
        self._sel_marker = None
        self._line_cache = {}
        self._label_cache = {}
        self._legend = None

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

        # NEW: callback for time-dependent edge weights
        self.edge_weight_fn = None

    #Just a placeholder, not used now
    def _blit_draw(self):
        if not self._blit_enabled:
            self.fig.canvas.draw_idle()
            return
        else:
            return
    
    #Draw or get line between u and v (and if it is both ways or only one way)
    def _draw_or_get_line(self, u, v, in_out: bool):
            key = (u, v, in_out)
            ln = self._line_cache.get(key)
            if ln is None:
                x1, y1 = self.graph.nodes[u]["pos"]; x2, y2 = self.graph.nodes[v]["pos"]
                #Calculating angle and normal vector so we can plot the line correctly
                # (Slightly off from the mid of the path between the nodes so we can plot two 
                # lines )
                ang = np.atan2(y2-y1, x2-x1)
                n_vec = (-np.sin(ang), np.cos(ang))
                mult = 120
                offset = mult if in_out == True else -mult
                color = "blue" if in_out == True else "red"
                ln, = self.ax.plot([x1+offset*n_vec[0], x2+offset*n_vec[0]], 
                                    [y1+offset*n_vec[1], y2+offset*n_vec[1]], 
                                    color=color, linestyle="dashed", zorder=1)
                #adding it to the dict
                self._line_cache[key] = ln
            return ln
    
    #Similar logic to the drawing of the lines
    def _set_or_get_label(self, u, v):
        key = (u, v)
        lb = self._label_cache.get(key)
        if lb is None:
            x1, y1 = self.graph.nodes[u]["pos"]; x2, y2 = self.graph.nodes[v]["pos"]
            mid_x, mid_y = 0.5*(x1+x2), 0.5*(y1+y2)
            ang = np.atan2(y2-y1, x2-x1)
            if ang > np.pi/2: ang -= np.pi
            if ang < -np.pi/2: ang += np.pi
            n_vec = (-np.sin(ang), np.cos(ang))
            mult = 250
            #The label is empty for now ("") - we will add text later
            lb = self.ax.text(mid_x+mult*n_vec[0], mid_y+mult*n_vec[1], 
                              "", ha="center", va="bottom",
                          rotation=0, rotation_mode="anchor", fontsize=8,
                          color="black", zorder=2)
            self._label_cache[key] = lb
        return lb

    #Norming for colorbars
    def _expand_norm(self, norm: mcolors.Normalize, w: float):
        if norm.vmin is None or norm.vmax is None:
            norm.vmin = w
            norm.vmax = w + 1e-12
        else:
            changed = False
            if w < norm.vmin:
                norm.vmin = w; changed = True
            if w > norm.vmax:
                norm.vmax = w; changed = True
            if changed and norm.vmax == norm.vmin:
                norm.vmax = norm.vmin + 1e-12

    #Coloring edges baset on weights
    def _set_weight_color(self, ln, w: float, in_out: bool):
        if in_out == True:
            self._expand_norm(self.in_norm, w)
            ln.set_color(self.in_cmap(self.in_norm(w)))
        else:
            self._expand_norm(self.out_norm, w)
            ln.set_color(self.out_cmap(self.out_norm(w)))

    #After click logic
    def _update_selection(self, node):
        # move the red selection ring
        self._sel_marker.center = self.graph.nodes[node]["pos"]

        # clear previous highlighted lines & legend labels
        for ln in self._line_cache.values():
            ln.set_alpha(0.0)
            ln.set_label(None)

        # collect in/out weights to all neighbours
        edges = []
        for idxn in self.graph.neighbors(node):
            w_out = self.graph.get_edge_data(node, idxn)["weight"]
            w_in  = self.graph.get_edge_data(idxn, node)["weight"]
            edges.append((node, idxn, w_in, w_out))

        # optional: sort neighbours by total flow (largest first)
        edges.sort(key=lambda e: e[2] + e[3], reverse=True)

        # remove any old numeric labels that were drawn on edges
        for lb in self._label_cache.values():
            lb.set_text("")

        # highlight edges + create clean legend entries
        for u, v, w_in, w_out in edges:
            # in-flow line (blue)
            ln1 = self._draw_or_get_line(u, v, True)
            ln1.set_linewidth(1.5)
            self._set_weight_color(ln1, w_in, True)
            ln1.set_alpha(1.0)

            # out-flow line (red)
            ln2 = self._draw_or_get_line(u, v, False)
            ln2.set_linewidth(1.5)
            self._set_weight_color(ln2, w_out, False)
            ln2.set_alpha(1.0)

            # ONE legend entry per neighbour, attached to ln1
            # self.code[v] is the SIG of the neighbour
            label = f"{self.code[v]}: in {w_in:.3g}, out {w_out:.3g}"
            ln1.set_label(label)
            ln2.set_label(None)  # don't duplicate in legend

        self.selected_node = node
        if self.values_ts is not None and self.ax_ts is not None:
            self._update_timeseries_lines()
   
    #On click logic
    def _on_click_graph(self, event):
        
        if event.inaxes is self.ax and event.xdata and event.ydata:
            nodes = list(self.graph)     
            point = (event.xdata, event.ydata)
            distances = np.array([])
            for node in nodes:
                difference = np.subtract(point, self.position[node])
                distance = difference[0]**2+difference[1]**2
                distances = np.append(distances, distance)

            idx = nodes[np.argmin(distances)]
            self._update_selection(idx)
            
        if self._legend is not None:
            self._legend.remove()

        self._legend = self.ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),   
            borderaxespad=0.0,
            frameon=False,
        )

        self.fig.canvas.draw_idle()  
   
    #Callable command that draws the graph and then connects clicking listening
    def interactive_graph(self):
        self.draw_graph()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click_graph)

    #Drawing of graph logic 
    def draw_graph(self):
        self._sel_marker =Circle((float('nan'), float('nan')),
                              radius=1000,   
                              fill=False,
                              edgecolor='crimson',
                              linewidth=2.5,
                              zorder=80)
        self.ax.add_patch(self._sel_marker)
        for node in self.graph:
            circ = Circle(
                    self.position[node], 
                    radius=1000,
                    color = "white",
                    zorder = 50) 
            self.ax.add_patch(circ)

            self.node_patches[node] = circ

            self.ax.text(self.position[node][0], self.position[node][1],
                         f"{self.code[node]}",
                         path_effects=[
                            pe.withStroke(linewidth=2, foreground="white")],
                         ha='center', va='center', zorder=100, fontsize=9)

        # --- build a list of line segments for edges ---
        segments = []
        weights = []

        for u, v in self.graph.edges():
            w = self.graph.get_edge_data(u, v)["weight"]
            if w <= 0:
                continue

            x1, y1 = self.graph.nodes[u]["pos"]
            x2, y2 = self.graph.nodes[v]["pos"]
            segments.append([(x1, y1), (x2, y2)])
            weights.append(w)

        # OPTIONAL: randomly subsample edges if there are too many to draw nicely
        MAX_EDGES_TO_DRAW = 3000
        if len(segments) > MAX_EDGES_TO_DRAW:
            idx = random.sample(range(len(segments)), MAX_EDGES_TO_DRAW)
            segments = [segments[i] for i in idx]
            weights = [weights[i] for i in idx]

        # create a single LineCollection instead of thousands of ax.plot calls
        edge_collection = LineCollection(
            segments,
            linewidths=0.4,
            alpha=0.0
        )
        self.ax.add_collection(edge_collection)
        self.ax.autoscale_view() 
        self.ax.set_aspect('equal')
        self.ax.set_autoscale_on(False)
        self.ax.set_xticks([]) 
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)

    #Attaching/drawing of the timeseries (plot with the evolution of S/I/R)
    def attach_timeseries(self, values_ts, dt_snapshot,
                          var_names=None):
        vals = np.asarray(values_ts)
        if vals.ndim == 2:
            vals = vals[:, np.newaxis, :]   # (T, 1, N)
        if vals.ndim != 3:
            raise ValueError("values_ts must have shape (T, K, N) or (T, N)")

        T, K, N = vals.shape
        if len(list(self.graph.nodes())) != N:
            raise ValueError("len(node_order) must equal N (third axis of values_ts)")

        if var_names is None:
            var_names = [f"var{k}" for k in range(K)]
        if len(var_names) != K:
            raise ValueError("len(var_names) must equal K (second axis of values_ts)")

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
        node_vals = vals[:, :, idx]   # (T, K)

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
        # explicitly show the full [0, T-1] (or 0..T*dt)
        self.ax_ts.set_xlim(self.t_vec[0], self.t_vec[-1])
        self.ax_ts.relim()
        self.ax_ts.autoscale_view()
        self.ts_fig.canvas.draw_idle()
   
    #Updating of the vertical line
    def _update_timeseries_lines(self):
        if self.values_ts is None or not self.lines_vars:
            return
        if self.selected_node not in self.node_idxs:
            return

        idx = self.node_idxs[self.selected_node]
        vals = self.values_ts[:, :, idx]  # (T, K)

        for k, line in enumerate(self.lines_vars):
            line.set_ydata(vals[:, k])

        self.ax_ts.set_title(f"values at node {self.selected_node}")
        if self.t_vec is not None:
            self.ax_ts.set_xlim(self.t_vec[0], self.t_vec[-1])
        self.ax_ts.relim()
        self.ax_ts.autoscale_view()
        self.ts_fig.canvas.draw_idle()
   
    #Animation inside the graph and the timeseries
    def start_animation(self, var_names, values_ts, dt_snapshot, var=0, interval=60,
                        cmap_name="plasma", fps=10, edge_weight_fn=None): 
        #Parameters to input:
        #   var - int or str
        #       Which variable to animate (only 1) index or name from var_names.
        #   interval - int
        #       Milliseconds between frames.
        #   cmap_name - str
        #       Matplotlib colormap name (I choose from default cmaps here...)
        #   save - str or None
        #       If given, path to save GIF.
        #   fps - int
        #       Frames per second for GIF.

        vals = np.asarray(values_ts)
        if vals.ndim == 2:
            vals = vals[:, np.newaxis, :]   # (T, 1, N)
        if vals.ndim != 3:
            raise ValueError("values_ts must have shape (T, K, N) or (T, N)")

        T, K, N = vals.shape
        self.attach_timeseries(values_ts, dt_snapshot, var_names)

        # store callback for dynamic edge weights
        self.edge_weight_fn = edge_weight_fn
        
        # Which variable is attached - raise error if unknown
        if isinstance(var, str):
            if self.var_names is None or var not in self.var_names:
                raise ValueError(f"Unknown variable name {var!r}")
            k = self.var_names.index(var)
        else:
            k = int(var)
            if not (0 <= k < K):
                raise ValueError(f"var index {k} out of range [0, {K})")

        field_ts = vals[:, k, :]  # (T, N)

        #Colour normalization (So the colorbars are useful)
        vmin = float(field_ts.min())
        vmax = float(field_ts.max())
        if vmax <= vmin:
            vmax = vmin + 1e-12

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        cmap = cm.get_cmap(cmap_name)

        #colorbar under the graph
        # ScalarMappable just to connect cmap + norm to a colorbar
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        self.cbar = self.fig.colorbar(
            sm,
            ax=self.ax,
            orientation="horizontal",
            fraction=0.05,   
            pad=0.10,        
        )

        #Label the colorbar with the variable name
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

            # update edge weights if callback provided
            if self.edge_weight_fn is not None:
                for u, v in self.graph.edges():
                    self.graph[u][v]["weight"] = self.edge_weight_fn(frame, u, v)

                # refresh selection highlighting with new weights
                if self.selected_node is not None:
                    self._update_selection(self.selected_node)

            #Move time marker in the time-series window
            if self.time_marker is not None and self.t_vec is not None:
                t = self.t_vec[frame]
                self.time_marker.set_xdata([t, t])
                if self.ax_ts is not None:
                    self.ax_ts.figure.canvas.draw_idle()
            #t_vec holds the physical time of each snapshot if attach_timeseries
            #was called with a valid dt_snapshot - so we use it here to know at what time
            #we are
            if self.t_vec is not None:
                t = self.t_vec[frame]
                self.ax.set_title(f"{self.var_names[k]}(t = {t:.3f}), frame {frame}")
            else:
                # fall back to frame index as "time"
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

