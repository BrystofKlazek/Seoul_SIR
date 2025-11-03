import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.patches import Circle

import networkx as nx

from shapely import set_precision
from shapely.geometry import Point

import nametocode as ntc

import numpy as np

def maptograph(map, MIN_BORDER = 25, label="SIG_KOR_NM", pairs={}, mode="all",
               name_dict = ntc.code_dict()): 


    GRID = 0.2   
    map["geometry"] = map.geometry.buffer(0).apply(
            lambda geom: set_precision(geom, GRID)
            )

    sindex = map.sindex
    G = nx.DiGraph()
    
    if mode=="neighbours":
        for i, row_i in map.iterrows():
            i_code = name_dict.lookup(row_i[label])
            G.add_node(i_code, name=row_i[label], pos=(row_i["x"], row_i["y"]),
                       code = i_code)
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
   
    if mode=="all":
        for i, row_i in map.iterrows():
            i_code = name_dict.lookup(row_i[label])
            G.add_node(i_code, name=row_i[label], pos=(row_i["x"], row_i["y"]),
                       code = i_code)
            for j, row_j in map.iterrows():
                if j <= i:
                    continue
                j_code = name_dict.lookup(row_j[label])
                pair_ij = (i_code, j_code)
                pair_ji = (j_code, i_code)
                G.add_edge(i_code, j_code, weight=pairs[pair_ij])
                G.add_edge(j_code, i_code, weight=pairs[pair_ji])


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


class graphDisplay:
    def __init__(self, graph, gdt = None, map_name_col = "SIG_KOR_NM",
                 name_dict = ntc.code_dict()):
        self.map = gdt
        self.graph = graph
        self.fig, self.ax = plt.subplots()
        self.cid = None
        self.name_dict = name_dict
        self.position = nx.get_node_attributes(self.graph, 'pos')
        self.code = nx.get_node_attributes(self.graph, 'code')
        self.map_name_col = map_name_col

        self.in_cmap  = mcolors.LinearSegmentedColormap.from_list(
            "inmap", ["#a1dab4", "#00a6b2", "#3366cc", "#5e3c99"])
        self.out_cmap = mcolors.LinearSegmentedColormap.from_list(
            "outmap", ["#ffb74d", "#ff7043", "#c62828"])

        self.in_norm  = mcolors.Normalize(vmin=None, vmax=None, clip=True)
        self.out_norm = mcolors.Normalize(vmin=None, vmax=None, clip=True)

        #For further speedup down the line it seems promising to add blitting if
        #it is needed to make this part of the code run faster. Depends
        #on further implementation. This is left as an opportunity. But I hope
        #it wont be needed 
        self._blit_enabled = False
        self._sel_marker = None
        self._line_cache = {}
        self._label_cache = {}
        self._legend = None

    #Just a placeholder, not used now
    def _blit_draw(self):
        if not self._blit_enabled:
            self.fig.canvas.draw_idle()
            return
        else:
            return

    def _draw_or_get_line(self, u, v, in_out: bool):
            key = (u, v, in_out)
            ln = self._line_cache.get(key)
            if ln is None:
                x1, y1 = self.graph.nodes[u]["pos"]; x2, y2 = self.graph.nodes[v]["pos"]
                angle = np.atan2(y2-y1, x2-x1)
                n_vec = (-np.sin(angle), np.cos(angle))
                mult = 120
                offset = mult if in_out == True else -mult
                color = "blue" if in_out == True else "red"
                ln, = self.ax.plot([x1+offset*n_vec[0], x2+offset*n_vec[0]], 
                                    [y1+offset*n_vec[1], y2+offset*n_vec[1]], 
                                    color=color, linestyle="dashed", zorder=1)
                self._line_cache[key] = ln
            return ln

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
            lb = self.ax.text(mid_x+mult*n_vec[0], mid_y+mult*n_vec[1], 
                              "", ha="center", va="bottom",
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
                norm.vmin = w; changed = True
            if w > norm.vmax:
                norm.vmax = w; changed = True
            if changed and norm.vmax == norm.vmin:
                norm.vmax = norm.vmin + 1e-12

    def _set_weight_color(self, ln, w: float, in_out: bool):
        if in_out == True:
            self._expand_norm(self.in_norm, w)
            ln.set_color(self.in_cmap(self.in_norm(w)))
        else:
            self._expand_norm(self.out_norm, w)
            ln.set_color(self.out_cmap(self.out_norm(w)))

    def _update_selection(self, node):
        self._sel_marker.center = self.graph.nodes[node]["pos"]
        edges = []

        for ln in self._line_cache.values():
            ln.set_alpha(0.00)
            ln.set_label(None)
            
        for idxn in self.graph.neighbors(node):
            weight_out = self.graph.get_edge_data(node, idxn)["weight"]
            weight_in = self.graph.get_edge_data(idxn, node)["weight"]
            edges.append((node, idxn, weight_in, weight_out))

            
        for lb in self._label_cache.values():
            lb.set_text("")

        for i, (u, v, w_in, w_out) in enumerate(edges, start=1):
            lb = self._set_or_get_label(u, v)
            lb.set_text(f"{i}")

            x1, y1 = self.graph.nodes[u]["pos"]; x2, y2 = self.graph.nodes[v]["pos"]
            ang = np.degrees(np.atan2(y2 - y1, x2 - x1))
            if ang > 90: ang -= 180
            if ang < -90: ang += 180
            lb.set_rotation(ang)
            
            ln1 = self._draw_or_get_line(u, v, True)
            ln1.set_linewidth(1)
            self._set_weight_color(ln1, w_in, True)
            ln1.set_alpha(1)
            ln1.set_label(f"{i}: flow in: {w_in}")

            ln2 = self._draw_or_get_line(u, v, False)
            ln2.set_linewidth(1)
            self._set_weight_color(ln2, w_out, False)
            ln2.set_alpha(1)
            ln2.set_label(f"{i}: flow out: {w_out}")


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
            
        self._legend = self.ax.legend(loc="upper left", frameon=False, 
                              fontsize=8, ncol=1)
        self.fig.canvas.draw_idle()  
    
    def interactive_graph(self):
        self.draw_graph()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click_graph)

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
                    color = "lightgreen",
                    zorder = 50) 
            self.ax.add_patch(circ)

            self.ax.text(self.position[node][0],self.position[node][1],
                    f"{self.code[node]}",
                    ha='center', va='center', zorder=100)
        
        for u, v in self.graph.edges():
            if u < v:
                continue
            x1, y1 = self.position[u][0], self.position[u][1]
            x2, y2 = self.position[v][0], self.position[v][1]
            self.ax.plot(
                    [x1, x2], [y1, y2], 
                    linestyle="dashed", color="gray", linewidth=0.02
                    )
        self.ax.autoscale_view() 
        self.ax.set_aspect('equal')
        self.ax.set_autoscale_on(False)
        self.ax.set_xticks([]) 
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)
