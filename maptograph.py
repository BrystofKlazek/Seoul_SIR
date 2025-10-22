import geopandas as gpd
import matplotlib.pyplot as plt
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
            G.add_node(i_code, name=row_i[label], pos=(row_i["x"], row_i["y"]))
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
    def __init__(self, graph, gdt = None, style = "double", map_name_col = "SIG_KOR_NM",
                 name_dict = ntc.code_dict()):
        self.style = style
        self.map = gdt
        self.graph = graph
        self.fig, self.ax = plt.subplots()
        self.cid = None
        self.name_dict = name_dict
        self.position = nx.get_node_attributes(self.graph, 'pos')
        self.code = nx.get_node_attributes(self.graph, 'code')

        self.map_name_col = map_name_col
        if self.map is None:
            self.sidx = None
        else:
            self.sidx = self.map.sindex
        self._last_collection = None
        self._last_circle = None
        self._lines = []
        self._annotations = []
        self._legend = None
        self._legend_proxies = []

    def _clear_lines(self):
        for ln in self._lines:
            ln.remove()
        self._lines.clear()

    def _clear_annotations(self):
        for annotation in self._annotations:
            annotation.remove()
        self._annotations.clear()

    def _clear_legend(self):
        if self._legend is not None:
            self._legend.remove()
            self._legend = None

    def _draw_lines(self, idx):
        for i_n, idxn in enumerate(self.graph.neighbors(idx)):
            weight_out = self.graph.get_edge_data(idx, idxn)["weight"]
            weight_in = self.graph.get_edge_data(idxn, idx)["weight"]
            
            if self.style == "double":

                line_normal_angle = np.arctan2(
                    self.position[idxn][1]-self.position[idx][1],
                    self.position[idxn][0]-self.position[idx][0]) + np.pi/2

                normal_vec = (np.cos(line_normal_angle), np.sin(line_normal_angle))

                ofst = 250
            
                ln1, = self.ax.plot([self.position[idx][0] + ofst*normal_vec[0], 
                                    self.position[idxn][0] + ofst*normal_vec[0]],
                                   [self.position[idx][1] + ofst*normal_vec[1], 
                                    self.position[idxn][1] + ofst*normal_vec[1]],
                                    linewidth = 0.5/300 * weight_out,
                                    color = "red", zorder = 55)
    
                ln2, = self.ax.plot([self.position[idxn][0] - ofst*normal_vec[0], 
                                    self.position[idx][0] - ofst*normal_vec[0]],
                                   [self.position[idxn][1] - ofst*normal_vec[1], 
                                    self.position[idx][1] - ofst*normal_vec[1]],
                                    linewidth = 0.5/300 * weight_in,
                                    color = "blue", zorder = 55,
                                    linestyle = "dashed")
 
                self._lines.append(ln1)
                self._lines.append(ln2)           

                #mid_x =  (self.position[idxn][0]+self.position[idx][0])/2
                #mid_y = (self.position[idxn][1]+self.position[idx][1])/2 
                #annot1 = self.ax.annotate(
                #        "label", xy=(mid_x, mid_y), xycoords = 'data',  
                #        ha='center', va='bottom',
                #        rotation = line_angle)

            else:
                ln, = self.ax.plot([self.position[idxn][0], self.position[idx][0]],
                                   [self.position[idxn][1], self.position[idx][1]],
                                    linewidth =(0.5/500 * np.abs(weight_in+weight_out)),
                                    color = "red" if weight_in > weight_out else
                                    "blue", zorder = 55, 
                                   label=f"{i_n+1}: flow in = {weight_in}, flow out = {weight_out}")
 
                self._lines.append(ln)

                mid_x =  (self.position[idxn][0]+self.position[idx][0])/2
                mid_y = (self.position[idxn][1]+self.position[idx][1])/2 

                line_angle = np.degrees(np.arctan2(
                    self.position[idxn][1]-self.position[idx][1],
                    self.position[idxn][0]-self.position[idx][0]))
                
                if line_angle > 90:
                    line_angle -= 180
                elif line_angle < -90:
                    line_angle += 180

                annot1 = self.ax.annotate(
                        f"{i_n+1}", xy=(mid_x, mid_y), xycoords = 'data',  
                        ha='center', va='bottom',
                        rotation = line_angle, rotation_mode="anchor")

                self._annotations.append(annot1)
        

            print(
            f"flow {idx} -> {idxn}: {weight_out}, flow {idxn} -> {idx}: {weight_in}"
            )

    def on_click_map(self, event): 
        if self.map is None:
            print("No map supplied")
            return
        if self._last_collection is not None:
            self._last_collection.remove()
            self._last_collection = None
        
        if self._last_circle is not None:
            self._last_circle.remove()
            self._last_circle = None

        self._clear_lines()
        self._clear_annotations()

        if event.inaxes is self.ax and event.xdata and event.ydata:
            
            idx = intercept_check(event.xdata, event.ydata, self.map)
            map_row_pd = self.map.iloc[idx]
            map_row = self.map.iloc[[idx]]
            dist_code = map_row_pd[self.map_name_col]
            SIG_code = self.name_dict.lookup(dist_code)
            map_row.plot(
                    ax=self.ax, 
                    facecolor="#f7f7f7", 
                    edgecolor="red", 
                    linewidth=0.8)
            
            self._last_collection = self.ax.collections[-1]

            self._last_circle = Circle(
                    self.position[SIG_code], 1000, 
                    zorder = 60, color = "orange") 

            self.ax.add_patch(self._last_circle)
            self._draw_lines(SIG_code)

        self._legend = self.ax.legend(loc="upper left", frameon=False, 
                                      fontsize=8, ncol=1)
        self.fig.canvas.draw_idle()  

    def on_click_graph(self, event):
        
        if self._last_circle is not None:
            self._last_circle.remove()
            self._last_circle = None

        self._clear_lines()
        self._clear_annotations()

        if event.inaxes is self.ax and event.xdata and event.ydata:
            nodes = list(self.graph)     
            point = (event.xdata, event.ydata)
            distances = np.array([])
            for node in nodes:
                difference = np.subtract(point, self.position[node])
                distance = difference[0]**2+difference[1]**2
                distances = np.append(distances, distance)

            idx = nodes[np.argmin(distances)]
            
            self._last_circle = Circle(
                    self.position[idx], 1000, 
                    zorder = 60, color = "orange") 
           
            annot_graph = self.ax.annotate(
                        f"{self.code[idx]}", xy=self.position[idx], xycoords = 'data',  
                        ha='center', va='center', zorder=100)

            self._annotations.append(annot_graph)
            
            self.ax.add_patch(self._last_circle)
            self._draw_lines(idx)
        
        self._legend = self.ax.legend(loc="upper left", frameon=False, 
                              fontsize=8, ncol=1)
        self.fig.canvas.draw_idle()  
    
    def interactive_graph(self, draw_map = True):
        if draw_map == True:
            if self.map is None:
                print("No map supplied")
                return
            else:
                self.map.plot(ax=self.ax, 
                            facecolor="#f7f7f7", 
                            edgecolor="black", 
                            linewidth=0.8)
                self.draw_graph()
                self.fig.canvas.mpl_connect('button_press_event', self.on_click_map) 

        else:
            self.draw_graph()
            self.fig.canvas.mpl_connect('button_press_event', self.on_click_graph)

    def draw_graph(self):
        for node in self.graph:
            circ = Circle(self.position[node], 1000, zorder = 50) 
            #Pak třeba zakoponovat lidnatost
            self.ax.add_patch(circ)

            self.ax.annotate(
                    f"{self.code[node]}", xy=self.position[node], xycoords = 'data',  
                    ha='center', va='center', zorder=100)
        
        for u, v in self.graph.edges():
            if u < v:
                continue
            x1, y1 = self.position[u][0], self.position[u][1]
            x2, y2 = self.position[v][0], self.position[v][1]
            self.ax.plot(
                    [x1, x2], [y1, y2], 
                    linestyle="dashed", color="gray", linewidth="0.3"
                    )
        self.ax.autoscale_view() 
        self.ax.set_aspect('equal')

            
