import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Circle
import networkx as nx
from shapely import set_precision
from shapely.geometry import Point

import numpy as np


def maptograph(map, MIN_BORDER = 25, label="SIG_KOR_NM", mode="all"):
# 2) Clean + snap to grid (no growth like buffer does)
#    - buffer(0) fixes invalid rings
#    - set_precision(..., GRID) snaps coords to GRID-meter lattice
    GRID = 0.2   # 20 cm grid; adjust 0.1–0.5 if needed
    map["geometry"] = map.geometry.buffer(0).apply(
            lambda geom: set_precision(geom, GRID)
            )

    sindex = map.sindex
    G = nx.DiGraph()

    if mode=="neighbours":
        for i, row in map.iterrows():
            G.add_node(i, name=row[label], pos=(row["x"], row["y"]))
            for j in sindex.intersection(row.geometry.bounds):
                if j <= i:
                    continue
                other = map.geometry.iloc[j]
                inter = row.geometry.boundary.intersection(other.boundary)
                L = float(inter.length)
                if L >= MIN_BORDER:
                    # add both directions; weights can diverge later in your model
                    G.add_edge(i, j, weight=L, border_m=L)
                    G.add_edge(j, i, weight=L, border_m=L)
    
    if mode=="all":
        for i, row in map.iterrows():
            G.add_node(i, name=row[label], pos=(row["x"], row["y"]))
            for j, row in map.iterrows():
                if j <= i:
                    continue
                G.add_edge(i, j, weight=25)
                G.add_edge(j, i, weight=25)


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
    def __init__(self, graph, gdt = None):
        self.map = gdt
        self.graph = graph
        self.fig, self.ax = plt.subplots()
        self.cid = None
        if self.map is None:
            self.sidx = None
        else:
            self.sidx = self.map.sindex
        self._last_collection = None
        self._last_circle = None
        self._lines = []
        self.position = nx.get_node_attributes(self.graph, 'pos')

    def _clear_lines(self):
        for ln in self._lines:
            ln.remove()
        self._lines.clear()

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

        if event.inaxes is self.ax and event.xdata and event.ydata:
            
            idx = intercept_check(event.xdata, event.ydata, self.map)
           
            self.map.iloc[[idx]].plot(
                    ax=self.ax, 
                    facecolor="#f7f7f7", 
                    edgecolor="red", 
                    linewidth=0.8)
            
            self._last_collection = self.ax.collections[-1]

            self._last_circle = Circle(
                    self.position[idx], 1000, 
                    zorder = 60, color = "orange") 
            #Pak třeba zakoponovat lidnatost

            self.ax.add_patch(self._last_circle)
            for node in self.graph.neighbors(idx):
                ln, = self.ax.plot([self.position[idx][0], self.position[node][0]],
                                   [self.position[idx][1], self.position[node][1]],
                                   linewidth = 0.5, color = "red", zorder = 55)
                self._lines.append(ln)

        self.fig.canvas.draw_idle()  

    def on_click_graph(self, event):
        
        if self._last_circle is not None:
            self._last_circle.remove()
            self._last_circle = None

        self._clear_lines()

        if event.inaxes is self.ax and event.xdata and event.ydata:
            
            point = (event.xdata, event.ydata)
            distances = np.array([])
            for node in self.graph:
                difference = np.subtract(point, self.position[node])
                distance = difference[0]**2+difference[1]**2
                distances = np.append(distances, distance)

            idx = np.argmin(distances)
            
            self._last_circle = Circle(
                    self.position[idx], 1000, 
                    zorder = 60, color = "orange") 
            #Pak třeba zakoponovat lidnatost

            self.ax.add_patch(self._last_circle)
            for node in self.graph.neighbors(idx):
                ln, = self.ax.plot([self.position[idx][0], self.position[node][0]],
                                   [self.position[idx][1], self.position[node][1]],
                                   linewidth = 0.5, color = "red", zorder = 55)
                self._lines.append(ln)

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
        for u, v in self.graph.edges():
            if u < v:
                continue
            x1, y1 = self.position[u][0], self.position[u][1]
            x2, y2 = self.position[v][0], self.position[v][1]
            self.ax.plot(
                    [x1, x2], [y1, y2], 
                    linestyle="dashed", color="gray", linewidth="0.1"
                    )
        self.ax.autoscale_view() 
        self.ax.set_aspect('equal')

            
