import json
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import geopandas as gpd
import networkx as nx
import numpy as np
from shapely import set_precision
from shapely.geometry import Point

import nametocode as ntc

# Plotly (interactive, browser-based)
import plotly.graph_objects as go
from plotly.subplots import make_subplots


FRAME_SKIP = 1


def maptograph(
    gdf: gpd.GeoDataFrame,
    MIN_BORDER: float = 25,
    label: str = "SIG_KOR_NM",
    pairs: Optional[Dict[Tuple[int, int], float]] = None,
    mode: str = "all",
    weight_fn=None,
    name_dict: Optional[ntc.code_dict] = None,
) -> nx.DiGraph:
    """Build a directed graph from a GeoDataFrame.

    This is unchanged logically from the matplotlib version.

    Node attributes:
      - pos: (x,y) in the GeoDataFrame CRS
      - index: the GeoDataFrame row index at creation time
    """

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
            for v in nodes[a + 1 :]:
                w_uv = get_w(u, v)
                w_vu = get_w(v, u)
                if w_uv:
                    G.add_edge(u, v, weight=w_uv)
                if w_vu:
                    G.add_edge(v, u, weight=w_vu)

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
                    w_uv = get_w(u, v)
                    w_vu = get_w(v, u)
                    if w_uv:
                        G.add_edge(u, v, weight=w_uv)
                    if w_vu:
                        G.add_edge(v, u, weight=w_vu)

    elif mode == "from_file":
        if pairs:
            for (u, v), w in pairs.items():
                if u in G and v in G and w:
                    G.add_edge(u, v, weight=float(w))
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    return G


def intercept_check(x_coord, y_coord, polygons: gpd.GeoDataFrame):
    pt = Point(x_coord, y_coord)
    cand = polygons.sindex.query(pt)
    deep_check = polygons.iloc[cand]

    mask = deep_check.covers(pt)
    containing_polygon = mask[mask]
    if containing_polygon.empty is False:
        return containing_polygon.index[0]
    return None


class graphDisplay:
    def __init__(
        self,
        graph: nx.DiGraph,
        gdt: Optional[gpd.GeoDataFrame] = None,
        max_flow: Optional[float] = None,
        map_name_col: str = "SIG_KOR_NM",
        gif_name: Optional[str] = None,
    ):
        self.map = gdt
        self.graph = graph
        self.map_name_col = map_name_col
        self.save_name = gif_name

        self.max_flow = float(max_flow) if (max_flow is not None and max_flow > 0.0) else 1.0

        # Node ordering consistent with your solver pipeline
        self._node_order: List[Any] = sorted(
            list(self.graph.nodes()),
            key=lambda n: self.graph.nodes[n].get("index", 0),
        )
        self._node_to_idx: Dict[Any, int] = {n: i for i, n in enumerate(self._node_order)}

        # Prepare map geojson in EPSG:4326 so Plotly geo traces work.
        self._gdf_plot: Optional[gpd.GeoDataFrame] = None
        self._geojson: Optional[Dict[str, Any]] = None
        self._codes_present: List[Any] = []
        self._codes_to_idx: List[int] = []

        if self.map is not None:
            gdf = self.map.copy()

            # Ensure a stable code column matching graph node ids
            if "code" not in gdf.columns:
                name_dict = ntc.code_dict()
                gdf["code"] = gdf[self.map_name_col].map(name_dict.lookup)

            gdf_plot = gdf.to_crs(4326)
            anchors = gdf_plot.representative_point()
            gdf_plot["lon"] = anchors.x
            gdf_plot["lat"] = anchors.y

            self._gdf_plot = gdf_plot
            self._geojson = json.loads(gdf_plot.to_json())

            gdf_codes = set(gdf_plot["code"].tolist())
            self._codes_present = [c for c in self._node_order if c in gdf_codes]
            self._codes_to_idx = [self._node_to_idx[c] for c in self._codes_present]

        self.fig: Optional[go.Figure] = None

    def interactive_graph(self):
        # retained for backwards compatibility
        return

    def start_animation(
        self,
        var_names: Sequence[str],
        values_ts: np.ndarray,
        dt_snapshot: Optional[float],
        var: Union[int, str] = 0,
        interval: int = 60,
        cmap_name: str = "Plasma",
        fps: int = 10,
        edge_weight_fn=None,
        output_html: Optional[str] = None,
        selected_node: Optional[Any] = None,
        auto_show: bool = True,
    ) -> go.Figure:
        vals = np.asarray(values_ts)
        if vals.ndim == 2:
            vals = vals[:, np.newaxis, :]
        if vals.ndim != 3:
            raise ValueError("values_ts must have shape (T, K, N) or (T, N)")

        T, K, N = vals.shape

        if isinstance(var, str):
            if var not in var_names:
                raise ValueError(f"Unknown variable name {var!r}")
            k = list(var_names).index(var)
        else:
            k = int(var)
            if not (0 <= k < K):
                raise ValueError(f"var index {k} out of range [0, {K})")

        # Time vector
        if dt_snapshot is None or dt_snapshot <= 0:
            t_vec = np.arange(T, dtype=float)
            t_label = "frame"
        else:
            t_vec = np.arange(T, dtype=float) * float(dt_snapshot)
            t_label = "time"

        if self._gdf_plot is None or self._geojson is None or not self._codes_present:
            raise ValueError(
                "graphDisplay requires a GeoDataFrame map (gdt) so it can build a choropleth."
            )

        field_ts = vals[:, k, :]  # (T, N)

        # Default selected node for time series
        if selected_node is None:
            selected_node = self._codes_present[0]
        if selected_node not in self._node_to_idx:
            raise ValueError(f"selected_node {selected_node!r} not in graph")

        sel_idx = self._node_to_idx[selected_node]

        # Choropleth values must be ordered the same as `locations`
        z0 = field_ts[0, self._codes_to_idx]

        # Build time series for selected node
        y_series = [vals[:, kk, sel_idx] for kk in range(K)]

        # Global min/max for the vertical time marker
        y_min = float(np.min(vals))
        y_max = float(np.max(vals))
        if y_max <= y_min:
            y_max = y_min + 1e-12

        # Map centroids (for hover) in the same order as codes_present
        gdf_idxed = self._gdf_plot.set_index("code")
        cent_lon = [float(gdf_idxed.loc[c, "lon"]) for c in self._codes_present]
        cent_lat = [float(gdf_idxed.loc[c, "lat"]) for c in self._codes_present]

        fig = make_subplots(
            rows=1,
            cols=2,
            specs=[[{"type": "geo"}, {"type": "xy"}]],
            column_widths=[0.62, 0.38],
            subplot_titles=(
                f"{var_names[k]} over districts",
                f"Time series (selected: {selected_node})",
            ),
            horizontal_spacing=0.03,
        )

        # Choropleth on the left
        fig.add_trace(
            go.Choropleth(
                geojson=self._geojson,
                featureidkey="properties.code",
                locations=self._codes_present,
                z=z0,
                colorscale=cmap_name,
                marker_line_width=0.3,
                colorbar=dict(title=str(var_names[k])),
                hovertemplate="code=%{location}<br>value=%{z:.4g}<extra></extra>",
            ),
            row=1,
            col=1,
        )

        # Centroid markers for hover
        fig.add_trace(
            go.Scattergeo(
                lon=cent_lon,
                lat=cent_lat,
                mode="markers",
                marker=dict(size=6, opacity=0.75),
                text=[str(c) for c in self._codes_present],
                hovertemplate="code=%{text}<extra></extra>",
                showlegend=False,
            ),
            row=1,
            col=1,
        )

        # Time series on the right
        for kk, nm in enumerate(var_names):
            fig.add_trace(
                go.Scatter(x=t_vec, y=y_series[kk], mode="lines", name=str(nm)),
                row=1,
                col=2,
            )

        # Vertical time marker (updated via frames)
        fig.add_trace(
            go.Scatter(
                x=[t_vec[0], t_vec[0]],
                y=[y_min, y_max],
                mode="lines",
                line=dict(dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=1,
            col=2,
        )

        # Frames: update choropleth z + vertical marker x
        frames = []
        for f in range(T):
            frames.append(
                go.Frame(
                    name=str(f),
                    data=[
                        go.Choropleth(z=field_ts[f, self._codes_to_idx], locations=self._codes_present),
                        go.Scatter(x=[t_vec[f], t_vec[f]]),
                    ],
                    traces=[0, 6],
                )
            )
        fig.frames = frames

        slider_steps = [
            dict(
                method="animate",
                args=[[str(f)], {"mode": "immediate", "frame": {"duration": interval, "redraw": True}, "transition": {"duration": 0}}],
                label=str(f),
            )
            for f in range(T)
        ]

        fig.update_layout(
            height=650,
            width=1100,
            margin=dict(l=10, r=10, t=60, b=10),
            hovermode="closest",
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    x=0.12,
                    y=1.08,
                    buttons=[
                        dict(
                            label="Play",
                            method="animate",
                            args=[None, {"fromcurrent": True, "frame": {"duration": interval, "redraw": True}, "transition": {"duration": 0}}],
                        ),
                        dict(
                            label="Pause",
                            method="animate",
                            args=[[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False}, "transition": {"duration": 0}}],
                        ),
                    ],
                )
            ],
            sliders=[
                dict(
                    active=0,
                    x=0.10,
                    y=1.02,
                    len=0.46,
                    pad=dict(t=10, b=10),
                    currentvalue=dict(prefix=f"{t_label}: "),
                    steps=slider_steps,
                )
            ],
        )

        fig.update_geos(fitbounds="locations", visible=False)
        fig.update_xaxes(title_text=t_label, row=1, col=2)
        fig.update_yaxes(title_text="value", row=1, col=2)

        # Dropdown to change selected node time series (traces 2..5)
        node_buttons = []
        for node in self._codes_present:
            idx = self._node_to_idx[node]
            y_new = [vals[:, kk, idx] for kk in range(K)]
            node_buttons.append(
                dict(
                    label=str(node),
                    method="restyle",
                    args=[{"y": y_new}, [2, 3, 4, 5]],
                )
            )

        if node_buttons:
            fig.update_layout(
                updatemenus=list(fig.layout.updatemenus)
                + [
                    dict(
                        type="dropdown",
                        x=0.70,
                        y=1.08,
                        showactive=True,
                        buttons=node_buttons,
                    )
                ]
            )

        self.fig = fig

        if output_html is None and self.save_name is not None:
            output_html = self.save_name

        if output_html is not None:
            if not output_html.lower().endswith(".html"):
                output_html = output_html + ".html"
            fig.write_html(output_html, include_plotlyjs="cdn")
            print(f"Saved Plotly animation to {output_html}")

        if auto_show:
            fig.show()

        return fig
