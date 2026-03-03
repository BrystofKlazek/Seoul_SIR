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

    def start_animation(
        self,
        var_names: Sequence[str],
        values_ts: np.ndarray,
        dt_snapshot: Optional[float],
        var: Union[int, str] = 0,
        interval: int = 60,
        cmap_name: str = "Plasma",
        fps: int = 20,
        edge_weight_fn=None,
        hourly_weights: Optional[Dict[Tuple[Any, Any], np.ndarray]] = None,
        edge_pairs: Optional[Sequence[Tuple[Any, Any]]] = None,
        output_html: Optional[str] = None,
        selected_u: Optional[Any] = None,
        selected_v: Optional[Any] = None,
        flow_top_k: Optional[int] = None,
        auto_show: bool = True,
        **_ignored_kwargs: Any,
    ) -> go.Figure:
        """
        Plotly animation:

        Left: animated choropleth (chosen `var`).
        Right-top: state S/E/I/R at selected source node (u).
        Right-bottom: selected directed edge flow:
            - outflow u -> v
            - inflow v -> u

        Selections:
          - source node u and target node v are controlled by HTML <select> elements
            when output_html is provided (static HTML with a small JS callback).
          - If output_html is None, the figure still renders, but u/v selection
            stays at the initial values.
        """
        max_var = float(values_ts[:, :, 2].max())

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
            raise ValueError("graphDisplay requires a GeoDataFrame map (gdt) so it can build a choropleth.")

        field_ts = vals[:, k, :]  # (T, N)

        # Defaults for u, v
        if selected_u is None:
            selected_u = self._codes_present[0]
        if selected_u not in self._node_to_idx:
            raise ValueError(f"selected_u {selected_u!r} not in graph")

        if selected_v is None:
            # Prefer a neighbour if possible, else fall back to itself
            outs = list(self.graph.successors(selected_u)) if selected_u in self.graph else []
            selected_v = outs[0] if outs else selected_u
        if selected_v not in self._node_to_idx:
            raise ValueError(f"selected_v {selected_v!r} not in graph")

        u_idx = self._node_to_idx[selected_u]

        # Choropleth values must be ordered the same as `locations`
        z0 = field_ts[0, self._codes_to_idx]

        # State time series for selected_u
        y_state = [vals[:, kk, u_idx] for kk in range(K)]

        # Flow series helper
        def _flow_series(u: Any, v: Any) -> np.ndarray:
            if edge_weight_fn is None:
                return np.zeros(T, dtype=float)
            return np.array([float(edge_weight_fn(f, u, v)) for f in range(T)], dtype=float)

        y_out = _flow_series(selected_u, selected_v)
        y_in = _flow_series(selected_v, selected_u)

        # y-limits
        y_state_min, y_state_max = -0.02, 1.02
        y_flow_max = float(self.max_flow) if self.max_flow is not None else float(max(np.max(y_out), np.max(y_in), 1e-12))
        y_flow_min, y_flow_max = 0.0, max(y_flow_max, 1e-12)

        # Map centroids (for hover) in the same order as codes_present
        gdf_idxed = self._gdf_plot.set_index("code")
        cent_lon = [float(gdf_idxed.loc[c, "lon"]) for c in self._codes_present]
        cent_lat = [float(gdf_idxed.loc[c, "lat"]) for c in self._codes_present]

        # ---- Layout: map (left, 2 rows) + state (top-right) + out/in (bottom-right split) ----
        fig = make_subplots(
            rows=2,
            cols=3,
            specs=[
                [{"type": "geo", "rowspan": 2}, {"type": "xy", "colspan": 2}, None],
                [None, {"type": "xy"}, {"type": "xy"}],
            ],
            column_widths=[0.52, 0.24, 0.24],
            row_heights=[0.58, 0.42],
            subplot_titles=(
                f"{var_names[k]} over districts",
                f"State at source node u = {selected_u}",
                "Selected outflow (u → v)",
                "Selected inflow (v → u)",
            ),
            horizontal_spacing=0.05,
            vertical_spacing=0.12,
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
                colorbar=dict(
                    title=str(var_names[k]),
                    # Keep the colorbar strictly within the left column.
                    x=1.15,
                    y=0.53,
                    zmax=float(max_var)
                    len=0.88,
                    thickness=14,
                    xanchor="right",
                    yanchor="middle",
                ),
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

        # State time series
        for kk, nm in enumerate(var_names):
            fig.add_trace(
                go.Scatter(x=t_vec, y=y_state[kk], mode="lines", name=str(nm)),
                row=1,
                col=2,
            )

        # Vertical time marker for state
        fig.add_trace(
            go.Scatter(
                x=[t_vec[0], t_vec[0]],
                y=[y_state_min, y_state_max],
                mode="lines",
                line=dict(dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=1,
            col=2,
        )

        # Outflow trace + vertical marker
        fig.add_trace(
            go.Scatter(x=t_vec, y=y_out, mode="lines", name=f"{selected_u} → {selected_v}", showlegend=False),
            row=2,
            col=2,
        )
        fig.add_trace(
            go.Scatter(
                x=[t_vec[0], t_vec[0]],
                y=[y_flow_min, y_flow_max],
                mode="lines",
                line=dict(dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=2,
            col=2,
        )

        # Inflow trace + vertical marker
        fig.add_trace(
            go.Scatter(x=t_vec, y=y_in, mode="lines", name=f"{selected_v} → {selected_u}", showlegend=False),
            row=2,
            col=3,
        )
        fig.add_trace(
            go.Scatter(
                x=[t_vec[0], t_vec[0]],
                y=[y_flow_min, y_flow_max],
                mode="lines",
                line=dict(dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=2,
            col=3,
        )

        # Trace indices for frame updates
        choropleth_idx = 0
        state_marker_idx = 2 + K
        out_marker_idx = 4 + K
        in_marker_idx = 6 + K

        # Frames: update choropleth z + vertical markers x (state/out/in)
        frames = []
        for f in range(T):
            frames.append(
                go.Frame(
                    name=str(f),
                    data=[
                        go.Choropleth(z=field_ts[f, self._codes_to_idx], locations=self._codes_present),
                        go.Scatter(x=[t_vec[f], t_vec[f]]),
                        go.Scatter(x=[t_vec[f], t_vec[f]]),
                        go.Scatter(x=[t_vec[f], t_vec[f]]),
                    ],
                    traces=[choropleth_idx, state_marker_idx, out_marker_idx, in_marker_idx],
                )
            )
        fig.frames = frames

        slider_steps = [
            dict(
                method="animate",
                args=[
                    [str(f)],
                    {"mode": "immediate", "frame": {"duration": interval, "redraw": True}, "transition": {"duration": 0}},
                ],
                label=str(f),
            )
            for f in range(T)
        ]

        # Make axes readable and non-overlapping
        fig.update_xaxes(title_text=t_label, row=1, col=2)
        fig.update_xaxes(title_text=t_label, row=2, col=2)
        fig.update_xaxes(title_text=t_label, row=2, col=3)

        fig.update_yaxes(title_text="fraction", range=[y_state_min, y_state_max], row=1, col=2)
        fig.update_yaxes(title_text="weight", range=[y_flow_min, y_flow_max], row=2, col=2)
        fig.update_yaxes(title_text="weight", range=[y_flow_min, y_flow_max], row=2, col=3)

        fig.update_geos(fitbounds="locations", visible=False)

        # Put animation controls below the plots (no overlap)
        fig.update_layout(
            height=720,
            width=1280,
            margin=dict(l=10, r=240, t=70, b=120),
            hovermode="x unified",
            legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top"),
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    x=0.05,
                    y=-0.14,
                    xanchor="left",
                    yanchor="top",
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
                    x=0.05,
                    y=-0.07,
                    len=0.65,
                    pad=dict(t=0, b=0),
                    currentvalue=dict(prefix=f"{t_label}: "),
                    steps=slider_steps,
                )
            ],
        )

        self.fig = fig

        # Optional HTML output: plain Plotly figure (no custom HTML/JS injected)
        if output_html is None and self.save_name is not None:
            output_html = self.save_name

        if output_html is not None:
            if not output_html.lower().endswith(".html"):
                output_html = output_html + ".html"

            fig.write_html(output_html, auto_open=auto_show, include_plotlyjs="cdn")
            print(f"Saved Plotly animation to {output_html}")

        if auto_show and output_html is None:
            fig.show()

        return fig
