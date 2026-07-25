#!/usr/bin/env python3
"""C. elegans nervous system explorer -- interactive 3D.

Every neuron at its real position inside the animal, with its real neurite
morphology. Click a cell body to see what it wires to and which genes it switches
on, or colour the whole nervous system by any one of ~13,700 genes.

    pip install -r requirements.txt
    python fetch_data.py        # once: downloads expression, wiring and geometry
    python app.py               # then open the URL it prints

Drag to rotate, scroll to zoom, click a cell body to inspect it.
"""

import sys
import webbrowser
from pathlib import Path

from dash import Dash, Input, Output, State, ctx, dcc, html

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from wormview import anatomy, connectome as cn, figure, panel   # noqa: E402
from wormview.data import MissingData, load_matrix              # noqa: E402
from wormview.theme import (BG, FONT, GROUP_COLOR, INK, INK_DIM,  # noqa: E402
                            LINE, PANEL)

PORT = 8050


def load_everything():
    """Expression, wiring and geometry, assembled into one payload."""
    tpm = cn.merge_ambiguous_classes(load_matrix(2)[0])
    edges, _ = cn.load_edges(list(tpm.columns), "chemical")
    payload = anatomy.build_payload(tpm, edges)

    # Real positions from the Virtual Worm model, in microns. Axes are remapped so
    # the body runs left-to-right on screen: x = anterior-posterior,
    # y = left-right, z = dorsal-ventral.
    positions, counts, missing = anatomy.load_positions(list(tpm.columns))
    for node in payload["nodes"]:
        pos = positions.get(node["name"])
        node["real"] = pos is not None
        node["n_cells"] = counts.get(node["name"], 1)
        if pos:
            node["x"], node["y"], node["z"] = pos

    payload["true_positions"] = positions
    payload["morphology"] = anatomy.load_morphology(list(tpm.columns))
    # the wall is fitted to the neurites, not just the cell bodies: nerve cords run
    # right along the body wall and reach much further out than the somata
    payload["surface"] = anatomy.worm_surface(positions, payload["morphology"])
    payload["missing_positions"] = missing
    return tpm, payload


CONTROL_BAR = {
    "display": "flex", "alignItems": "flex-end", "gap": "26px",
    "flexWrap": "wrap", "padding": "12px 20px",
    "background": "#f7f7f5", "borderTop": f"1px solid {LINE}",
    "borderBottom": f"1px solid {LINE}",
}


def _field(label_text, control, width=None):
    """A labelled control, so every label sits on the same baseline."""
    style = {"display": "flex", "flexDirection": "column", "gap": "4px"}
    if width:
        style["width"] = width
    return html.Div([
        html.Label(label_text, style={"color": INK_DIM, "fontSize": "11px",
                                      "letterSpacing": "0.02em"}),
        control,
    ], style=style)


GROUPS = ("sensory", "interneuron", "motor", "other")


def _colour_key():
    """Colour key that is also the on/off control.

    Plotly's own legend was click-to-toggle, and removing it took that away. This
    puts the toggles back as real checkboxes: same function, but in a fixed place
    instead of floating over the animal.
    """
    swatch = {"display": "inline-block", "width": "16px", "height": "3px",
              "borderRadius": "2px", "marginRight": "6px",
              "verticalAlign": "middle"}
    options = [
        {"value": g,
         "label": html.Span([html.Span(style={**swatch,
                                              "background": GROUP_COLOR[g]}), g],
                            style={"fontSize": "11.5px", "color": INK,
                                   "marginRight": "13px"})}
        for g in GROUPS
    ]
    return _field(
        "show neuron groups  (click to toggle)",
        dcc.Checklist(id="groups", options=options, value=list(GROUPS),
                      inline=True, style={"paddingBottom": "3px"}))


def build_layout(names, gene_options):
    return html.Div([
        # ---- header ------------------------------------------------------
        html.Div([
            html.H1("C. elegans nervous system explorer",
                    style={"margin": 0, "fontSize": "19px", "fontWeight": 600,
                           "color": INK, "letterSpacing": "-0.01em"}),
            html.Div("302 neurons at their real positions, with real neurite "
                     "morphology  ·  wiring from Cook et al. 2019  ·  expression "
                     "from CeNGEN",
                     style={"color": INK_DIM, "fontSize": "12px",
                            "marginTop": "3px"}),
        ], style={"padding": "16px 20px 12px"}),

        # ---- control bar -------------------------------------------------
        html.Div([
            _field("colour by", dcc.RadioItems(
                id="mode", value="type", inline=True,
                options=[{"label": " neuron type", "value": "type"},
                         {"label": " one gene", "value": "gene"}],
                labelStyle={"marginRight": "12px"},
                style={"color": INK, "fontSize": "12px",
                       "paddingBottom": "3px"})),

            _field("gene", dcc.Dropdown(
                id="gene", options=gene_options, value="fmi-1", clearable=False,
                style={"fontSize": "13px"}), width="180px"),

            _field("jump to neuron", dcc.Dropdown(
                id="pick", options=names, value=None, placeholder="e.g. AWA",
                style={"fontSize": "13px"}), width="160px"),

            _field("show", dcc.Checklist(
                id="layers", inline=True,
                options=[{"label": " body wall", "value": "body"},
                         {"label": " labels", "value": "labels"},
                         {"label": " synaptic links", "value": "syn"}],
                value=["body"], labelStyle={"marginRight": "12px"},
                style={"color": INK, "fontSize": "12px",
                       "paddingBottom": "3px"})),

            _field("spread crowded cell bodies", html.Div(
                dcc.Slider(id="spread", min=0, max=3, step=0.5, value=0,
                           marks={0: {"label": "true"}, 1.5: {"label": "1.5x"},
                                  3: {"label": "3x"}},
                           tooltip=None),
                # the slider draws its end labels outside its own box, so without
                # this padding the first and last are clipped
                style={"width": "190px", "padding": "0 14px"})),

            _colour_key(),
        ], style=CONTROL_BAR),

        dcc.Store(id="sel", data=None),
        # trace index -> neuron names, rebuilt with every figure.
        # Needed because plotly >= 6.0 drops customdata from clickData.
        dcc.Store(id="clickmap", data={}),

        # ---- plot + panel ------------------------------------------------
        html.Div([
            html.Div(
                dcc.Graph(
                    id="net", style={"height": "74vh", "width": "100%"},
                    # responsive=True makes Plotly resize its canvas to the
                    # container. Without it the canvas keeps its initial size, the
                    # container crops it, and every pointer coordinate is offset --
                    # which silently breaks rotate, zoom AND click at once.
                    responsive=True,
                    config={"displaylogo": False, "scrollZoom": True,
                            "responsive": True, "doubleClick": "reset",
                            "modeBarButtonsToRemove": ["toImage"]}),
                # no overflow:hidden here -- it would crop the canvas
                style={"minWidth": 0}),
            html.Div(id="panel", style={
                "background": PANEL, "borderLeft": f"1px solid {LINE}",
                "padding": "20px 22px", "overflowY": "auto", "height": "74vh",
                "boxSizing": "border-box", "fontSize": "13px"}),
        ], style={"display": "grid", "gridTemplateColumns": "minmax(0,1fr) 330px",
                  "alignItems": "start"}),
    ], style={"background": BG, "minHeight": "100vh", "fontFamily": FONT,
              "WebkitFontSmoothing": "antialiased"})


def main():
    try:
        print("loading expression, wiring and geometry ...")
        tpm, payload = load_everything()
    except (MissingData, cn.MissingConnectome) as exc:
        print(f"error: {exc}\n\nRun:  python fetch_data.py", file=sys.stderr)
        return 1

    node_index = {n["name"]: n for n in payload["nodes"]}
    names = sorted(node_index)
    watchlist = payload["watchlist"]
    # watch-list genes first so they are easy to find, then everything else
    gene_options = watchlist + [g for g in tpm.index if g not in set(watchlist)]

    print(f"  {len(names)} neuron types, {len(payload['edges'])} connections, "
          f"{payload['n_genes_total']:,} genes, "
          f"{sum(len(v) for v in payload['morphology'].values()):,} neurite segments")

    app = Dash(__name__)
    app.title = "C. elegans nervous system explorer"
    app.layout = build_layout(names, gene_options)

    @app.callback(Output("sel", "data"),
                  Input("net", "clickData"), Input("pick", "value"),
                  State("sel", "data"), State("clickmap", "data"))
    def choose(click, picked, current, clickmap):
        """Resolve what was clicked.

        Since plotly 6.0, clickData no longer carries `customdata`, so the neuron
        is looked up by curveNumber + pointNumber against the map the figure
        emitted. customdata is still tried first, for older plotly versions.

        Only a click that lands on a cell body changes the selection: clicks on a
        neurite, the body wall or empty space leave it alone, so a near miss does
        not wipe the panel.
        """
        if ctx.triggered_id == "pick":
            return picked or current
        if not (click and click.get("points")):
            return current

        point = click["points"][0]

        name = point.get("customdata")
        if isinstance(name, list) and name:
            name = name[0]
        if isinstance(name, str) and name in node_index:
            return name

        curve, index = point.get("curveNumber"), point.get("pointNumber")
        if curve is not None and index is not None:
            # a Store round-trips through JSON, so keys arrive as strings
            trace = (clickmap or {}).get(str(curve)) or (clickmap or {}).get(curve)
            if trace and 0 <= index < len(trace) and trace[index] in node_index:
                return trace[index]
        return current

    @app.callback(Output("net", "figure"), Output("panel", "children"),
                  Output("clickmap", "data"),
                  Input("sel", "data"), Input("mode", "value"),
                  Input("gene", "value"), Input("layers", "value"),
                  Input("spread", "value"), Input("groups", "value"))
    def update(selected, mode, gene, layers, spread, groups):
        fig, click_map = figure.make_figure(
            payload, node_index, selected=selected,
            gene=gene if mode == "gene" else None, tpm=tpm,
            show_labels="labels" in (layers or []),
            show_synapses="syn" in (layers or []),
            show_body="body" in (layers or []),
            groups=groups if groups is not None else list(GROUPS),
            spread=float(spread or 0))
        return (fig, panel.panel_for(node_index.get(selected), node_index),
                click_map)

    url = f"http://127.0.0.1:{PORT}"
    print(f"\n  open {url}\n  (ctrl-c to stop)\n")
    try:
        webbrowser.open(url)
    except Exception:                                  # noqa: BLE001
        pass
    app.run(debug=False, port=PORT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
