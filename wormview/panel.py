"""The side panel: what a neuron connects to, and what it expresses."""

from dash import html

from . import genes as gene_meta
from .theme import CHIP_BG, CHIP_INK, CHIP_WATCH, GROUP_COLOR, INK, INK_DIM


# --- side panel ----------------------------------------------------------



def _chip(text, color=None, bold=False, light_text=False):
    return html.Span(text, style={
        "display": "inline-block", "padding": "3px 9px", "marginRight": "6px",
        "marginBottom": "6px", "borderRadius": "11px", "fontSize": "11px",
        "background": color or CHIP_BG,
        "color": "#ffffff" if light_text else CHIP_INK,
        "border": "1px solid rgba(0,0,0,0.06)",
        "fontWeight": "600" if bold else "400"})


def panel_for(node, node_index):
    if node is None:
        return html.Div([
            html.H3("Click a coloured dot", style={"marginTop": 0, "color": INK}),
            html.P("Each dot is a neuron type, at its real position in the animal "
                   "— head on the left, tail on the right. Lines are real "
                   "connections. Drag to rotate, scroll to zoom.",
                   style={"color": INK_DIM, "lineHeight": "1.65"}),
            html.P("Click directly on a dot (not a line, not empty space) and its "
                   "connections light up gold and its details appear here. If the "
                   "dots are too crowded to hit, raise the spread slider or use "
                   "\u201cjump to neuron\u201d.",
                   style={"color": INK_DIM, "lineHeight": "1.65"}),
            html.P("Switch \u201ccolour by\u201d to one gene to paint the whole "
                   "nervous system by that gene\u2019s expression.",
                   style={"color": INK_DIM, "lineHeight": "1.65"}),
        ])

    partners = node["partners"]
    out_ = [p for p in partners if p["dir"] == "out"]
    in_ = [p for p in partners if p["dir"] == "in"]

    def partner_rows(items, arrow):
        return [html.Div([
            html.Span(f"{arrow} ", style={"color": INK_DIM}),
            html.Span(p["name"], style={"color": INK, "fontWeight": "600"}),
            html.Span(f"  {p['w']} synapses",
                      style={"color": INK_DIM, "fontSize": "11px"}),
        ], style={"padding": "1px 0"}) for p in items[:12]]

    watch = node["watch_genes"]
    return html.Div([
        html.Div([
            html.H2(node["name"], style={"margin": "0 6px 0 0", "color": INK,
                                         "display": "inline-block"}),
            _chip(node["group"], GROUP_COLOR[node["group"]], bold=True,
                  light_text=True),
        ]),
        html.Div([
            _chip(f"{node['degree']} connections"),
            _chip(f"{node['n_genes_on']:,} genes on"),
            _chip(f"{node['n_cells']} cell"
                  f"{'s' if node['n_cells'] > 1 else ''} in this class"),
        ], style={"marginTop": "8px"}),

        html.H4(f"Sends to  ({len(out_)})",
                style={"color": INK, "marginBottom": "4px"}),
        html.Div(partner_rows(out_, "→") or
                 [html.Div("none", style={"color": INK_DIM})]),

        html.H4(f"Receives from  ({len(in_)})",
                style={"color": INK, "marginBottom": "4px"}),
        html.Div(partner_rows(in_, "←") or
                 [html.Div("none", style={"color": INK_DIM})]),

        html.H4("Most strongly expressed genes",
                style={"color": INK, "marginBottom": "6px"}),
        html.Div([_chip(f"{g['g']}  {g['tpm']:,}",
                        CHIP_WATCH if g["watch"] else None, bold=g["watch"])
                  for g in node["top_genes"]]),

        html.H4("Watch-list genes switched on here",
                style={"color": INK, "marginBottom": "6px"}),
        html.Div([_chip(f"{g['g']}  {g['tpm']:,}", CHIP_WATCH, bold=True)
                  for g in watch] or
                 [html.Span("none of the watch-list detected here",
                            style={"color": INK_DIM, "fontSize": "12px"})]),

        *( [html.Div([
              html.B("Read the watch-list genes with care. "),
              html.Span("; ".join(
                  f"{g['g']}: {gene_meta.autonomy_caveat(g['g'])}"
                  for g in watch if gene_meta.autonomy_caveat(g["g"]))),
              html.Span(". CeNGEN sequenced neurons only, so a low signal for "
                        "these does not mean the gene is not involved."),
           ], style={"color": INK_DIM, "fontSize": "10.5px", "marginTop": "12px",
                     "lineHeight": "1.55"})]
           if any(gene_meta.autonomy_caveat(g["g"]) for g in watch) else [] ),

        html.P("Expression from CeNGEN (L4 animals); wiring from Cook et al. "
               "2019. A gene being switched on does not mean it acts at these "
               "connections.",
               style={"color": INK_DIM, "fontSize": "10.5px",
                      "marginTop": "16px", "lineHeight": "1.5"}),
    ])
