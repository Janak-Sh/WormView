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
        row = {"color": INK_DIM, "fontSize": "12px", "lineHeight": "1.5",
               "marginBottom": "9px"}
        return html.Div([
            html.Div("Nothing selected", style={
                "fontSize": "14px", "fontWeight": 600, "color": INK,
                "marginBottom": "12px"}),
            html.Div([html.B("Click a dot"),
                      " to see what it connects to and what it expresses."],
                     style=row),
            html.Div([html.B("Drag"), " to rotate, ", html.B("scroll"),
                      " to zoom, ", html.B("double-click"), " to reset."],
                     style=row),
            html.Div("Dots are cell bodies at their real positions. Lines are "
                     "their real neurites.", style=row),
            html.Div("Hard to hit a dot? Turn off “body wall” — the cell "
                     "bodies sit inside it, so it can absorb the click. Or use "
                     "“jump to neuron”, or raise the spread slider.", style=row),
        ])

    partners = node["partners"]
    out_ = [p for p in partners if p["dir"] == "out"]
    in_ = [p for p in partners if p["dir"] == "in"]

    def partner_rows(items, arrow):
        return [html.Div([
            html.Span(f"{arrow} ", style={"color": INK_DIM}),
            html.Span(p["name"], style={"color": INK, "fontWeight": "600"}),
            html.Span(f"  {p['w']}",
                      style={"color": INK_DIM, "fontSize": "11px"}),
        ], style={"padding": "1px 0", "fontSize": "12.5px"})
            for p in items[:10]]

    watch = node["watch_genes"]
    return html.Div([
        html.Div([
            html.Span(node["name"], style={
                "fontSize": "22px", "fontWeight": 600, "color": INK,
                "marginRight": "8px", "letterSpacing": "-0.01em"}),
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
                style={"color": INK, "marginBottom": "4px", "fontSize": "12px",
                       "marginTop": "18px", "fontWeight": 600,
                       "textTransform": "uppercase",
                       "letterSpacing": "0.04em"}),
        html.Div(partner_rows(out_, "→") or
                 [html.Div("none", style={"color": INK_DIM})]),

        html.H4(f"Receives from  ({len(in_)})",
                style={"color": INK, "marginBottom": "4px", "fontSize": "12px",
                       "marginTop": "18px", "fontWeight": 600,
                       "textTransform": "uppercase",
                       "letterSpacing": "0.04em"}),
        html.Div(partner_rows(in_, "←") or
                 [html.Div("none", style={"color": INK_DIM})]),

        html.H4("Most strongly expressed genes",
                style={"color": INK, "marginBottom": "6px", "fontSize": "12px",
                       "marginTop": "18px", "fontWeight": 600,
                       "textTransform": "uppercase",
                       "letterSpacing": "0.04em"}),
        html.Div([_chip(f"{g['g']}  {g['tpm']:,}",
                        CHIP_WATCH if g["watch"] else None, bold=g["watch"])
                  for g in node["top_genes"]]),

        html.H4("Watch-list genes switched on here",
                style={"color": INK, "marginBottom": "6px", "fontSize": "12px",
                       "marginTop": "18px", "fontWeight": 600,
                       "textTransform": "uppercase",
                       "letterSpacing": "0.04em"}),
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
