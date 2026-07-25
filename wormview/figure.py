"""The 3D scene.

Layers, back to front: translucent body wall, every neurite, optional straight
synaptic links, the selected neuron's neurites in gold, then the cell bodies --
which are the only clickable layer.
"""

import numpy as np
import plotly.graph_objects as go

from . import anatomy
from .theme import (BG, BODY, EDGE_DIM, FONT, GENE_SCALE, GROUP_COLOR, HILITE,
                    INK, INK_DIM, LINE)


def apply_spread(payload, node_index, spread):
    """Set node coordinates for this render: true anatomy, optionally spread."""
    positions = payload["true_positions"]
    coords = (anatomy.declutter(positions, strength=spread) if spread > 0
              else dict(positions))
    for name, xyz in coords.items():
        if name in node_index:
            node_index[name]["x"], node_index[name]["y"], node_index[name]["z"] = xyz


def edge_lines(payload, node_index, only=None):
    """Flattened 3D line segments. `only` restricts to one neuron's connections."""
    xs, ys, zs = [], [], []
    for e in payload["edges"]:
        if only and only not in (e["a"], e["b"]):
            continue
        a, b = node_index[e["a"]], node_index[e["b"]]
        xs += [a["x"], b["x"], None]
        ys += [a["y"], b["y"], None]
        zs += [a["z"], b["z"], None]
    return xs, ys, zs


def make_figure(payload, node_index, selected=None, gene=None, tpm=None,
                show_labels=False, spread=0.0, show_synapses=False,
                groups=None, show_body=True):
    """The 3D scene: translucent body, real neurite morphology, clickable somata.

    Neurites are the 13,566 real segments from the Virtual Worm model, not straight
    lines between cell bodies. The straight synaptic lines are off by default --
    with 1,764 of them they formed the grey sheet that used to bury the anatomy.

    Returns (figure, click_map). click_map maps a trace index to the neuron names
    of that trace's points, in order, so a click can be resolved from
    curveNumber + pointNumber.

    That indirection is necessary, not stylistic: since plotly 6.0, `customdata`
    is dropped from Dash's clickData (plotly/plotly.py#5119), so resolving a click
    through customdata alone silently fails on every click. curveNumber and
    pointNumber are always present.
    """
    apply_spread(payload, node_index, spread)
    click_map = {}
    all_nodes = payload["nodes"]
    visible = set(payload["groups"] if groups is None else groups)
    nodes = [n for n in all_nodes if n["group"] in visible]
    names = [n["name"] for n in nodes]
    morph = payload.get("morphology", {})
    if not nodes:                      # every group switched off
        nodes, names = all_nodes[:1], [all_nodes[0]["name"]]

    fig = go.Figure()

    # 1. the body wall.
    #
    # Toggleable, and it matters: the cell bodies sit INSIDE this surface, so a
    # click aimed at a neuron travels through it. If the surface captures the
    # click instead, the callback sees a trace that is not in the click map and
    # keeps the current selection -- which reads as "clicking does nothing".
    # Turning the wall off is the reliable way to hit a crowded neuron.
    if show_body and payload.get("surface"):
        X, Y, Z = payload["surface"]
        fig.add_trace(go.Surface(
            x=X, y=Y, z=Z, showscale=False, hoverinfo="skip", opacity=0.13,
            colorscale=BODY, showlegend=False,
            lighting=dict(ambient=0.62, diffuse=0.75, specular=0.12,
                          roughness=0.9, fresnel=0.2),
            lightposition=dict(x=0, y=-400, z=400),
            contours=dict(x=dict(highlight=False), y=dict(highlight=False),
                          z=dict(highlight=False)), name="body"))

    # 2. every neurite, thin, coloured by functional group
    for group in payload["groups"]:
        if group not in visible:
            continue
        members = [n["name"] for n in nodes if n["group"] == group]
        segs = [seg for m in members for seg in morph.get(m, [])]
        if not segs:
            continue
        xs, ys, zs = anatomy.segments_to_lines(segs)
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip", name=group,
            showlegend=False, opacity=0.45,
            line=dict(color=GROUP_COLOR[group], width=1.3)))

    # 3. optional straight synaptic connections
    if show_synapses:
        ex_, ey, ez = edge_lines(payload, node_index)
        fig.add_trace(go.Scatter3d(
            x=ex_, y=ey, z=ez, mode="lines", hoverinfo="skip", showlegend=False,
            name="synaptic links", line=dict(color=EDGE_DIM, width=0.8)))

    # 4. the selected neuron: its own neurites, bright and thick
    if selected and selected in morph:
        xs, ys, zs = anatomy.segments_to_lines(morph[selected])
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip", showlegend=False,
            name=f"{selected} neurites", line=dict(color=HILITE, width=5)))

    # 5. cell bodies -- the clickable layer
    sizes = np.array([n["degree"] for n in nodes], dtype=float)
    # Generous sizes: these are the click targets, and in 3D the hit area is the
    # marker itself. Small dots are precise but nearly unclickable.
    sizes = 8.0 + 10.0 * (sizes / max(sizes.max(), 1)) ** 0.6

    if gene and tpm is not None and gene in tpm.index:
        vals = np.array([float(tpm.loc[gene, n["name"]])
                         if n["name"] in tpm.columns else 0.0 for n in nodes])
        logged = np.log10(vals + 1)
        click_map[len(fig.data)] = list(names)
        fig.add_trace(go.Scatter3d(
            x=[n["x"] for n in nodes], y=[n["y"] for n in nodes],
            z=[n["z"] for n in nodes], mode="markers", customdata=names,
            showlegend=False, name="cell bodies",
            hovertext=[f"<b>{n['name']}</b><br>{gene}: "
                       f"{format(v, ',.0f') + ' TPM' if v > 0 else 'not detected'}"
                       f"<br>{n['degree']} connections" for n, v in zip(nodes, vals)],
            hoverinfo="text",
            marker=dict(size=sizes, color=logged, colorscale=GENE_SCALE,
                        cmin=0, cmax=max(logged.max(), 1e-6), opacity=1.0,
                        line=dict(color="rgba(255,255,255,0.9)", width=1.0),
                        colorbar=dict(
                            title=dict(text=f"{gene}<br>log10 TPM", side="right",
                                       font=dict(size=11, color=INK)),
                            tickfont=dict(size=10, color=INK_DIM),
                            thickness=12, len=0.5, x=1.0, outlinewidth=0))))
    else:
        for group in payload["groups"]:
            if group not in visible:
                continue
            members = [n for n in nodes if n["group"] == group]
            if not members:
                continue
            idx = [names.index(m["name"]) for m in members]
            click_map[len(fig.data)] = [m["name"] for m in members]
            fig.add_trace(go.Scatter3d(
                x=[m["x"] for m in members], y=[m["y"] for m in members],
                z=[m["z"] for m in members], mode="markers",
                name=group, showlegend=False,
                customdata=[m["name"] for m in members],
                hovertext=[
                    f"<b>{m['name']}</b><br>{group}"
                    f"<br>{m['n_cells']} neuron{'s' if m['n_cells'] > 1 else ''} "
                    f"in this class<br>{m['degree']} connections"
                    f"<br>{m['n_genes_on']:,} genes switched on" for m in members],
                hoverinfo="text",
                marker=dict(size=[sizes[i] for i in idx],
                            color=GROUP_COLOR[group], opacity=1.0,
                            line=dict(color="rgba(255,255,255,0.9)", width=1.0))))

    # 6. ring the selection so it stays findable while rotating
    if selected and selected in node_index:
        sel = node_index[selected]
        fig.add_trace(go.Scatter3d(
            x=[sel["x"]], y=[sel["y"]], z=[sel["z"]], mode="markers+text",
            text=[selected], textposition="top center", showlegend=False,
            textfont=dict(color="#8a6400", size=13), hoverinfo="skip",
            marker=dict(size=20, color="rgba(0,0,0,0)",
                        line=dict(color=HILITE, width=3))))

    if show_labels:
        fig.add_trace(go.Scatter3d(
            x=[n["x"] for n in nodes], y=[n["y"] for n in nodes],
            z=[n["z"] for n in nodes], mode="text", text=names,
            textfont=dict(color=INK_DIM, size=8), hoverinfo="skip",
            showlegend=False))

    pts = np.array([[n["x"], n["y"], n["z"]] for n in nodes])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    pad = 0.17 * (hi - lo)
    pad[pad == 0] = 1.0
    axis = dict(showbackground=False, showgrid=False, zeroline=False,
                showticklabels=False, title="", showspikes=False, color=INK_DIM)
    fig.update_layout(
        scene=dict(
            annotations=_end_labels(pts, lo, hi, pad),
            xaxis={**axis, "range": [lo[0] - pad[0], hi[0] + pad[0]],
                   "title": dict(text="anterior  \u2192  posterior",
                                 font=dict(size=10, color=INK_DIM))},
            yaxis={**axis, "range": [lo[1] - pad[1], hi[1] + pad[1]]},
            zaxis={**axis, "range": [lo[2] - pad[2], hi[2] + pad[2]]},
            # The animal is ~13:1. At true proportions it is an unreadable
            # thread, so the cross-section is exaggerated; position ALONG the
            # body stays to scale.
            aspectmode="manual", aspectratio=dict(x=4.6, y=1, z=1), bgcolor=BG,
            camera=dict(eye=dict(x=0.0, y=-2.1, z=0.8), center=dict(x=0, y=0, z=0),
                        up=dict(x=0, y=0, z=1),
                        projection=dict(type="perspective"))),
        paper_bgcolor=BG, plot_bgcolor=BG,
        font=dict(family=FONT, color=INK, size=12),
        # bottom-left: the top-left corner is where the HEAD label sits, since
        # the head is the low end of the body axis
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0), uirevision="keep")
    return fig, click_map


# --------------------------------------------------------------------------
# which end is which
# --------------------------------------------------------------------------

# In the Virtual Worm coordinate frame the anterior-posterior axis is plotted as
# x, and anterior (the head) is the LOW end -- the pharyngeal neurons I1-I6 sit at
# about -300 um and the posterior touch neuron PLM at about +410 um.
HEAD_IS_LOW_X = True


def _end_labels(pts, lo, hi, pad):
    """3D annotations marking the head and the tail.

    Scene annotations rather than a text trace, for two reasons: they stay
    legible at any zoom, and they add no trace, so the click map's trace indices
    are unaffected.

    Each label is parked just beyond its end of the body and vertically level
    with the neurons actually there, so it reads as belonging to that end rather
    than floating in space.
    """
    span = hi[0] - lo[0]
    near_head = pts[pts[:, 0] < lo[0] + 0.08 * span]
    near_tail = pts[pts[:, 0] > hi[0] - 0.08 * span]
    head_z = float(near_head[:, 2].mean()) if len(near_head) else float(pts[:, 2].mean())
    tail_z = float(near_tail[:, 2].mean()) if len(near_tail) else float(pts[:, 2].mean())

    def label(text, sub, x, z, xanchor):
        return dict(
            x=x, y=0.0, z=z,
            text=f"<b>{text}</b><br><span style='font-size:10px'>{sub}</span>",
            showarrow=True, arrowhead=2, arrowsize=1.1, arrowwidth=1.4,
            arrowcolor=INK_DIM, ax=0, ay=-60,
            font=dict(size=16, color=INK), xanchor="center",
            bgcolor="rgba(255,255,255,0.82)", bordercolor=LINE, borderwidth=1,
            borderpad=5, opacity=0.97)

    # Just INSIDE each end, not beyond it: with the body filling the frame,
    # anything placed outside the data range lands off-canvas and is never seen.
    inset = 0.03 * span
    head_x = lo[0] + inset
    tail_x = hi[0] - inset
    if not HEAD_IS_LOW_X:
        head_x, tail_x = tail_x, head_x
        head_z, tail_z = tail_z, head_z

    return [label("HEAD", "anterior", head_x, head_z, "right"),
            label("TAIL", "posterior", tail_x, tail_z, "left")]
