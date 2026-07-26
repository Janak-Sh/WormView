"""Tests for the parts where a silent error would be costly.

Priorities, in order:
  1. the neuron-name mapping between three datasets that name things differently
  2. the neurite morphology parser, which had two silent bugs
  3. the figure builder, across every combination of options
  4. the settings that keep the 3D scene interactive
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from wormview import anatomy, connectome as cn, genes as gene_meta   # noqa: E402
from wormview.data import MissingData, load_matrix                   # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def tpm():
    try:
        matrix, _ = load_matrix(2)
    except MissingData:
        pytest.skip("expression data not downloaded; run fetch_data.py")
    return cn.merge_ambiguous_classes(matrix)


@pytest.fixture(scope="module")
def edges(tpm):
    try:
        table, _ = cn.load_edges(list(tpm.columns), "chemical")
    except cn.MissingConnectome:
        pytest.skip("connectome not downloaded; run fetch_data.py")
    return table


@pytest.fixture(scope="module")
def loaded():
    import app
    try:
        return app.load_everything()
    except Exception as exc:                                   # noqa: BLE001
        pytest.skip(f"data not available: {exc}")


# --------------------------------------------------------------------------
# 1. neuron-name mapping
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,expected", [
    # CeNGEN merges all VD and DD GABAergic motor neurons into one class
    ("VD13", "VD_DD"), ("DD01", "VD_DD"), ("VD01", "VD_DD"),
    # but pulls a few individual ventral cord neurons out of their class
    ("DA09", "DA9"), ("DA05", "DA"), ("VA12", "VA12"), ("VA07", "VA"),
    ("VB01", "VB01"), ("VB02", "VB02"), ("VB07", "VB"),
    ("DB01", "DB01"), ("DB04", "DB"), ("AS11", "AS"),
    ("VC04", "VC_4_5"), ("VC01", "VC"),
    # plain left/right pairs
    ("ADAL", "ADA"), ("PVQL", "PVQ"), ("ASEL", "ASEL"), ("ASER", "ASER"),
    # dorsoventral-plus-side names
    ("CEPDL", "CEP"), ("IL1VR", "IL1"), ("I1L", "I1"), ("M3L", "M3"),
    ("RMDDL", "RMD_DV"), ("RMDL", "RMD_LR"),
    ("RMED", "RME_DV"), ("RMEL", "RME_LR"),
    ("IL2DL", "IL2_DV"), ("IL2L", "IL2_LR"),
    # singletons
    ("AVL", "AVL"), ("DVB", "DVB"), ("PVT", "PVT"), ("MC", "MC"),
])
def test_neuron_name_maps_to_expression_class(name, expected, tpm):
    assert cn._map_one(name, set(tpm.columns)) == expected


@pytest.mark.parametrize("name", [
    "BWM-DL01", "dBWML7", "vBWMR12",     # body wall muscle
    "hyp", "GLRDL", "CEPshVL", "mu_int",  # hypodermis, glia, other muscle
    "mc1DL", "mc2v", "mc3R",              # pharyngeal marginal cells (epithelial)
])
def test_non_neurons_are_rejected(name, tpm):
    assert cn._map_one(name, set(tpm.columns)) is None


def test_uppercase_mc_is_a_neuron_but_lowercase_mc_is_not(tpm):
    """MC is a pharyngeal motor neuron; mc1/mc2/mc3 are marginal cells. Case is
    the only thing distinguishing them, so a case-insensitive filter conflates
    them and silently drops a real neuron."""
    valid = set(tpm.columns)
    assert cn._map_one("MC", valid) == "MC"
    assert cn._map_one("mc2dl", valid) is None


def test_all_neuron_to_neuron_edges_map(tpm):
    """Regression: an earlier mapping silently dropped 44% of the wiring."""
    report = cn.mapping_report(list(tpm.columns), "chemical")
    assert report["fraction_of_neuronal_kept"] == 1.0
    assert report["unmapped_names"] == []


def test_mapping_never_invents_a_class(tpm, edges):
    valid = set(tpm.columns)
    assert set(edges.pre) <= valid
    assert set(edges.post) <= valid


def test_ventral_cord_motor_neurons_are_present(edges):
    """The Witvliet reconstructions have none of these, which is why the Cook
    connectome is used instead. Losing them would gut the picture."""
    classes = set(edges.pre) | set(edges.post)
    assert {"VD_DD", "DA", "DB", "VA", "VB"} <= classes


def test_awc_is_merged_rather_than_guessed(tpm):
    """CeNGEN splits AWC by functional state (ON/OFF), which anatomy cannot
    assign. The merge must exist so neither state is picked arbitrarily."""
    assert "AWC" in tpm.columns
    assert cn._map_one("AWCL", set(tpm.columns)) == "AWC"


# --------------------------------------------------------------------------
# 2. geometry and morphology
# --------------------------------------------------------------------------

def test_every_class_in_the_wiring_has_a_position(tpm, edges):
    positions, _, _ = anatomy.load_positions(list(tpm.columns))
    needed = set(edges.pre) | set(edges.post)
    assert not (needed - set(positions)), f"no position for {needed - set(positions)}"


@pytest.mark.parametrize("neuron,lo,hi,why", [
    ("I1", -320, -280, "pharyngeal, extreme anterior"),
    ("AWA", -290, -240, "head sensory"),
    ("HSN", 40, 80, "mid-body, at the vulva"),
    ("PLM", 380, 430, "posterior touch neuron"),
])
def test_positions_are_anatomically_correct(neuron, lo, hi, why, tpm):
    """Spot-checks along the body axis. If these drift, the geometry is wrong."""
    positions, _, _ = anatomy.load_positions(list(tpm.columns))
    if neuron not in positions:
        pytest.skip(f"{neuron} absent")
    ap = positions[neuron][0]
    assert lo <= ap <= hi, f"{neuron} at AP {ap:.0f} um, expected {lo}..{hi} ({why})"


def test_morphology_is_a_connected_tree_not_a_zigzag():
    """Regression on two parser bugs. First version treated all proximal and
    distal points as one polyline in file order, so neurites zigzagged across the
    animal. Second required an explicit <proximal> per segment, keeping only
    branch roots -- 3 segments for AVAL instead of 56."""
    import json
    path = ROOT / "data" / "neuron_positions.json"
    if not path.exists():
        pytest.skip("run fetch_data.py")
    raw = json.loads(path.read_text())

    assert len(raw) == 302, f"expected 302 neurons, got {len(raw)}"

    segs = raw["AVAL"]["segments"]
    assert len(segs) > 40, f"AVAL collapsed to {len(segs)} segments"
    length = sum(float(np.linalg.norm(np.array(b) - np.array(a))) for a, b in segs)
    assert length > 500, f"AVAL is only {length:.0f} um; it spans the whole cord"

    ends = {tuple(b) for a, b in segs}
    joined = sum(1 for a, b in segs if tuple(a) in ends)
    assert joined >= 0.8 * len(segs), "segments do not form a connected tree"


def test_short_and_long_neurons_differ_as_expected():
    """A head sensory neuron must not come out as long as a command interneuron."""
    import json
    path = ROOT / "data" / "neuron_positions.json"
    if not path.exists():
        pytest.skip("run fetch_data.py")
    raw = json.loads(path.read_text())

    def ap_extent(name):
        pts = np.array([p for seg in raw[name]["segments"] for p in seg])
        return float(np.ptp(pts[:, 1]))

    assert ap_extent("AVAL") > 500     # runs the whole ventral cord
    assert ap_extent("AWAL") < 200     # head sensory, local


def test_body_wall_contains_the_neurites(loaded):
    """The wall must be fitted to the neurites, not just the cell bodies:
    nerve cords run right along the body wall and reach much further out."""
    _, payload = loaded
    if not payload.get("surface"):
        pytest.skip("no surface")
    _, Y, _ = payload["surface"]
    pts = np.array([p for segs in payload["morphology"].values()
                    for seg in segs for p in seg])
    assert np.abs(pts[:, 1]).max() <= np.abs(Y).max() + 1e-6


def test_spread_never_moves_a_neuron_along_the_body(loaded):
    """The declutter step may only act in cross-section. Position along the body
    is the axis that carries the anatomy."""
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    true_ap = {k: v[0] for k, v in payload["true_positions"].items()}

    from wormview import figure
    figure.apply_spread(payload, node_index, 3.0)
    for name, ap in true_ap.items():
        if name in node_index:
            assert np.isclose(node_index[name]["x"], ap), f"{name} moved along body"


def test_spread_actually_separates_crowded_cell_bodies(loaded):
    _, payload = loaded
    positions = payload["true_positions"]

    def min_sep(coords):
        pts = np.array(list(coords.values()))
        d = np.linalg.norm(pts[:, None] - pts[None, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        return d.min()

    assert min_sep(anatomy.declutter(positions, 1.5)) > min_sep(positions)


# --------------------------------------------------------------------------
# 3. the figure
# --------------------------------------------------------------------------

def test_figure_builds_for_every_option_combination(loaded):
    """A NameError inside a Dash callback only shows up at runtime in a browser,
    so exercise the whole option matrix here instead."""
    from wormview import figure
    tpm, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    for kwargs in (
        dict(),
        dict(selected="AVA"),
        dict(gene="unc-25", tpm=tpm),
        dict(selected="VD_DD", show_labels=True, show_synapses=True, spread=1.5),
        dict(selected="PLM", gene="cle-1", tpm=tpm, show_synapses=True),
        dict(spread=3.0),
    ):
        fig, click_map = figure.make_figure(payload, node_index, **kwargs)
        assert len(fig.data) > 0, kwargs
        assert click_map, f"no clickable traces for {kwargs}"


def test_every_clickable_trace_carries_a_neuron_name(loaded):
    """Selection resolves through customdata. A marker trace without it is a dead
    click, which is what made the app feel broken."""
    from wormview import figure
    tpm, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    names = set(node_index)

    for kwargs in (dict(), dict(gene="fmi-1", tpm=tpm)):
        fig, _ = figure.make_figure(payload, node_index, **kwargs)
        markers = [t for t in fig.data if getattr(t, "mode", None) == "markers"]
        assert markers, "no clickable markers at all"
        for trace in markers:
            assert trace.customdata is not None, f"{trace.name} is unclickable"
            assert set(trace.customdata) <= names


def test_selecting_a_neuron_highlights_its_own_neurites(loaded):
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig, _ = figure.make_figure(payload, node_index, selected="AVA")
    assert any("AVA" in (t.name or "") for t in fig.data)


# --------------------------------------------------------------------------
# 4. interaction settings
# --------------------------------------------------------------------------

def test_graph_is_configured_for_interaction():
    """Rotate, zoom and click all break together if the canvas is cropped or not
    responsive: the container sizes the canvas wrongly and pointer coordinates
    stop matching what is on screen."""
    src = (ROOT / "app.py").read_text()
    assert "responsive=True" in src, "dcc.Graph must be responsive"
    assert '"scrollZoom": True' in src, "scroll zoom must be enabled"
    assert '"overflow": "hidden"' not in src, \
        "overflow:hidden on the graph wrapper crops the canvas and offsets clicks"


# --------------------------------------------------------------------------
# 5. the gene watch-list
# --------------------------------------------------------------------------

def test_watchlist_genes_all_exist_in_the_expression_data(tpm):
    missing = [g for g in gene_meta.WATCHLIST if g not in tpm.index]
    assert not missing, f"watch-list genes absent from CeNGEN: {missing}"


def test_secreted_genes_carry_an_autonomy_caveat():
    """CeNGEN sequenced neurons only. For a secreted protein made by muscle,
    absence in neurons proves nothing, and the UI must say so."""
    for gene in ("nid-1", "cle-1", "unc-52", "lin-44"):
        assert gene_meta.autonomy_caveat(gene), f"{gene} needs a caveat"


def test_click_resolves_without_customdata(loaded):
    """The real bug behind "clicking does nothing".

    Since plotly 6.0, Dash's clickData no longer carries `customdata`
    (plotly/plotly.py#5119), so a callback that reads only customdata silently
    fails on every click. The click map must resolve a neuron from
    curveNumber + pointNumber alone.

    An earlier HTTP test passed only because it hand-crafted clickData WITH
    customdata -- a payload the browser never sends. This asserts the path the
    browser actually takes.
    """
    from wormview import figure
    tpm, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    for kwargs in (dict(), dict(gene="fmi-1", tpm=tpm)):
        fig, click_map = figure.make_figure(payload, node_index, **kwargs)
        assert click_map, "no clickable traces"

        for curve, names in click_map.items():
            trace = fig.data[curve]
            assert getattr(trace, "mode", None) == "markers", \
                f"trace {curve} is in the click map but is not a marker trace"
            # one name per plotted point, or an index lookup goes to the wrong neuron
            assert len(names) == len(trace.x), \
                f"trace {curve}: {len(names)} names for {len(trace.x)} points"
            for name in names:
                assert name in node_index

        # simulate what the browser actually sends: no customdata
        curve = sorted(click_map)[0]
        point = {"curveNumber": curve, "pointNumber": 0}
        resolved = click_map[curve][point["pointNumber"]]
        assert resolved in node_index


def test_click_map_point_order_matches_the_plotted_coordinates(loaded):
    """A name is looked up by position, so the order must match exactly. If the
    two drift apart, clicks silently select the wrong neuron -- worse than not
    working at all, because it looks like it works."""
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig, click_map = figure.make_figure(payload, node_index)

    for curve, names in click_map.items():
        trace = fig.data[curve]
        for i, name in enumerate(names):
            node = node_index[name]
            assert np.isclose(trace.x[i], node["x"]), f"{name} x mismatch at {i}"
            assert np.isclose(trace.y[i], node["y"]), f"{name} y mismatch at {i}"
            assert np.isclose(trace.z[i], node["z"]), f"{name} z mismatch at {i}"


# --------------------------------------------------------------------------
# 6. group toggles and the body wall
# --------------------------------------------------------------------------

def test_group_toggle_hides_that_group(loaded):
    """The colour key doubles as on/off controls, replacing what Plotly's
    click-to-toggle legend used to do."""
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    _, all_map = figure.make_figure(payload, node_index)
    _, motor_map = figure.make_figure(payload, node_index, groups=["motor"])

    all_names = {n for names in all_map.values() for n in names}
    motor_names = {n for names in motor_map.values() for n in names}
    assert motor_names < all_names
    assert all(node_index[n]["group"] == "motor" for n in motor_names)


def test_body_wall_can_be_turned_off(loaded):
    """Cell bodies sit inside the wall, so a click aimed at one passes through it.
    Turning the wall off has to be possible or crowded neurons stay unclickable."""
    import plotly.graph_objects as go
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    def has_surface(**kw):
        fig, _ = figure.make_figure(payload, node_index, **kw)
        return any(isinstance(t, go.Surface) for t in fig.data)

    assert has_surface(show_body=True)
    assert not has_surface(show_body=False)


def test_click_map_stays_consistent_under_every_filter(loaded):
    """Filtering changes how many traces exist, so a stale index would resolve a
    click to the wrong neuron -- worse than not working, because it looks fine."""
    from wormview import figure
    tpm, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    for kwargs in (
        dict(),
        dict(groups=["motor"]),
        dict(groups=["sensory", "motor"], show_body=False),
        dict(groups=["interneuron"], gene="fmi-1", tpm=tpm),
        dict(groups=["motor"], show_synapses=True, show_labels=True),
        dict(groups=[]),
    ):
        fig, click_map = figure.make_figure(payload, node_index, **kwargs)
        for curve, names in click_map.items():
            trace = fig.data[curve]
            assert getattr(trace, "mode", None) == "markers", \
                f"{kwargs}: trace {curve} in the click map is not markers"
            assert len(names) == len(trace.x), \
                f"{kwargs}: trace {curve} has {len(names)} names, {len(trace.x)} points"
            for i, name in enumerate(names):
                assert np.isclose(trace.x[i], node_index[name]["x"])


def test_markers_are_large_enough_to_click(loaded):
    """In a 3D scene the hit area is the marker, so tiny dots are unclickable."""
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig, click_map = figure.make_figure(payload, node_index)
    for curve in click_map:
        sizes = fig.data[curve].marker.size
        assert min(sizes) >= 7, f"smallest marker is {min(sizes)}px"


# --------------------------------------------------------------------------
# 7. expression on the whole cell, and partner highlighting
# --------------------------------------------------------------------------

def test_gene_mode_colours_the_neurites_not_only_the_cell_bodies(loaded):
    """"Show expression on the neuron" means the whole cell. Colouring only the
    cell bodies leaves a gene like mec-7 as four dots instead of four long
    processes running down the body."""
    from wormview import figure
    tpm, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    fig, _ = figure.make_figure(payload, node_index, gene="mec-7", tpm=tpm)
    neurites = [t for t in fig.data if t.name == "neurites"]
    assert len(neurites) == 1, "gene mode should draw one gradient neurite trace"

    colours = neurites[0].line.color
    assert not isinstance(colours, str), "neurites must carry a per-point colour array"
    # Plotly rejects None inside a colour array, so separators need a real number
    assert all(c is not None for c in colours)
    assert len(colours) == len(neurites[0].x), "colour array must match coordinates"
    assert max(colours) > min(colours), "no variation, so nothing is being shown"


def test_type_mode_still_colours_neurites_by_group(loaded):
    """The original view must survive the addition."""
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig, _ = figure.make_figure(payload, node_index)
    names = {t.name for t in fig.data}
    assert {"sensory", "interneuron", "motor"} <= names


def test_selecting_a_neuron_rings_its_partners(loaded):
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    fig, _ = figure.make_figure(payload, node_index, selected="ALM")
    rings = [t for t in fig.data if t.name == "partners"]
    assert rings, "no partner rings"

    expected = {q["name"] for q in node_index["ALM"]["partners"]
                if q["name"] != "ALM"}
    assert set(rings[0].customdata) == expected


def test_partner_rings_are_clickable(loaded):
    """Rings sit on top of the cell bodies, so if they were not in the click map a
    click landing on one would resolve to nothing and appear to do nothing. Being
    clickable also makes them useful: click a partner to jump to it."""
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}

    fig, click_map = figure.make_figure(payload, node_index, selected="ALM")
    ring_curves = [c for c in click_map if fig.data[c].name == "partners"]
    assert ring_curves, "partner rings are not in the click map"
    for curve in ring_curves:
        names = click_map[curve]
        assert len(names) == len(fig.data[curve].x)
        for i, name in enumerate(names):
            assert np.isclose(fig.data[curve].x[i], node_index[name]["x"])


def test_no_partner_rings_without_a_selection(loaded):
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig, _ = figure.make_figure(payload, node_index)
    assert not any(t.name == "partners" for t in fig.data)


def test_partner_rings_respect_the_group_filter(loaded):
    """A hidden group must not get rings, or the filter silently leaks."""
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig, _ = figure.make_figure(payload, node_index, selected="ALM",
                                groups=["motor"])
    rings = [t for t in fig.data if t.name == "partners"]
    if rings:
        for name in rings[0].customdata:
            assert node_index[name]["group"] == "motor"


def test_partner_ring_colour_is_not_a_group_colour(loaded):
    """The rings were invisible at first because they used the motor colour."""
    from wormview import theme
    assert theme.PARTNER not in theme.GROUP_COLOR.values()
    assert theme.PARTNER != theme.HILITE
