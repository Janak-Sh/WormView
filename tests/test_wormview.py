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
        fig = figure.make_figure(payload, node_index, **kwargs)
        assert len(fig.data) > 0, kwargs


def test_every_clickable_trace_carries_a_neuron_name(loaded):
    """Selection resolves through customdata. A marker trace without it is a dead
    click, which is what made the app feel broken."""
    from wormview import figure
    tpm, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    names = set(node_index)

    for kwargs in (dict(), dict(gene="fmi-1", tpm=tpm)):
        fig = figure.make_figure(payload, node_index, **kwargs)
        markers = [t for t in fig.data if getattr(t, "mode", None) == "markers"]
        assert markers, "no clickable markers at all"
        for trace in markers:
            assert trace.customdata is not None, f"{trace.name} is unclickable"
            assert set(trace.customdata) <= names


def test_selecting_a_neuron_highlights_its_own_neurites(loaded):
    from wormview import figure
    _, payload = loaded
    node_index = {n["name"]: n for n in payload["nodes"]}
    fig = figure.make_figure(payload, node_index, selected="AVA")
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
