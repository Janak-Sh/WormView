"""Where the neurons are, and what shape they are.

Reads the real 3D geometry of all 302 C. elegans neurons -- soma position and full
neurite morphology, in microns -- from the WormBase Virtual Worm model, and
collapses individual neurons onto the CeNGEN classes the expression atlas uses.

Also builds the translucent body wall the nervous system sits inside, and a
declutter step that separates crowded cell bodies without moving them along the
body axis.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import connectome as cn
from . import genes as gene_meta

CELL_INFO = cn.DATA_DIR / "all_cell_info.csv"
CELL_INFO_URL = ("https://raw.githubusercontent.com/openworm/ConnectomeToolbox/"
                 "main/cect/data/all_cell_info.csv")

# Broad functional groups, in the order they appear in the legend.
GROUPS = ["sensory", "interneuron", "motor", "other"]


def load_cell_types():
    """CeNGEN class -> broad functional group, from the Cook cell annotations.

    The file has a couple of malformed rows, hence the tolerant read.
    """
    if not CELL_INFO.exists():
        return {}
    try:
        info = pd.read_csv(CELL_INFO, on_bad_lines="skip")
    except TypeError:                          # older pandas
        info = pd.read_csv(CELL_INFO, error_bad_lines=False)

    out = {}
    for row in info.itertuples():
        name = str(getattr(row, "_1", "") or "").strip()   # "Cell name"
        text = " ".join(str(getattr(row, c, "")) for c in
                        ("Type", "Classification")).lower()
        if not name:
            continue
        if "sensory" in text or "receptor" in text:
            group = "sensory"
        elif "motor" in text or "muscle" in text:
            group = "motor"
        elif "interneuron" in text or "ring" in text:
            group = "interneuron"
        else:
            group = "other"
        out[name] = group
    return out


def group_for(cengen_class, cell_types):
    """Best functional group for a CeNGEN class name."""
    if cengen_class in cell_types:
        return cell_types[cengen_class]
    # try the individual neurons that collapse into this class
    stem = re.sub(r"_(DV|LR|\d+_\d+)$", "", cengen_class)
    for suffix in ("", "L", "R", "01", "1", "DL", "VL"):
        if stem + suffix in cell_types:
            return cell_types[stem + suffix]
    # ventral cord motor neuron families
    if re.match(r"^(VD_DD|VA|VB|DA|DB|AS|VC)", cengen_class):
        return "motor"
    return "other"


def build_payload(tpm_all, edges, top_genes=24, watchlist=None):
    """Everything the page needs, as plain Python ready for JSON.

    tpm_all is the FULL CeNGEN matrix (all ~13,700 genes) so a neuron's gene list
    is genuinely its transcriptome, not just this lab's shortlist.
    """
    watchlist = set(watchlist or gene_meta.WATCHLIST)
    edges = edges[~edges.self_loop]
    classes = sorted(set(edges.pre) | set(edges.post))

    cell_types = load_cell_types()
    edge_list = [(r.pre, r.post, int(r.weight)) for r in edges.itertuples()]

    # per-neuron connection lists, with direction and strength
    partners = {c: [] for c in classes}
    degree = {c: 0 for c in classes}
    for a, b, weight in edge_list:
        partners[a].append({"name": b, "w": weight, "dir": "out"})
        partners[b].append({"name": a, "w": weight, "dir": "in"})
        degree[a] += 1
        degree[b] += 1
    for c in partners:
        partners[c].sort(key=lambda p: -p["w"])

    nodes = []
    for c in classes:
        col = tpm_all[c] if c in tpm_all.columns else None
        genes = []
        if col is not None:
            on = col[col > 0].sort_values(ascending=False)
            for g in on.head(top_genes).index:
                genes.append({"g": g, "tpm": round(float(on[g])),
                              "watch": g in watchlist})
            n_on = int((col > 0).sum())
            watch_here = [{"g": g, "tpm": round(float(tpm_all.loc[g, c]))}
                          for g in sorted(watchlist)
                          if g in tpm_all.index and tpm_all.loc[g, c] > 0]
        else:
            n_on, watch_here = 0, []

        nodes.append({
            "name": c,
            "group": group_for(c, cell_types),
            "degree": degree[c],
            "n_genes_on": n_on,
            "top_genes": genes,
            "watch_genes": watch_here,
            "partners": partners[c],
        })

    return {
        "nodes": nodes,
        "edges": [{"a": a, "b": b, "w": weight} for a, b, weight in edge_list],
        "groups": GROUPS,
        "watchlist": sorted(watchlist),
        "n_genes_total": int(tpm_all.shape[0]),
    }


# --------------------------------------------------------------------------
# Real anatomical positions
# --------------------------------------------------------------------------

POSITIONS_FILE = Path(__file__).resolve().parent.parent / "data" / "neuron_positions.json"

# The NeuroML files name ventral-cord neurons without zero padding (DB1), while
# CeNGEN pads them (DB01). Normalise before mapping or those classes lose their
# position and silently fall back to the graph layout.
_PAD = re.compile(r"^([A-Z]+)(\d)$")


def _padded(name):
    m = _PAD.match(name)
    return f"{m.group(1)}0{m.group(2)}" if m else name


def load_positions(valid_classes):
    """CeNGEN class -> (ap, lr, dv) centroid in microns, from the Virtual Worm.

    Individual neurons are averaged into their CeNGEN class, so VD_DD sits at the
    mean position of all nineteen VD and DD neurons -- the middle of the ventral
    cord, which is honest for a merged class but worth saying out loud.

    Returns (positions, per_class_counts, missing).
    """
    import json
    from collections import defaultdict

    if not POSITIONS_FILE.exists():
        return {}, {}, sorted(valid_classes)

    raw = json.loads(POSITIONS_FILE.read_text())
    valid = set(valid_classes)
    grouped = defaultdict(list)
    for name, entry in raw.items():
        padded = _padded(name)
        cls = (cn._map_one(padded, valid) if padded != name else None) \
            or cn._map_one(name, valid)
        if cls:
            grouped[cls].append(entry["soma"])

    positions, counts = {}, {}
    for cls, somas in grouped.items():
        arr = np.array(somas, dtype=float)
        lr, ap, dv = arr.mean(axis=0)
        positions[cls] = (float(ap), float(lr), float(dv))   # AP first
        counts[cls] = len(somas)
    missing = sorted(valid - set(positions))
    return positions, counts, missing


def worm_surface(positions, morphology=None, n_rings=140, n_theta=44):
    """A smooth translucent body wall for the nervous system to sit inside.

    A tapered tube with an elliptical cross-section, following the animal's
    dorsoventral curve. Deliberately schematic: the silhouette of a worm, not a
    scan of one.

    It must be fitted to the NEURITES, not just the cell bodies. Nerve cords and
    sensory dendrites run right along the body wall and reach much further out
    than the somata do, so a tube sized to cell bodies leaves most of the
    morphology hanging outside it -- which is what made the body and the wireframe
    look like two unrelated objects.
    """
    if not positions:
        return None
    cloud = list(positions.values())
    if morphology:
        for segs in morphology.values():
            for a, b in segs:
                cloud.append(a)
                cloud.append(b)
    pts = np.array(cloud, dtype=float)
    ap, lr, dv = pts[:, 0], pts[:, 1], pts[:, 2]

    ap_lo, ap_hi = ap.min(), ap.max()
    length = ap_hi - ap_lo
    centres = np.linspace(ap_lo - 0.02 * length, ap_hi + 0.02 * length, n_rings)

    # the body drifts dorsally in the head; follow it with a smooth fit
    fit = np.polyfit(ap, dv, 3)
    mid_dv = np.polyval(fit, centres)

    # generous: the wall should sit outside everything it contains
    half_lr = np.percentile(np.abs(lr), 99.5) + 8.0
    half_dv = np.percentile(np.abs(dv - np.polyval(fit, ap)), 99.5) + 10.0

    # taper: blunt nose, fuller middle, pointed tail
    t = np.linspace(0.0, 1.0, n_rings)
    taper = np.sin(np.pi * np.clip(t, 0, 1)) ** 0.42
    taper = np.clip(taper * (1.0 - 0.30 * t), 0.10, 1.0)

    theta = np.linspace(0, 2 * np.pi, n_theta)
    X = np.repeat(centres[:, None], n_theta, axis=1)
    Y = (half_lr * taper)[:, None] * np.cos(theta)[None, :]
    Z = mid_dv[:, None] + (half_dv * taper)[:, None] * np.sin(theta)[None, :]
    return X, Y, Z


def declutter(positions, strength=1.0, min_dist=9.0, iterations=200):
    """Spread overlapping neurons in cross-section only, never along the body.

    Most of the 302 neurons sit in the head, so at true coordinates the head is an
    unclickable clump. This pushes neurons apart in the left-right and
    dorsoventral plane while leaving the anterior-posterior coordinate exactly as
    measured -- so position along the body, the axis that carries the anatomy,
    stays truthful.
    """
    if not positions or strength <= 0:
        return dict(positions)
    names = list(positions)
    arr = np.array([positions[n] for n in names], dtype=float)
    ap = arr[:, 0].copy()
    cross = arr[:, 1:]

    rng = np.random.default_rng(5)
    cross = cross + rng.normal(scale=0.4, size=cross.shape)
    target = min_dist * strength

    for _ in range(iterations):
        # only neurons at a similar position along the body can visually overlap
        d_ap = np.abs(ap[:, None] - ap[None, :])
        delta = cross[:, None, :] - cross[None, :, :]
        dist = np.linalg.norm(delta, axis=-1)
        np.fill_diagonal(dist, np.inf)
        near = (dist < target) & (d_ap < target * 1.6)
        if not near.any():
            break
        unit = delta / (dist[:, :, None] + 1e-9)
        push = (unit * np.where(near, (target - dist) / target, 0.0)[:, :, None]
                ).sum(axis=1)
        cross += push * 0.35 * target

    return {n: (float(ap[i]), float(cross[i, 0]), float(cross[i, 1]))
            for i, n in enumerate(names)}


def load_morphology(valid_classes):
    """CeNGEN class -> the neurite segments of every neuron in that class.

    Each entry is a list of [[x0,y0,z0],[x1,y1,z1]] pairs in microns, with the
    axes reordered to match load_positions(): (anterior-posterior, left-right,
    dorsal-ventral).
    """
    import json
    from collections import defaultdict

    if not POSITIONS_FILE.exists():
        return {}

    raw = json.loads(POSITIONS_FILE.read_text())
    valid = set(valid_classes)
    out = defaultdict(list)
    for name, entry in raw.items():
        padded = _padded(name)
        cls = (cn._map_one(padded, valid) if padded != name else None) \
            or cn._map_one(name, valid)
        if not cls:
            continue
        for a, b in entry.get("segments", []):
            # (lr, ap, dv) in the file -> (ap, lr, dv) on screen
            out[cls].append([[a[1], a[0], a[2]], [b[1], b[0], b[2]]])
    return dict(out)


def segments_to_lines(segments):
    """Flatten [[p0,p1], ...] into x/y/z lists with None between pieces."""
    xs, ys, zs = [], [], []
    for a, b in segments:
        xs += [a[0], b[0], None]
        ys += [a[1], b[1], None]
        zs += [a[2], b[2], None]
    return xs, ys, zs


# --------------------------------------------------------------------------
# where synapses actually are
# --------------------------------------------------------------------------

CONTACTS_FILE = POSITIONS_FILE.parent / "contact_points.json"


def contact_points(edges, morphology, cache=True):
    """For each connection, the point where the two neurons' neurites come closest.

    A synapse forms where two *neurites* touch, not at the cell bodies. PVR sits in
    the tail but runs its axon the length of the body and synapses onto IL1 in the
    head -- so a line drawn between those two somata is 690 um long and points at
    the wrong end of the animal, while the actual contact is 0.0 um wide.

    Measured over all 1,764 connections: the neurites come within 5 um in 99% of
    cases, but soma-to-soma distance exceeds 300 um for 22% of them. Drawing links
    between cell bodies was therefore an abstraction pretending to be anatomy.

    The connectome gives partner pairs and synapse counts, not synapse coordinates,
    so closest approach is the best available estimate of where the contact is.

    Returns a list of {pre, post, weight, x, y, z, gap} and caches to disk, since
    the all-pairs distance work takes a few seconds.
    """
    import json

    if cache and CONTACTS_FILE.exists():
        try:
            return json.loads(CONTACTS_FILE.read_text())
        except (ValueError, OSError):
            pass

    clouds = {}
    for cls, segs in morphology.items():
        if segs:
            clouds[cls] = np.array([p for seg in segs for p in seg], dtype=float)

    out = []
    for row in edges.itertuples():
        if row.pre == row.post:
            continue
        a, b = clouds.get(row.pre), clouds.get(row.post)
        if a is None or b is None:
            continue
        d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
        i, j = np.unravel_index(np.argmin(d), d.shape)
        mid = (a[i] + b[j]) / 2.0
        out.append({"pre": row.pre, "post": row.post, "weight": int(row.weight),
                    "x": round(float(mid[0]), 2), "y": round(float(mid[1]), 2),
                    "z": round(float(mid[2]), 2), "gap": round(float(d[i, j]), 2)})

    if cache:
        try:
            CONTACTS_FILE.write_text(json.dumps(out))
        except OSError:
            pass
    return out
