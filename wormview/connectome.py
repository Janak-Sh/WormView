"""The wiring diagram, and mapping it onto CeNGEN's neuron classes.

Connectivity comes from Cook et al. 2019 (Nature), the whole-animal hermaphrodite
connectome -- chosen over the Witvliet 2021 reconstructions because Witvliet covers
the brain/nerve ring only and contains *no* ventral cord motor neurons. Cook
includes all 75 of them.

The awkward part is that the connectome names individual neurons (VD13, ADAL,
DA09) while CeNGEN reports classes (VD_DD, ADA, DA9). `_map_one()` does that
translation, and `mapping_report()` states how much of the connectome survives it --
a mapping that silently dropped a third of the wiring would invalidate everything
downstream.
"""

import re
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "connectome"

EDGE_FILE = "herm_full_edgelist.csv"
ADJACENCY_FILE = "Dataset8_adjacency.csv"

EDGE_URL = ("https://raw.githubusercontent.com/openworm/ConnectomeToolbox/"
            "main/cect/data/herm_full_edgelist.csv")
ADJACENCY_URL = ("https://raw.githubusercontent.com/dwitvliet/nature2021/"
                 "master/data/physical_contact/Dataset8_adjacency.csv")

# Anything matching these is not a neuron: body wall muscle (BWM, and the
# dBWML/vBWMR variants Cook uses), other muscle, hypodermis, glia, gonad.
# These are not failures of the mapping -- CeNGEN sequenced neurons only, so a
# neuron-to-muscle junction genuinely has no expression counterpart and must be
# excluded rather than counted as lost.
_NOT_A_NEURON = re.compile(
    r"^([dv]?BWM|mu_|hyp|exc|int|GLR|CEPsh|vm|um|sph|g1|g2|hmc|HSNL_|obl|"
    r"defecation|anal|intestine|BAG_|PM|MI_)", re.I)

# Case matters here, so this cannot fold into the pattern above: lowercase
# mc1/mc2/mc3 are pharyngeal *marginal cells* (epithelial), while uppercase MC is
# a genuine pharyngeal motor neuron that CeNGEN sequenced.
_MARGINAL_CELL = re.compile(r"^mc\d")


def _rstrip_lr(name):
    """ADAL -> ADA, RIMR -> RIM. Leaves ASEL/ASER alone (handled by exact map)."""
    return name[:-1] if len(name) > 3 and name[-1] in "LR" else name


def _map_one(name, valid):
    """Map an individual neuron name to a CeNGEN class, or None.

    `valid` is the set of CeNGEN class names, so the mapping can never invent a
    class that does not exist in the expression matrix.
    """
    n = name.strip()
    if not n or _NOT_A_NEURON.match(n) or _MARGINAL_CELL.match(n):
        return None

    # 1. exact hit (ASEL, ASER, AVL, DVB, PVT, VA12, VB01, VB02, DB01 ...)
    if n in valid:
        return n

    # 2. CeNGEN pulls a few individual neurons out of their class using a
    #    different zero-padding convention than the connectome does.
    for pad, plain in (("DA09", "DA9"), ("VA12", "VA12"),
                       ("VB01", "VB01"), ("VB02", "VB02"), ("DB01", "DB01")):
        if n == pad and plain in valid:
            return plain

    # 3. GABAergic ventral cord motor neurons: CeNGEN merges VD and DD.
    if re.match(r"^(VD|DD)\d+$", n) and "VD_DD" in valid:
        return "VD_DD"

    # 4. VC4 and VC5 are their own CeNGEN class; the rest are plain VC.
    if re.match(r"^VC0?[45]$", n) and "VC_4_5" in valid:
        return "VC_4_5"

    # 5. classes CeNGEN splits by dorsoventral vs left-right position
    for stem in ("IL2", "RMD", "RME"):
        if n.startswith(stem):
            suffix = n[len(stem):]
            if re.match(r"^[DV][LR]?$", suffix) and f"{stem}_DV" in valid:
                return f"{stem}_DV"
            if re.match(r"^[LR]$", suffix) and f"{stem}_LR" in valid:
                return f"{stem}_LR"

    # 6. Strip positional suffixes, longest name first so we never over-strip
    #    into a shorter class that happens to exist. This covers numbered members
    #    (DA05 -> DA), left/right pairs (ADAL -> ADA) and dorsoventral-plus-side
    #    names (CEPDL -> CEP, IL1VR -> IL1, I1L -> I1, M3L -> M3).
    for cut in range(1, 4):
        if len(n) - cut < 2:
            break
        stem, tail = n[:-cut], n[-cut:]
        # only strip characters that encode position or an index
        if not re.fullmatch(r"[DVLR0-9]+", tail):
            continue
        if stem in valid:
            return stem
        bare = re.sub(r"\d+$", "", stem)
        if bare and bare != stem and bare in valid:
            return bare

    return None


# AWC is a special case worth naming rather than hiding. CeNGEN splits it into
# AWC_ON and AWC_OFF, which are functional states, not anatomical sides, so the
# connectome's AWCL/AWCR cannot be assigned to one or the other. Rather than pick
# arbitrarily, merge_ambiguous_classes() adds a combined "AWC" column (a gene
# counts as present if it is present in either state) and the mapper targets that.
AMBIGUOUS = {"AWCL", "AWCR"}
MERGED_CLASSES = {"AWC": ("AWC_ON", "AWC_OFF")}


def merge_ambiguous_classes(tpm):
    """Add combined columns for classes the connectome cannot disambiguate.

    Call this on the expression matrix before using it with the connectome.
    Returns a copy; the original split columns are left in place.
    """
    out = tpm.copy()
    for merged, parts in MERGED_CLASSES.items():
        have = [p for p in parts if p in out.columns]
        if have and merged not in out.columns:
            out[merged] = out[have].max(axis=1)
    return out


class MissingConnectome(RuntimeError):
    pass


def load_edges(valid_classes, kind="chemical"):
    """Load the connectome and collapse it onto CeNGEN classes.

    kind: 'chemical', 'electrical', or 'both'.

    Returns (edges, unmapped) where edges has one row per class-pair with
    columns pre, post, weight (summed synapses), n_pairs (how many individual
    neuron pairs were merged into it), and self_loop.
    """
    path = DATA_DIR / EDGE_FILE
    if not path.exists():
        raise MissingConnectome(
            f"{path} not found.\nRun:  python run_wiring.py --fetch-data")

    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    for col in ("Source", "Target", "Type"):
        df[col] = df[col].astype(str).str.strip()

    if kind != "both":
        df = df[df.Type == kind]
    if df.empty:
        raise ValueError(f"no edges of type {kind!r}")

    valid = set(valid_classes)
    cache = {}

    def mp(name):
        if name not in cache:
            cache[name] = _map_one(name, valid)
        return cache[name]

    df = df.assign(pre=df.Source.map(mp), post=df.Target.map(mp))

    seen = set(df.Source) | set(df.Target)
    unmapped = sorted(n for n in seen if mp(n) is None
                      and not _NOT_A_NEURON.match(n.strip())
                      and not _MARGINAL_CELL.match(n.strip()))

    df = df.dropna(subset=["pre", "post"])
    grouped = (df.groupby(["pre", "post"], as_index=False)
                 .agg(weight=("Weight", "sum"), n_pairs=("Weight", "size")))
    grouped["self_loop"] = grouped.pre == grouped.post
    return grouped, unmapped


def mapping_report(valid_classes, kind="chemical"):
    """How much of the connectome survives the class mapping?

    Reported against neuron-to-neuron edges, which is the honest denominator:
    CeNGEN sequenced neurons only, so a neuron-to-muscle junction genuinely has no
    expression counterpart and is excluded rather than counted as lost.
    """
    path = DATA_DIR / EDGE_FILE
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    for col in ("Source", "Target", "Type"):
        df[col] = df[col].astype(str).str.strip()
    if kind != "both":
        df = df[df.Type == kind]

    valid = set(valid_classes)

    def is_neuron(name):
        n = name.strip()
        return not (_NOT_A_NEURON.match(n) or _MARGINAL_CELL.match(n))

    nn = df[[is_neuron(s) and is_neuron(t) for s, t in zip(df.Source, df.Target)]]
    ok = sum(1 for s, t in zip(nn.Source, nn.Target)
             if _map_one(s, valid) and _map_one(t, valid))
    edges, unmapped = load_edges(valid_classes, kind)
    return {
        "raw_edges": len(df),
        "neuron_neuron_edges": len(nn),
        "non_neuronal_edges": len(df) - len(nn),
        "mapped_edges": ok,
        "fraction_of_neuronal_kept": ok / len(nn) if len(nn) else 0.0,
        "class_pairs": len(edges),
        "self_loops": int(edges.self_loop.sum()),
        "unmapped_names": unmapped,
        "ambiguous": sorted(AMBIGUOUS),
    }
