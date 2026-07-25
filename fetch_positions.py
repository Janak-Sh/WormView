#!/usr/bin/env python3
"""Download real 3D neuron positions and neurite shapes for C. elegans.

Called by fetch_data.py; can also be run on its own.

Source: openworm/CElegansNeuroML -- one NeuroML file per neuron, generated from
the WormBase Virtual Worm model. Each file carries the soma's xyz position and
the full neurite path, in microns, in the coordinate frame of an actual worm.

That means the explorer can place neurons where they really sit in the animal
instead of in an abstract graph layout.

    .venv/bin/python fetch_positions.py        # ~300 small files, cached
"""

import concurrent.futures as futures
import json
import re
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data" / "neuroml"
OUT = ROOT / "data" / "neuron_positions.json"

TREE_URL = ("https://api.github.com/repos/openworm/CElegansNeuroML/git/trees/"
            "master?recursive=1")
RAW = ("https://raw.githubusercontent.com/openworm/CElegansNeuroML/master/"
       "CElegans/generatedNeuroML2")

# Files in that directory that are not neurons.
NOT_NEURON = re.compile(r"^(MDL|MDR|MVL|MVR|MVU|MDU|BWM|GLR|CEPsh|exc|hyp|int|"
                        r"mu_|vm|um|sph|g1|g2|hmc|.*Muscle)", re.I)

SEG = re.compile(
    r'<segment id="(\d+)"[^>]*name="([^"]*)"[^>]*>(.*?)</segment>', re.S)
PROX = re.compile(r'<proximal x="([-\d.eE]+)" y="([-\d.eE]+)" z="([-\d.eE]+)"')
DIST = re.compile(r'<distal x="([-\d.eE]+)" y="([-\d.eE]+)" z="([-\d.eE]+)"')
PARENT = re.compile(r'<parent segment="(\d+)"')


def list_cells():
    with urllib.request.urlopen(TREE_URL, timeout=60) as fh:
        tree = json.load(fh)
    names = []
    for entry in tree.get("tree", []):
        path = entry["path"]
        if not path.endswith(".cell.nml"):
            continue
        name = Path(path).name[: -len(".cell.nml")]
        if NOT_NEURON.match(name):
            continue
        names.append(name)
    return sorted(set(names))


def fetch_one(name):
    dest = CACHE / f"{name}.cell.nml"
    if dest.exists() and dest.stat().st_size > 200:
        return name, dest.read_text()
    try:
        with urllib.request.urlopen(f"{RAW}/{name}.cell.nml", timeout=60) as fh:
            text = fh.read().decode("utf-8", "replace")
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return name, None
    dest.write_text(text)
    return name, text


def parse(text):
    """Soma xyz plus the neurite as a list of straight SEGMENTS, in microns.

    Each NeuroML segment has its own proximal and distal point, and the morphology
    is a branching tree, not one path. An earlier version appended every proximal
    and distal point in document order and called it a polyline -- drawn as a
    connected line that zigzags across the animal, because consecutive segments in
    the file are not consecutive in space. Storing (proximal, distal) pairs and
    drawing each independently reproduces the real shape.
    """
    soma, segments, distal_of = None, [], {}

    for seg_id, seg_name, body in SEG.findall(text):
        prox = PROX.search(body)
        dist = DIST.search(body)
        parent = PARENT.search(body)
        if not dist:
            continue

        end = [round(float(v), 2) for v in dist.groups()]
        distal_of[seg_id] = end

        if prox:
            start = [round(float(v), 2) for v in prox.groups()]
        elif parent and parent.group(1) in distal_of:
            # NeuroML convention: a segment with no <proximal> starts where its
            # parent ended. Requiring an explicit proximal (an earlier version)
            # kept only the handful of branch roots -- three segments for AVAL
            # instead of sixty, so the neurite collapsed to a stub.
            start = distal_of[parent.group(1)]
        else:
            continue

        if soma is None and "soma" in seg_name.lower():
            soma = start
        if start != end:
            segments.append([start, end])

    if soma is None and segments:
        soma = segments[0][0]
    if soma is None:
        return None
    return {"soma": [round(v, 2) for v in soma], "segments": segments}


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    print("listing neuron files ...")
    names = list_cells()
    print(f"  {len(names)} neurons")

    print("downloading (cached after the first run) ...")
    out, failed = {}, []
    with futures.ThreadPoolExecutor(max_workers=16) as pool:
        for i, (name, text) in enumerate(pool.map(fetch_one, names), 1):
            if text is None:
                failed.append(name)
                continue
            parsed = parse(text)
            if parsed:
                out[name] = parsed
            else:
                failed.append(name)
            if i % 50 == 0:
                print(f"  {i}/{len(names)}")

    OUT.write_text(json.dumps(out))
    print(f"\nwrote {OUT}  ({OUT.stat().st_size / 1e6:.1f} MB)")
    print(f"  {len(out)} neurons with 3D positions")
    if failed:
        print(f"  {len(failed)} without usable geometry: {failed[:8]}")

    xs = [v["soma"][0] for v in out.values()]
    ys = [v["soma"][1] for v in out.values()]
    zs = [v["soma"][2] for v in out.values()]
    print(f"\n  extent (microns): x {min(xs):.0f}..{max(xs):.0f}  "
          f"y {min(ys):.0f}..{max(ys):.0f}  z {min(zs):.0f}..{max(zs):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
