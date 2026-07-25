#!/usr/bin/env python3
"""Download everything the explorer needs. Run once.

Three public datasets, about 50 MB in total:

  1. CeNGEN single-cell expression -- which genes are on in which neuron class
  2. The whole-animal connectome (Cook et al. 2019) -- who wires to whom
  3. 3D geometry from the WormBase Virtual Worm model -- where each of the 302
     neurons sits, and the shape of its neurites

    python fetch_data.py           # skips anything already downloaded
    python fetch_data.py --force   # re-download everything
"""

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from wormview import connectome as cn                      # noqa: E402
from wormview.data import DATA_DIR, THRESHOLD_FILES        # noqa: E402

CENGEN_BASE = "https://cengen.org/storage"

EXTRA = {
    # functional group annotations (sensory / interneuron / motor)
    "all_cell_info.csv":
        "https://raw.githubusercontent.com/openworm/ConnectomeToolbox/"
        "main/cect/data/all_cell_info.csv",
}


def get(url, dest, force=False):
    if dest.exists() and dest.stat().st_size > 200 and not force:
        print(f"  have     {dest.name}")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  getting  {dest.name} ...", end=" ", flush=True)
    urllib.request.urlretrieve(url, dest)
    size = dest.stat().st_size
    print(f"{size / 1e6:.1f} MB" if size > 1e6 else f"{size / 1e3:.0f} KB")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="re-download everything")
    args = ap.parse_args()

    print("\n1. gene expression  (CeNGEN, Taylor et al. 2021)")
    for name in THRESHOLD_FILES.values():
        get(f"{CENGEN_BASE}/{name}", DATA_DIR / name, args.force)

    print("\n2. wiring  (Cook et al. 2019 connectome; Witvliet et al. 2021 contact)")
    get(cn.EDGE_URL, cn.DATA_DIR / cn.EDGE_FILE, args.force)
    for name, url in EXTRA.items():
        get(url, cn.DATA_DIR / name, args.force)

    print("\n3. 3D geometry  (WormBase Virtual Worm, via openworm/CElegansNeuroML)")
    positions = ROOT / "data" / "neuron_positions.json"
    if positions.exists() and not args.force:
        print(f"  have     {positions.name}")
    else:
        import fetch_positions
        fetch_positions.main()

    print("\nDone. Start the explorer with:  python app.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
