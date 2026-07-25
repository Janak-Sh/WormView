"""Loading CeNGEN thresholded expression matrices.

CeNGEN published four matrices at increasing stringency. In each, a value of 0
means "not detected in this neuron class at this threshold"; any value above 0 is
the average TPM for that gene in that class. So these files carry *both* a
presence call and an expression level -- which matters, because several of the
lab's genes are detected almost everywhere and only their level is informative.

Source: https://www.cengen.org/downloads/
Reference: Taylor et al. 2021, Cell -- the CeNGEN L4 hermaphrodite atlas.
"""

from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

THRESHOLD_FILES = {
    1: "021821_liberal_threshold1.csv",
    2: "021821_medium_threshold2.csv",
    3: "021821_conservative_threshold3.csv",
    4: "021821_stringent_threshold4.csv",
}

THRESHOLD_LABELS = {
    1: "1 - liberal",
    2: "2 - medium",
    3: "3 - conservative",
    4: "4 - stringent",
}

_META_COLS = ("gene_name", "Wormbase_ID")


class MissingData(RuntimeError):
    pass


def load_matrix(threshold=2):
    """Return (tpm, neuron_classes) for one threshold.

    tpm is a DataFrame indexed by gene_name with one column per neuron class.
    Values are average TPM; 0 means not detected at this threshold.
    """
    if threshold not in THRESHOLD_FILES:
        raise ValueError(f"threshold must be one of {sorted(THRESHOLD_FILES)}")

    path = DATA_DIR / THRESHOLD_FILES[threshold]
    if not path.exists():
        raise MissingData(
            f"{path} not found.\nRun:  python run_atlas.py --fetch-data"
        )

    raw = pd.read_csv(path, index_col=0)
    neurons = [c for c in raw.columns if c not in _META_COLS]

    # A handful of gene names appear more than once (different WormBase IDs
    # mapping to the same name). Keep the row with the strongest signal so a
    # duplicate can never silently mask real expression.
    raw = raw.assign(_total=raw[neurons].sum(axis=1))
    raw = raw.sort_values("_total", ascending=False).drop_duplicates("gene_name")

    tpm = raw.set_index("gene_name")[neurons].astype(float)
    return tpm, neurons


def load_all_thresholds():
    """{threshold: DataFrame} for every threshold file present on disk."""
    out = {}
    for t in sorted(THRESHOLD_FILES):
        try:
            out[t], _ = load_matrix(t)
        except MissingData:
            continue
    if not out:
        raise MissingData(
            "No CeNGEN files found.\nRun:  python run_atlas.py --fetch-data"
        )
    return out


def subset_genes(tpm, gene_names):
    """Rows for the requested genes, preserving the requested order.

    Returns (subset, missing) so a typo or renamed gene is reported rather than
    silently dropped.
    """
    present = [g for g in gene_names if g in tpm.index]
    missing = [g for g in gene_names if g not in tpm.index]
    return tpm.loc[present], missing
