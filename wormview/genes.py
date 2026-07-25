"""Genes to highlight in the explorer.

A short watch-list of cell adhesion, extracellular matrix and synapse-formation
genes -- the molecules that guide a growing neurite and help decide which partner
it connects to. They make interesting things to paint across a nervous system.

Nothing here restricts what the explorer can show: the gene dropdown covers all
~13,700 genes in the CeNGEN atlas. This list only decides which genes get called
out in the side panel, and it is meant to be edited.
"""

# gene -> what it is, in one line
WATCHLIST = {
    # --- cell adhesion molecules ------------------------------------------
    "fmi-1":  "Flamingo/CELSR cadherin - follower axons navigate along pioneers",
    "cdh-4":  "Fat-like cadherin - neuroblast migration, acts non-cell-autonomously",
    "ptp-3":  "LAR receptor phosphatase - axon guidance and synapse formation",
    "syg-1":  "Immunoglobulin superfamily - marks where a synapse should form",
    "syg-2":  "SYG-1's binding partner (a matched pair, not self-sticking)",
    "efn-4":  "Ephrin - axon outgrowth and branching",

    # --- basement membrane / extracellular matrix -------------------------
    "nid-1":  "Nidogen - basement membrane, organises neuromuscular junctions",
    "cle-1":  "Type XVIII collagen / endostatin - migration and axon guidance",
    "unc-52": "Perlecan - works with nidogen in dendrite morphogenesis",

    # --- calcium-dependent synapse formation ------------------------------
    "unc-2":  "Voltage-gated calcium channel, alpha subunit",
    "unc-36": "Voltage-gated calcium channel, alpha-2/delta subunit",
    "unc-43": "CaMKII - calcium-dependent kinase",

    # --- axon termination and patterning ----------------------------------
    "rpm-1":  "Ubiquitin ligase - axon termination and synaptogenesis",
    "egl-5":  "Hox transcription factor - selector for posterior neurons",
    "mab-5":  "Hox transcription factor - posterior patterning",
    "lin-44": "Wnt ligand - secreted from the tail; expect little neuronal signal",
}

# Genes whose product is secreted, or which act from a different cell than the one
# they affect. CeNGEN sequenced NEURONS ONLY, so for these a low or absent neuronal
# signal is uninformative: the source tissue was never sampled. Reading absence as
# "not involved" would be a real error, so flag them wherever they appear.
NON_AUTONOMOUS = {
    "nid-1":  "secreted basement membrane protein, made largely by muscle and "
              "hypodermis",
    "cle-1":  "secreted collagen XVIII, with substantial non-neuronal sources",
    "unc-52": "secreted perlecan, made mainly by muscle",
    "lin-44": "secreted Wnt, from tail hypodermis rather than neurons",
    "cdh-4":  "known to act cell-non-autonomously",
    "syg-2":  "in the canonical model, acts from a non-neuronal guidepost cell",
    "fmi-1":  "acts partly from the pioneer neuron, so the relevant cell may not "
              "be the affected one",
}


def describe(gene):
    """One-line description, or None if the gene is not on the watch-list."""
    return WATCHLIST.get(gene)


def autonomy_caveat(gene):
    """Why absence-in-neurons cannot be trusted for this gene, or None."""
    return NON_AUTONOMOUS.get(gene)
