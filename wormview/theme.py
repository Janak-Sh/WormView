"""Visual theme for the explorer -- light, print-safe, colourblind-checked.

The four group colours are the first three slots of a validated categorical
palette plus a neutral grey, so no two of them collapse together under common
colour-vision deficiencies.
"""

BG = "#fcfcfb"            # page and scene background
PANEL = "#f5f5f3"         # side panel
INK = "#16181d"
INK_DIM = "#6b6f7a"
LINE = "#e2e2de"

GROUP_COLOR = {
    "sensory": "#2a78d6",       # blue
    "interneuron": "#eb6834",   # orange
    "motor": "#1baf7a",         # aqua
    "other": "#9aa0ac",         # muted grey
}

HILITE = "#d4a017"                        # the selected neuron
SYNAPSE = "#7d3cc4"                       # synapse contact points
BODY = [[0, "#9fb6d8"], [1, "#c8d6e8"]]   # body wall

# Sequential ramp for gene expression: one light-to-dark hue, then hot at the top.
# The low end is deliberately darker than the background so "not detected" neurons
# stay visible instead of vanishing into the page.
GENE_SCALE = [[0.0, "#c9cbc8"], [0.18, "#9ec5f4"], [0.42, "#3987e5"],
              [0.68, "#1c5cab"], [0.86, "#eb6834"], [1.0, "#a8340c"]]

# side-panel chips
CHIP_BG = "#ececea"
CHIP_WATCH = "#fdeccf"
CHIP_INK = "#24262c"

FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
