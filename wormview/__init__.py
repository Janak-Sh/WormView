"""wormview -- an interactive 3D explorer for the C. elegans nervous system.

Joins three public datasets into one browsable picture:

  * real 3D geometry of all 302 neurons (soma positions and full neurite
    morphology) from the WormBase Virtual Worm model, via NeuroML;
  * the whole-animal connectome (Cook et al. 2019);
  * single-cell gene expression for every neuron class (CeNGEN).

Modules:
  data        loading the CeNGEN expression matrix
  connectome  the wiring diagram, and mapping neuron names onto CeNGEN classes
  anatomy     3D positions, neurite morphology, body wall, declutter
  genes       the watch-list of genes highlighted in the panel
  figure      the Plotly 3D scene
  panel       the side panel
  theme       colours and fonts
"""

__version__ = "1.0.0"
