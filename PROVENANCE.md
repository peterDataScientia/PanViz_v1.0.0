# PanViz provenance

## Historical lineage

Early PanViz prototypes incorporated and modified portions of the open-source
`cmwoodley/plip_to_2D` project by Christopher Woodley. That upstream project
is distributed under GPL-3.0 and uses PLIP, RDKit, Open Babel, PyMOL and Cairo
to generate static 2D protein-ligand interaction figures.

PanViz 6.0.0 replaces the active legacy rendering implementation with a new
`panviz_engine.py` architecture. The former `utils.py` and
`interactive_engine.py` implementations were removed from the active branch.

The PanViz 6 engine uses PLIP as the scientific interaction detector, but its
canonical interaction model, ligand-graph adapter, radial constraint layout,
scene construction, interaction routing, browser editor integration and static
SVG/PNG serialization are implemented in PanViz.

## Scientific and software dependencies

PanViz 6 depends on and should cite the relevant upstream projects:

- PLIP — protein-ligand interaction detection and report data.
- RDKit — molecular graph handling and deterministic 2D depiction.
- Open Babel — molecular format interpretation/conversion.
- CairoSVG / Cairo — SVG rasterization for static PNG export.
- Streamlit — web application framework.

## Historical acknowledgement

Removal of the old implementation does not erase the project's history.
The conceptual influence of `plip_to_2D` on early PanViz development should
remain acknowledged in publications and software documentation.

This file records provenance; it is not a substitute for the license terms of
PanViz or any dependency.
