# PanViz Publication Renderer Contract

The current molecule, residue, interaction, distance-label and legend appearance is a protected PanViz asset.

## Protected visual behavior

The following must remain visually unchanged unless a deliberate publication-design revision is approved:

- ligand/molecule bond depiction and stroke hierarchy;
- heteroatom typography and element colors;
- residue label typography, white knockout and placement;
- residue leader behavior;
- non-covalent interaction routing, endpoint offsets, colors, dash pattern and width;
- interaction-distance placement and magenta typography;
- legend ordering, sample style, typography, spacing and placement;
- canvas fitting, whitespace balance and publication margins.

Refactoring, renaming, moving code between modules, improving scientific-data handling, adding tests, fixing PDBQT processing, improving exports and changing non-visual architecture are allowed **only if the approved rendered appearance remains equivalent**.

## Baseline

The protected baseline is Git commit:

`237db30af8639379ed3976657f9bbcce28110725`

A safety branch was created:

`publication-renderer-baseline-2026-10-01`

At that baseline:

- `utils.py::_draw_mol` defines the static publication renderer.
- `interactive_engine.py::_exact_panviz_layout` reproduces its established residue, interaction and distance geometry for the editor.

## Regression strategy

Automated tests use two complementary safeguards:

1. **Golden geometry** — a representative synthetic scene is run through the layout engine and residue positions, interaction endpoints, distance labels and font sizes are compared to an approved snapshot.
2. **Visual style contract** — publication stroke widths, interaction offsets, collision-search geometry, font scaling, colors and distance-label styling are checked explicitly.

These checks protect output behavior rather than the names of implementation functions. If the code is reorganized or renamed, the tests may be adapted to the new API, but the approved visual output must continue to satisfy the same contract.

## Development principle

> Improve PanViz around the renderer; do not casually redesign the renderer.

A visual change should be treated as a scientific/publication design decision, not as a side effect of code cleanup.
