# v6.0.0 — Architecture Upgrade, Publication Renderer Preserved

- Promoted the application to PanViz 6.0.0 while retaining the approved publication renderer and editor geometry unchanged from baseline commit `237db30af8639379ed3976657f9bbcce28110725`.
- Added a canonical scientific-record layer that preserves all eight PLIP interaction classes and exports unified CSV/JSON records.
- Added a presentation-independent SHA-256 scientific-record signature to project manifests.
- Fixed PDBQT complex assembly so ligand atom serials begin after the receptor maximum serial.
- Reused one temporary workspace per Streamlit session instead of creating a new directory on every rerun.
- Centralized version loading through `VERSION.txt` / `panviz_version.py`.
- Made PyMOL opt-in for the command-line workflow.
- Added automated v6 architecture tests alongside the golden publication-renderer regression guard.
- The renderer scene schema remains separately recorded for compatibility; application version changes do not imply visual-renderer changes.

# v5.8.7 — Fully Editable Legend

- Legend entries are first-class editable content: add from a dropdown, remove, reorder, rename, recolor, restyle line width/dash, toggle visibility, and edit legend/background/text styling.
- Adding a standalone graphical bond automatically ensures a **Bond** legend entry. Adding a custom interaction automatically ensures the corresponding interaction-type legend entry.
- Removing a legend entry affects only presentation; it never deletes or changes the underlying interaction or molecular structure.
- Legend sizing recalculates dynamically as entries and labels change.

# Changelog

## 5.8.7
- Professional Word-style object selection and keyboard control layer.
- Original molecular structure remains strictly immutable.
- Molecular covalent bonds can be clicked only for copy-as-graphic; copied bonds are independent graphical objects.
- Added standalone graphical covalent bonds with endpoint drag editing.
- Added multi-select, group copy/paste/cut/duplicate, nudging, rotation, resize handles, text shortcuts, and z-order controls.
- Added dynamic bond/interaction expand and compress controls.
- Added contextual bond order/color/width/dash/direction controls.

# PanViz Interactive v5.8.7

## v5.8.7 — molecular structure hard lock
- Restored the original covalent molecular renderer directly from immutable atom coordinates and bond topology.
- Removed active rendering and interaction of legacy graphical bond-edit objects.
- Molecular atoms and covalent bonds cannot be selected, moved, resized, rotated, restyled, copied, pasted, duplicated, deleted, or newly created.
- Disabled canvas-resize controls so molecular coordinates cannot be changed through scene scaling.
- Legacy `bondEdits` are discarded when migrating or loading layouts.
- Presentation annotations remain independently editable.

# PanViz Interactive v5.8.2

## Independent figure editing
- Existing covalent bonds get independent graphical twins derived from immutable chemical topology.
- Bonds and non-covalent interactions support direct endpoint dragging, click-based expand/compress, copy/paste/duplicate, rotation, restyling, and deletion.
- Added creation of new bond, interaction, and distance annotations.
- Existing scientific PLIP measurements remain preserved separately from graphical edits.

## v5.8.1 — Editor control completion

- Added keyboard copy/duplicate/paste controls for editable annotations (`Ctrl+C`, `Ctrl+D`, `Ctrl+V`).
- Added non-covalent interaction `Expand` and `Min` length controls centered on the interaction midpoint.
- Added an absolute direction-angle control for selected non-covalent interactions while retaining 15° rotation buttons.
- Fixed rotated non-covalent interaction endpoint geometry so non-zero direction/rotation preserves both endpoints correctly.
- Scientific PLIP interaction and distance records remain protected from duplication.
- Updated application/version metadata to v5.8.1.

## v5.8.0 — Complete workflow release

- Reuses one completed PLIP analysis across original PNG/SVG rendering and the interactive editor.
- Added complete project ZIP export with source inputs, prepared PDB, interaction CSVs, figures, initial layout, manifest, and project README.
- Added Vina-score-aware docking pose selection when `REMARK VINA RESULT` records are present.
- Added layout validation and reload while preserving locked molecular topology and scientific interaction data.
- Added editable free-text and arrow annotations.
- Added 1×/2×/4×/6× PNG export scaling and keyboard shortcuts for save/undo/redo.
- Corrected editor text-annotation Apply wiring and strengthened live-style undo behavior.
- Standardized release versioning and documentation to v5.8.0.

## v5.7.13

- Redesigned Docking PDBQT inputs into exactly two cards.
- Each card now combines the role/instructions and its uploader in the same visual container.
- Removed the previous separate four-card presentation (two information cards + two uploader areas).
- Protein / Receptor remains the left card; Ligand / Docking Poses remains the right card.

## v5.7.12

- Matched edited PNG export to the original PanViz SVG/PNG canvas dimensions.
- Removed the 4× edited PNG raster path that could make exported stroke hierarchy differ from the original figure.
- Kept high-quality browser image smoothing for SVG-to-PNG rasterization.
- SVG export remains vector.

## v5.7.11
- Automatic duplicate-interaction weighting for shared residue–ligand-atom connections.

## v5.7.11
- Residue/3D bubble diameter is global across all residues.
- Bubble size control remains in the normal Styles tab and uses 1 px steps.
- Removed per-residue bubble-size editing from the contextual ribbon.
- Existing per-residue bubble-size overrides are normalized back to the global value.


## v5.7.9
- Added vector-preserving PDF export via a print-optimized SVG document.
- Added print CSS to remove browser page chrome and fit the current canvas to the PDF page.

## v5.7.11 — Adaptive Canvas Fit

- Increasing/decreasing canvas dimensions proportionally scales the complete figure scene and then fits it inside the page with a publication-safe margin.
- Canvas viewport uses responsive width with no internal horizontal/vertical scrollbars.
- Streamlit editor iframe is non-scrolling.


## v5.7.7
- Changed the global default residue presentation from Plain Text to 3D glass Bubble.
- The Styles tab still offers Plain as an alternative; it is no longer the initial/default shape.
- The 3D bubble remains independently size-adjustable and uses RGB (10, 255, 239) / #0AFFEF by default.
## v5.7.4 — Windows PDBQT temp-file lock fix

- Fixed `WinError 32` during docking PDBQT pose preparation by explicitly closing the Windows file descriptor returned by `tempfile.mkstemp()` before Open Babel access and cleanup.
- Temporary pose files are now safely removable after conversion.

# PanViz Interactive v5.7.3

## v5.7.3 export-quality correction

- Fixed faint/light edited PNG exports caused by SVG `vector-effect=non-scaling-stroke` on locked molecular bonds during 4× rasterization. Exported SVG now removes non-scaling stroke behavior so bond widths scale correctly with the exported artwork.
- Export SVG uses `shape-rendering: geometricPrecision` for cleaner vector rendering.
- Removed the global non-covalent color control from the toolbar; the publication cyan `RGB(10,255,239) / #0AFFEF` is reserved for the 3D bubble style.
- Original covalent bond color remains independent from the 3D bubble color.


## v5.7.3 glass-bubble refinement

- Upgraded the residue **3D bubble** to a glossy, translucent, STRING-style spherical treatment with radial shading, specular highlight, rim highlight, and soft shadow.
- Bubble diameter is adjustable globally (**24–140 px**) from Styles and individually on a selected bubble.
- Bubble size is preserved in the editable SVG/layout state.
- The underlying residue identity and PLIP interaction data remain unchanged.


- Added dual PDB / Docking PDBQT input workflow with receptor + ligand pose selection.
- Made 2.0 px covalent and 3.0 px non-covalent widths the global base defaults shared by the original renderer and editable scene.
- Reserved RGB(10,255,239) / #0AFFEF as the 3D glass-bubble default; covalent bond color remains independent.
- Removed the separate Publication preset control; the Styles tab remains the normal override interface.
- Added 4× PNG export while preserving SVG as the primary vector export.
- Fixed ligand atom labels to correctly follow the covalent bond color when the linked mode is selected.

# PanViz Interactive v5.6

- Compact STRING-inspired 3D bubble residue nodes; node diameter is independent of font size.
- Added polished editor/app styling.
- Added dynamic covalent bond color and ligand color mode (element colors or follow bond color).
- Fixed global non-covalent interaction width precedence; type-specific width is explicitly custom only when applied.
- Added a reset-to-global-width control for interaction types.
- Selection bounds now match circle/bubble residue geometry.
- Preserved immutable molecular topology and PLIP scientific data.

## v5.7.3 functional PDBQT patch
- Hardened Docking PDBQT preparation for receptor + docked-ligand workflows.
- Validates that converted receptor and ligand structures contain ATOM/HETATM records.
- Remaps Open Babel `CONECT` records to rewritten ligand atom serials.
- Writes remapped `CONECT` records before `TER/END` so ligand connectivity is retained.
- Verified multi-pose MODEL parsing, pose selection/bounds, complex construction, empty-structure validation, Python syntax, and editor JavaScript syntax.

## v5.8.0 — Complete reusable analysis/editor workflow

- Reused one completed PLIP analysis for the original PNG render, original SVG render, and interactive editor construction.
- Added persistent Streamlit result state so generated results remain available until the selected input/pose/site/canvas parameters change.
- Added Vina pose-score display when `REMARK VINA RESULT` records are present.
- Added complete project ZIP download containing original inputs, prepared structure, figures, interaction CSVs, initial layout, and manifest.
- Added portable `.panviz.json` layout loading with site and ligand atom-fingerprint validation.
- Added free text annotations with move, rotate, style, delete, undo/redo, save, load, and export support.
- Added arrow annotations with move, rotate, color, width, dash style, delete, save, load, and export support.
- Added 1×/2×/4×/6× PNG export scale controls while retaining vector SVG export.
- Added Ctrl/Cmd+S, Ctrl/Cmd+Z, and Ctrl/Cmd+Y editor shortcuts.
- Wired the text-edit popup Apply button to the underlying text mutation function.
- Fixed live style-control undo behavior for opacity controls.
- Updated version identifiers and documentation to v5.8.0.
