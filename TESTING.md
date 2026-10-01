# PanViz Interactive v6.0.0 — Testing and validation

## v6.0.0 independent-engine gates

PanViz 6 removes the active legacy `utils.py` and `interactive_engine.py`
implementations and routes both the Streamlit app and CLI through
`panviz_engine.py`.

Automated CI now checks:

- Python compilation for `app.py`, `PanViz.py`, `panviz_engine.py`, and `cloud_diagnostic.py`.
- all eight PLIP interaction classes are present in the canonical model:
  HPI, HB, WB, SB, PS, PC, XB, and MC.
- the removed legacy function names are absent from the active engine.
- `DetermineBonds(..., charge=0)` and `set_to_neutral_pH()` are absent from the depiction pipeline.
- standalone SVG serialization contains molecule, interaction, residue, distance, and legend content.
- CairoSVG produces a non-empty PNG from the same canonical scene.

### Independent-implementation overlap check

A direct normalized-line comparison of PanViz 6 `panviz_engine.py` against
`cmwoodley/plip_to_2D/utils.py` found:

- shared upstream function names: **0**
- exact nonblank/non-comment lines: 104 / 1108 PanViz normalized lines; these are overwhelmingly imports, punctuation, `try`, `continue`, and other boilerplate.
- after excluding imports/boilerplate and retaining substantive lines of at least 20 characters: **1 exact line out of 637 PanViz lines (~0.16%)**.
- the one substantive exact line was the generic Python/RDKit loop `for atom in mol.GetAtoms():`.

This comparison is a development check, not a legal test of copyright status.
Historical lineage remains documented in `PROVENANCE.md`.

## Historical v5 validation notes

## v5.8.7 molecular hard-lock checks
- Molecular covalent bonds are rendered directly from immutable `state.bonds`.
- No active bond hit targets or endpoint-edit handles are created.
- The molecule layer uses `pointer-events:none` and cannot be selected or dragged.
- `bondEdits` is discarded during scene migration and layout loading.
- The Add Bond control is removed.
- Molecular style controls no longer mutate molecular state.
- Annotation selection/editing remains available independently.

## v5.7.9 PDF export check
- `Export PDF` button is present.
- Export uses the current `sceneSvg()` output rather than a raster screenshot.
- Print document uses `@page { margin: 0 }` and responsive SVG fitting.
- Browser print dialog is used to save the vector figure as PDF.

## v5.7.11 — Adaptive Canvas Fit

- Increasing/decreasing canvas dimensions proportionally scales the complete figure scene and then fits it inside the page with a publication-safe margin.
- Canvas viewport uses responsive width with no internal horizontal/vertical scrollbars.
- Streamlit editor iframe is non-scrolling.

## General validation notes

Historical static checks completed for the v5 codebase:
- Python syntax compilation was performed for the then-active v5 modules.
- Embedded editor JavaScript syntax check with Node.js.

Interactive editor smoke test completed with Chromium + Playwright using a representative PanViz scene:
- molecule layer present once
- residue selection and contextual ribbon
- residue font-size change
- residue dragging without creating extra annotation nodes
- residue rotation
- canvas width adjustment
- covalent bond visual-width adjustment
- legend selection and dragging
- distance-label selection and rotation
- SVG export contains the legend and excludes editor hit targets/UI
- PNG export produces a non-empty image
- no browser console/page errors during the test

The complete Streamlit + PLIP runtime was not executed in this build environment. The editor interaction layer was tested independently with a representative scene.


v5.6 additional checks:
- Compact bubble sizing is independent of residue font size.
- Global interaction-width precedence and type-specific custom width behavior.
- Dynamic molecule bond color and follow-bond ligand color mode.
- Shape-aware residue selection geometry.
- Editor chrome/style panel visual structure.

## v5.6 controlled editor smoke test
A representative scene was exercised in Chromium/Playwright. Verified:
- 3D bubble node diameter resolves to a compact ~49 px at a 46 px node-size setting.
- Global non-covalent width changes visible interaction strokes (2.0 → 4.0 px).
- Per-type custom H-bond width changes to 5.0 px.
- Resetting a type to Global returns it to the global 4.0 px width.
- Bond color changes dynamically; with "Follow bond color (linked)", displayed ligand atom labels follow the chosen molecule color.
- Residue selection box remains shape-aware and stays around the compact bubble.
- Residue drag leaves a single residue object; no duplicate label node is created.
- Legend remains present and is included in exported SVG.
- No browser console/page errors in this controlled editor test.

The complete Streamlit + PLIP runtime was not executed in this build environment because Streamlit and PLIP are not installed here. The Python source was syntax-checked successfully.

### v5.7.3 style-mapping checks
- Verified `#0AFFEF` is reserved for the 3D bubble style and is exposed as `bubbleColor`.
- Verified the global base default preserves the original covalent bond color while the 3D bubble uses #0AFFEF.
- Verified original semantic interaction palette is restored/migrated for non-custom interactions.
- Verified legend defaults to bold, larger typography with increased padding.
- Verified embedded editor JavaScript passes `node --check`.
- Verified Python source files pass `py_compile`.
- Full PLIP runtime import was not available in this validation environment because `plip` is not installed.

### v5.7.6 global-default checks
- Confirmed there is no Publication preset button or preset mode.
- Confirmed 2.0 px covalent and 3.0 px non-covalent widths are base defaults in both the original renderer and editable scene.
- Confirmed per-interaction semantic colors remain the default interaction colors unless explicitly customized.
- Confirmed RGB(10,255,239) / #0AFFEF is reserved for the editable 3D glass-bubble residue style.
- Confirmed the Styles tab remains the normal style override interface.
- Confirmed the original renderer uses bold, enlarged legend typography for legibility.

### v5.7.11 duplicate-interaction test
- Synthetic interaction rows with the same residue and ligand atom are assigned `multiplicity >= 2`.
- Editor width calculation uses `baseWidth × multiplicity`.
- Existing scientific interaction and distance records remain separate.

### Backward-compatible multiplicity hydration
- Existing interactions are regrouped by `residueId + anchorAtom` during editor load.
- This restores automatic thickening even when an older saved layout has no explicit `multiplicity` field.

### v5.8.0 completion checks
- Python syntax check passes for `app.py`, `interactive_engine.py`, `utils.py`, and `PanViz.py`.
- Embedded editor JavaScript passes `node --check`.
- Browser smoke test verifies one locked molecule layer and one interaction layer.
- Browser smoke test verifies creation and editing of free text annotations.
- Browser smoke test verifies creation of arrow annotations.
- Browser smoke test verifies `.panviz.json` export contains notes, arrows, and scientific data.
- Browser smoke test verifies `.panviz.json` reload restores saved annotation state after a later modification.
- Browser smoke test verifies SVG export excludes editor-only hit targets/UI.
- Browser smoke test verifies default 4× PNG export produces a 3600×2600 image for a 900×650 test scene.
- Browser smoke test reports no page or console errors.
- Full Streamlit + PLIP runtime was not executed in this build environment because Streamlit, PLIP, and Open Babel are not installed here; the source has been syntax-checked and the browser editor has been exercised independently with a representative scene.


### v5.8.1 control checks
- Embedded editor JavaScript passes `node --check`.
- Python source files pass `py_compile`.
- Keyboard handler supports `Ctrl+C`, `Ctrl+D`, and `Ctrl+V` without intercepting active text/form fields.
- Editable annotation duplication preserves the locked molecular and scientific PLIP data model.
- Non-covalent `Expand` and `Min` operations scale endpoints about the interaction midpoint.
- Direction control changes the interaction angle without adding arrowheads or changing PLIP interaction semantics.


## v5.8.2 independent-editing smoke test
Validated in Chromium with a synthetic scene: existing bond rendering, new bond creation, draggable endpoint handles, click-based expand/compress, Ctrl+D bond duplication, Ctrl+D interaction-to-independent-custom duplication, Add interaction, Ctrl+C/Ctrl+V interaction copy/paste, Add distance and text editing, and layout save trigger. Python modules compile with `py_compile`; editor JavaScript passes `node --check`.


## v5.8.7 professional object-editor checks
- Python source files pass `py_compile` for `app.py`, `interactive_engine.py`, `utils.py`, and `PanViz.py`.
- Embedded editor JavaScript passes `node --check`.
- Original molecular atoms and covalent bonds remain byte-for-byte unchanged after attempted source-bond drag and all presentation editing operations.
- Original covalent bonds expose zero endpoint-edit handles and can only be selected as copy sources.
- `Ctrl+C` + `Ctrl+V` on an original covalent bond creates a standalone graphical bond with no `a`, `b`, `sourceBondId`, or `scientificSourceId` linkage.
- Standalone bond endpoint dragging changes its graphical coordinates and does not change the molecular topology.
- Bond order, styling, dynamic expand/compress, and 15° keyboard rotation are functional on standalone graphical bonds.
- New graphical bonds and new non-covalent interaction annotations can be created independently.
- Copied residue labels have `sourceResidue=null`, `residueId=null`, and `leaderVisible=false`, so they are detached from the original scientific residue annotation.
- Multi-selection, group copy/paste, keyboard nudge, cut/paste, duplicate, delete, undo, and redo were exercised in Chromium.
- Cut removes independent objects; first paste restores the cut object and subsequent paste creates a copy.
- SVG export excludes editor hit targets and UI.
- No browser console/page errors were reported in the professional object-editor smoke test.
- Full Streamlit + PLIP runtime remains dependent on the packages listed in `requirements.txt`; those runtime packages are not installed in this build environment.
