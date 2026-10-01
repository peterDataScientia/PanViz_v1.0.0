import inspect
import json
from pathlib import Path

import pytest

import interactive_engine
import utils


GOLDEN = json.loads(
    (Path(__file__).parent / "fixtures" / "publication_layout_golden.json").read_text(
        encoding="utf-8"
    )
)


DATA_POINTS = [
    (-1.5, 0.0, "C", "C1"),
    (0.0, 0.0, "O", "O1"),
    (1.5, 0.0, "N", "N1"),
    (0.0, 1.5, "centroid", "ring1"),
    (-3.0, 2.0, "residue", "ASP34_A"),
    (3.0, 2.5, "residue", "TYR192_A"),
    (0.0, -3.0, "residue", "GLY216_A"),
]

CONNECTIONS = [
    (0, 1, "SINGLE"),
    (1, 2, "DOUBLE"),
    (1, 4, "HB", 2.80, 0),
    (3, 5, "PS", 4.72, 1),
    (2, 6, "HPI", 3.61, 2),
]


def approx(value, expected, tolerance=0.5):
    assert float(value) == pytest.approx(float(expected), abs=tolerance)


def test_publication_geometry_matches_approved_golden_layout():
    layout = interactive_engine._exact_panviz_layout(
        DATA_POINTS,
        CONNECTIONS,
        GOLDEN["canvas"]["width"],
        GOLDEN["canvas"]["height"],
    )

    assert layout["residue_font_size"] == GOLDEN["font_sizes"]["residue"]
    assert layout["distance_font_size"] == GOLDEN["font_sizes"]["distance"]

    actual_residues = layout["residue_draw"]
    assert [x["text"] for x in actual_residues] == [
        x["text"] for x in GOLDEN["residues"]
    ]
    for actual, expected in zip(actual_residues, GOLDEN["residues"]):
        for key in ("x", "y", "anchorX", "anchorY"):
            approx(actual[key], expected[key])
        assert bool(actual["leader"]) is bool(expected["leader"])

    actual_interactions = layout["interaction_lines"]
    assert len(actual_interactions) == len(GOLDEN["interactions"])
    for actual, expected in zip(actual_interactions, GOLDEN["interactions"]):
        assert actual["index"] == expected["index"]
        assert actual["row_idx"] == expected["row_idx"]
        assert actual["type"] == expected["type"]
        for key in ("x1", "y1", "x2", "y2", "distance"):
            approx(actual[key], expected[key])

    actual_distances = layout["distance_lines"]
    for expected in GOLDEN["distance_labels"]:
        actual = actual_distances[expected["index"]]
        assert actual["type"] == expected["type"]
        assert actual["distanceText"] == expected["text"]
        approx(actual["distanceX"], expected["x"])
        approx(actual["distanceY"], expected["y"])


def test_static_publication_renderer_style_contract_is_unchanged():
    """Guard the approved publication styling while allowing code refactors.

    These are visual behavior constants, not function-name checks. If the renderer
    is reorganized, this test may move with it, but the approved values must remain.
    """
    source = inspect.getsource(utils._draw_mol)

    required_visual_rules = [
        "residue_font_size = max(17, min(23, int(base * 0.025)))",
        "atom_font_size = max(14, min(19, int(base * 0.019)))",
        "distance_font_size = max(14, min(18, int(base * 0.018)))",
        "single_width = 2.0",
        "double_width = 2.0",
        "triple_width = 2.0",
        "interaction_width = 3.0",
        'ctx.set_source_rgb(1, 1, 1)',
        'for radial in (20, 28, 38, 50, 64, 80, 98):',
        'for tangential in (0, 8, -8, 16, -16):',
        'for frac in (0.50, 0.43, 0.57, 0.35, 0.65, 0.25, 0.75):',
        'for off in (0, 4, -4, 8, -8, 14, -14, 20, -20):',
        '"HB": 6',
        '"PS": 10',
        '"PC": -10',
        '"SB": -5',
        '"WB": 4',
        '"XB": -8',
        '"MC": 8',
        'ctx.set_source_rgb(0.82, 0.0, 0.58)',
    ]
    for rule in required_visual_rules:
        assert rule in source, f"Publication renderer visual rule changed: {rule}"


def test_editor_geometry_uses_same_approved_visual_vocabulary():
    assert interactive_engine.INTERACTION_COLORS == {
        "HPI": "#595959",
        "HB": "#0000E0",
        "PS": "#008500",
        "PC": "#E08000",
        "SB": "#D900B0",
        "WB": "#1596B8",
        "XB": "#7A5CC7",
        "MC": "#A45700",
    }
    assert interactive_engine.BUBBLE_COLOR == "#0AFFEF"
    assert interactive_engine.NONCOVALENT_COLOR == "#595959"


def test_visual_renderer_files_exist_and_legacy_geometry_is_active():
    assert Path("utils.py").exists()
    assert Path("interactive_engine.py").exists()
    assert hasattr(utils, "_draw_mol")
    assert hasattr(interactive_engine, "_exact_panviz_layout")


def test_all_eight_plip_interaction_classes_are_supported():
    assert interactive_engine.INTERACTION_TYPES == {
        "HPI", "HB", "WB", "SB", "PS", "PC", "XB", "MC"
    }
    assert interactive_engine.INTERACTION_LABELS["WB"] == "Water bridge"
    assert interactive_engine.INTERACTION_LABELS["XB"] == "Halogen bond"
    assert interactive_engine.INTERACTION_LABELS["MC"] == "Metal coordination"


def test_publication_legend_has_fixed_scientific_order_and_formal_labels():
    assert interactive_engine.INTERACTION_ORDER == (
        "HPI", "HB", "WB", "SB", "PS", "PC", "XB", "MC"
    )
    assert interactive_engine.INTERACTION_LABELS["HPI"] == "Hydrophobic contact"
    assert interactive_engine.INTERACTION_LABELS["HB"] == "Hydrogen bond"


def test_editor_legend_wraps_within_canvas():
    editor = Path("editor.html").read_text(encoding="utf-8")
    assert "const LEGEND_ORDER=['HPI','HB','WB','SB','PS','PC','XB','MC'];" in editor
    assert "function legendLayout(l)" in editor
    assert "Math.min(Number(state.width||1200)-80,900)" in editor
    assert "for(const row of layout.rows)" in editor
    assert "let cursor=ox-row.width/2;" in editor


def test_static_legend_uses_same_order_and_wrap_limit():
    source = Path("utils.py").read_text(encoding="utf-8")
    expected = [
        '("HPI", "Hydrophobic contact"',
        '("HB",  "Hydrogen bond"',
        '("WB",  "Water bridge"',
        '("SB",  "Salt bridge"',
        '("PS",  "π-Stacking"',
        '("PC",  "π-Cation"',
        '("XB",  "Halogen bond"',
        '("MC",  "Metal coordination"',
    ]
    positions = [source.index(token) for token in expected]
    assert positions == sorted(positions)
    assert "max_legend_width = max(260.0, min(float(canvas_width) - 80.0, 900.0))" in source
    assert "for row, row_width in rows:" in source


def test_legend_samples_match_rendered_interaction_styles():
    editor = Path("editor.html").read_text(encoding="utf-8")
    static = Path("utils.py").read_text(encoding="utf-8")

    assert "HPI:{label:'Hydrophobic contact',color:'#595959',lineWidth:3,dash:'9 5'}" in editor
    assert "HB:{label:'Hydrogen bond',color:'#0000E0',lineWidth:3,dash:'9 5'}" in editor
    assert "SB:{label:'Salt bridge',color:'#D900B0',lineWidth:3,dash:'9 5'}" in editor
    assert '("HPI", "Hydrophobic contact", (0.35, 0.35, 0.35), [9, 5])' in static
    assert '("HB",  "Hydrogen bond",       (0.0, 0.0, 0.88), [9, 5])' in static
    assert '("SB",  "Salt bridge",         (0.85, 0.0, 0.70), [9, 5])' in static
