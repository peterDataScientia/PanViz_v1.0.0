from pathlib import Path

from panviz_engine import (
    INTERACTION_SPECS,
    scene_to_svg,
    write_static_exports,
)


def _scene():
    return {
        "width": 640,
        "height": 420,
        "style": {
            "bondWidth": 2.0,
            "moleculeLabelSize": 16,
            "interactions": {
                code: {
                    "color": spec["color"],
                    "width": 3.0,
                    "dash": "9 5",
                    "custom": False,
                }
                for code, spec in INTERACTION_SPECS.items()
            },
        },
        "atoms": [
            {
                "id": 0,
                "name": "C1",
                "element": "C",
                "label": "C",
                "x": 280.0,
                "y": 210.0,
                "showLabel": False,
            },
            {
                "id": 1,
                "name": "O1",
                "element": "O",
                "label": "O",
                "x": 340.0,
                "y": 210.0,
                "showLabel": True,
            },
        ],
        "bonds": [{"id": "b0_1", "a": 0, "b": 1, "order": 2, "aromatic": False}],
        "labels": [
            {
                "id": "res_0",
                "text": "ASP34:A",
                "x": 470.0,
                "y": 150.0,
                "w": 92.0,
                "h": 34.0,
                "fontSize": 18,
                "visible": True,
            }
        ],
        "interactions": [
            {
                "id": "int_0",
                "type": "HB",
                "x1": 340.0,
                "y1": 210.0,
                "x2": 430.0,
                "y2": 160.0,
                "color": INTERACTION_SPECS["HB"]["color"],
                "lineWidth": 3.0,
                "multiplicity": 1,
                "dash": "9 5",
                "visible": True,
            }
        ],
        "distances": [
            {
                "id": "dist_0",
                "displayText": "2.80 Å",
                "x": 385.0,
                "y": 175.0,
                "fontSize": 15,
                "visible": True,
            }
        ],
        "legend": {
            "x": 320.0,
            "y": 390.0,
            "fontSize": 15,
            "visible": True,
            "items": [
                {
                    "type": "HB",
                    "label": "H-bond",
                    "color": INTERACTION_SPECS["HB"]["color"],
                }
            ],
        },
    }


def test_all_plip_interaction_classes_are_canonical():
    assert set(INTERACTION_SPECS) == {
        "HPI",
        "HB",
        "WB",
        "SB",
        "PS",
        "PC",
        "XB",
        "MC",
    }


def test_legacy_renderer_functions_and_chemistry_mutations_do_not_return():
    source = Path("panviz_engine.py").read_text(encoding="utf-8")
    forbidden = [
        "def _save_pymol",
        "def _get_interactions",
        "def _get_res_info",
        "def _draw_mol",
        "def set_to_neutral_pH",
        "def convert_and_write_pdb",
        "def plip_2d_interactions",
        "DetermineBonds(",
        "charge=0",
        "set_to_neutral_pH(",
    ]
    for token in forbidden:
        assert token not in source


def test_scene_svg_is_standalone_and_contains_scientific_labels():
    svg = scene_to_svg(_scene())
    assert 'viewBox="0 0 640 420"' in svg
    assert "ASP34:A" in svg
    assert "2.80 Å" in svg
    assert "H-bond" in svg


def test_static_export_writes_svg_and_png(tmp_path):
    svg_path = tmp_path / "figure.svg"
    png_path = tmp_path / "figure.png"
    write_static_exports(_scene(), svg_path, png_path, png_scale=1)
    assert svg_path.exists() and svg_path.stat().st_size > 200
    assert png_path.exists() and png_path.stat().st_size > 200
