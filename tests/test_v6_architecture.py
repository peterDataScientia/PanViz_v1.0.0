from pathlib import Path

import pandas as pd

import scientific_records
from panviz_version import PANVIZ_VERSION


class FakeBindingSiteReport:
    hydrophobic_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST"]
    hydrophobic_info = [["ALA", 10, "A", 3.80]]

    hbond_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST_D-A"]
    hbond_info = [["ASP", 34, "A", 2.80]]

    waterbridge_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST_A-W", "DIST_D-W"]
    waterbridge_info = [["SER", 79, "A", 2.70, 2.95]]

    saltbridge_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST"]
    saltbridge_info = [["GLU", 120, "A", 3.10]]

    pistacking_features = ["RESTYPE", "RESNR", "RESCHAIN", "CENTDIST", "LIG_IDX_LIST"]
    pistacking_info = [["TYR", 192, "A", 4.72, [1, 2, 3, 4, 5, 6]]]

    pication_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST"]
    pication_info = [["PHE", 111, "A", 4.20]]

    halogen_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST"]
    halogen_info = [["ASN", 76, "A", 3.30]]

    metal_features = ["RESTYPE", "RESNR", "RESCHAIN", "DIST", "COORDINATION"]
    metal_info = [["HIS", 88, "A", 2.15, "tetrahedral"]]


def test_v6_version():
    assert PANVIZ_VERSION == "6.0.0"


def test_scientific_layer_normalizes_all_eight_plip_classes(monkeypatch):
    monkeypatch.setattr(
        scientific_records,
        "BindingSiteReport",
        lambda _interactions: FakeBindingSiteReport(),
    )

    records, tables = scientific_records.build_scientific_records(object())

    assert set(records["Code"]) == {"HPI", "HB", "WB", "SB", "PS", "PC", "XB", "MC"}
    assert set(tables) == {"HPI", "HB", "WB", "SB", "PS", "PC", "XB", "MC"}
    assert len(records) == 8

    rendered = set(records.loc[records["Renderer-supported class"], "Code"])
    preserved_only = set(records.loc[~records["Renderer-supported class"], "Code"])
    assert rendered == {"HPI", "HB", "WB", "SB", "PS", "PC", "XB", "MC"}
    assert preserved_only == set()

    ps = records.loc[records["Code"] == "PS"].iloc[0]
    assert ps["Residue"] == "TYR192:A"
    assert ps["Distance (Å)"] == 4.72
    assert "LIG_IDX_LIST" in ps["Raw PLIP record"]


def test_scientific_signature_is_stable_and_presentation_independent(monkeypatch):
    monkeypatch.setattr(
        scientific_records,
        "BindingSiteReport",
        lambda _interactions: FakeBindingSiteReport(),
    )
    records, _ = scientific_records.build_scientific_records(object())

    sig1 = scientific_records.scientific_signature(records)
    sig2 = scientific_records.scientific_signature(records.copy())
    assert sig1 == sig2
    assert len(sig1) == 64

    # Columns outside the signed scientific subset must not affect the signature.
    altered = records.copy()
    altered["Presentation note"] = "moved label"
    assert scientific_records.scientific_signature(altered) == sig1


def test_scientific_exports_include_unified_and_per_class_files(monkeypatch, tmp_path):
    monkeypatch.setattr(
        scientific_records,
        "BindingSiteReport",
        lambda _interactions: FakeBindingSiteReport(),
    )
    records, tables = scientific_records.build_scientific_records(object())
    result = scientific_records.write_scientific_exports(records, tables, tmp_path)

    assert Path(result["csv"]).exists()
    assert Path(result["json"]).exists()
    assert result["signature"] == scientific_records.scientific_signature(records)
    for code in scientific_records.INTERACTION_CLASSES:
        assert (tmp_path / f"PanViz_{code}.csv").exists()


def test_v6_app_integrity_guards_are_present():
    source = Path("app.py").read_text(encoding="utf-8")

    # One temporary workspace per Streamlit session, not one per rerun.
    assert '"panviz_work_root" not in st.session_state' in source

    # Docked ligand atom serials start after the receptor maximum serial.
    assert "ligand_start=(max(receptor_serials)+1) if receptor_serials else 1" in source
    assert "start_serial=ligand_start" in source

    # v6 scientific records are built from the completed PLIP analysis while the
    # established renderer calls remain active.
    assert 'build_scientific_records(analysis_obj["my_interactions"])' in source
    assert "plip_2d_interactions(" in source
    assert "build_editor_scene(" in source


def test_main_figure_records_are_derived_from_scene_exactly():
    scene = {
        "labels": [
            {"id": "res_0", "text": "ASP34:A", "sourceResidue": "ASP34_A"},
            {"id": "res_1", "text": "TYR192:A", "sourceResidue": "TYR192_A"},
        ],
        "scientificData": {
            "interactions": [
                {
                    "id": "int_0",
                    "type": "HB",
                    "residueId": "res_0",
                    "anchorAtom": "O1",
                    "originalDistance": 2.80,
                    "multiplicity": 1,
                },
                {
                    "id": "int_1",
                    "type": "PS",
                    "residueId": "res_1",
                    "anchorAtom": "ring1",
                    "originalDistance": 4.72,
                    "multiplicity": 1,
                },
            ]
        },
    }

    figure = scientific_records.build_figure_records(scene)
    assert len(figure) == 2
    assert figure["Residue"].tolist() == ["ASP34:A", "TYR192:A"]
    assert figure["Code"].tolist() == ["HB", "PS"]
    assert figure["Distance (Å)"].tolist() == [2.80, 4.72]


def test_figure_record_labels_cover_new_interaction_classes():
    scene = {
        "labels": [{"id": "r", "text": "SER79:A"}],
        "scientificData": {
            "interactions": [
                {"id": "w", "type": "WB", "residueId": "r", "anchorAtom": "O2", "originalDistance": 2.7},
                {"id": "x", "type": "XB", "residueId": "r", "anchorAtom": "CL1", "originalDistance": 3.3},
                {"id": "m", "type": "MC", "residueId": "r", "anchorAtom": "ZN1", "originalDistance": 2.15},
            ]
        },
    }
    figure = scientific_records.build_figure_records(scene)
    assert figure["Interaction"].tolist() == [
        "Water bridge",
        "Halogen bond",
        "Metal coordination",
    ]
