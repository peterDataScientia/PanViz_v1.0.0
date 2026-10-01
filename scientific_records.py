from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from plip.exchange.report import BindingSiteReport


INTERACTION_CLASSES = {
    "HPI": ("hydrophobic", "Hydrophobic"),
    "HB": ("hbond", "Hydrogen bond"),
    "WB": ("waterbridge", "Water bridge"),
    "SB": ("saltbridge", "Salt bridge"),
    "PS": ("pistacking", "Pi-stacking"),
    "PC": ("pication", "Pi-cation"),
    "XB": ("halogen", "Halogen bond"),
    "MC": ("metal", "Metal coordination"),
}

# The protected publication renderer currently visualizes these five classes.
RENDERED_INTERACTION_CLASSES = {"HPI", "HB", "SB", "PS", "PC"}


def _is_missing(value) -> bool:
    if value is None:
        return True
    if isinstance(value, float):
        return pd.isna(value)
    return False


def _first_present(row: pd.Series, names):
    for name in names:
        if name in row.index and not _is_missing(row[name]):
            return row[name]
    return None


def _residue_label(row: pd.Series) -> str:
    restype = str(_first_present(row, ["RESTYPE"]) or "").strip()
    resnr = str(_first_present(row, ["RESNR"]) or "").strip()
    chain = str(_first_present(row, ["RESCHAIN"]) or "").strip()
    base = f"{restype}{resnr}".strip()
    return f"{base}:{chain}" if chain else (base or "—")


def _distance(row: pd.Series):
    for name in (
        "DIST",
        "DIST_D-A",
        "DIST_H-A",
        "DIST_A-W",
        "DIST_D-W",
        "CENTDIST",
    ):
        value = _first_present(row, [name])
        if value is None:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            pass
    return None


def build_scientific_records(my_interactions) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Normalize PLIP BindingSiteReport output without changing the renderer.

    This is a scientific-data layer only. It intentionally does not alter
    molecule depiction, residue placement, interaction routing, or the
    protected publication renderer.
    """
    report = BindingSiteReport(my_interactions)
    tables: dict[str, pd.DataFrame] = {}
    rows = []

    for code, (stem, label) in INTERACTION_CLASSES.items():
        features = list(getattr(report, stem + "_features", ()) or ())
        info = list(getattr(report, stem + "_info", ()) or ())
        table = pd.DataFrame(info, columns=features)
        tables[code] = table

        for idx, row in table.iterrows():
            raw = {
                str(key): (None if _is_missing(value) else value)
                for key, value in row.items()
            }
            # Ensure JSON-safe scalar representations in the exported record.
            raw_json = json.dumps(raw, default=str, sort_keys=True, ensure_ascii=False)
            rows.append(
                {
                    "Record ID": f"{code}_{idx + 1}",
                    "Residue": _residue_label(row),
                    "Interaction": label,
                    "Code": code,
                    "Distance (Å)": _distance(row),
                    "Rendered in figure": code in RENDERED_INTERACTION_CLASSES,
                    "Raw PLIP record": raw_json,
                }
            )

    df = pd.DataFrame(
        rows,
        columns=[
            "Record ID",
            "Residue",
            "Interaction",
            "Code",
            "Distance (Å)",
            "Rendered in figure",
            "Raw PLIP record",
        ],
    )
    return df, tables


def scientific_signature(records: pd.DataFrame) -> str:
    """Stable presentation-independent signature for scientific records."""
    import hashlib

    if records.empty:
        payload = "[]"
    else:
        stable = records[
            ["Record ID", "Residue", "Interaction", "Code", "Distance (Å)", "Raw PLIP record"]
        ].copy()
        payload = stable.to_json(
            orient="records",
            force_ascii=False,
            double_precision=10,
        )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_scientific_exports(
    records: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
    output_dir: str | Path,
) -> dict[str, str]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    unified_csv = output_dir / "PanViz_scientific_records.csv"
    unified_json = output_dir / "PanViz_scientific_records.json"

    records.to_csv(unified_csv, index=False)
    unified_json.write_text(
        records.to_json(
            orient="records",
            indent=2,
            force_ascii=False,
            double_precision=10,
        ),
        encoding="utf-8",
    )

    for code, table in tables.items():
        if not table.empty:
            table.to_csv(output_dir / f"PanViz_{code}.csv", index=False)

    return {
        "csv": str(unified_csv),
        "json": str(unified_json),
        "signature": scientific_signature(records),
    }


FIGURE_INTERACTION_LABELS = {
    "HPI": "Hydrophobic",
    "HB": "Hydrogen bond",
    "SB": "Salt bridge",
    "PS": "Pi-stacking",
    "PC": "Pi-cation",
}


def build_figure_records(scene: dict) -> pd.DataFrame:
    """Return the exact scientific records represented by the protected figure.

    The count is derived from scene["scientificData"]["interactions"], which is
    produced by the same protected interaction/layout pipeline used by the
    publication figure. This deliberately avoids estimating figure counts from
    the broader eight-class PLIP record table.
    """
    labels = {
        str(item.get("id")): item
        for item in scene.get("labels", [])
    }
    rows = []
    scientific = scene.get("scientificData", {}).get("interactions", [])

    for index, item in enumerate(scientific, start=1):
        residue_obj = labels.get(str(item.get("residueId")), {})
        residue = (
            residue_obj.get("text")
            or residue_obj.get("sourceResidue")
            or item.get("residueId")
            or "—"
        )
        code = str(item.get("type") or "")
        rows.append(
            {
                "Figure record ID": item.get("id") or f"figure_{index}",
                "Residue": residue,
                "Interaction": FIGURE_INTERACTION_LABELS.get(code, code),
                "Code": code,
                "Ligand atom": item.get("anchorAtom"),
                "Distance (Å)": item.get("originalDistance"),
                "Multiplicity": int(item.get("multiplicity", 1) or 1),
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "Figure record ID",
            "Residue",
            "Interaction",
            "Code",
            "Ligand atom",
            "Distance (Å)",
            "Multiplicity",
        ],
    )
