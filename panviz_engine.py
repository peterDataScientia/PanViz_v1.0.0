from __future__ import annotations

import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

# Native-library load order matters on Linux/Streamlit Cloud:
# RDKit must be imported before Open Babel.
from rdkit import Chem
from rdkit.Chem import rdDepictor

import numpy as np
import pandas as pd

from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex
from openbabel import pybel


PANVIZ_ENGINE_VERSION = "6.0.0"

INTERACTION_SPECS = {
    "HPI": {
        "report": "hydrophobic",
        "label": "Hydrophobic",
        "color": "#595959",
        "distance": ("DIST",),
    },
    "HB": {
        "report": "hbond",
        "label": "H-bond",
        "color": "#2455D6",
        "distance": ("DIST_D-A", "DIST_H-A"),
    },
    "WB": {
        "report": "waterbridge",
        "label": "Water bridge",
        "color": "#1596B8",
        "distance": ("DIST_A-W", "DIST_D-W"),
    },
    "SB": {
        "report": "saltbridge",
        "label": "Salt bridge",
        "color": "#D900B0",
        "distance": ("DIST",),
    },
    "PS": {
        "report": "pistacking",
        "label": "π-Stacking",
        "color": "#008500",
        "distance": ("CENTDIST",),
    },
    "PC": {
        "report": "pication",
        "label": "π-Cation",
        "color": "#E08000",
        "distance": ("DIST",),
    },
    "XB": {
        "report": "halogen",
        "label": "Halogen bond",
        "color": "#7A5CC7",
        "distance": ("DIST",),
    },
    "MC": {
        "report": "metal",
        "label": "Metal coordination",
        "color": "#A45700",
        "distance": ("DIST",),
    },
}

REPORT_TO_CODE = {v["report"]: k for k, v in INTERACTION_SPECS.items()}


def _site_parts(binding_site: str) -> tuple[str, str, str]:
    parts = str(binding_site).split(":")
    if len(parts) < 3:
        raise ValueError(
            f"Binding-site identifier {binding_site!r} is not in RES:CHAIN:NUMBER form."
        )
    return parts[0].strip(), parts[1].strip(), parts[2].strip()


def _element_from_pdb_line(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.title()
    raw = "".join(ch for ch in line[12:16] if ch.isalpha()).strip()
    if len(raw) >= 2 and raw[:2].lower() in {"cl", "br"}:
        return raw[:2].title()
    return raw[:1].upper() or "C"


def _extract_ligand_block(pdb_path: str | Path, binding_site: str):
    resname, chain, resnum = _site_parts(binding_site)
    lines = Path(pdb_path).read_text(encoding="utf-8", errors="replace").splitlines()

    atom_lines = []
    atoms = []
    serials = set()

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        padded = line.ljust(80)
        if (
            padded[17:20].strip() != resname
            or padded[21:22].strip() != chain
            or padded[22:26].strip() != resnum
        ):
            continue

        try:
            serial = int(padded[6:11])
            x = float(padded[30:38])
            y = float(padded[38:46])
            z = float(padded[46:54])
        except ValueError as exc:
            raise ValueError(
                f"Invalid ligand atom record encountered for {binding_site}: {line}"
            ) from exc

        name = padded[12:16].strip() or f"A{len(atoms) + 1}"
        atom = {
            "serial": serial,
            "name": name,
            "element": _element_from_pdb_line(padded),
            "coord": (x, y, z),
        }
        atom_lines.append(padded.rstrip())
        atoms.append(atom)
        serials.add(serial)

    if not atoms:
        raise ValueError(
            f"No ATOM/HETATM records for binding site {binding_site} were found in {pdb_path}."
        )

    conect_lines = []
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        values = []
        for token in line.split()[1:]:
            try:
                values.append(int(token))
            except ValueError:
                pass
        if not values or values[0] not in serials:
            continue
        kept = [v for v in values if v in serials]
        if len(kept) >= 2:
            conect_lines.append("CONECT" + "".join(f"{v:5d}" for v in kept))

    block = "\n".join(atom_lines + conect_lines + ["END"]) + "\n"
    return block, atoms


def _nearest_source_atom(
    element: str,
    xyz: tuple[float, float, float],
    source_atoms: list[dict],
    used: set[int],
):
    candidates = [
        a
        for a in source_atoms
        if a["serial"] not in used and a["element"].upper() == element.upper()
    ]
    if not candidates:
        candidates = [a for a in source_atoms if a["serial"] not in used]
    if not candidates:
        return None

    x, y, z = xyz
    return min(
        candidates,
        key=lambda a: (
            (a["coord"][0] - x) ** 2
            + (a["coord"][1] - y) ** 2
            + (a["coord"][2] - z) ** 2
        ),
    )


def _build_ligand_graph(pdb_path: str | Path, binding_site: str):
    """Build the PanViz ligand graph without changing protonation or total charge.

    Open Babel is used only to interpret the source ligand records and serialize
    the perceived molecular graph to MOL2. PanViz does not call AddHydrogens,
    does not force a total charge, and does not run a neutral-pH transformation.
    RDKit is then used only for sanitization and deterministic 2D depiction.
    """

    block, source_atoms = _extract_ligand_block(pdb_path, binding_site)

    obmol = pybel.readstring("pdb", block)
    mol2 = obmol.write("mol2")
    mol = Chem.MolFromMol2Block(mol2, sanitize=True, removeHs=False)

    if mol is None:
        # Fallback uses only PDB connectivity/proximity perception. It deliberately
        # avoids rdDetermineBonds and therefore never imposes charge=0.
        mol = Chem.MolFromPDBBlock(
            block,
            sanitize=True,
            removeHs=False,
            proximityBonding=True,
        )
    if mol is None:
        raise ValueError(
            f"PanViz could not construct a ligand graph for binding site {binding_site}."
        )

    if mol.GetNumConformers() == 0:
        raise ValueError("The ligand graph has no coordinates.")

    conf = mol.GetConformer()
    used_serials: set[int] = set()

    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        source = _nearest_source_atom(
            atom.GetSymbol(),
            (float(p.x), float(p.y), float(p.z)),
            source_atoms,
            used_serials,
        )
        if source is not None:
            used_serials.add(source["serial"])
            atom.SetIntProp("_PanVizSourceSerial", int(source["serial"]))
            atom.SetProp("_PanVizSourceName", str(source["name"]))

    try:
        mol = Chem.RemoveHs(mol, sanitize=True)
    except Exception:
        # Keeping explicit hydrogens is scientifically safer than mutating the
        # ligand to make removal succeed.
        pass

    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        pass

    # Kekulization is a depiction operation here: it converts aromatic bonds to
    # an alternating single/double drawing representation without changing
    # protonation, atom identity, or formal charge.
    draw_mol = Chem.Mol(mol)
    try:
        Chem.Kekulize(draw_mol, clearAromaticFlags=True)
    except Exception:
        pass

    rdDepictor.Compute2DCoords(draw_mol, canonOrient=True)
    conf2d = draw_mol.GetConformer()

    atoms = []
    serial_to_atom_id = {}
    for atom in draw_mol.GetAtoms():
        p = conf2d.GetAtomPosition(atom.GetIdx())
        serial = (
            atom.GetIntProp("_PanVizSourceSerial")
            if atom.HasProp("_PanVizSourceSerial")
            else None
        )
        source_name = (
            atom.GetProp("_PanVizSourceName")
            if atom.HasProp("_PanVizSourceName")
            else f"A{atom.GetIdx() + 1}"
        )
        element = atom.GetSymbol()
        charge = int(atom.GetFormalCharge())
        label = element
        if charge:
            magnitude = abs(charge)
            sign = "+" if charge > 0 else "−"
            label = element + (str(magnitude) if magnitude > 1 else "") + sign

        row = {
            "id": atom.GetIdx(),
            "name": source_name,
            "sourceSerial": serial,
            "element": element,
            "formalCharge": charge,
            "label": label,
            "rawX": float(p.x),
            "rawY": float(p.y),
            "showLabel": element != "C" or charge != 0,
        }
        atoms.append(row)
        if serial is not None:
            serial_to_atom_id[int(serial)] = atom.GetIdx()

    bonds = []
    for bond in draw_mol.GetBonds():
        value = float(bond.GetBondTypeAsDouble())
        if value >= 2.5:
            order = 3
        elif value >= 1.5:
            order = 2
        else:
            order = 1
        bonds.append(
            {
                "id": f"b{bond.GetBeginAtomIdx()}_{bond.GetEndAtomIdx()}",
                "a": bond.GetBeginAtomIdx(),
                "b": bond.GetEndAtomIdx(),
                "order": order,
                "aromatic": bool(bond.GetIsAromatic()),
            }
        )

    return {
        "atoms": atoms,
        "bonds": bonds,
        "sourceAtoms": source_atoms,
        "serialToAtomId": serial_to_atom_id,
    }


def _safe_float(value):
    try:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return None
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _safe_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def _serials(value) -> list[int]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, np.ndarray)):
        tokens = value
    else:
        tokens = str(value).replace(";", ",").split(",")
    out = []
    for token in tokens:
        try:
            out.append(int(str(token).strip()))
        except (TypeError, ValueError):
            pass
    return out


def _coord(value):
    if isinstance(value, (tuple, list, np.ndarray)) and len(value) >= 3:
        try:
            return tuple(float(value[i]) for i in range(3))
        except (TypeError, ValueError):
            return None
    if isinstance(value, str):
        cleaned = value.strip().strip("()[]")
        parts = cleaned.replace(",", " ").split()
        if len(parts) >= 3:
            try:
                return tuple(float(parts[i]) for i in range(3))
            except ValueError:
                return None
    return None


def _report_tables(binding_site_interactions):
    report = BindingSiteReport(binding_site_interactions)
    tables = {}
    for code, spec in INTERACTION_SPECS.items():
        stem = spec["report"]
        features = list(getattr(report, stem + "_features", ()))
        info = list(getattr(report, stem + "_info", ()))
        tables[code] = pd.DataFrame(info, columns=features)
    return report, tables


def _distance_values(code: str, row: pd.Series):
    fields = INTERACTION_SPECS[code]["distance"]
    values = [_safe_float(row.get(field)) for field in fields]
    values = [v for v in values if v is not None]
    primary = values[0] if values else None
    secondary = values[1] if len(values) > 1 else None
    return primary, secondary


def _record_ligand_serials(code: str, row: pd.Series, ligand_serial_set: set[int]):
    values: list[int] = []

    if code == "HPI":
        values = _serials(row.get("LIGCARBONIDX"))
    elif code == "HB":
        values = _serials(
            row.get("ACCEPTORIDX") if _safe_bool(row.get("PROTISDON")) else row.get("DONORIDX")
        )
    elif code == "WB":
        values = _serials(
            row.get("ACCEPTOR_IDX")
            if _safe_bool(row.get("PROTISDON"))
            else row.get("DONOR_IDX")
        )
    elif code in {"SB", "PS", "PC"}:
        values = _serials(row.get("LIG_IDX_LIST"))
    elif code == "XB":
        donor = _serials(row.get("DON_IDX"))
        acceptor = _serials(row.get("ACC_IDX"))
        values = donor if any(v in ligand_serial_set for v in donor) else acceptor
    elif code == "MC":
        metal = _serials(row.get("METAL_IDX"))
        target = _serials(row.get("TARGET_IDX"))
        if any(v in ligand_serial_set for v in metal):
            values = metal
        elif any(v in ligand_serial_set for v in target):
            values = target

    return [v for v in values if v in ligand_serial_set]


def _nearest_serial_to_coord(coord, source_atoms):
    if coord is None or not source_atoms:
        return None
    x, y, z = coord
    atom = min(
        source_atoms,
        key=lambda a: (
            (a["coord"][0] - x) ** 2
            + (a["coord"][1] - y) ** 2
            + (a["coord"][2] - z) ** 2
        ),
    )
    return int(atom["serial"])


def _canonical_records(tables: dict[str, pd.DataFrame], source_atoms: list[dict]):
    ligand_serial_set = {int(a["serial"]) for a in source_atoms}
    source_by_serial = {int(a["serial"]): a for a in source_atoms}
    records = []

    for code, table in tables.items():
        for row_index, row in table.iterrows():
            restype = str(row.get("RESTYPE", "")).strip()
            resnr = str(row.get("RESNR", "")).strip()
            reschain = str(row.get("RESCHAIN", "")).strip()
            residue = f"{restype}{resnr}_{reschain}".strip("_")

            serial_list = _record_ligand_serials(code, row, ligand_serial_set)
            if not serial_list:
                fallback_coord = _coord(row.get("LIGCOO"))
                if fallback_coord is None and code == "MC":
                    for key in ("METALCOO", "TARGETCOO"):
                        candidate = _coord(row.get(key))
                        serial = _nearest_serial_to_coord(candidate, source_atoms)
                        if serial in ligand_serial_set:
                            serial_list = [serial]
                            break
                else:
                    serial = _nearest_serial_to_coord(fallback_coord, source_atoms)
                    if serial is not None:
                        serial_list = [serial]

            primary, secondary = _distance_values(code, row)
            names = [
                source_by_serial[s]["name"]
                for s in serial_list
                if s in source_by_serial
            ]

            details = {}
            for key in (
                "TYPE",
                "ANGLE",
                "OFFSET",
                "DON_ANGLE",
                "ACC_ANGLE",
                "WATER_ANGLE",
                "COORDINATION",
                "GEOMETRY",
                "LOCATION",
                "LIG_GROUP",
                "PROTISDON",
                "PROTISPOS",
                "PROTCHARGED",
            ):
                if key in row.index and pd.notna(row.get(key)):
                    details[key] = str(row.get(key))

            records.append(
                {
                    "scientificId": f"{code}_{row_index}",
                    "type": code,
                    "label": INTERACTION_SPECS[code]["label"],
                    "residue": residue,
                    "ligandSerials": serial_list,
                    "ligandAtomNames": names,
                    "distance": primary,
                    "secondaryDistance": secondary,
                    "details": details,
                }
            )

    return records


def _interaction_summary(records: list[dict]):
    rows = []
    for record in records:
        distance = record["distance"]
        secondary = record.get("secondaryDistance")
        detail = ""
        if distance is not None and secondary is not None:
            detail = f"{distance:.2f} / {secondary:.2f}"
        elif distance is not None:
            detail = f"{distance:.2f}"
        rows.append(
            {
                "Residue": record["residue"].replace("_", ":"),
                "Interaction": record["label"],
                "Distance (Å)": distance,
                "Distance details": detail,
            }
        )
    return pd.DataFrame(
        rows,
        columns=["Residue", "Interaction", "Distance (Å)", "Distance details"],
    )


def run_panviz_analysis(
    pdb_file: str | Path,
    binding_site: str,
    output_root: str | Path,
):
    """Run PLIP once and normalize its official report into PanViz records."""

    pdb_file = Path(pdb_file)
    output_root = Path(output_root)
    site_dir = output_root / binding_site.replace(":", "_")
    structures_dir = site_dir / "structures"
    interaction_dir = site_dir / "interactions"
    figures_dir = site_dir / "figures"

    for directory in (structures_dir, interaction_dir, figures_dir):
        directory.mkdir(parents=True, exist_ok=True)

    complex_obj = PDBComplex()
    complex_obj.load_pdb(str(pdb_file))
    complex_obj.analyze()

    if binding_site not in complex_obj.interaction_sets:
        raise RuntimeError(
            f"PLIP did not return interaction data for binding site {binding_site}."
        )

    binding_interactions = complex_obj.interaction_sets[binding_site]
    report, tables = _report_tables(binding_interactions)

    for code, table in tables.items():
        if table.empty:
            continue
        table.to_csv(
            interaction_dir / f"{binding_site.replace(':', '_')}_{code}.csv",
            index=False,
        )

    prepared_pdb = structures_dir / "PanViz_PLIP_structure.pdb"
    try:
        complex_obj.protcomplex.write("pdb", str(prepared_pdb), overwrite=True)
    except Exception:
        shutil.copy2(pdb_file, prepared_pdb)

    ligand_block, source_atoms = _extract_ligand_block(pdb_file, binding_site)
    records = _canonical_records(tables, source_atoms)

    return {
        "engineVersion": PANVIZ_ENGINE_VERSION,
        "bindingSite": binding_site,
        "file_prot": str(prepared_pdb),
        "site_dir": str(site_dir),
        "figures_dir": str(figures_dir),
        "interaction_dir": str(interaction_dir),
        "plip_complex": complex_obj,
        "plip_interactions": binding_interactions,
        "plip_report": report,
        "tables": tables,
        "records": records,
        "interaction_summary": _interaction_summary(records),
        "ligand_block": ligand_block,
        "source_atoms": source_atoms,
    }


def _display_residue(residue: str) -> str:
    text = str(residue)
    if "_" in text:
        base, chain = text.rsplit("_", 1)
        if base and chain:
            return f"{base}:{chain}"
    return text


def _rect_overlap(a, b, gap=8.0):
    return not (
        a[0] + a[2] / 2 + gap < b[0] - b[2] / 2
        or a[0] - a[2] / 2 - gap > b[0] + b[2] / 2
        or a[1] + a[3] / 2 + gap < b[1] - b[3] / 2
        or a[1] - a[3] / 2 - gap > b[1] + b[3] / 2
    )


def _scale_ligand(graph, width: int, height: int):
    xs = [a["rawX"] for a in graph["atoms"]]
    ys = [a["rawY"] for a in graph["atoms"]]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    target_w = max(220.0, width * 0.46)
    target_h = max(180.0, height * 0.48)
    scale = min(target_w / span_x, target_h / span_y)

    cx_raw = (min_x + max_x) / 2
    cy_raw = (min_y + max_y) / 2
    cx = width / 2
    cy = height / 2 - 18

    atoms = []
    for atom in graph["atoms"]:
        x = cx + (atom["rawX"] - cx_raw) * scale
        y = cy - (atom["rawY"] - cy_raw) * scale
        row = dict(atom)
        row["x"] = float(x)
        row["y"] = float(y)
        atoms.append(row)

    return atoms, (cx, cy)


def _anchor_for_record(record, atom_by_serial, ligand_center):
    points = [
        (atom_by_serial[s]["x"], atom_by_serial[s]["y"])
        for s in record.get("ligandSerials", [])
        if s in atom_by_serial
    ]
    if not points:
        return ligand_center
    return (
        float(sum(p[0] for p in points) / len(points)),
        float(sum(p[1] for p in points) / len(points)),
    )


def _label_dimensions(text: str, font_size: int):
    width = max(72.0, min(170.0, len(text) * font_size * 0.61 + 28.0))
    height = max(32.0, font_size * 1.55)
    return width, height


def _place_residues(
    records,
    atom_by_serial,
    ligand_center,
    width,
    height,
    residue_font_size,
):
    grouped = defaultdict(list)
    for record in records:
        grouped[record["residue"]].append(record)

    cx, cy = ligand_center
    specs = []
    for residue, items in grouped.items():
        anchors = [_anchor_for_record(r, atom_by_serial, ligand_center) for r in items]
        ax = sum(p[0] for p in anchors) / len(anchors)
        ay = sum(p[1] for p in anchors) / len(anchors)
        vx, vy = ax - cx, ay - cy
        angle = math.atan2(vy, vx) if abs(vx) + abs(vy) > 1e-6 else 0.0
        text = _display_residue(residue)
        w, h = _label_dimensions(text, residue_font_size)
        specs.append(
            {
                "residue": residue,
                "text": text,
                "anchor": (ax, ay),
                "preferredAngle": angle,
                "w": w,
                "h": h,
            }
        )

    specs.sort(key=lambda x: x["preferredAngle"])
    count = max(1, len(specs))
    base_radius = min(width, height) * (0.31 if count <= 10 else 0.35)
    base_radius = max(base_radius, 190.0)

    placed = []
    by_residue = {}

    angle_offsets = [0]
    for step in range(1, 10):
        angle_offsets.extend([math.radians(9 * step), -math.radians(9 * step)])
    radial_offsets = [0, 48, 92, 138]

    for index, spec in enumerate(specs):
        preferred = spec["preferredAngle"]
        if len(specs) > 1 and all(
            abs(s["preferredAngle"]) < 1e-8 for s in specs
        ):
            preferred = -math.pi + 2 * math.pi * index / len(specs)

        best = None
        best_score = float("inf")
        for radial in radial_offsets:
            radius = base_radius + radial
            for delta in angle_offsets:
                angle = preferred + delta
                px = cx + math.cos(angle) * radius
                py = cy + math.sin(angle) * radius
                box = (px, py, spec["w"], spec["h"])

                if (
                    px - spec["w"] / 2 < 14
                    or px + spec["w"] / 2 > width - 14
                    or py - spec["h"] / 2 < 14
                    or py + spec["h"] / 2 > height - 62
                ):
                    continue

                overlap_count = sum(_rect_overlap(box, prev, gap=10) for prev in placed)
                if overlap_count:
                    continue

                angle_cost = abs(delta) * 95.0
                radial_cost = radial * 0.45
                anchor = spec["anchor"]
                tether = math.hypot(px - anchor[0], py - anchor[1]) * 0.035
                score = angle_cost + radial_cost + tether
                if score < best_score:
                    best_score = score
                    best = (px, py)

        if best is None:
            angle = preferred
            radius = base_radius + 138
            best = (
                min(width - spec["w"] / 2 - 14, max(spec["w"] / 2 + 14, cx + math.cos(angle) * radius)),
                min(height - spec["h"] / 2 - 62, max(spec["h"] / 2 + 14, cy + math.sin(angle) * radius)),
            )

        px, py = best
        box = (px, py, spec["w"], spec["h"])
        placed.append(box)
        by_residue[spec["residue"]] = {
            "text": spec["text"],
            "x": float(px),
            "y": float(py),
            "anchorX": float(spec["anchor"][0]),
            "anchorY": float(spec["anchor"][1]),
            "w": float(spec["w"]),
            "h": float(spec["h"]),
        }

    return by_residue


def _label_boundary(label, from_x, from_y):
    dx = from_x - label["x"]
    dy = from_y - label["y"]
    if abs(dx) < 1e-8 and abs(dy) < 1e-8:
        return label["x"], label["y"]

    half_w = label["w"] / 2
    half_h = label["h"] / 2
    tx = half_w / abs(dx) if abs(dx) > 1e-8 else float("inf")
    ty = half_h / abs(dy) if abs(dy) > 1e-8 else float("inf")
    t = min(tx, ty)
    return label["x"] + dx * t, label["y"] + dy * t


def build_editor_scene(
    pdb_file: str | Path,
    binding_site: str,
    width: int = 1200,
    height: int = 850,
    analysis: dict | None = None,
):
    """Build the PanViz 6 scene from a canonical PLIP record model."""

    if analysis is None:
        raise ValueError(
            "PanViz 6 requires the reusable analysis object from run_panviz_analysis()."
        )

    graph = _build_ligand_graph(pdb_file, binding_site)
    atoms, ligand_center = _scale_ligand(graph, int(width), int(height))
    atom_by_serial = {
        int(atom["sourceSerial"]): atom
        for atom in atoms
        if atom.get("sourceSerial") is not None
    }

    residue_font_size = max(16, min(22, int(min(width, height) * 0.024)))
    distance_font_size = max(13, min(17, int(min(width, height) * 0.017)))

    residue_layout = _place_residues(
        analysis["records"],
        atom_by_serial,
        ligand_center,
        int(width),
        int(height),
        residue_font_size,
    )

    labels = []
    residue_id = {}
    for idx, (source, item) in enumerate(residue_layout.items()):
        rid = f"res_{idx}"
        residue_id[source] = rid
        labels.append(
            {
                "id": rid,
                "text": item["text"],
                "sourceResidue": source,
                "x": item["x"],
                "y": item["y"],
                "anchorX": item["anchorX"],
                "anchorY": item["anchorY"],
                "w": item["w"],
                "h": item["h"],
                "fontSize": residue_font_size,
                "fontFamily": "Arial",
                "textColor": "#000000",
                "nodeShape": "bubble",
                "styleOverrides": {},
                "bold": True,
                "italic": False,
                "underline": False,
                "rotation": 0,
                "backgroundEnabled": True,
                "backgroundColor": "#ffffff",
                "backgroundOpacity": 0.98,
                "bubbleColor": "#0AFFEF",
                "leaderVisible": True,
                "visible": True,
            }
        )

    pair_counts = Counter()
    for record in analysis["records"]:
        anchor_name = (
            record["ligandAtomNames"][0]
            if record.get("ligandAtomNames")
            else "ligand-group"
        )
        pair_counts[(record["residue"], anchor_name)] += 1

    interactions = []
    distances = []
    residue_record_index = defaultdict(int)

    for idx, record in enumerate(analysis["records"]):
        label = residue_layout.get(record["residue"])
        rid = residue_id.get(record["residue"])
        if label is None or rid is None:
            continue

        ax, ay = _anchor_for_record(record, atom_by_serial, ligand_center)
        bx, by = _label_boundary(label, ax, ay)
        dx, dy = bx - ax, by - ay
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length

        local_index = residue_record_index[record["residue"]]
        residue_record_index[record["residue"]] += 1
        parallel_offset = ((local_index % 5) - 2) * 2.4

        x1 = ax + nx * parallel_offset
        y1 = ay + ny * parallel_offset
        x2 = bx + nx * parallel_offset
        y2 = by + ny * parallel_offset

        anchor_name = (
            record["ligandAtomNames"][0]
            if record.get("ligandAtomNames")
            else "ligand-group"
        )
        multiplicity = max(
            1,
            pair_counts[(record["residue"], anchor_name)],
        )

        iid = f"int_{idx}"
        spec = INTERACTION_SPECS[record["type"]]
        interactions.append(
            {
                "id": iid,
                "scientificSourceId": record["scientificId"],
                "type": record["type"],
                "residueId": rid,
                "anchorAtom": anchor_name,
                "anchorSerials": record.get("ligandSerials", []),
                "x1": float(x1),
                "y1": float(y1),
                "x2": float(x2),
                "y2": float(y2),
                "originalDistance": record.get("distance"),
                "secondaryDistance": record.get("secondaryDistance"),
                "multiplicity": multiplicity,
                "color": spec["color"],
                "visible": True,
                "lineWidth": 3.0,
                "customLineWidth": False,
                "customColor": False,
                "customDash": False,
                "dash": "9 5",
                "rotation": 0,
            }
        )

        primary = record.get("distance")
        secondary = record.get("secondaryDistance")
        if primary is not None:
            if secondary is not None:
                distance_text = f"{primary:.2f}/{secondary:.2f} Å"
            else:
                distance_text = f"{primary:.2f} Å"

            midpoint_x = (x1 + x2) / 2
            midpoint_y = (y1 + y2) / 2
            side = -1 if local_index % 2 else 1
            offset = 13 + 7 * (local_index // 2)
            tx = midpoint_x + nx * side * offset
            ty = midpoint_y + ny * side * offset
            text_w = max(48.0, len(distance_text) * distance_font_size * 0.56)

            distances.append(
                {
                    "id": f"dist_{idx}",
                    "interactionId": iid,
                    "originalDistance": primary,
                    "secondaryDistance": secondary,
                    "displayText": distance_text,
                    "x": float(tx),
                    "y": float(ty),
                    "w": float(text_w),
                    "h": float(distance_font_size * 1.3),
                    "fontSize": distance_font_size,
                    "fontFamily": "Arial",
                    "textColor": "#D100A0",
                    "backgroundColor": "#ffffff",
                    "backgroundOpacity": 0.96,
                    "styleOverrides": {},
                    "bold": True,
                    "italic": False,
                    "underline": False,
                    "rotation": 0,
                    "visible": True,
                }
            )

    present = []
    for item in interactions:
        if item["type"] not in present:
            present.append(item["type"])

    legend_items = [
        {
            "type": code,
            "label": INTERACTION_SPECS[code]["label"],
            "color": INTERACTION_SPECS[code]["color"],
        }
        for code in present
    ]
    legend_font_size = 15 if height <= 850 else 17
    legend_w = max(
        140,
        sum(66 + len(x["label"]) * 7 for x in legend_items)
        + max(0, len(legend_items) - 1) * 16,
    )

    scene = {
        "version": PANVIZ_ENGINE_VERSION,
        "width": int(width),
        "height": int(height),
        "site": binding_site,
        "style": {
            "bondWidth": 2.0,
            "interactionWidth": 3.0,
            "noncovalentColor": "#595959",
            "moleculeLabelSize": max(14, min(19, int(min(width, height) * 0.019))),
            "residue": {
                "labelFormat": "nameNumChain",
                "fontFamily": "Arial",
                "fontSize": residue_font_size,
                "fontWeight": 700,
                "textColor": "#000000",
                "shape": "bubble",
                "backgroundColor": "#ffffff",
                "backgroundOpacity": 0.98,
                "bubbleColor": "#0AFFEF",
                "colorMode": "neutral",
            },
            "distance": {
                "fontFamily": "Arial",
                "fontSize": distance_font_size,
                "fontWeight": 700,
                "textColor": "#D100A0",
                "backgroundColor": "#ffffff",
                "backgroundOpacity": 0.96,
            },
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
        "atoms": atoms,
        "bonds": graph["bonds"],
        "labels": labels,
        "interactions": interactions,
        "distances": distances,
        "legend": {
            "x": float(width / 2),
            "y": float(height - 28),
            "dx": 0.0,
            "dy": 0.0,
            "w": float(legend_w),
            "h": float(legend_font_size + 24),
            "fontSize": legend_font_size,
            "fontFamily": "Arial",
            "bold": True,
            "italic": False,
            "underline": False,
            "rotation": 0,
            "items": legend_items,
            "visible": True,
            "movable": True,
        },
        "metadata": {
            "editorMode": "panviz-6-canonical-scene",
            "chemistryEditable": False,
            "annotationLayerFullyEditable": True,
            "plipDataImmutable": True,
            "moleculeLocked": True,
            "canonicalInteractionTypes": list(INTERACTION_SPECS.keys()),
            "chemistryPolicy": (
                "Ligand protonation/formal charge is not neutralized or forced to zero for depiction."
            ),
            "analysisReuse": True,
            "layoutEngine": "deterministic-radial-constraint-layout",
        },
        "scientificData": {
            "interactions": [
                {
                    "id": item["id"],
                    "scientificSourceId": item.get("scientificSourceId"),
                    "type": item["type"],
                    "residueId": item["residueId"],
                    "anchorAtom": item["anchorAtom"],
                    "anchorSerials": item.get("anchorSerials", []),
                    "originalDistance": item.get("originalDistance"),
                    "secondaryDistance": item.get("secondaryDistance"),
                    "multiplicity": item.get("multiplicity", 1),
                }
                for item in interactions
            ]
        },
    }

    return scene, None


def _svg_line(x1, y1, x2, y2, color, width, dash=""):
    dash_attr = f' stroke-dasharray="{escape(str(dash))}"' if dash else ""
    return (
        f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" '
        f'stroke="{escape(str(color))}" stroke-width="{width:.3f}" '
        f'stroke-linecap="round"{dash_attr}/>'
    )


def scene_to_svg(scene: dict) -> str:
    """Serialize the canonical PanViz scene to a standalone SVG."""

    width = int(scene["width"])
    height = int(scene["height"])
    atom_by_id = {a["id"]: a for a in scene.get("atoms", [])}
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">'
        ),
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
        '<g id="interactions">',
    ]

    for item in scene.get("interactions", []):
        if item.get("visible", True) is False:
            continue
        style = scene["style"]["interactions"].get(item["type"], {})
        color = item.get("color") or style.get("color", "#595959")
        base_width = float(item.get("lineWidth", 3.0))
        width_value = base_width * max(1, int(item.get("multiplicity", 1)))
        dash = item.get("dash", style.get("dash", "9 5"))
        parts.append(
            _svg_line(
                float(item["x1"]),
                float(item["y1"]),
                float(item["x2"]),
                float(item["y2"]),
                color,
                width_value,
                dash,
            )
        )
    parts.append("</g><g id=\"molecule\">")

    bond_color = "#000000"
    bond_width = float(scene["style"].get("bondWidth", 2.0))
    for bond in scene.get("bonds", []):
        a = atom_by_id.get(bond["a"])
        b = atom_by_id.get(bond["b"])
        if not a or not b:
            continue
        x1, y1 = float(a["x"]), float(a["y"])
        x2, y2 = float(b["x"]), float(b["y"])
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length
        order = int(bond.get("order", 1))
        if order == 1:
            parts.append(_svg_line(x1, y1, x2, y2, bond_color, bond_width))
        elif order == 2:
            off = max(3.8, bond_width * 2.8) / 2
            parts.append(_svg_line(x1 + nx * off, y1 + ny * off, x2 + nx * off, y2 + ny * off, bond_color, bond_width * 0.92))
            parts.append(_svg_line(x1 - nx * off, y1 - ny * off, x2 - nx * off, y2 - ny * off, bond_color, bond_width * 0.92))
        else:
            off = max(4.2, bond_width * 3.0)
            for shift in (-off, 0, off):
                parts.append(_svg_line(x1 + nx * shift, y1 + ny * shift, x2 + nx * shift, y2 + ny * shift, bond_color, bond_width * 0.88))

    atom_size = int(scene["style"].get("moleculeLabelSize", 16))
    element_colors = {
        "O": "#D90000",
        "N": "#0000CC",
        "S": "#B78A00",
        "P": "#D96600",
        "B": "#E57373",
        "F": "#178A17",
        "Cl": "#178A17",
        "Br": "#8A3B12",
        "I": "#6B3FA0",
    }
    for atom in scene.get("atoms", []):
        if not atom.get("showLabel", False):
            continue
        x, y = float(atom["x"]), float(atom["y"])
        color = element_colors.get(atom.get("element"), "#000000")
        label = escape(str(atom.get("label") or atom.get("element") or ""))
        parts.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="8.5" fill="#ffffff" fill-opacity="0.95"/>')
        parts.append(
            f'<text x="{x:.3f}" y="{y + atom_size * 0.34:.3f}" text-anchor="middle" '
            f'font-family="Arial,Helvetica,sans-serif" font-size="{atom_size}" font-weight="700" '
            f'fill="{color}">{label}</text>'
        )
    parts.append("</g><g id=\"residues\">")

    for label in scene.get("labels", []):
        if label.get("visible", True) is False:
            continue
        x, y = float(label["x"]), float(label["y"])
        w, h = float(label["w"]), float(label["h"])
        parts.append(
            f'<rect x="{x - w / 2:.3f}" y="{y - h / 2:.3f}" width="{w:.3f}" height="{h:.3f}" '
            f'rx="{h / 2:.3f}" fill="#ffffff" stroke="#0AFFEF" stroke-width="2"/>'
        )
        text = escape(str(label["text"]))
        fs = int(label.get("fontSize", 18))
        parts.append(
            f'<text x="{x:.3f}" y="{y + fs * 0.34:.3f}" text-anchor="middle" '
            f'font-family="Arial,Helvetica,sans-serif" font-size="{fs}" font-weight="700" '
            f'fill="#000000">{text}</text>'
        )
    parts.append("</g><g id=\"distances\">")

    for dist in scene.get("distances", []):
        if dist.get("visible", True) is False:
            continue
        text = escape(str(dist.get("displayText", "")))
        x, y = float(dist["x"]), float(dist["y"])
        fs = int(dist.get("fontSize", 15))
        parts.append(
            f'<text x="{x:.3f}" y="{y + fs * 0.34:.3f}" text-anchor="middle" '
            f'font-family="Arial,Helvetica,sans-serif" font-size="{fs}" font-weight="700" '
            f'fill="#D100A0" stroke="#ffffff" stroke-width="3" paint-order="stroke">{text}</text>'
        )
    parts.append("</g>")

    legend = scene.get("legend", {})
    if legend.get("visible", True) and legend.get("items"):
        items = legend["items"]
        fs = int(legend.get("fontSize", 15))
        total = sum(66 + len(str(i["label"])) * 7 for i in items)
        x = float(legend.get("x", width / 2)) - total / 2
        y = float(legend.get("y", height - 28))
        parts.append('<g id="legend">')
        for item in items:
            color = item.get("color", "#595959")
            parts.append(_svg_line(x, y, x + 30, y, color, 3.0, "9 5"))
            x += 38
            label = escape(str(item.get("label", item.get("type", ""))))
            parts.append(
                f'<text x="{x:.3f}" y="{y + fs * 0.32:.3f}" font-family="Arial,Helvetica,sans-serif" '
                f'font-size="{fs}" font-weight="700" fill="#111111">{label}</text>'
            )
            x += 28 + len(str(item.get("label", ""))) * 7
        parts.append("</g>")

    parts.append("</svg>")
    return "\n".join(parts)


def write_static_exports(
    scene: dict,
    svg_path: str | Path,
    png_path: str | Path,
    png_scale: int = 2,
):
    svg_path = Path(svg_path)
    png_path = Path(png_path)
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    svg_text = scene_to_svg(scene)
    svg_path.write_text(svg_text, encoding="utf-8")

    try:
        import cairosvg
    except ImportError as exc:
        raise RuntimeError(
            "CairoSVG is required for PanViz static PNG export. Install the project requirements."
        ) from exc

    cairosvg.svg2png(
        bytestring=svg_text.encode("utf-8"),
        write_to=str(png_path),
        output_width=max(1, int(scene["width"]) * max(1, int(png_scale))),
        output_height=max(1, int(scene["height"]) * max(1, int(png_scale))),
        background_color="#ffffff",
    )
    return svg_path, png_path
