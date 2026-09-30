from __future__ import annotations

import os
import tempfile
import math
from collections import Counter
import numpy as np
import pandas as pd
import cairo
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

from utils import (
    convert_and_write_pdb,
    _get_interactions,
    _get_res_info,
    hbkeys,
    keys,
    pybel,
    set_to_neutral_pH,
    _distance_text,
)
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport

INTERACTION_TYPES = {"HPI", "HB", "PS", "PC", "SB"}
# Interaction colors retain the original PanViz semantic mapping.
# RGB(10, 255, 239) is reserved for the optional 3D glass-bubble residue node.
BUBBLE_COLOR = "#0AFFEF"  # RGB(10, 255, 239)
NONCOVALENT_COLOR = "#595959"  # neutral global fallback
INTERACTION_COLORS = {
    "HPI": "#595959",  # Hydrophobic
    "HB": "#0000E0",   # H-bond
    "PS": "#008500",   # π-Stacking
    "PC": "#E08000",   # π-Cation
    "SB": "#D900B0",   # Salt bridge
}
INTERACTION_LABELS = {
    "HPI": "Hydrophobic",
    "HB": "H-bond",
    "PS": "π-Stacking",
    "PC": "π-Cation",
    "SB": "Salt bridge",
}


def format_residue_display(residue_id: str) -> str:
    """Format an internal residue identifier for human-facing PanViz labels.

    PLIP/PanViz may internally represent a residue as ``TYR192_A``. The editor
    should display this as ``TYR192:A`` while preserving the original identifier
    internally for interaction matching and scientific-data integrity.
    """
    residue_id = str(residue_id)
    if "_" in residue_id:
        head, chain = residue_id.rsplit("_", 1)
        if head and chain and len(chain) <= 2:
            return f"{head}:{chain}"
    return residue_id


def _prepare_interaction_tables(my_interactions):
    bsr = BindingSiteReport(my_interactions)
    interactions = {
        k: [getattr(bsr, k + "_features")] + getattr(bsr, k + "_info")
        for k in keys
    }
    hydrophobic_df = pd.DataFrame(
        interactions["hydrophobic"][1:],
        columns=interactions["hydrophobic"][0],
    )
    hbond_rows = []
    for hb in my_interactions.all_hbonds_pdon + my_interactions.all_hbonds_ldon:
        hbond_rows.append([getattr(hb, k) for k in hbkeys])
    if hbond_rows:
        hbond_df = pd.DataFrame(np.stack(hbond_rows), columns=hbkeys)
        hbond_df["h"] = [x.idx for x in hbond_df["h"]]
    else:
        hbond_df = pd.DataFrame()
    pi_stacking_df = pd.DataFrame(
        interactions["pistacking"][1:],
        columns=interactions["pistacking"][0],
    )
    pi_cation_df = pd.DataFrame(
        interactions["pication"][1:],
        columns=interactions["pication"][0],
    )
    saltbridge_df = pd.DataFrame(
        interactions["saltbridge"][1:],
        columns=interactions["saltbridge"][0],
    )
    return hydrophobic_df, hbond_df, pi_stacking_df, pi_cation_df, saltbridge_df


def _extract_ligand_mol(file_prot, bsid):
    with open(file_prot, "r", encoding="utf-8") as handle:
        pdb_lines = [line for line in handle.readlines() if line.startswith(("ATOM", "HETATM"))]
    parts = bsid.split(":")
    lig = "".join(
        line
        for line in pdb_lines
        if line[17:20] == parts[0]
        and line[21] == parts[1]
        and line[22:26].strip() == parts[2]
    )
    with tempfile.TemporaryDirectory() as td:
        lig_path = os.path.join(td, "lig.pdb")
        try:
            mol = Chem.MolFromPDBBlock(lig, removeHs=False)
            if mol is None:
                raise ValueError("RDKit could not read ligand PDB block")
            rdDetermineBonds.DetermineBonds(mol, charge=0)
        except Exception:
            py_mol = pybel.readstring("pdb", lig)
            py_mol.write("pdb", lig_path, overwrite=True)
            mol = Chem.MolFromPDBFile(lig_path, removeHs=False)
    if mol is None:
        raise ValueError("Unable to construct ligand molecule from the PDB complex.")
    AllChem.EmbedMolecule(mol)
    set_to_neutral_pH(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetBonds():
            bound_atom = atom.GetBonds()[0].GetOtherAtom(atom)
            if bound_atom.GetSymbol() in {"O", "N", "S"}:
                atom.SetAtomicNum(100)
    mol = Chem.RemoveHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 100:
            atom.SetAtomicNum(1)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
    except Exception:
        pass
    Chem.rdDepictor.Compute2DCoords(mol)
    return mol


def _boxes_overlap(cx, cy, w, h, bx, by, bw, bh, gap=4):
    return not (
        cx + w / 2 + gap < bx - bw / 2
        or cx - w / 2 - gap > bx + bw / 2
        or cy + h / 2 + gap < by - bh / 2
        or cy - h / 2 - gap > by + bh / 2
    )


def _point_line_distance(px, py, ax, ay, bx, by):
    vx, vy = bx - ax, by - ay
    wx, wy = px - ax, py - ay
    vv = vx * vx + vy * vy
    if vv <= 1e-12:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, (wx * vx + wy * vy) / vv))
    qx, qy = ax + t * vx, ay + t * vy
    return math.hypot(px - qx, py - qy)


def _exact_panviz_layout(data_points, connections, canvas_width, canvas_height, padding=35.0):
    """Replicate the established PanViz annotation-placement geometry."""
    base = min(canvas_width, canvas_height)
    residue_font_size = max(17, min(23, int(base * 0.025)))
    distance_font_size = max(14, min(18, int(base * 0.018)))

    min_x = min(p[0] for p in data_points)
    min_y = min(p[1] for p in data_points)
    max_x = max(p[0] for p in data_points)
    max_y = max(p[1] for p in data_points)
    x_range = max(max_x - min_x, 1e-6)
    y_range = max(max_y - min_y, 1e-6)
    usable_w = max(canvas_width - 2.0 * padding, 1.0)
    usable_h = max(canvas_height - 2.0 * padding, 1.0)
    scale = min(usable_w / x_range, usable_h / y_range)
    cx0, cy0 = canvas_width / 2.0, canvas_height / 2.0
    dcx, dcy = (min_x + max_x) / 2.0, (min_y + max_y) / 2.0

    def screen_xy(x, y):
        return cx0 + (x - dcx) * scale, cy0 + (y - dcy) * scale

    atom_screen, residue_points = [], []
    for x, y, label, res in data_points:
        sx, sy = screen_xy(x, y)
        if label == "residue":
            residue_points.append((sx, sy, res))
        elif label != "centroid":
            atom_screen.append((sx, sy, 9 if label != "C" else 7))
    lig_cx = np.mean([x for x, _, _ in atom_screen]) if atom_screen else cx0
    lig_cy = np.mean([y for _, y, _ in atom_screen]) if atom_screen else cy0

    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, 10, 10)
    ctx = cairo.Context(surf)
    ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(residue_font_size)

    placed_residue_boxes, residue_draw = [], []
    residue_points.sort(key=lambda p: math.hypot(p[0] - lig_cx, p[1] - lig_cy))
    for anchor_x, anchor_y, res in residue_points:
        ext = ctx.text_extents(res)
        tw, th = ext[2], ext[3]
        vx, vy = anchor_x - lig_cx, anchor_y - lig_cy
        vlen = math.hypot(vx, vy)
        if vlen < 1e-6:
            vx, vy, vlen = 1.0, 0.0, 1.0
        ux, uy = vx / vlen, vy / vlen
        nx, ny = -uy, ux
        candidates = []
        for radial in (20, 28, 38, 50, 64, 80, 98):
            for tangential in (0, 8, -8, 16, -16):
                candidates.append(
                    (anchor_x + ux * radial + nx * tangential,
                     anchor_y + uy * radial + ny * tangential)
                )
        for radial in (24, 36, 52):
            candidates.append((anchor_x - ux * radial, anchor_y - uy * radial))
        chosen, best = None, float("inf")
        for px, py in candidates:
            if (
                px - tw / 2 < 10
                or px + tw / 2 > canvas_width - 10
                or py - th / 2 < 10
                or py + th / 2 > canvas_height - 10
            ):
                continue
            score = math.hypot(px - anchor_x, py - anchor_y) * 0.30
            for ax, ay, radius in atom_screen:
                d = math.hypot(px - ax, py - ay)
                if d < max(13, th * 0.90):
                    score += 1200
                elif d < 28:
                    score += (28 - d) * 12
            if any(_boxes_overlap(px, py, tw, th, bx, by, bw, bh, gap=8)
                   for bx, by, bw, bh in placed_residue_boxes):
                continue
            score -= ((px - lig_cx) * ux + (py - lig_cy) * uy) * 0.06
            if score < best:
                best, chosen = score, (px, py)
        if chosen is None:
            chosen = (
                max(tw / 2 + 10, min(canvas_width - tw / 2 - 10, anchor_x)),
                max(th / 2 + 10, min(canvas_height - th / 2 - 10, anchor_y)),
            )
        px, py = chosen
        placed_residue_boxes.append((px, py, tw, th))
        dx, dy = px - anchor_x, py - anchor_y
        length = math.hypot(dx, dy)
        leader = length > 26
        residue_draw.append({
            "text": res,
            "x": px,
            "y": py,
            "anchorX": anchor_x,
            "anchorY": anchor_y,
            "w": tw,
            "h": th,
            "leader": leader,
        })

    interaction_lines = []
    interaction_distance_lines = []
    shifts = {"HPI": 0, "HB": 6, "PS": 10, "PC": -10, "SB": -5}
    for idx, connection in enumerate(connections):
        a, b, kind = connection[:3]
        distance = connection[3] if len(connection) >= 4 else None
        if kind not in INTERACTION_TYPES:
            continue
        x0, y0, _, _ = data_points[a]
        x1, y1, _, _ = data_points[b]
        sx0, sy0 = screen_xy(x0, y0)
        sx1, sy1 = screen_xy(x1, y1)
        length = math.hypot(sx1 - sx0, sy1 - sy0)
        if length < 1e-9:
            continue
        dx, dy = sx1 - sx0, sy1 - sy0
        shift = shifts[kind]
        ox, oy = (-dy / length) * shift, (dx / length) * shift
        row_idx = connection[4] if len(connection) >= 5 else None
        interaction_lines.append({
            "index": idx,
            "row_idx": row_idx,
            "type": kind,
            "x1": sx0 + ox,
            "y1": sy0 + oy,
            "x2": sx1,
            "y2": sy1,
            "distance": None if distance is None else float(distance),
        })
        interaction_distance_lines.append({
            "index": idx,
            "row_idx": row_idx,
            "type": kind,
            "x1": sx0,
            "y1": sy0,
            "x2": sx1,
            "y2": sy1,
            "distance": None if distance is None else float(distance),
        })

    ctx.set_font_size(distance_font_size)
    placed_distance_boxes = []
    for line in sorted(
        interaction_distance_lines,
        key=lambda z: math.hypot(z["x2"] - z["x1"], z["y2"] - z["y1"]),
    ):
        if line["distance"] is None:
            continue
        x0, y0, x1, y1 = line["x1"], line["y1"], line["x2"], line["y2"]
        dx, dy = x1 - x0, y1 - y0
        length = math.hypot(dx, dy)
        if length < 1e-9:
            continue
        tx, ty = dx / length, dy / length
        nx, ny = -ty, tx
        text = _distance_text(line["distance"])
        ext = ctx.text_extents(text)
        tw, th = ext[2], ext[3]
        candidates = []
        for frac in (0.50, 0.43, 0.57, 0.35, 0.65, 0.25, 0.75):
            for off in (0, 4, -4, 8, -8, 14, -14, 20, -20):
                mx, my = x0 + dx * frac, y0 + dy * frac
                candidates.append((mx + nx * off, my + ny * off))
        chosen, best = None, float("inf")
        for px, py in candidates:
            if (
                px - tw / 2 < 7
                or px + tw / 2 > canvas_width - 7
                or py - th / 2 < 7
                or py + th / 2 > canvas_height - 7
            ):
                continue
            score = _point_line_distance(px, py, x0, y0, x1, y1) * 1.2
            for bx, by, bw, bh in placed_residue_boxes:
                if _boxes_overlap(px, py, tw, th, bx, by, bw, bh, gap=4):
                    score += 1000
            if any(_boxes_overlap(px, py, tw, th, bx, by, bw, bh, gap=4)
                   for bx, by, bw, bh in placed_distance_boxes):
                score += 700
            if score < best:
                best, chosen = score, (px, py)
        if chosen is None:
            continue
        px, py = chosen
        placed_distance_boxes.append((px, py, tw, th))
        line.update({
            "distanceX": px,
            "distanceY": py,
            "distanceW": tw,
            "distanceH": th,
            "distanceText": text,
        })

    return {
        "screen_xy": screen_xy,
        "residue_draw": residue_draw,
        "interaction_lines": interaction_lines,
        "distance_lines": {x["index"]: x for x in interaction_distance_lines if "distanceX" in x},
        "residue_font_size": residue_font_size,
        "distance_font_size": distance_font_size,
    }


def build_editor_scene(pdb_file, bsid, width=1200, height=850, base_svg=None, analysis=None):
    """Create the v5 scene.

    v5 deliberately does NOT edit or transform the original annotation nodes.
    The molecule is rendered once into a locked SVG layer. Editable annotations are
    separate SVG objects, one object per residue, leader line, interaction and
    distance label. This removes the stale/ghost "fingerprint" caused by moving
    children of a static original SVG.
    """
    root = tempfile.mkdtemp(prefix="panviz_scene_")
    if analysis is not None:
        file_prot = os.path.abspath(str(analysis.get("file_prot")))
        if not os.path.exists(file_prot):
            raise RuntimeError("The reusable PanViz analysis structure is no longer available.")
        input_pdb = analysis["input_pdb"]
        my_mol = analysis["my_mol"]
        my_interactions = analysis["my_interactions"]
        hpi = analysis["hydrophobic_df"]
        hb = analysis["hbond_df"]
        ps = analysis["pi_stacking_df"]
        pc = analysis["pi_cation_df"]
        sb = analysis["saltbridge_df"]
    else:
        structures = os.path.join(root, "structures")
        os.makedirs(structures, exist_ok=True)
        file_prot = os.path.join(structures, "PanViz_prot.pdb")
        input_pdb = convert_and_write_pdb(pdb_file, file_prot, bsid)
        my_mol = PDBComplex()
        my_mol.load_pdb(file_prot)
        my_mol.analyze()
        if bsid not in my_mol.interaction_sets:
            raise RuntimeError(f"Binding site {bsid} was not found after PLIP analysis.")
        my_interactions = my_mol.interaction_sets[bsid]
        hpi, hb, ps, pc, sb = _prepare_interaction_tables(my_interactions)
    mol = _extract_ligand_mol(file_prot, bsid)

    raw_atoms, bonds = [], []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        name = info.GetName().strip() if info else f"A{atom.GetIdx() + 1}"
        p = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        raw_atoms.append((float(p.x), float(p.y), atom.GetSymbol(), name))
    for bond in mol.GetBonds():
        bt = str(bond.GetBondType()).split(".")[-1]
        bonds.append({
            "id": f"b{bond.GetBeginAtomIdx()}_{bond.GetEndAtomIdx()}",
            "a": bond.GetBeginAtomIdx(),
            "b": bond.GetEndAtomIdx(),
            "order": 2 if bt in {"DOUBLE", "AROMATIC"} else 3 if bt == "TRIPLE" else 1,
            "aromatic": bt == "AROMATIC",
        })

    coord_dict = {name: (x, y) for x, y, _, name in raw_atoms}
    interaction_rows, centroids, used_res = _get_interactions(
        input_pdb, hpi, hb, ps, pc, sb, coord_dict
    )
    res_info = _get_res_info(used_res, coord_dict, interaction_rows)
    data_points = [(x, y, symbol, name) for x, y, symbol, name in raw_atoms]
    data_points += [(c[0], c[1], c[2], c[3]) for c in centroids]
    data_points += list(res_info)
    atom_index_by_name = {p[3]: i for i, p in enumerate(data_points)}
    connections = [(b["a"], b["b"], "DOUBLE" if b["order"] == 2 else "TRIPLE" if b["order"] == 3 else "SINGLE") for b in bonds]
    for row_idx, (atom_name, residue, kind, distance) in enumerate(interaction_rows):
        aidx = atom_index_by_name.get(atom_name)
        ridx = atom_index_by_name.get(residue)
        if aidx is not None and ridx is not None:
            connections.append((aidx, ridx, kind, distance, row_idx))

    layout = _exact_panviz_layout(data_points, connections, width, height)
    sx = layout["screen_xy"]

    atoms = []
    for i, (x, y, symbol, name) in enumerate(raw_atoms):
        px, py = sx(x, y)
        atoms.append({
            "id": i,
            "name": name,
            "element": symbol,
            "x": float(px),
            "y": float(py),
            "showLabel": symbol != "C",
        })

    residues = []
    residue_id_by_source = {}
    for i, r in enumerate(layout["residue_draw"]):
        rid = f"res_{i}"
        source_residue = str(r["text"])
        display_residue = format_residue_display(source_residue)
        residue_id_by_source.setdefault(source_residue, rid)
        residues.append({
            "id": rid,
            "text": display_residue,
            "sourceResidue": source_residue,
            "x": float(r["x"]),
            "y": float(r["y"]),
            "anchorX": float(r["anchorX"]),
            "anchorY": float(r["anchorY"]),
            "w": float(r["w"]),
            "h": float(r["h"]),
            "fontSize": int(layout["residue_font_size"]),
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
            "bubbleColor": BUBBLE_COLOR,
            "leaderVisible": bool(r["leader"]),
            "visible": True,
        })

    interactions, distances = [], []
    # Count how many PLIP interactions connect the same residue to the same
    # ligand atom.  Multiple interaction records for one endpoint pair should
    # remain scientifically distinct (and retain their separate distances),
    # but their rendered connector must become visibly thicker.
    pair_counts = Counter(
        (str(residue), str(atom_name))
        for atom_name, residue, _kind, _distance in interaction_rows
    )
    for line in layout["interaction_lines"]:
        row_idx = line.get("row_idx")
        if row_idx is None or row_idx >= len(interaction_rows):
            continue
        atom_name, residue, kind, distance = interaction_rows[row_idx]
        rid = residue_id_by_source.get(str(residue))
        if rid is None:
            continue
        iid = f"int_{row_idx}"
        multiplicity = max(1, int(pair_counts.get((str(residue), str(atom_name)), 1)))
        interactions.append({
            "id": iid,
            "type": kind,
            "residueId": rid,
            "anchorAtom": atom_name,
            "x1": float(line["x1"]),
            "y1": float(line["y1"]),
            "x2": float(line["x2"]),
            "y2": float(line["y2"]),
            "originalDistance": None if distance is None else float(distance),
            "multiplicity": multiplicity,
            "color": INTERACTION_COLORS[kind],
            "visible": True,
            "lineWidth": 3.0,
            "customLineWidth": False,
            "customColor": False,
            "customDash": False,
            "dash": "9 5",
            "rotation": 0,
        })
        dl = layout["distance_lines"].get(line["index"])
        if dl is not None and distance is not None:
            distances.append({
                "id": f"dist_{row_idx}",
                "interactionId": iid,
                "originalDistance": float(distance),
                "displayText": dl["distanceText"],
                "x": float(dl["distanceX"]),
                "y": float(dl["distanceY"]),
                "w": float(dl["distanceW"]),
                "h": float(dl["distanceH"]),
                "fontSize": int(layout["distance_font_size"]),
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
            })

    present_types = []
    for it in interactions:
        if it["type"] not in present_types:
            present_types.append(it["type"])
    legend_items = [
        {"type": k, "label": INTERACTION_LABELS[k], "color": INTERACTION_COLORS[k]}
        for k in present_types
        if k in INTERACTION_TYPES
    ]
    # Publication-friendly legend: larger, bold, and padded for legibility.
    legend_font_size = 15 if height <= 800 else 17
    item_gap = 24
    sample_width = 32
    text_gap = 8
    horizontal_padding = 12
    vertical_padding = 8
    estimated_item_widths = [sample_width + text_gap + max(54, 8.4 * len(x["label"])) for x in legend_items]
    legend_w = sum(estimated_item_widths) + item_gap * max(0, len(legend_items) - 1) + horizontal_padding * 2
    legend_h = legend_font_size + vertical_padding * 2 + 4
    legend = {
        "x": float(width / 2.0),
        "y": float(height - 28.0),
        "dx": 0.0,
        "dy": 0.0,
        "w": float(max(legend_w, 100)),
        "h": float(legend_h),
        "fontSize": int(legend_font_size),
        "fontFamily": "Arial",
        "bold": True,
        "italic": False,
        "underline": False,
        "rotation": 0,
        "items": legend_items,
        "visible": True,
        "movable": True,
    }

    return {
        "version": "5.8.6",
        "width": int(width),
        "height": int(height),
        "style": {
            "bondWidth": 2.0,
            "interactionWidth": 3.0,
            "noncovalentColor": NONCOVALENT_COLOR,
            "moleculeLabelSize": max(14, min(19, int(min(width, height) * 0.019))),
            "residue": {
                "labelFormat": "nameNumChain",
                "fontFamily": "Arial",
                "fontSize": int(layout["residue_font_size"]),
                "fontWeight": 700,
                "textColor": "#000000",
                "shape": "bubble",
                "backgroundColor": "#ffffff",
                "backgroundOpacity": 0.98,
                "bubbleColor": BUBBLE_COLOR,
                "colorMode": "neutral",
            },
            "distance": {
                "fontFamily": "Arial",
                "fontSize": int(layout["distance_font_size"]),
                "fontWeight": 700,
                "textColor": "#D100A0",
                "backgroundColor": "#ffffff",
                "backgroundOpacity": 0.96,
            },
            "interactions": {
                k: {"color": INTERACTION_COLORS[k], "width": 3.0, "dash": "9 5", "custom": False}
                for k in INTERACTION_TYPES
            },
        },
        "site": bsid,
        "atoms": atoms,
        "bonds": bonds,
        "labels": residues,
        "interactions": interactions,
        "distances": distances,
        "legend": legend,
        "metadata": {
            "editorMode": "layer-separated-scene",
            "chemistryEditable": False,
            "bondGraphicsEditable": True,
            "annotationLayerFullyEditable": True,
            "plipDataImmutable": True,
            "initialView": "panviz-layout-scene",
            "canvasResize": "proportional-content-scale-fit",
            "canvasScrollbars": False,
            "contextualRibbon": "external-overlay",
            "ribbonFrozenDuringDrag": True,
            "separateInteractionWidthControl": True,
            "semanticResidueNodeStyles": ["plain", "circle", "plate", "rounded", "bubble"],
            "globalLocalStyleHierarchy": True,
            "globalDefaults": {"covalentWidth": 2.0, "noncovalentWidth": 3.0, "bubbleColor": BUBBLE_COLOR},
            "duplicateInteractionWeight": "baseWidth × interactionMultiplicity",
            "smartResidueLeaderAttachment": True,
            "moleculeLocked": True,
            "annotationOnly": True,
            "ghostFree": True,
            "reusesCompletedPLIPAnalysis": analysis is not None,
        },
        "scientificData": {
            "interactions": [
                {
                    "id": x["id"],
                    "type": x["type"],
                    "residueId": x["residueId"],
                    "anchorAtom": x["anchorAtom"],
                    "originalDistance": x["originalDistance"],
                    "multiplicity": x.get("multiplicity", 1),
                }
                for x in interactions
            ]
        },
    }, root
