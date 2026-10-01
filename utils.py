"""PanViz visualization and PLIP interaction-processing utilities.

PanViz uses PLIP for interaction detection and provides the customized
2D renderer plus optional PyMOL session generation. The visualization
logic is retained from the supplied original PLIPViz utility.
"""

print("UTILCHK 00 — utils.py entered", flush=True)
print("UTILCHK 00R — preloading RDKit before PLIP/Open Babel", flush=True)
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
print("UTILCHK 00S — RDKit preload OK", flush=True)

from plip.basic import config
print("UTILCHK 01 — plip.basic.config OK", flush=True)
from plip.structure.preparation import PDBComplex
print("UTILCHK 02 — PDBComplex OK", flush=True)
from plip.exchange.report import BindingSiteReport
print("UTILCHK 03 — BindingSiteReport OK", flush=True)

print("UTILCHK 04 — before optional PLIP visualization imports", flush=True)
try:
    from plip.visualization.visualize import PyMOLVisualizer
    print("UTILCHK 05a — PyMOLVisualizer OK", flush=True)
    from plip.basic.remote import VisualizerData
    print("UTILCHK 05b — VisualizerData OK", flush=True)
    from plip.basic.supplemental import start_pymol
    print("UTILCHK 05c — start_pymol OK", flush=True)
except Exception as exc:
    print(f"UTILCHK 05x — optional PLIP visualization unavailable: {type(exc).__name__}: {exc}", flush=True)
    PyMOLVisualizer = None
    VisualizerData = None
    start_pymol = None

print("UTILCHK 06 — before pandas", flush=True)
import pandas as pd
print("UTILCHK 07 — pandas OK", flush=True)
import numpy as np
print("UTILCHK 08 — numpy OK", flush=True)
import math
print("UTILCHK 09 — math OK", flush=True)

print("UTILCHK 10 — before openbabel.pybel", flush=True)
from openbabel import pybel
print("UTILCHK 11 — openbabel.pybel OK", flush=True)

print("UTILCHK 12 — RDKit already preloaded safely", flush=True)
print("UTILCHK 13 — RDKit Chem OK", flush=True)
print("UTILCHK 14 — RDKit AllChem OK", flush=True)

print("UTILCHK 15 — before cairo", flush=True)
import cairo
print("UTILCHK 16 — cairo OK", flush=True)

print("UTILCHK 17 — before optional pymol.cmd", flush=True)
try:
    from pymol import cmd
    print("UTILCHK 18a — pymol.cmd OK", flush=True)
except Exception as exc:
    print(f"UTILCHK 18x — optional pymol unavailable: {type(exc).__name__}: {exc}", flush=True)
    cmd = None

import os
import tempfile
print("UTILCHK 19 — os/tempfile OK", flush=True)
print("UTILCHK 20 — rdDetermineBonds already preloaded OK", flush=True)

config.NOHYDRO = True
ob = pybel.ob
print("UTILCHK 21 — utils.py imports complete", flush=True)

def _save_pymol(my_mol, my_id, outdir):
    '''Save a PyMOL session when PyMOL support is available.'''
    if cmd is None or PyMOLVisualizer is None or VisualizerData is None or start_pymol is None:
        raise RuntimeError(
            "PyMOL is not available in this installation. Run the web app with PyMOL disabled."
        )
    complex = VisualizerData(my_mol, my_id)
    vis = PyMOLVisualizer(complex)

    lig_members = complex.lig_members
    chain = complex.chain

    ligname = vis.ligname
    hetid = complex.hetid

    metal_ids = complex.metal_ids
    metal_ids_str = '+'.join([str(i) for i in metal_ids])

    start_pymol(run=True, options='-pcq', quiet=True)
    vis.set_initial_representations()

    pdb_path = os.path.join(outdir, complex.pdbid + "_temp.pdb"); open(pdb_path, "w").write(my_mol.corrected_pdb); cmd.load(pdb_path)

    current_name = cmd.get_object_list(selection='(all)')[0]

    current_name = cmd.get_object_list(selection='(all)')[0]
    cmd.set_name(current_name, complex.pdbid)

    cmd.hide('everything', 'all')
    cmd.select(ligname, 'resn %s and chain %s and resi %s*' % (hetid, chain, complex.position))


    # Visualize and color metal ions if there are any
    if not len(metal_ids) == 0:
        vis.select_by_ids(ligname, metal_ids, selection_exists=True)
        cmd.show('spheres', 'id %s and %s' % (metal_ids_str))

    # Additionally, select all members of composite ligands
    if len(lig_members) > 1:
        for member in lig_members:
            resid, chain, resnr = member[0], member[1], str(member[2])
            cmd.select(ligname, '%s or (resn %s and chain %s and resi %s)' % (ligname, resid, chain, resnr))

    cmd.show('sticks', ligname)
    cmd.color('myblue')
    cmd.color('myorange', ligname)
    cmd.util.cnc('all')
    if not len(metal_ids) == 0:
        cmd.color('hotpink', 'id %s' % metal_ids_str)
        cmd.hide('sticks', 'id %s' % metal_ids_str)
        cmd.set('sphere_scale', 0.3, ligname)
    cmd.deselect()

    vis.make_initial_selections()
    vis.show_hydrophobic()  # Hydrophobic Contacts
    vis.show_hbonds()  # Hydrogen Bonds
    vis.show_halogen()  # Halogen Bonds
    vis.show_stacking()  # pi-Stacking Interactions
    vis.show_cationpi()  # pi-Cation Interactions
    vis.show_sbridges()  # Salt Bridges
    vis.show_wbridges()  # Water Bridges
    vis.show_metal()  # Metal Coordination
    vis.refinements()
    vis.zoom_to_ligand()
    vis.selections_cleanup()
    vis.selections_group()
    vis.additional_cleanup()
    vis.save_session(outdir)


def _distance_from_row(row, *names):
    """Return the first available numeric interaction distance in Angstrom."""
    for name in names:
        try:
            if name in row.index:
                value = float(row[name])
                if np.isfinite(value):
                    return value
        except (TypeError, ValueError):
            pass
    return None


def _distance_text(distance):
    """Format a PLIP interaction distance for the figure."""
    if distance is None:
        return None
    return f"{float(distance):.2f} Å"

def _get_interactions(
        input_pdb,
        hydrophobic_df,
        hbond_df,
        pi_stacking_df,
        pi_cation_df,
        saltbridge_df,
        waterbridge_df,
        halogen_df,
        metal_df,
        coord_dict):
    interactions = []
    centroids = []
    centroid_counter = 0

    # Hydrogen bonds: use PLIP donor-acceptor distance (distance_ad).
    for i, row in hbond_df.iterrows():
        if row.protisdon:
            atom = input_pdb.atoms[row.a_orig_idx-1].OBAtom
            int_atom = atom.GetResidue().GetAtomID(atom).strip()
        else:
            atom = input_pdb.atoms[row.d_orig_idx-1].OBAtom
            int_atom = atom.GetResidue().GetAtomID(atom).strip()

        distance = _distance_from_row(row, "distance_ad", "distance_ah")
        interactions.append((
            int_atom,
            row["restype"] + str(row["resnr"]) + "_" + row["reschain"],
            "HB",
            distance
        ))

    # Hydrophobic contacts: PLIP normally reports DIST.
    for i, row in hydrophobic_df.drop_duplicates(
            subset=["LIGCARBONIDX", "RESTYPE", "RESNR", "RESCHAIN"]).iterrows():
        atom = input_pdb.atoms[row.LIGCARBONIDX-1].OBAtom
        int_atom = atom.GetResidue().GetAtomID(atom).strip()
        distance = _distance_from_row(row, "DIST", "distance", "dist")
        interactions.append((
            int_atom,
            row["RESTYPE"] + str(row["RESNR"]) + "_" + row["RESCHAIN"],
            "HPI",
            distance
        ))

    # Pi-stacking, pi-cation and salt-bridge interactions.
    for df, int_type in zip(
            [pi_stacking_df, pi_cation_df, saltbridge_df],
            ["PS", "PC", "SB"]):
        for i, row in df.drop_duplicates(
                subset=["LIG_IDX_LIST", "RESTYPE", "RESNR", "RESCHAIN"]).iterrows():

            atoms = [
                input_pdb.atoms[x-1].OBAtom
                for x in np.array(row["LIG_IDX_LIST"].split(","), dtype=int)
            ]
            com = np.stack([
                coord_dict[atom.GetResidue().GetAtomID(atom).strip()]
                for atom in atoms
            ]).mean(axis=0)

            centroid_name = "centroid_{}".format(centroid_counter)
            centroid_key = (com[0], com[1], "centroid", centroid_name)

            if centroid_key not in centroids:
                centroids.append(centroid_key)
                coord_dict[centroid_name] = (com[0], com[1])
                centroid_counter += 1
                current_centroid = "centroid_{}".format(centroid_counter - 1)
            else:
                current_centroid = centroid_name

            distance = _distance_from_row(row, "DIST", "distance", "dist")
            interactions.append((
                current_centroid,
                row["RESTYPE"] + str(row["RESNR"]) + "_" + row["RESCHAIN"],
                int_type,
                distance
            ))

    # Water bridges: anchor at the ligand donor/acceptor atom reported by PLIP.
    for _, row in waterbridge_df.iterrows():
        protisdon = bool(row["PROTISDON"])
        ligand_idx = int(row["ACCEPTOR_IDX"] if protisdon else row["DONOR_IDX"])
        atom = input_pdb.atoms[ligand_idx - 1].OBAtom
        int_atom = atom.GetResidue().GetAtomID(atom).strip()
        distance = _distance_from_row(
            row,
            "DIST_A-W" if protisdon else "DIST_D-W",
            "DIST_D-W",
            "DIST_A-W",
        )
        interactions.append((
            int_atom,
            row["RESTYPE"] + str(row["RESNR"]) + "_" + row["RESCHAIN"],
            "WB",
            distance
        ))

    # Halogen bonds: PLIP defines the halogen donor as the ligand-side atom.
    for _, row in halogen_df.iterrows():
        atom = input_pdb.atoms[int(row["DON_IDX"]) - 1].OBAtom
        int_atom = atom.GetResidue().GetAtomID(atom).strip()
        distance = _distance_from_row(row, "DIST")
        interactions.append((
            int_atom,
            row["RESTYPE"] + str(row["RESNR"]) + "_" + row["RESCHAIN"],
            "XB",
            distance
        ))

    # Metal coordination: use the PLIP metal atom as the ligand-side anchor.
    for _, row in metal_df.iterrows():
        atom = input_pdb.atoms[int(row["METAL_IDX"]) - 1].OBAtom
        int_atom = atom.GetResidue().GetAtomID(atom).strip()
        distance = _distance_from_row(row, "DIST")
        interactions.append((
            int_atom,
            row["RESTYPE"] + str(row["RESNR"]) + "_" + row["RESCHAIN"],
            "MC",
            distance
        ))

    used_res = np.unique(np.array([x[1] for x in interactions]))
    return interactions, centroids, used_res

def _get_res_info(used_res, coord_dict, interactions):
    '''
    Greedy label placement algorithm.
    Calculates number of overlaps in a grid of candidate positions.
    Returns a solution with fewest overlaps.
    '''
    res_info = []
    _coord_temp = np.array([np.array(coord_dict[key]) for key in coord_dict.keys()])
    for i,res in enumerate(used_res):
        res_ints = [x[0] for x in interactions if x[1] == res]
        initial_lab_coords = np.array([coord_dict[x] for x in res_ints]).mean(axis=0)
        candidates = (np.mgrid[-3.5:3.6:0.1, -3.5:3.6:0.1].reshape(2,-1).T) + initial_lab_coords

        overlaps = []
        for cand in candidates:
            overlap_count = 0
            for coord in _coord_temp:
                if np.abs(cand-coord)[0] < 2 and np.abs(cand-coord)[1] < 2:
                    overlap_count += 1
            for coord in [(x[0],x[1]) for x in res_info]:
                if np.abs(cand-coord)[0] < 2.5 and np.abs(cand-coord)[1] < 2:
                    overlap_count += 1

            overlaps.append(overlap_count)

        lab_coords = candidates[np.argmin(overlaps)]

        coord_dict[res] = (lab_coords[0],lab_coords[1])
        res_info.append((lab_coords[0],lab_coords[1],"residue",res))

    return res_info


def _draw_mol(atom_info, connections, padding, canvas_height, canvas_width, out_name):
    """
    Publication-quality PLIP 2D renderer.

    Improvements over the original renderer:
      * larger default working area
      * compact but controlled whitespace
      * stronger residue typography
      * collision-aware residue labels
      * collision-aware interaction-distance labels
      * white knockout behind residue labels
      * all labels constrained inside the canvas
      * short leader lines for displaced residue labels
      * darker, higher-contrast interaction graphics
      * stronger chemical bond stroke hierarchy
      * bold magenta interaction-distance labels
      * true parallel double/triple bond rendering
      * Kekule-style aromatic bond representation
    """
    padding = max(float(padding), 35.0)

    # Typography and stroke weights tuned for a publication-style canvas.
    base = min(canvas_width, canvas_height)
    residue_font_size = max(17, min(23, int(base * 0.025)))
    atom_font_size = max(14, min(19, int(base * 0.019)))
    distance_font_size = max(14, min(18, int(base * 0.018)))

    # Global PanViz defaults shared by the original renderer and editor.
    # Covalent chemistry is rendered at 2.0 px; non-covalent interactions at 3.0 px.
    single_width = 2.0
    double_width = 2.0
    double_offset = max(5.0, base * 0.0054)
    triple_width = 2.0
    triple_offset = max(5.0, base * 0.0052)
    interaction_width = 3.0

    color_dict = {
        "O": (0.85, 0.0, 0.0),
        "N": (0.0, 0.0, 0.80),
        "S": (0.72, 0.58, 0.02),
        "P": (0.85, 0.38, 0.0),
        "B": (0.90, 0.45, 0.45),
    }

    if out_name is None:
        surface = cairo.ImageSurface(
            cairo.FORMAT_ARGB32, canvas_width, canvas_height
        )
    elif out_name.lower().endswith(".svg"):
        surface = cairo.SVGSurface(
            out_name, canvas_width, canvas_height
        )
    else:
        surface = cairo.ImageSurface(
            cairo.FORMAT_ARGB32, canvas_width, canvas_height
        )

    ctx = cairo.Context(surface)

    # White background for conventional publication-style PLIP diagrams.
    ctx.set_source_rgb(1, 1, 1)
    ctx.paint()
    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    ctx.set_line_join(cairo.LINE_JOIN_ROUND)

    data_points = atom_info

    min_x = min(p[0] for p in data_points)
    min_y = min(p[1] for p in data_points)
    max_x = max(p[0] for p in data_points)
    max_y = max(p[1] for p in data_points)

    x_range = max(max_x - min_x, 1e-6)
    y_range = max(max_y - min_y, 1e-6)

    # Fit the complete diagram while retaining a controlled margin.
    usable_w = max(canvas_width - 2.0 * padding, 1.0)
    usable_h = max(canvas_height - 2.0 * padding, 1.0)
    scale = min(usable_w / x_range, usable_h / y_range)

    center_screen_x = canvas_width / 2.0
    center_screen_y = canvas_height / 2.0
    data_center_x = (min_x + max_x) / 2.0
    data_center_y = (min_y + max_y) / 2.0

    def screen_xy(x, y):
        return (
            center_screen_x + (x - data_center_x) * scale,
            center_screen_y + (y - data_center_y) * scale,
        )

    # ---------------------------------------------------------
    # Draw ligand bonds and interaction lines.
    # ---------------------------------------------------------
    for connection in connections:
        start, end, bond_type = connection[:3]

        x0, y0, _, _ = data_points[start]
        x1, y1, _, _ = data_points[end]
        x0, y0 = screen_xy(x0, y0)
        x1, y1 = screen_xy(x1, y1)

        if bond_type == "SINGLE":
            ctx.set_source_rgb(0, 0, 0)
            ctx.set_line_width(single_width)
            ctx.set_dash([], 0)
            ctx.move_to(x0, y0)
            ctx.line_to(x1, y1)
            ctx.stroke()

        elif bond_type in {"DOUBLE", "AROMATIC"}:
            # Draw true parallel bond strokes instead of a thick line + white knockout.
            # Aromatic bonds are kekulized before rendering, so valid aromatic systems
            # appear as alternating single/double bonds (e.g. benzene = 3 double bonds).
            dx, dy = x1 - x0, y1 - y0
            length = math.hypot(dx, dy)
            if length < 1e-9:
                length = 1.0
            nx, ny = -dy / length, dx / length
            offset = double_offset / 2.0

            ctx.set_source_rgb(0, 0, 0)
            ctx.set_line_width(double_width)
            ctx.set_dash([], 0)

            ctx.move_to(x0 + nx * offset, y0 + ny * offset)
            ctx.line_to(x1 + nx * offset, y1 + ny * offset)
            ctx.stroke()

            ctx.move_to(x0 - nx * offset, y0 - ny * offset)
            ctx.line_to(x1 - nx * offset, y1 - ny * offset)
            ctx.stroke()

        elif bond_type == "TRIPLE":
            # Three clean, parallel strokes for triple bonds.
            dx, dy = x1 - x0, y1 - y0
            length = math.hypot(dx, dy)
            if length < 1e-9:
                length = 1.0
            nx, ny = -dy / length, dx / length

            ctx.set_source_rgb(0, 0, 0)
            ctx.set_line_width(triple_width)
            ctx.set_dash([], 0)

            for offset in (-triple_offset, 0.0, triple_offset):
                ctx.move_to(x0 + nx * offset, y0 + ny * offset)
                ctx.line_to(x1 + nx * offset, y1 + ny * offset)
                ctx.stroke()

        elif bond_type in {"HPI", "HB", "PS", "PC", "SB", "WB", "XB", "MC"}:
            colors = {
                "HPI": (0.35, 0.35, 0.35),
                "HB": (0.0, 0.0, 0.88),
                "PS": (0.0, 0.52, 0.0),
                "PC": (0.88, 0.52, 0.0),
                "SB": (0.85, 0.0, 0.70),
                "WB": (0.08, 0.58, 0.72),
                "XB": (0.48, 0.36, 0.78),
                "MC": (0.64, 0.34, 0.00),
            }
            color = colors[bond_type]

            ctx.set_source_rgba(color[0], color[1], color[2], 0.92)
            ctx.set_line_width(interaction_width)
            ctx.set_dash({
                "WB": [3, 3],
                "XB": [6, 4],
                "MC": [2, 3],
            }.get(bond_type, [9, 5]), 0)

            dx, dy = x1 - x0, y1 - y0
            length = math.hypot(dx, dy)
            if length < 1e-9:
                length = 1.0

            # Small endpoint displacement prevents the dashed line
            # from sitting directly on top of the atom graphic.
            shift = {
                "HPI": 0,
                "HB": 6,
                "PS": 10,
                "PC": -10,
                "SB": -5,
                "WB": 4,
                "XB": -8,
                "MC": 8,
            }[bond_type]

            ox = (-dy / length) * shift
            oy = (dx / length) * shift

            ctx.move_to(x0 + ox, y0 + oy)
            ctx.line_to(x1, y1)
            ctx.stroke()

    # ---------------------------------------------------------
    # Text utility functions.
    # ---------------------------------------------------------
    def text_extents(text, size, bold=False):
        ctx.select_font_face(
            "Sans",
            cairo.FONT_SLANT_NORMAL,
            cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL,
        )
        ctx.set_font_size(size)
        return ctx.text_extents(text)

    def boxes_overlap(cx, cy, w, h, bx, by, bw, bh, gap=4):
        return not (
            cx + w / 2 + gap < bx - bw / 2
            or cx - w / 2 - gap > bx + bw / 2
            or cy + h / 2 + gap < by - bh / 2
            or cy - h / 2 - gap > by + bh / 2
        )

    def point_line_distance(px, py, ax, ay, bx, by):
        vx, vy = bx - ax, by - ay
        wx, wy = px - ax, py - ay
        vv = vx * vx + vy * vy
        if vv <= 1e-12:
            return math.hypot(px - ax, py - ay)
        t = max(0.0, min(1.0, (wx * vx + wy * vy) / vv))
        qx, qy = ax + t * vx, ay + t * vy
        return math.hypot(px - qx, py - qy)

    # ---------------------------------------------------------
    # Determine ligand/atom screen positions.
    # ---------------------------------------------------------
    atom_screen = []
    residue_points = []

    for x, y, label, res in data_points:
        sx, sy = screen_xy(x, y)

        if label == "residue":
            residue_points.append((sx, sy, res))
        elif label not in {"centroid"}:
            radius = 9 if label != "C" else 7
            atom_screen.append((sx, sy, radius))

    if atom_screen:
        lig_cx = np.mean([x for x, y, r in atom_screen])
        lig_cy = np.mean([y for x, y, r in atom_screen])
    else:
        lig_cx = center_screen_x
        lig_cy = center_screen_y

    # ---------------------------------------------------------
    # Residue labels: place outward from ligand with collision
    # avoidance. This is the main readability improvement.
    # ---------------------------------------------------------
    ctx.select_font_face(
        "Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD
    )
    ctx.set_font_size(residue_font_size)

    placed_residue_boxes = []
    residue_draw = []

    # Give crowded/central labels priority.
    residue_points = sorted(
        residue_points,
        key=lambda p: math.hypot(p[0] - lig_cx, p[1] - lig_cy)
    )

    for anchor_x, anchor_y, res in residue_points:
        ext = text_extents(res, residue_font_size, bold=True)
        tw, th = ext[2], ext[3]

        vx, vy = anchor_x - lig_cx, anchor_y - lig_cy
        vlen = math.hypot(vx, vy)
        if vlen < 1e-6:
            vx, vy, vlen = 1.0, 0.0, 1.0

        ux, uy = vx / vlen, vy / vlen
        nx, ny = -uy, ux

        candidates = []

        # Compact radial positions first.
        for radial in (20, 28, 38, 50, 64, 80, 98):
            for tangential in (0, 8, -8, 16, -16):
                candidates.append((
                    anchor_x + ux * radial + nx * tangential,
                    anchor_y + uy * radial + ny * tangential
                ))

        # In exceptional crowded cases test the opposite side.
        for radial in (24, 36, 52):
            candidates.append((
                anchor_x - ux * radial,
                anchor_y - uy * radial
            ))

        chosen = None
        best_score = float("inf")

        for cx, cy in candidates:
            # Strict canvas boundary.
            if (
                cx - tw / 2 < 10
                or cx + tw / 2 > canvas_width - 10
                or cy - th / 2 < 10
                or cy + th / 2 > canvas_height - 10
            ):
                continue

            score = math.hypot(cx - anchor_x, cy - anchor_y) * 0.30

            # Avoid ligand atom labels/centres.
            for ax, ay, radius in atom_screen:
                d = math.hypot(cx - ax, cy - ay)
                if d < max(13, th * 0.90):
                    score += 1200
                elif d < 28:
                    score += (28 - d) * 12

            # Avoid existing residue labels.
            collision = False
            for bx, by, bw, bh in placed_residue_boxes:
                if boxes_overlap(cx, cy, tw, th, bx, by, bw, bh, gap=8):
                    collision = True
                    break
            if collision:
                continue

            # Prefer positions pointing away from ligand centre.
            radial_projection = (cx - lig_cx) * ux + (cy - lig_cy) * uy
            score -= radial_projection * 0.06

            if score < best_score:
                best_score = score
                chosen = (cx, cy)

        if chosen is None:
            chosen = (
                max(tw / 2 + 10, min(canvas_width - tw / 2 - 10, anchor_x)),
                max(th / 2 + 10, min(canvas_height - th / 2 - 10, anchor_y)),
            )

        cx, cy = chosen
        placed_residue_boxes.append((cx, cy, tw, th))
        residue_draw.append(
            (cx, cy, tw, th, res, anchor_x, anchor_y)
        )

    # Leader lines, behind text.
    for cx, cy, tw, th, res, ax, ay in residue_draw:
        dx, dy = cx - ax, cy - ay
        length = math.hypot(dx, dy)

        if length > 26:
            ux, uy = dx / length, dy / length
            sx, sy = ax + ux * 7, ay + uy * 7
            ex, ey = cx - ux * (tw / 2 + 5), cy - uy * (th / 2 + 3)

            ctx.set_source_rgba(0, 0, 0, 0.30)
            ctx.set_line_width(1.1)
            ctx.set_dash([], 0)
            ctx.move_to(sx, sy)
            ctx.line_to(ex, ey)
            ctx.stroke()

    # ---------------------------------------------------------
    # Interaction-distance labels.
    # ---------------------------------------------------------
    interaction_lines = []

    for idx, connection in enumerate(connections):
        a, b, bond_type = connection[:3]
        distance = connection[3] if len(connection) >= 4 else None

        if bond_type not in {"HPI", "HB", "PS", "PC", "SB", "WB", "XB", "MC"}:
            continue
        if distance is None:
            continue

        x0, y0, _, _ = data_points[a]
        x1, y1, _, _ = data_points[b]
        x0, y0 = screen_xy(x0, y0)
        x1, y1 = screen_xy(x1, y1)

        length = math.hypot(x1 - x0, y1 - y0)
        if length < 1e-9:
            continue

        interaction_lines.append({
            "x0": x0, "y0": y0,
            "x1": x1, "y1": y1,
            "length": length,
            "distance": distance,
        })

    # Shorter lines get priority in crowded regions.
    interaction_lines.sort(key=lambda z: z["length"])

    ctx.select_font_face(
        "Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD
    )
    ctx.set_font_size(distance_font_size)

    placed_distance_boxes = []

    for line in interaction_lines:
        x0, y0 = line["x0"], line["y0"]
        x1, y1 = line["x1"], line["y1"]
        dx, dy = x1 - x0, y1 - y0
        length = line["length"]

        tx, ty = dx / length, dy / length
        nx, ny = -ty, tx

        label = _distance_text(line["distance"])
        ext = ctx.text_extents(label)
        tw, th = ext[2], ext[3]

        candidates = []
        for frac in (0.50, 0.43, 0.57, 0.35, 0.65, 0.25, 0.75):
            for off in (0, 4, -4, 8, -8, 14, -14, 20, -20):
                mx = x0 + dx * frac
                my = y0 + dy * frac
                candidates.append((
                    mx + nx * off,
                    my + ny * off
                ))

        chosen = None
        best_score = float("inf")

        for cx, cy in candidates:
            if (
                cx - tw / 2 < 7
                or cx + tw / 2 > canvas_width - 7
                or cy - th / 2 < 7
                or cy + th / 2 > canvas_height - 7
            ):
                continue

            score = point_line_distance(
                cx, cy, x0, y0, x1, y1
            ) * 1.2

            # Avoid residue labels.
            for bx, by, bw, bh in placed_residue_boxes:
                if boxes_overlap(cx, cy, tw, th, bx, by, bw, bh, gap=4):
                    score += 1000

            # Avoid previous distance labels.
            bad = False
            for bx, by, bw, bh in placed_distance_boxes:
                if boxes_overlap(cx, cy, tw, th, bx, by, bw, bh, gap=4):
                    bad = True
                    break

            if bad:
                score += 700

            if score < best_score:
                best_score = score
                chosen = (cx, cy)

        if chosen is None:
            continue

        cx, cy = chosen
        placed_distance_boxes.append((cx, cy, tw, th))

        # White knockout keeps numbers readable.
        ctx.set_source_rgba(1, 1, 1, 0.94)
        ctx.rectangle(
            cx - tw / 2 - 3,
            cy - th / 2 - 2,
            tw + 6,
            th + 4
        )
        ctx.fill()

        # Bold magenta distance labels for high-contrast readability.
        ctx.set_source_rgb(0.82, 0.0, 0.58)
        ctx.move_to(
            cx - tw / 2,
            cy + th / 2 - 1
        )
        ctx.show_text(label)

    # ---------------------------------------------------------
    # Ligand atom labels.
    # ---------------------------------------------------------
    for x, y, label, res in data_points:
        if label in {"C", "centroid", "residue"}:
            continue

        sx, sy = screen_xy(x, y)

        # Small white knockout around heteroatom labels.
        ctx.set_source_rgba(1, 1, 1, 0.95)
        ctx.arc(sx, sy, 8.5, 0, 2 * math.pi)
        ctx.fill()

        bold = label in {"O", "N", "S", "P", "B"}

        ctx.select_font_face(
            "Sans",
            cairo.FONT_SLANT_NORMAL,
            cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL
        )
        ctx.set_font_size(atom_font_size)

        text = label
        ext = ctx.text_extents(text)

        if label in color_dict:
            color = color_dict[label]
            ctx.set_source_rgba(color[0], color[1], color[2], 1)
        else:
            ctx.set_source_rgba(0, 0, 0, 1)

        ctx.move_to(
            sx - ext[2] / 2,
            sy + ext[3] / 2
        )
        ctx.show_text(text)

    # ---------------------------------------------------------
    # Residue labels last: strongest visual hierarchy.
    # ---------------------------------------------------------
    for cx, cy, tw, th, res, ax, ay in residue_draw:
        # Opaque knockout.
        ctx.set_source_rgba(1, 1, 1, 0.98)
        ctx.rectangle(
            cx - tw / 2 - 6,
            cy - th / 2 - 5,
            tw + 12,
            th + 10
        )
        ctx.fill()

        ctx.select_font_face(
            "Sans",
            cairo.FONT_SLANT_NORMAL,
            cairo.FONT_WEIGHT_BOLD
        )
        ctx.set_font_size(residue_font_size)
        ctx.set_source_rgba(0, 0, 0, 1)

        ext = ctx.text_extents(res)
        ctx.move_to(
            cx - ext[2] / 2,
            cy + ext[3] / 2
        )
        ctx.show_text(res)

    return ctx, surface

###### Function from  matteoferla on Github https://gist.github.com/matteoferla/94eb8e4f8441ddfb458bfc45722469b8 ######

def set_to_neutral_pH(mol: Chem):
    """
    Not great, but does the job.
    
    * Protonates amines, but not aromatic bound amines.
    * Deprotonates carboxylic acid, phosphoric acid and sulfuric acid, without ruining esters.
    """
    protons_added = 0
    protons_removed = 0
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('[N;D1]')):
        atom = mol.GetAtomWithIdx(indices[0])
        if atom.GetNeighbors()[0].GetIsAromatic():
            continue # aniline
        atom.SetFormalCharge(1)
        protons_added += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[O;D1]')):
        atom = mol.GetAtomWithIdx(indices[2])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        protons_removed += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('P(=O)[Oh1]')):
        atom = mol.GetAtomWithIdx(indices[2])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        protons_removed += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('S(=O)(=O)[Oh1]')):
        atom = mol.GetAtomWithIdx(indices[3])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        protons_removed += 1
    return (protons_added, protons_removed)

# Define periodic table dictionary
periodic_table = {
    1: 'H',  2: 'He', 3: 'Li', 4: 'Be', 5: 'B',  6: 'C',  7: 'N',  8: 'O',  9: 'F',  10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U',  93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
    110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og',
}

def convert_and_write_pdb(input_pdb, output_pdb, bsid):
    # Load PDB from string
    pdb = [x for x in pybel.readfile("pdb", input_pdb)][0]

    for res in ob.OBResidueIter(pdb.OBMol):
        if not res.GetResidueProperty(0):
            if not res.GetName().isalnum():
                res.SetName("LIG")
            if res.GetChain().isspace():
                res.SetChain(bsid.split(":")[1])    
        for i,atom in enumerate(ob.OBResidueAtomIter(res)):
            if not res.GetAtomID(atom).isalnum():
                atom_name = periodic_table[atom.GetAtomicNum()]+str(i+1)
                res.SetAtomID(atom, atom_name)

    pdb.addh()
    pdb.write("pdb", output_pdb, overwrite=True)

    return pdb

keys = (
    "hydrophobic",
    "hbond",
    "waterbridge",
    "saltbridge",
    "pistacking",
    "pication",
    "halogen",
    "metal",
)

hbkeys = [
    "resnr",
    "restype",
    "reschain",
    "resnr_l",
    "restype_l",
    "reschain_l",
    "sidechain",
    "distance_ah", 
    "distance_ad", 
    "angle", 
    "type",
    "protisdon",
    "d_orig_idx",
    "a_orig_idx",
    "h"
]

def plip_2d_interactions(file, bsid, padding=35, canvas_height=700, canvas_width=1000,
                         save_files=True, save_pymol=True, out_name="PanViz_interactions.png",
                         output_dir=None, analysis=None):
    """Run PLIP analysis for one binding site and render a PanViz interaction diagram.

    Parameters
    ----------
    file : str
        Input PDB complex.
    bsid : str
        PLIP binding-site identifier, e.g. ``LIG:X:999``.
    output_dir : str or None
        Root directory for organized PanViz outputs. If None, preserve the original
        native output location beside the input PDB.
    analysis : dict or None
        Optional previously completed PanViz/PLIP analysis. When supplied, PLIP
        is not re-run; the existing normalized structure and interaction tables are
        reused for an additional figure render (for example, PNG and SVG parity).
    """

    if not save_files:
        out_name = None
        save_pymol = False

    input_dir = os.path.dirname(os.path.abspath(file))
    input_stem = os.path.splitext(os.path.basename(file))[0]

    # Preserve the original native output behavior when no explicit output directory is given.
    if output_dir is None:
        root_outdir = os.path.join(input_dir, input_stem + "_output")
    else:
        root_outdir = os.path.abspath(output_dir)

    binding_site_dir = os.path.join(root_outdir, "_".join(bsid.split(":")))
    figures_dir = os.path.join(binding_site_dir, "figures")
    interactions_dir = os.path.join(binding_site_dir, "interactions")
    structures_dir = os.path.join(binding_site_dir, "structures")
    pymol_dir = os.path.join(binding_site_dir, "pymol")

    for directory in [figures_dir, interactions_dir, structures_dir]:
        os.makedirs(directory, exist_ok=True)
    if save_pymol:
        os.makedirs(pymol_dir, exist_ok=True)

    # The corrected PDB used for PLIP analysis belongs with structures.
    # A completed analysis can be reused so one scientific analysis feeds all
    # downstream renders and the interactive editor.
    file_prot = os.path.join(structures_dir, f"{input_stem}_prot.pdb")
    if analysis is not None:
        file_prot = os.path.abspath(str(analysis.get("file_prot", file_prot)))
        input_pdb = analysis["input_pdb"]
        my_mol = analysis["my_mol"]
        my_interactions = analysis["my_interactions"]
        hydrophobic_df = analysis["hydrophobic_df"]
        hbond_df = analysis["hbond_df"]
        pi_stacking_df = analysis["pi_stacking_df"]
        pi_cation_df = analysis["pi_cation_df"]
        saltbridge_df = analysis["saltbridge_df"]
        waterbridge_df = analysis.get("waterbridge_df", pd.DataFrame())
        halogen_df = analysis.get("halogen_df", pd.DataFrame())
        metal_df = analysis.get("metal_df", pd.DataFrame())
    else:
        input_pdb = convert_and_write_pdb(file, file_prot, bsid)

        my_mol = PDBComplex()
        my_mol.load_pdb(file_prot)
        my_mol.analyze()
        my_interactions = my_mol.interaction_sets[bsid]

        bsr = BindingSiteReport(my_interactions)
        interactions = {
            k: [getattr(bsr, k + "_features")] + getattr(bsr, k + "_info")
            for k in keys
        }

        hydrophobic_df = pd.DataFrame(interactions["hydrophobic"][1:], columns=interactions["hydrophobic"][0])
        hbond_df = []
        for hb in my_interactions.all_hbonds_pdon + my_interactions.all_hbonds_ldon:
            hb_interactions = []
            for k in hbkeys:
                hb_interactions.append(getattr(hb, k))
            hbond_df.append(np.array(hb_interactions))
        if len(hbond_df) != 0:
            hbond_df = pd.DataFrame(np.stack(hbond_df), columns=hbkeys)
            hbond_df["h"] = [x.idx for x in hbond_df["h"]]
        else:
            hbond_df = pd.DataFrame()
        pi_stacking_df = pd.DataFrame(interactions["pistacking"][1:], columns=interactions["pistacking"][0])
        pi_cation_df = pd.DataFrame(interactions["pication"][1:], columns=interactions["pication"][0])
        saltbridge_df = pd.DataFrame(interactions["saltbridge"][1:], columns=interactions["saltbridge"][0])
        waterbridge_df = pd.DataFrame(interactions["waterbridge"][1:], columns=interactions["waterbridge"][0])
        halogen_df = pd.DataFrame(interactions["halogen"][1:], columns=interactions["halogen"][0])
        metal_df = pd.DataFrame(interactions["metal"][1:], columns=interactions["metal"][0])

        if save_pymol:
            _save_pymol(my_mol, bsid, pymol_dir)

    if save_files:
        if len(hydrophobic_df) > 0:
            hydrophobic_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_HPI.csv"), index=False)
        if len(hbond_df) > 0:
            hbond_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_HB.csv"), index=False)
        if len(pi_stacking_df) > 0:
            pi_stacking_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_PS.csv"), index=False)
        if len(pi_cation_df) > 0:
            pi_cation_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_PC.csv"), index=False)
        if len(saltbridge_df) > 0:
            saltbridge_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_SB.csv"), index=False)
        if len(waterbridge_df) > 0:
            waterbridge_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_WB.csv"), index=False)
        if len(halogen_df) > 0:
            halogen_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_XB.csv"), index=False)
        if len(metal_df) > 0:
            metal_df.to_csv(os.path.join(interactions_dir, f"{input_stem}_MC.csv"), index=False)

    with open(file_prot,"r") as f:
        pdb = f.readlines()

    pdb = [line for line in pdb if line.startswith(("ATOM","HETATM"))]

    with tempfile.TemporaryDirectory() as temp_dir:
        lig_path = os.path.join(temp_dir, "lig.pdb")

        lig = "".join([line for line in pdb if ((line[17:20] == bsid.split(":")[0])&(line[21] == bsid.split(":")[1])&(line[22:26].strip() == bsid.split(":")[2]))])

        ## Pybel seems to protonate differently when running from a jupyter notebook versus command line
        ## No idea why this is, but this is a workaround

        try:
            mol = Chem.MolFromPDBBlock(lig, removeHs=False)
            rdDetermineBonds.DetermineBonds(mol, charge=0)
        except:
            mol = pybel.readstring("pdb",lig)
            mol.write("pdb",lig_path, overwrite=True)
            mol = Chem.MolFromPDBFile(lig_path, removeHs=False)

    AllChem.EmbedMolecule(mol)
    set_to_neutral_pH(mol) #### Helper function to protonate/deprotonate groups. ####

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            bound_atom = atom.GetBonds()[0].GetEndAtom()
            if bound_atom.GetIdx == atom.GetIdx():
                bound_atom = atom.GetBonds()[0].GetBeginAtom()

            if bound_atom.GetSymbol() in ["O","N","S"]: #### Workaround so RDKit keeps explicit polar H's ####
                atom.SetAtomicNum(100) 

    mol=Chem.RemoveHs(mol)

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 100:
            atom.SetAtomicNum(1)

    # Convert aromatic systems to a valid Kekule form before rendering.
    # This preserves chemical valence while allowing the Cairo renderer to draw
    # conventional alternating double bonds (benzene -> three double bonds).
    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
    except Exception:
        # Keep the original aromatic representation if kekulization is not possible.
        pass

    Chem.rdDepictor.Compute2DCoords(mol)

    atom_info = []
    charge_info = []
    bonds = []
    for i,atom in enumerate(mol.GetAtoms()):
        coords = mol.GetConformer().GetAtomPosition(i)
        atom_info.append((coords.x, coords.y, atom.GetSymbol(), atom.GetPDBResidueInfo().GetName().strip()))

        if atom.GetFormalCharge() < 0:
            charge_info.append((coords.x+0.3, coords.y-0.2, "–", "Charge"))
        if atom.GetFormalCharge() > 0:
            charge_info.append((coords.x+0.3, coords.y-0.2, "+", "Charge"))
            
        startatoms = [bond.GetBeginAtomIdx() for bond in atom.GetBonds()]
        endatoms = [bond.GetEndAtomIdx() for bond in atom.GetBonds()]
        bond_type = [str(bond.GetBondType()).split(".")[-1] for bond in atom.GetBonds()]

        for a,b,c in zip(startatoms,endatoms,bond_type):
            if (a,b,c) not in bonds and (b,a,c) not in bonds:
                bonds.append((a,b,c))

    coord_dict = {}
    for entry in atom_info:
        coord_dict[entry[3]] = (entry[0],entry[1])

    interactions, centroids, used_res = _get_interactions(
        input_pdb,
        hydrophobic_df,
        hbond_df,
        pi_stacking_df,
        pi_cation_df,
        saltbridge_df,
        waterbridge_df,
        halogen_df,
        metal_df,
        coord_dict,
    )

    res_info = _get_res_info(used_res, coord_dict, interactions)

    atom_info = atom_info + centroids
    lines = []
    for i in range(len(res_info)):
        res = res_info[i][3]
        for j in range(len(interactions)):
            if interactions[j][1] == res:
                atom = interactions[j][0]
                atom_index = [x[3] for x in atom_info].index(atom)
                lines.append((
                    i + len(atom_info),
                    atom_index,
                    interactions[j][2],
                    interactions[j][3]
                ))

    atom_info = atom_info + res_info + charge_info

    connections = bonds + lines

    if out_name is not None:
        safe_out_name = os.path.basename(out_name)
        if not safe_out_name.lower().endswith((".png", ".svg")):
            raise ValueError("Unsupported file format. Please provide either PNG or SVG extension.")
        out_name = os.path.join(figures_dir, safe_out_name)

    ctx, surface = _draw_mol(atom_info, connections, 40, canvas_height, canvas_width, out_name)


    # -----------------------------------------------------------------------
    # Compact Discovery Studio-inspired interaction legend
    # -----------------------------------------------------------------------
    # The legend is drawn manually so that line samples, typography and
    # spacing remain consistent with the interaction graphics above.
    # Show ONLY interaction types that are actually present in this
    # binding site.  Keep the order consistent with the interaction
    # vocabulary used by the renderer.
    legend_definitions = [
        ("HPI", "Hydrophobic", (0.5, 0.5, 0.5), [10, 5]),
        ("HB",  "H-bond",      (0.0, 0.0, 1.0), [10, 5]),
        ("PS",  "π-Stacking",  (0.0, 0.6, 0.0), [10, 5]),
        ("PC",  "π-Cation",    (1.0, 0.7, 0.0), [10, 5]),
        ("SB",  "Salt bridge", (1.0, 0.0, 1.0), [10, 5]),
        ("WB",  "Water bridge", (0.08, 0.58, 0.72), [3, 3]),
        ("XB",  "Halogen bond", (0.48, 0.36, 0.78), [6, 4]),
        ("MC",  "Metal coordination", (0.64, 0.34, 0.00), [2, 3]),
    ]

    present_interactions = {
        connection[2]
        for connection in connections
        if len(connection) >= 3
    }

    legend_items = [
        (label, color, dash)
        for interaction_code, label, color, dash in legend_definitions
        if interaction_code in present_interactions
    ]

    # If no displayed interaction type is present, don't draw an empty box.
    if not legend_items:
        legend_items = []

    # Legend defaults are intentionally bold and legible at publication scale.
    legend_font_size = 15 if canvas_height <= 800 else 17
    legend_line_width = 3

    ctx.select_font_face(
        "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD
    )
    ctx.set_font_size(legend_font_size)

    if legend_items:
        # Measure each item so the box is always correctly sized.
        item_metrics = []
        line_sample_width = 32
        line_to_text_gap = 8
        item_gap = 24
        horizontal_padding = 12
        vertical_padding = 8

        total_width = horizontal_padding * 2

        for label, color, dash in legend_items:
            ext = ctx.text_extents(label)
            text_width = ext[2]
            item_width = line_sample_width + line_to_text_gap + text_width
            item_metrics.append((label, color, dash, item_width, text_width))
            total_width += item_width

        total_width += item_gap * (len(legend_items) - 1)

        legend_height = legend_font_size + vertical_padding * 2 + 2
        legend_x = (canvas_width - total_width) / 2
        legend_y = canvas_height - padding / 3

        # White background with a restrained black border.
        box_x = legend_x
        box_y = legend_y - legend_height + 4

        ctx.set_dash([], 0)
        ctx.set_line_width(1.5)
        ctx.set_source_rgb(1, 1, 1)
        ctx.rectangle(box_x, box_y, total_width, legend_height)
        ctx.fill_preserve()
        ctx.set_source_rgb(0, 0, 0)
        ctx.stroke()

        # Draw compact, uniform legend entries.
        cursor_x = legend_x + horizontal_padding
        baseline_y = box_y + vertical_padding + legend_font_size - 1

        for label, color, dash, item_width, text_width in item_metrics:
            # Interaction sample.
            ctx.set_source_rgba(color[0], color[1], color[2], 0.8)
            ctx.set_line_width(legend_line_width)
            ctx.set_dash(dash, 0)

            line_y = baseline_y - legend_font_size * 0.38
            ctx.move_to(cursor_x, line_y)
            ctx.line_to(cursor_x + line_sample_width, line_y)
            ctx.stroke()

            # Label.
            ctx.set_dash([], 0)
            ctx.set_source_rgb(0, 0, 0)
            ctx.move_to(
                cursor_x + line_sample_width + line_to_text_gap,
                baseline_y
            )
            ctx.show_text(label)

            cursor_x += item_width + item_gap


    # Save the image to a file
    if save_files:
        if out_name.lower().endswith(".png"):
            surface.write_to_png(out_name)
        elif out_name.lower().endswith(".svg"):
            surface.finish()
        else:
            raise ValueError("Unsupported file format. Please provide either PNG or SVG extension.")
    else:
        try:
            from IPython.display import display, Image
            import io
            image_stream = io.BytesIO()
            image_surface = surface
            if isinstance(image_surface, cairo.ImageSurface):
                image_surface.write_to_png(image_stream)
                display(Image(data=image_stream.getvalue(), format="png"))
        except Exception:
            pass

    # Keep the completed scientific analysis reusable by downstream renderers.
    return {
        "file": os.path.abspath(file),
        "file_prot": os.path.abspath(file_prot),
        "bsid": str(bsid),
        "input_pdb": input_pdb,
        "my_mol": my_mol,
        "my_interactions": my_interactions,
        "hydrophobic_df": hydrophobic_df,
        "hbond_df": hbond_df,
        "pi_stacking_df": pi_stacking_df,
        "pi_cation_df": pi_cation_df,
        "saltbridge_df": saltbridge_df,
        "waterbridge_df": waterbridge_df,
        "halogen_df": halogen_df,
        "metal_df": metal_df,
        "binding_site_dir": os.path.abspath(binding_site_dir),
        "figures_dir": os.path.abspath(figures_dir),
        "interactions_dir": os.path.abspath(interactions_dir),
        "structures_dir": os.path.abspath(structures_dir),
    }
