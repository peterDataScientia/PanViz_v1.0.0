"""PanViz Streamlit Community Cloud native-dependency diagnostic.

Run with:
    streamlit run cloud_diagnostic.py

This file intentionally imports PanViz dependencies one at a time and flushes a
checkpoint to stdout before and after each import. If the interpreter is killed
by a native segmentation fault, the Streamlit Cloud log will show the last
checkpoint reached.

It does not analyze structures and does not modify PanViz scientific data.
"""

from __future__ import annotations

import importlib
import platform
import sys
import traceback
from pathlib import Path

import streamlit as st


def log(message: str) -> None:
    """Write the same checkpoint to Streamlit UI and process stdout."""
    print(message, flush=True)
    st.write(message)


def module_path(module) -> str:
    return str(getattr(module, "__file__", "built-in / unavailable"))


def module_version(module, *attrs: str) -> str:
    for attr in attrs:
        value = getattr(module, attr, None)
        if value is not None:
            return str(value)
    return "unknown"


st.set_page_config(
    page_title="PanViz Cloud Diagnostic",
    page_icon="🧪",
    layout="centered",
)

st.title("PanViz Cloud Diagnostic")
st.caption(
    "Sequential native/scientific dependency test. "
    "If Streamlit Cloud segfaults, inspect the log for the last completed checkpoint."
)

log("CHECKPOINT 00 — diagnostic script started")
log(f"Python: {sys.version.replace(chr(10), ' ')}")
log(f"Executable: {sys.executable}")
log(f"Platform: {platform.platform()}")
log(f"Working directory: {Path.cwd()}")
log(f"Streamlit: {st.__version__}")


def run_checkpoint(number: int, label: str, importer):
    prefix = f"CHECKPOINT {number:02d} — {label}"
    log(prefix + " — START")
    try:
        result = importer()
    except Exception as exc:
        log(prefix + f" — PYTHON EXCEPTION: {type(exc).__name__}: {exc}")
        st.code(traceback.format_exc())
        raise
    log(prefix + " — OK")
    return result


np = run_checkpoint(1, "NumPy", lambda: importlib.import_module("numpy"))
log(f"  NumPy version={module_version(np, '__version__')} path={module_path(np)}")

pd = run_checkpoint(2, "Pandas", lambda: importlib.import_module("pandas"))
log(f"  Pandas version={module_version(pd, '__version__')} path={module_path(pd)}")

pa = run_checkpoint(3, "PyArrow", lambda: importlib.import_module("pyarrow"))
log(f"  PyArrow version={module_version(pa, '__version__')} path={module_path(pa)}")

cairo = run_checkpoint(4, "PyCairo", lambda: importlib.import_module("cairo"))
log(
    "  PyCairo version="
    + module_version(cairo, "version", "__version__")
    + f" path={module_path(cairo)}"
)

rdkit = run_checkpoint(5, "RDKit", lambda: importlib.import_module("rdkit"))
log(f"  RDKit version={module_version(rdkit, '__version__')} path={module_path(rdkit)}")


def import_rdkit_chemistry():
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    # Exercise compiled bindings, not merely the top-level package import.
    mol = Chem.MolFromSmiles("CCO")
    if mol is None or mol.GetNumAtoms() != 3:
        raise RuntimeError("RDKit smoke molecule construction failed.")
    return Chem, AllChem, rdDetermineBonds


run_checkpoint(6, "RDKit compiled chemistry bindings", import_rdkit_chemistry)


def import_openbabel():
    from openbabel import openbabel as ob

    version = ob.OBReleaseVersion()
    return ob, version


ob, ob_version = run_checkpoint(7, "Open Babel core", import_openbabel)
log(f"  Open Babel release={ob_version} path={module_path(ob)}")


def import_pybel():
    from openbabel import pybel

    # Exercise the shared native Open Babel layer.
    mol = pybel.readstring("smi", "CCO")
    if mol is None:
        raise RuntimeError("Pybel smoke molecule construction failed.")
    return pybel


pybel = run_checkpoint(8, "Open Babel pybel", import_pybel)
log(f"  pybel path={module_path(pybel)}")


plip = run_checkpoint(9, "PLIP package", lambda: importlib.import_module("plip"))
log(f"  PLIP version={module_version(plip, '__version__')} path={module_path(plip)}")


def import_plip_core():
    from plip.structure.preparation import PDBComplex
    from plip.exchange.report import BindingSiteReport

    return PDBComplex, BindingSiteReport


run_checkpoint(10, "PLIP core preparation/report", import_plip_core)


def import_plip_visualization():
    from plip.visualization.visualize import PyMOLVisualizer
    from plip.basic.remote import VisualizerData
    from plip.basic.supplemental import start_pymol

    return PyMOLVisualizer, VisualizerData, start_pymol


try:
    run_checkpoint(11, "PLIP visualization / optional PyMOL bridge", import_plip_visualization)
except Exception:
    # PanViz treats this path as optional. A normal Python import exception here
    # should be visible, but should not prevent the remaining startup tests.
    log("CHECKPOINT 11 — optional visualization bridge unavailable; continuing")


run_checkpoint(
    12,
    "PanViz independent canonical engine",
    lambda: importlib.import_module("panviz_engine"),
)
run_checkpoint(
    13,
    "CairoSVG static export dependency",
    lambda: importlib.import_module("cairosvg"),
)


def read_editor():
    path = Path(__file__).with_name("editor.html")
    text = path.read_text(encoding="utf-8")
    if "__PANVIZ_SCENE__" not in text:
        raise RuntimeError("editor.html scene placeholder was not found.")
    return len(text)


editor_size = run_checkpoint(14, "editor.html read", read_editor)
log(f"  editor.html characters={editor_size}")

log("CHECKPOINT 15 — ALL PANVIZ STARTUP CHECKS PASSED")
st.success("All startup/import checkpoints passed in this environment.")
st.info(
    "If the normal PanViz app still segfaults while this diagnostic passes, "
    "the failure is after startup and should be isolated by testing specific runtime actions."
)


# ---------------------------------------------------------------------------
# Phase 2: reproduce app.py's pre-upload Streamlit execution path.
# ---------------------------------------------------------------------------
import tempfile

log("CHECKPOINT 16 — PHASE 2 pre-upload Streamlit path — START")

log("CHECKPOINT 17 — create temporary PanViz workspace — START")
_phase2_work_root = Path(tempfile.mkdtemp(prefix="panviz_diag_"))
log(f"CHECKPOINT 17 — create temporary PanViz workspace — OK: {_phase2_work_root}")

log("CHECKPOINT 18 — render PanViz-style HTML/CSS — START")
st.markdown(
    """
    <style>
    .pv-diag-shell{
        border:1px solid #dbe3ee;
        border-radius:14px;
        padding:12px 14px;
        background:#ffffff;
        margin:10px 0;
    }
    </style>
    <div class="pv-diag-shell">
      <strong>PanViz pre-upload UI smoke test</strong><br>
      This block intentionally exercises the same unsafe-HTML rendering path used by app.py.
    </div>
    """,
    unsafe_allow_html=True,
)
log("CHECKPOINT 18 — render PanViz-style HTML/CSS — OK")

log("CHECKPOINT 19 — Streamlit radio widget — START")
_diag_mode = st.radio(
    "Diagnostic input mode",
    ["PDB complex", "Docking PDBQT"],
    horizontal=True,
    key="panviz_diag_input_mode",
)
log(f"CHECKPOINT 19 — Streamlit radio widget — OK: {_diag_mode}")

log("CHECKPOINT 20 — Streamlit file uploader — START")
_diag_upload = st.file_uploader(
    "Diagnostic PDB uploader — no file is processed",
    type=["pdb"],
    key="panviz_diag_pdb_upload",
)
log("CHECKPOINT 20 — Streamlit file uploader — OK")

log("CHECKPOINT 21 — Streamlit info/markdown path — START")
if _diag_upload is None:
    st.info("No diagnostic file uploaded. This intentionally matches PanViz's initial no-upload state.")
    st.markdown("</div>", unsafe_allow_html=True)
log("CHECKPOINT 21 — Streamlit info/markdown path — OK")

log("CHECKPOINT 22 — st.stop() path — ABOUT TO EXECUTE")
st.caption(
    "Checkpoint 22 is the expected final line for this diagnostic. "
    "If the app remains online after this, Streamlit's normal st.stop() path is healthy."
)
st.stop()
