import os
import shutil
import tempfile
import subprocess
import re
from pathlib import Path
import json
import hashlib
import zipfile
from datetime import datetime, timezone
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from plip.structure.preparation import PDBComplex
from utils import plip_2d_interactions
from interactive_engine import build_editor_scene

st.set_page_config(page_title="PanViz v5.8.6", page_icon="🧬", layout="wide", initial_sidebar_state="expanded")
st.markdown("""
<style>
:root{--pv-navy:#143761;--pv-blue:#1f5aa6;--pv-ink:#16243a;--pv-muted:#65748b;--pv-line:#dce4ef}
.block-container{padding-top:1.2rem;padding-bottom:2.2rem;max-width:1520px}
.panviz-shell{border:1px solid #dbe3ee;border-radius:18px;padding:18px 20px 16px;background:linear-gradient(135deg,#f7faff 0%,#ffffff 52%,#f4f7fb 100%);box-shadow:0 8px 28px rgba(24,54,90,.07);margin-bottom:16px}
.panviz-brand{display:flex;align-items:center;gap:12px}
.panviz-mark{width:42px;height:42px;border-radius:12px;display:grid;place-items:center;color:white;font-size:22px;font-weight:800;background:linear-gradient(135deg,#1d5da9,#143761);box-shadow:0 5px 14px rgba(20,55,97,.24)}
.panviz-title{font-size:2.0rem;font-weight:800;letter-spacing:-.6px;color:var(--pv-ink);line-height:1.05}
.panviz-subtitle{color:var(--pv-muted);font-size:.94rem;margin-top:3px}
.panviz-badges{display:flex;gap:7px;flex-wrap:wrap;margin-top:12px}
.panviz-badge{font-size:.72rem;font-weight:700;color:#34506f;background:#eef4fb;border:1px solid #dce7f3;border-radius:999px;padding:5px 9px}
.panviz-section{border:1px solid var(--pv-line);border-radius:14px;padding:12px 14px;background:#fff;box-shadow:0 4px 14px rgba(31,55,88,.05);margin:10px 0}
.panviz-section h4{margin:0 0 8px;color:var(--pv-navy);font-size:13px}
div[data-testid="stFileUploader"]{border:1px dashed #b8c8de;border-radius:14px;background:#fbfdff;padding:4px}
.pdbqt-step{font-size:11px;font-weight:800;letter-spacing:.08em;color:#6d7d92;margin-bottom:4px}
.pdbqt-card-title{font-size:14px;font-weight:800;color:#143761;letter-spacing:.02em}
.pdbqt-card-desc{font-size:12px;color:#536579;margin-top:4px;line-height:1.45}
.pdbqt-example{font-size:11px;color:#68788d;margin-top:8px;line-height:1.45}
.pdbqt-example code{font-size:10.5px;background:#eef3f8;padding:2px 4px;border-radius:4px}
.stButton>button{border-radius:9px;font-weight:700}
.stDownloadButton>button{border-radius:9px}
.panviz-foot{color:#7a8798;font-size:.73rem;margin-top:10px}
</style>
""", unsafe_allow_html=True)
st.markdown("""<div class="panviz-shell"><div class="panviz-brand"><div class="panviz-mark">🧬</div><div><div class="panviz-title">PanViz</div><div class="panviz-subtitle">PLIP-based protein–ligand interaction visualization &amp; publication figure editor</div></div></div><div class="panviz-badges"><span class="panviz-badge">PLIP interaction analysis</span><span class="panviz-badge">Editable presentation layer</span><span class="panviz-badge">Molecular topology locked</span><span class="panviz-badge">v5.8.6</span></div></div>""", unsafe_allow_html=True)

EDITOR_HTML = (Path(__file__).with_name("editor.html")).read_text(encoding="utf-8")

def render_editor(scene):
    html = EDITOR_HTML.replace("__PANVIZ_SCENE__", json.dumps(scene, ensure_ascii=False))
    # Keep the iframe tall enough for the selected canvas and editor chrome.
    canvas_h = int(scene.get("height", 850) or 850)
    iframe_h = max(1200, min(5000, canvas_h + 520))
    components.html(html, height=iframe_h, scrolling=False)

def _eligible_binding_sites(pdb_path):
    mol = PDBComplex(); mol.load_pdb(str(pdb_path)); excluded={"ARN","ASH","GLH","LYN","HIE","HIP"}
    return [x for x in str(mol).split("\n")[1:] if x.strip() and x.split(":")[0] not in excluded]

def _find_obabel():
    for name in ("obabel", "obabel.exe"):
        path = shutil.which(name)
        if path:
            return path
    return None

def _extract_pdbqt_models(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    blocks = re.findall(r"MODEL\s+\d+.*?ENDMDL", text, flags=re.S|re.I)
    if blocks:
        return blocks
    return [text]

def _safe_remove(path, attempts=20, delay=0.20):
    """Best-effort Windows-safe removal of a temporary file.

    Antivirus/indexer/Open Babel child processes can briefly retain a handle
    even after subprocess.run() returns. Cleanup must never convert a
    successful PDBQT conversion into a fatal WinError 32.
    """
    if path is None:
        return
    target = Path(path)
    for _ in range(max(1, attempts)):
        try:
            target.unlink(missing_ok=True)
            return
        except PermissionError:
            import time
            time.sleep(delay)
        except OSError:
            import time
            time.sleep(delay)
    # Leave the file in place rather than failing the scientific workflow.

def _convert_with_obabel(input_path, output_path, input_format, selected_block=None):
    obabel = _find_obabel()
    if not obabel:
        raise RuntimeError("Open Babel executable 'obabel' was not found. Install openbabel-wheel in the PanViz environment.")
    source = Path(input_path)
    cleanup = None
    if selected_block is not None:
        # Create, close, and then write the selected pose so no Python handle
        # remains open while Open Babel reads the file on Windows.
        temp = tempfile.NamedTemporaryFile(
            mode="w", prefix="panviz_pose_", suffix="." + input_format,
            encoding="utf-8", newline="", delete=False
        )
        try:
            temp.write(selected_block + "\n")
            temp.flush()
        finally:
            temp.close()
        cleanup = Path(temp.name)
        source = cleanup
    try:
        cmd = [obabel, "-i", input_format, str(source), "-o", "pdb", "-O", str(output_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0 or not Path(output_path).exists():
            msg = (result.stderr or result.stdout or "Open Babel conversion failed.").strip()
            raise RuntimeError(msg)
    finally:
        # Cleanup failure must not mask a successful Open Babel conversion.
        _safe_remove(cleanup)

def _normalize_docked_ligand_pdb(pdb_path, out_path, chain="Z", residue_number=900, residue_name="LIG"):
    """Normalize a converted docking pose into a PLIP-friendly ligand residue.

    Atom serials are rewritten consistently and any CONECT records emitted by Open Babel
    are remapped to the new serials and placed before END, preserving ligand connectivity.
    """
    lines = Path(pdb_path).read_text(encoding="utf-8", errors="replace").splitlines()
    atom_lines=[]; conect=[]; serial_map={}
    serial=1
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            s=line.ljust(80)
            record="HETATM"
            atom_name=s[12:16]
            element=s[76:78].strip()
            if not element:
                raw = re.sub(r"[^A-Za-z]", "", atom_name).strip()
                element = raw[:2].title() if raw[:2].lower() in {"cl","br"} else raw[:1].upper()
            x=s[30:38]; y=s[38:46]; z=s[46:54]
            charge=s[78:80] if len(s)>=80 else "  "
            old_serial = s[6:11].strip()
            if old_serial.isdigit():
                serial_map[int(old_serial)] = serial
            new=(f"{record:<6}{serial:5d} {atom_name:>4} {residue_name:>3} {chain:1}{residue_number:4d}    "
                 f"{x:>8}{y:>8}{z:>8}  1.00  0.00          {element:>2}{charge:>2}")
            atom_lines.append(new)
            serial+=1
        elif line.startswith("CONECT"):
            fields=line.split()[1:]
            nums=[]
            for token in fields:
                try:
                    nums.append(serial_map.get(int(token)))
                except ValueError:
                    nums.append(None)
            nums=[n for n in nums if n is not None]
            if len(nums)>=2:
                conect.append("CONECT" + "".join(f"{n:5d}" for n in nums))
    if not atom_lines:
        raise RuntimeError("The selected PDBQT pose did not contain any atom records after conversion.")
    Path(out_path).write_text("\n".join(atom_lines + conect + ["TER", "END"]) + "\n", encoding="utf-8")

def _validate_pdb_has_atoms(path, role):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if not any(line.startswith(("ATOM", "HETATM")) for line in lines):
        raise RuntimeError(f"The uploaded {role} does not contain any ATOM/HETATM records.")

def _build_pdbqt_complex(receptor_path, ligand_path, work_root, pose_index=0):
    receptor_pdb = work_root / "receptor.pdb"
    ligand_pdb_raw = work_root / "ligand_raw.pdb"
    ligand_pdb = work_root / "ligand_normalized.pdb"

    receptor_suffix = Path(receptor_path).suffix.lower()
    if receptor_suffix == ".pdbqt":
        _convert_with_obabel(receptor_path, receptor_pdb, "pdbqt")
    else:
        receptor_pdb.write_bytes(Path(receptor_path).read_bytes())
    _validate_pdb_has_atoms(receptor_pdb, "protein / receptor")

    blocks = _extract_pdbqt_models(ligand_path)
    if pose_index < 0 or pose_index >= len(blocks):
        raise ValueError(f"Docking pose {pose_index + 1} is outside the available range (1–{len(blocks)}).")
    _convert_with_obabel(ligand_path, ligand_pdb_raw, "pdbqt", selected_block=blocks[pose_index])
    _validate_pdb_has_atoms(ligand_pdb_raw, "ligand / docking pose")
    _normalize_docked_ligand_pdb(ligand_pdb_raw, ligand_pdb)

    combined = work_root / "panviz_pdbqt_complex.pdb"
    receptor_lines=[x for x in receptor_pdb.read_text(encoding="utf-8", errors="replace").splitlines() if x[:6].strip() in {"ATOM", "HETATM", "TER"}]
    ligand_lines=[x for x in ligand_pdb.read_text(encoding="utf-8", errors="replace").splitlines() if x[:6].strip() in {"ATOM", "HETATM", "TER"}]
    combined.write_text("\n".join(receptor_lines + ["TER"] + ligand_lines + ["END"]) + "\n", encoding="utf-8")
    return combined, len(blocks)

def _first_col(df, names):
    lookup={str(c).lower():c for c in df.columns}
    for name in names:
        if name.lower() in lookup:return lookup[name.lower()]
    return None

def _read_interactions(interaction_dir):
    mapping={"HPI":"Hydrophobic interaction","HB":"Hydrogen bond","PS":"π-Stacking","PC":"π-Cation","SB":"Salt bridge"}
    rows=[]
    for path in sorted(Path(interaction_dir).glob("*.csv")):
        code=path.stem.rsplit("_",1)[-1]
        if code not in mapping: continue
        df=pd.read_csv(path)
        if df.empty: continue
        rt=_first_col(df,["RESTYPE","restype"]);rn=_first_col(df,["RESNR","resnr"]);rc=_first_col(df,["RESCHAIN","reschain"]);dist=_first_col(df,["DIST","distance","distance_ad","distance_ah","dist"])
        for _,r in df.iterrows():
            residue=""
            if rt is not None and pd.notna(r[rt]):residue+=str(r[rt]).strip()
            if rn is not None and pd.notna(r[rn]):residue+=str(r[rn]).strip()
            if rc is not None and pd.notna(r[rc]):residue+=str(r[rc]).strip()
            d=None
            if dist is not None:
                try:d=float(r[dist])
                except (TypeError,ValueError):pass
            rows.append({"Residue":residue or "—","Interaction":mapping[code],"Distance (Å)":d})
    return pd.DataFrame(rows)

def _sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _pose_scores(path):
    """Return Vina REMARK scores in MODEL order when available."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    blocks = _extract_pdbqt_models(path)
    scores = []
    for block in blocks:
        score = None
        for line in block.splitlines():
            m = re.search(r"REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)", line, flags=re.I)
            if m:
                try:
                    score = float(m.group(1))
                except ValueError:
                    pass
                break
        scores.append(score)
    return scores


def _zip_tree(source_root, output_zip):
    source_root = Path(source_root)
    output_zip = Path(output_zip)
    with zipfile.ZipFile(output_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in source_root.rglob("*"):
            if path.is_file() and path != output_zip:
                zf.write(path, path.relative_to(source_root))


def _write_project_manifest(result, manifest_path):
    manifest = {
        "panviz_version": "5.8.6",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_mode": result["input_mode"],
        "source_files": result["source_files"],
        "selected_site": result["selected_site"],
        "figure_width": result["figure_width"],
        "figure_height": result["figure_height"],
        "pose_index": result.get("pose_index"),
        "pose_score_kcal_mol": result.get("pose_score"),
        "binding_site_count": result["binding_site_count"],
        "interaction_count": int(len(result["interaction_df"])),
        "residue_count": int(result["interaction_df"]["Residue"].nunique()) if not result["interaction_df"].empty else 0,
        "interaction_types": sorted(result["interaction_df"]["Interaction"].dropna().unique().tolist()) if not result["interaction_df"].empty else [],
        "scientific_data_policy": "PLIP interaction measurements are not altered by editor styling.",
    }
    Path(manifest_path).write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


st.markdown('<div class="panviz-section"><h4>1 · Input structure</h4>', unsafe_allow_html=True)
input_mode=st.radio("Input mode",["PDB complex","Docking PDBQT"],horizontal=True)
work_root=Path(tempfile.mkdtemp(prefix="panviz_"))
source_files=[]
source_payloads=[]
pose_index=None
pose_score=None

if input_mode=="PDB complex":
    uploaded=st.file_uploader("Upload protein–ligand PDB complex",type=["pdb"],help="Upload a complete PDB complex containing the protein and ligand.")
    if not uploaded:
        st.info("Upload a PDB complex to begin.")
        st.markdown("</div>", unsafe_allow_html=True)
        st.stop()
    pdb_path=work_root/uploaded.name
    file_bytes=uploaded.getvalue()
    pdb_path.write_bytes(file_bytes)
    source_stem=Path(uploaded.name).stem
    source_files=[uploaded.name]
    source_payloads=[(uploaded.name,file_bytes)]
else:
    left_card, right_card = st.columns(2, gap="large")
    with left_card:
        with st.container(border=True):
            st.markdown('<div class="pdbqt-step">01 · PROTEIN INPUT</div>', unsafe_allow_html=True)
            st.markdown('<div class="pdbqt-card-title">PROTEIN / RECEPTOR</div>', unsafe_allow_html=True)
            st.markdown('<div class="pdbqt-card-desc">Upload the macromolecular receptor structure.</div>', unsafe_allow_html=True)
            st.markdown('<div class="pdbqt-example">Accepted: <b>.pdb</b> or <b>.pdbqt</b> · Example: <code>1LF2_receptor.pdbqt</code></div>', unsafe_allow_html=True)
            receptor_upload=st.file_uploader("Upload PROTEIN / RECEPTOR",type=["pdb","pdbqt"],help="This is the protein/macromolecular receptor. Upload a receptor PDB or PDBQT file.",key="pdbqt_receptor_upload")
    with right_card:
        with st.container(border=True):
            st.markdown('<div class="pdbqt-step">02 · LIGAND INPUT</div>', unsafe_allow_html=True)
            st.markdown('<div class="pdbqt-card-title">LIGAND / DOCKING POSES</div>', unsafe_allow_html=True)
            st.markdown('<div class="pdbqt-card-desc">Upload the docked ligand or Vina output containing one or more poses.</div>', unsafe_allow_html=True)
            st.markdown('<div class="pdbqt-example">Accepted: <b>.pdbqt</b> · Example: <code>compound_01_out.pdbqt</code></div>', unsafe_allow_html=True)
            ligand_upload=st.file_uploader("Upload LIGAND / DOCKING POSES",type=["pdbqt"],help="This is the docked ligand/Vina output. Upload a ligand PDBQT containing one or more MODEL poses.",key="pdbqt_ligand_upload")
    if not receptor_upload or not ligand_upload:
        st.info("Upload both files: first the **PROTEIN / RECEPTOR**, then the **LIGAND / DOCKING POSES**.")
        st.markdown("</div>", unsafe_allow_html=True)
        st.stop()
    receptor_path=work_root/receptor_upload.name; receptor_bytes=receptor_upload.getvalue(); receptor_path.write_bytes(receptor_bytes)
    ligand_path=work_root/ligand_upload.name; ligand_bytes=ligand_upload.getvalue(); ligand_path.write_bytes(ligand_bytes)
    source_files=[receptor_upload.name, ligand_upload.name]
    source_payloads=[(receptor_upload.name,receptor_bytes),(ligand_upload.name,ligand_bytes)]
    pose_blocks=_extract_pdbqt_models(ligand_path)
    scores=_pose_scores(ligand_path)
    pose_options=list(range(1,len(pose_blocks)+1))
    def _pose_label(x):
        score=scores[x-1] if x-1 < len(scores) else None
        return f"Pose {x}" + (f"  ·  Vina {score:.2f} kcal/mol" if score is not None else "")
    pose_index=st.selectbox("Docking pose",pose_options,format_func=_pose_label)-1
    pose_score=scores[pose_index] if pose_index < len(scores) else None
    try:
        pdb_path,nposes=_build_pdbqt_complex(receptor_path,ligand_path,work_root,pose_index=pose_index)
    except Exception as exc:
        st.error(f"PDBQT preparation failed: {exc}")
        st.markdown("</div>", unsafe_allow_html=True)
        st.stop()
    source_stem=Path(ligand_upload.name).stem
    st.success(f"Prepared docking pose {pose_index+1} of {nposes} as a PLIP-ready PDB complex." + (f" Vina score: {pose_score:.2f} kcal/mol." if pose_score is not None else ""))

st.markdown("</div>", unsafe_allow_html=True)

try:binding_sites=_eligible_binding_sites(pdb_path)
except Exception as exc:st.error(f"The structure could not be read by PLIP: {exc}");st.stop()
if not binding_sites:st.error("No eligible small-molecule binding site was detected in the prepared structure.");st.stop()
st.success(f"Detected {len(binding_sites)} eligible binding site(s).")
st.markdown('<div class="panviz-section"><h4>2 · Analysis setup</h4>', unsafe_allow_html=True)
left,right=st.columns([1,1])
with left:selected_site=st.selectbox("Ligand / binding site",binding_sites)
with right:out_width=st.number_input("Figure width",min_value=700,max_value=3000,value=1200,step=100)
out_height=st.number_input("Figure height",min_value=500,max_value=3000,value=850,step=50)
analyze=st.button("Generate PanViz interaction diagram",type="primary",use_container_width=True)
st.caption("One PLIP analysis is reused for the original figures and the interactive editor; alternate exports do not trigger a second scientific analysis.")
st.markdown("</div>", unsafe_allow_html=True)

result_key_payload = {
    "mode": input_mode,
    "source_hashes": [_sha256_bytes(x[1]) for x in source_payloads],
    "pose": pose_index,
    "site": selected_site,
    "width": int(out_width),
    "height": int(out_height),
}
result_key=_sha256_bytes(json.dumps(result_key_payload,sort_keys=True).encode())

if analyze:
    results_root=work_root/"PanViz_results";results_root.mkdir(parents=True,exist_ok=True)
    with st.spinner("Running PLIP once and generating PanViz outputs…"):
        try:
            analysis_obj=plip_2d_interactions(str(pdb_path),selected_site,save_files=True,save_pymol=False,canvas_height=int(out_height),canvas_width=int(out_width),out_name="PanViz_interactions.png",output_dir=str(results_root))
            plip_2d_interactions(str(pdb_path),selected_site,save_files=True,save_pymol=False,canvas_height=int(out_height),canvas_width=int(out_width),out_name="PanViz_interactions.svg",output_dir=str(results_root),analysis=analysis_obj)
            site_dir=results_root/selected_site.replace(":","_")
            png_path=site_dir/"figures"/"PanViz_interactions.png";svg_path=site_dir/"figures"/"PanViz_interactions.svg";interaction_dir=site_dir/"interactions"
            interaction_df=_read_interactions(interaction_dir)
            scene,scene_root=build_editor_scene(str(pdb_path),selected_site,width=int(out_width),height=int(out_height),base_svg=svg_path.read_text(encoding="utf-8"),analysis=analysis_obj)
            (site_dir/"PanViz_initial_layout.json").write_text(json.dumps(scene,indent=2,ensure_ascii=False),encoding="utf-8")
            result={
                "key": result_key,
                "input_mode": input_mode,
                "source_files": [
                    {"name":name,"sha256":_sha256_bytes(data)} for name,data in source_payloads
                ],
                "source_payloads": source_payloads,
                "selected_site": selected_site,
                "figure_width": int(out_width),
                "figure_height": int(out_height),
                "pose_index": pose_index,
                "pose_score": pose_score,
                "binding_site_count": len(binding_sites),
                "interaction_df": interaction_df,
                "png_path": str(png_path),
                "svg_path": str(svg_path),
                "site_dir": str(site_dir),
                "prepared_pdb": str(analysis_obj["file_prot"]),
                "results_root": str(results_root),
                "scene": scene,
                "analysis_obj": analysis_obj,
            }
            inputs_dir=results_root/"inputs";inputs_dir.mkdir(parents=True,exist_ok=True)
            for name,data in source_payloads:
                target=inputs_dir/Path(name).name
                target.write_bytes(data)
            project_readme = results_root/"PROJECT_README.md"
            project_readme.write_text(
                "# PanViz 5.8.6 project bundle\n\n"
                "This package contains the original uploaded input file(s), the PLIP-prepared complex, "
                "original PanViz PNG/SVG figures, PLIP interaction CSV tables, the initial editable layout, "
                "and a machine-readable manifest. Presentation styling in PanViz does not modify the underlying "
                "PLIP scientific interaction records. Use the editor's **Save layout** and **Load layout** controls "
                "to carry edited presentation state between sessions.\n",
                encoding="utf-8"
            )
            _write_project_manifest(result, site_dir/"PanViz_manifest.json")
            zip_path=results_root/f"{source_stem}_PanViz_project.zip"
            _zip_tree(results_root,zip_path)
            result["project_zip"]=str(zip_path)
            st.session_state.panviz_result=result
        except Exception as exc:
            st.error(f"PanViz analysis failed: {exc}")
            st.stop()

result=st.session_state.get("panviz_result")
if not result or result.get("key")!=result_key:
    st.info("Configure the analysis above, then click **Generate PanViz interaction diagram**. Generated results remain available until you change the input, pose, binding site, or canvas size.")
    st.stop()

png_path=Path(result["png_path"]); svg_path=Path(result["svg_path"]); site_dir=Path(result["site_dir"]); interaction_df=result["interaction_df"]

st.markdown('<div class="panviz-section"><h4>3 · Interactive figure editor</h4>', unsafe_allow_html=True)
render_editor(result["scene"])
st.caption("v5.8.6: the imported molecular structure and PLIP scientific records are immutable; only the separate presentation annotation layer can be edited, saved, reloaded, and exported.")
st.markdown("</div>", unsafe_allow_html=True)

st.markdown('<div class="panviz-section"><h4>4 · Source figures and complete project download</h4>', unsafe_allow_html=True)
st.image(str(png_path),use_container_width=True)
b1,b2,b3=st.columns(3)
with b1:st.download_button("Download original PNG",data=png_path.read_bytes(),file_name=f"{source_stem}_PanViz.png",mime="image/png",use_container_width=True)
with b2:st.download_button("Download original SVG",data=svg_path.read_bytes(),file_name=f"{source_stem}_PanViz.svg",mime="image/svg+xml",use_container_width=True)
with b3:
    zpath=Path(result["project_zip"])
    st.download_button("Download complete project ZIP",data=zpath.read_bytes(),file_name=zpath.name,mime="application/zip",use_container_width=True)
st.markdown("</div>", unsafe_allow_html=True)

st.markdown('<div class="panviz-section"><h4>5 · Scientific records and prepared structure</h4>', unsafe_allow_html=True)
rec1,rec2,rec3,rec4=st.columns(4)
rec1.metric("Interactions",len(interaction_df));rec2.metric("Residues",interaction_df["Residue"].nunique() if not interaction_df.empty else 0);rec3.metric("Interaction types",interaction_df["Interaction"].nunique() if not interaction_df.empty else 0);rec4.metric("Binding sites",result["binding_site_count"])
if not interaction_df.empty:
    st.dataframe(interaction_df,use_container_width=True,hide_index=True)
    st.download_button("Download interaction table (CSV)",data=interaction_df.to_csv(index=False).encode("utf-8"),file_name=f"{source_stem}_{selected_site.replace(':','_')}_interactions.csv",mime="text/csv",use_container_width=True)
else:
    st.info("No PLIP interaction records were returned for this binding site.")
prepared=Path(result["prepared_pdb"])
if prepared.exists():
    st.download_button("Download PLIP-prepared PDB complex",data=prepared.read_bytes(),file_name=f"{source_stem}_{selected_site.replace(':','_')}_prepared.pdb",mime="chemical/x-pdb",use_container_width=True)
st.download_button("Download initial editor layout JSON",data=json.dumps(result["scene"],indent=2,ensure_ascii=False).encode("utf-8"),file_name=f"{source_stem}_{selected_site.replace(':','_')}_initial_layout.json",mime="application/json",use_container_width=True)
st.markdown('</div><div class="panviz-foot">PanViz v5.8.6 · one reusable PLIP analysis → original figures + locked molecular scene + editable presentation annotations + complete project package.</div>', unsafe_allow_html=True)

