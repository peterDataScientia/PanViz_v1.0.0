from __future__ import annotations

import argparse
import json
from pathlib import Path

# Preserve the Linux native-library import order used by the web app.
from rdkit import Chem  # noqa: F401

from plip.structure.preparation import PDBComplex

from panviz_engine import (
    PANVIZ_ENGINE_VERSION,
    build_editor_scene,
    run_panviz_analysis,
    write_static_exports,
)


EXCLUDED_RESIDUES = {"ARN", "ASH", "GLH", "LYN", "HIE", "HIP"}


def discover_binding_sites(pdb_file: str | Path) -> list[str]:
    complex_obj = PDBComplex()
    complex_obj.load_pdb(str(pdb_file))
    return [
        site
        for site in str(complex_obj).splitlines()[1:]
        if site.strip() and site.split(":", 1)[0] not in EXCLUDED_RESIDUES
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="PanViz",
        description=(
            "PanViz 6: PLIP interaction analysis with an independent canonical "
            "2D scene model and publication-oriented SVG/PNG exports."
        ),
    )
    parser.add_argument(
        "-f",
        "--file",
        required=True,
        help="Protein-ligand PDB complex.",
    )
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument(
        "--site",
        help="Analyze one PLIP binding-site identifier, for example LIG:A:401.",
    )
    scope.add_argument(
        "--all",
        action="store_true",
        help="Analyze every eligible small-molecule binding site.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1200,
        help="Scene width in pixels (default: 1200).",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=850,
        help="Scene height in pixels (default: 850).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="PanViz_results",
        help="Output root directory (default: PanViz_results).",
    )
    return parser


def analyze_site(
    pdb_file: Path,
    site: str,
    output_root: Path,
    width: int,
    height: int,
) -> None:
    analysis = run_panviz_analysis(
        pdb_file,
        site,
        output_root=output_root,
    )
    scene, _ = build_editor_scene(
        pdb_file,
        site,
        width=width,
        height=height,
        analysis=analysis,
    )

    site_dir = Path(analysis["site_dir"])
    figures = site_dir / "figures"
    svg_path = figures / "PanViz_interactions.svg"
    png_path = figures / "PanViz_interactions.png"
    write_static_exports(scene, svg_path, png_path, png_scale=2)

    layout_path = site_dir / "PanViz_initial_layout.json"
    layout_path.write_text(
        json.dumps(scene, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(
        f"[PanViz {PANVIZ_ENGINE_VERSION}] {site}: "
        f"{len(analysis['records'])} interactions -> {site_dir}"
    )


def main() -> None:
    args = build_parser().parse_args()
    pdb_file = Path(args.file).expanduser().resolve()
    if not pdb_file.exists():
        raise SystemExit(f"Input file does not exist: {pdb_file}")

    sites = discover_binding_sites(pdb_file)
    if not sites:
        raise SystemExit("No eligible PLIP small-molecule binding sites were detected.")

    if args.site:
        if args.site not in sites:
            available = "\n  ".join(sites)
            raise SystemExit(
                f"Binding site {args.site!r} was not detected.\nAvailable sites:\n  {available}"
            )
        selected = [args.site]
    elif args.all:
        selected = sites
    else:
        available = "\n  ".join(sites)
        raise SystemExit(
            "Choose --site SITE or --all.\nDetected binding sites:\n  " + available
        )

    output_root = Path(args.output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    for site in selected:
        analyze_site(
            pdb_file,
            site,
            output_root,
            width=max(700, int(args.width)),
            height=max(500, int(args.height)),
        )


if __name__ == "__main__":
    main()
