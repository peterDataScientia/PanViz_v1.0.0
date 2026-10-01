from utils import plip_2d_interactions
from plip.structure.preparation import PDBComplex
from panviz_version import PANVIZ_VERSION
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        prog="PanViz",
        description=(
            "PanViz: PLIP-based protein-ligand interaction visualization. "
            "Generates publication-quality 2D interaction diagrams and, "
            "optionally, 3D PyMOL sessions from X-ray or simulated binding poses."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=r"""
WORKFLOW
--------
1. Load the supplied PDB complex.
2. Detect PLIP binding-site objects.
3. In interactive mode, prompt for each eligible binding-site object and enter
   y to analyze it or n to skip it; use -y/--analyse_all to analyze every
   eligible small-molecule binding site automatically.
4. Re-analyze each selected binding site with PLIP.
5. Export interaction tables and generate the 2D interaction diagram.
6. Optionally create a PyMOL session.

OUTPUT STRUCTURE
----------------
If --output_dir is supplied, PanViz writes each analyzed binding site as:

  OUTPUT_DIR/
  └── <binding_site_id>/
      ├── figures/
      │   └── <figure>.png or .svg
      ├── interactions/
      │   ├── *_HPI.csv    hydrophobic interactions
      │   ├── *_HB.csv     hydrogen bonds
      │   ├── *_PS.csv     pi-stacking interactions
      │   ├── *_PC.csv     pi-cation interactions
      │   └── *_SB.csv     salt bridges
      ├── structures/
      │   └── *_prot.pdb   corrected/protonated PDB used for analysis
      └── pymol/
          ├── *_temp.pdb   temporary visualization structure
          └── *.pse         PyMOL session

Only interaction CSV files for interaction types actually detected by PLIP
are written. The output filename controls the 2D figure name.

If --output_dir is omitted, PanViz preserves the original PLIPViz behavior
and creates <input_stem>_output beside the input PDB.

CUSTOMIZATION
-------------
--canvas_height / --canvas_width : 2D canvas dimensions in pixels
-o / --out_file                  : 2D figure filename (.png or .svg)
--output_dir                     : root directory for organized results
--pymol / --no-pymol             : enable/disable PyMOL session generation
-y / --analyse_all               : analyze all eligible small molecules

INTERACTION CODES
-----------------
HPI = Hydrophobic interaction
HB  = Hydrogen bond
PS  = Pi-stacking
PC  = Pi-cation
SB  = Salt bridge

EXAMPLE
-------
PanViz -f complex.pdb -y -o ligand_interactions.png \
       --output_dir PanViz_results --canvas_width 1400 --canvas_height 1000 \
       --no-pymol
"""
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'PanViz {PANVIZ_VERSION}',
    )

    parser.add_argument(
        '-f', '--file',
        type=str,
        required=True,
        metavar='FILE',
        help='Input PDB complex for PLIP/PanViz analysis.'
    )

    pymol_group = parser.add_mutually_exclusive_group()
    pymol_group.add_argument(
        '--pymol',
        dest='pymol',
        action='store_true',
        help='Generate and save a PyMOL session.'
    )
    pymol_group.add_argument(
        '--no-pymol',
        dest='pymol',
        action='store_false',
        help='Do not generate a PyMOL session.'
    )
    parser.set_defaults(pymol=False)

    parser.add_argument(
        '--canvas_height',
        type=int,
        default=700,
        metavar='HEIGHT',
        help='2D output canvas height in pixels (default: 700).'
    )

    parser.add_argument(
        '--canvas_width',
        type=int,
        default=1000,
        metavar='WIDTH',
        help='2D output canvas width in pixels (default: 1000).'
    )

    parser.add_argument(
        '-o', '--out_file',
        type=str,
        default='PanViz_interactions.png',
        metavar='FILE',
        help='2D figure filename; must end in .png or .svg (default: PanViz_interactions.png).'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        default=None,
        metavar='DIR',
        help='Root directory for organized PanViz results. If omitted, use the original native output location.'
    )

    parser.add_argument(
        '-y', '--analyse_all',
        dest='analyse_all',
        action='store_true',
        help='Analyze all eligible small-molecule binding sites automatically; otherwise prompt for each object.'
    )

    return parser.parse_args()


def main(args):
    my_mol = PDBComplex()
    my_mol.load_pdb(args.file)

    binding_sites = [
        x for x in str(my_mol).split("\n")[1:]
        if x.split(":")[0] not in ["ARN", "ASH", "GLH", "LYN", "HIE", "HIP"]
    ]

    if not binding_sites:
        raise RuntimeError("No eligible small-molecule binding sites were detected in the input PDB.")

    for my_id in binding_sites:
        user_input = ''

        if not args.analyse_all:
            while True:
                user_input = input('Do you want to analyze object {}? y/n: '.format(my_id))

                if user_input.lower() in ['n', 'y']:
                    break
                else:
                    print('Do you want to analyze object {}? y/n: '.format(my_id))
        else:
            user_input = 'y'

        if user_input.lower() == 'n':
            continue

        plip_2d_interactions(
            args.file,
            my_id,
            save_files=True,
            save_pymol=args.pymol,
            canvas_height=args.canvas_height,
            canvas_width=args.canvas_width,
            out_name=args.out_file,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    args = parse_args()
    main(args)
