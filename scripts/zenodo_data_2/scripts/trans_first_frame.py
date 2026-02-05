#!/usr/bin/env python3
import argparse
from pymol import cmd

def translate_to_fe_origin(pdb_file, output_pdb):
    """Translate so that FE atom of HEM/HEMO* residue is at the origin."""
    cmd.reinitialize()
    cmd.load(pdb_file, "structure")

    # Match any residue starting with HEM (HEM, HEMO, HEMOA, etc.)
    cmd.select("fe_atom", "resn HEM* and name FE")
    if cmd.count_atoms("fe_atom") == 0:
        print(f"⚠️ No FE atom found in residues HEM* in {pdb_file}")
        return False

    fe_coords = cmd.get_coords("fe_atom")[0]
    translation_vector = [-fe_coords[0], -fe_coords[1], -fe_coords[2]]

    cmd.translate(translation_vector, "all")
    cmd.save(output_pdb, "structure")
    print(f"✅ Saved translated structure to {output_pdb}")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate structure so FE in HEM*/HEMO* is at origin.")
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output PDB file")
    args = parser.parse_args()

    translate_to_fe_origin(args.input_pdb, args.output_pdb)

