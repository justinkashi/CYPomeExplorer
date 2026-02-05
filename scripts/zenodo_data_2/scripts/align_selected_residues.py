#!/usr/bin/env python3

import sys
from pymol import cmd
from pymol import finish_launching

def parse_ranges(ranges_str):
    """Parses input like '150-167,242-273' into a list of (start, end) tuples"""
    return [tuple(map(int, part.strip().split('-'))) for part in ranges_str.split(',')]

def align_custom(ref_pdb, mobile_pdb, output_pdb, residue_ranges, ref_name="ref", mobile_name="mobile"):
    # Start PyMOL without GUI
    finish_launching(['pymol', '-cq'])

    # Load the structures
    cmd.load(ref_pdb, ref_name)
    cmd.load(mobile_pdb, mobile_name)

    # Format the selection string
    range_str = "+".join([f"{start}-{end}" for start, end in residue_ranges])
    backbone = "n. n+ca+c+o"

    ref_sel = f"{ref_name} & (i. {range_str}) & {backbone}"
    mob_sel = f"{mobile_name} & (i. {range_str}) & {backbone}"

    # Perform alignment
    cmd.align(mob_sel, ref_sel)

    # Save the aligned structure
    cmd.save(output_pdb, mobile_name)

    # Finish PyMOL session
    cmd.quit()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 align_structures.py ref.pdb mobile.pdb '150-167,242-273,401-412,304-317' output_aligned.pdb")
        sys.exit(1)

    ref_pdb = sys.argv[1]
    mobile_pdb = sys.argv[2]
    ranges_str = sys.argv[3]
    output_pdb = sys.argv[4]

    ranges = parse_ranges(ranges_str)
    align_custom(ref_pdb, mobile_pdb, output_pdb, ranges)

