#!/bin/bash

# Set paths
REF_PDB="/pool/teakuvek/new_vectors/CYP3A4_4I3Q/reference.pdb"  # <-- Replace this with your actual reference PDB file
TRAJ_DIR="/pool/teakuvek/new_vectors/CYP3A4_4I3Q/snapshots/"
OUTPUT_DIR="/pool/teakuvek/new_vectors/CYP3A4_4I3Q/aligned_snapshots"
PYTHON_SCRIPT="align_selected_residues.py"
RANGES="5-8,32-40,45-49,52-56,62-66,97-105,115-136,146-158,193-197,219-230,265-277,283-295,301-309,320-324,330-339,346-349,366-369,371-375,418-433,464-468"

# Create output dir if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over frames
for pdb in "$TRAJ_DIR"/FRAME_*.pdb; do
    base=$(basename "$pdb" .pdb)
    frame_num=$(echo "$base" | sed 's/FRAME_//')
    
    # Treat frame_num as base 10 to avoid octal error
    output_file=$(printf "%s/%05d_bb.pdb" "$OUTPUT_DIR" "$((10#$frame_num))")

    echo "Aligning $pdb -> $output_file"
    python3 "$PYTHON_SCRIPT" "$REF_PDB" "$pdb" "$RANGES" "$output_file"
done
