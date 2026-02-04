"""
USAGE:
    python 03_prepare_reference.py data/4I3Q.cif

DEPENDENCIES:
    pip install biopython numpy matplotlib

OUTPUTS:
    - 4I3Q_reference.pdb: Centered and oriented master template.
    - reference_qc_plot.png: Visual confirmation of Heme alignment.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO, MMCIFParser

def prepare_reference(input_file, output_file="4I3Q_reference.pdb"):
    # Check if file exists
    if not os.path.exists(input_file):
        print(f"ERROR: File not found at {input_file}")
        sys.exit(1)

    parser = MMCIFParser(QUIET=True) if input_file.endswith('.cif') else PDBParser(QUIET=True)
    structure = parser.get_structure("REF", input_file)
    model = structure[0]

    # 1. Locate Heme Iron (FE)
    fe_atom = None
    for atom in model.get_atoms():
        res_name = atom.get_parent().get_resname().strip().upper()
        atom_name = atom.get_name().strip().upper()
        if res_name == "HEM" and atom_name == "FE":
            fe_atom = atom
            break

    if not fe_atom:
        raise ValueError("Reference Error: Could not find Heme Iron (FE).")

    # 2. Find Nitrogens and Proximal Cys
    fe_coord = fe_atom.get_coord()
    n_atoms_dict = {}
    prox_cys_sg = None
    min_dist_sg = 999.0

    for atom in model.get_atoms():
        res_name = atom.get_parent().get_resname().strip().upper()
        atom_name = atom.get_name().strip().upper()
        if res_name == "HEM" and atom_name in {"NA", "NB", "NC", "ND"}:
            n_atoms_dict[atom_name] = atom
        if res_name == "CYS" and atom_name == "SG":
            dist = np.linalg.norm(atom.get_coord() - fe_coord)
            if dist < 4.5 and dist < min_dist_sg:
                prox_cys_sg, min_dist_sg = atom, dist

    if len(n_atoms_dict) != 4:
        raise ValueError(f"Found {len(n_atoms_dict)} nitrogens. Need exactly 4.")

    # 3. Translation to Origin
    for atom in model.get_atoms():
        atom.set_coord(atom.get_coord() - fe_coord)

    # 4. Plane Fitting (SVD)
    n_coords = np.array([n.get_coord() for n in n_atoms_dict.values()])
    _, _, vh = np.linalg.svd(n_coords)
    normal_vector = vh[2, :]
    normal_vector /= np.linalg.norm(normal_vector)

    # 5. Flip for Distal side (+z)
    if prox_cys_sg:
        if np.dot(normal_vector, prox_cys_sg.get_coord()) > 0:
            normal_vector = -normal_vector

    # 6. Rotation Matrix
    target_z = np.array([0, 0, 1])
    v = np.cross(normal_vector, target_z)
    s, c = np.linalg.norm(v), np.dot(normal_vector, target_z)

    if s < 1e-6:
        rot = np.eye(3) if c > 0 else -np.eye(3)
        if c < 0: rot[2, 2] = 1
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rot = np.eye(3) + vx + (vx @ vx) * ((1 - c) / (s**2))

    for atom in model.get_atoms():
        atom.set_coord(np.dot(rot, atom.get_coord()))

    # --- QC REPORT ---
    final_fe = fe_atom.get_coord()
    final_ns = np.array([n.get_coord() for n in n_atoms_dict.values()])
    
    print("-" * 30)
    print(f"QC REPORT FOR: {input_file}")
    print(f"Iron Position: {final_fe} (Expected: [0,0,0])")
    print(f"Avg Nitrogen Z: {np.mean(final_ns[:, 2]):.4f} (Expected: ~0.0)")
    if prox_cys_sg:
        print(f"Proximal Cys Z: {prox_cys_sg.get_coord()[2]:.2f} (Expected: Negative)")
    print("-" * 30)

    # --- QC VISUALIZATION ---
    plt.figure(figsize=(6, 6))
    plt.scatter(final_ns[:, 0], final_ns[:, 1], color='blue', label='Heme Nitrogens')
    plt.scatter(final_fe[0], final_fe[1], color='red', marker='X', s=100, label='Iron (Origin)')
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(0, color='black', lw=0.5)
    plt.title("Heme Alignment QC (XY Plane View)")
    plt.xlabel("X (Å)")
    plt.ylabel("Y (Å)")
    plt.legend()
    plt.grid(True, linestyle='--')
    plt.savefig("reference_qc_plot.png")
    print("QC Plot saved: reference_qc_plot.png")

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python 03_prepare_reference.py <path_to_cif>")
    else:
        prepare_reference(sys.argv[1])