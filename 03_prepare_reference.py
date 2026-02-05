"""
USAGE:
    python 03_prepare_reference.py data/4I3Q.cif

DEPENDENCIES:
    pip install biopython numpy matplotlib

OUTPUTS:
    - 4I3Q_reference.pdb
    - reference_qc_plot.png
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO, MMCIFParser


def prepare_reference(input_file, output_file="4I3Q_reference.pdb"):

    if not os.path.exists(input_file):
        print(f"ERROR: File not found at {input_file}")
        sys.exit(1)

    parser = MMCIFParser(QUIET=True) if input_file.endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("REF", input_file)
    model = structure[0]

    # -------------------------------------------------
    # 1. Locate FE
    # -------------------------------------------------
    fe_atom = None
    for atom in model.get_atoms():
        if atom.get_parent().get_resname().upper() == "HEM" and atom.get_name().upper() == "FE":
            fe_atom = atom
            break

    if fe_atom is None:
        raise ValueError("Could not find heme FE atom.")

    fe_coord = fe_atom.get_coord()

    # -------------------------------------------------
    # 2. Find nitrogens + proximal cysteine
    # -------------------------------------------------
    n_atoms = {}
    prox_cys_sg = None
    min_dist = 999

    for atom in model.get_atoms():
        res = atom.get_parent().get_resname().upper()
        name = atom.get_name().upper()

        if res == "HEM" and name in {"NA", "NB", "NC", "ND"}:
            n_atoms[name] = atom

        if res == "CYS" and name == "SG":
            d = np.linalg.norm(atom.get_coord() - fe_coord)
            if d < 4.5 and d < min_dist:
                prox_cys_sg = atom
                min_dist = d

    if len(n_atoms) != 4:
        raise ValueError("Did not find 4 heme nitrogens.")

    # -------------------------------------------------
    # 3. Translate FE → origin
    # -------------------------------------------------
    for atom in model.get_atoms():
        atom.set_coord(atom.get_coord() - fe_coord)

    # -------------------------------------------------
    # 4. Fit heme plane via SVD
    # -------------------------------------------------
    coords = np.array([a.get_coord() for a in n_atoms.values()])
    _, _, vh = np.linalg.svd(coords)
    normal = vh[-1]
    normal /= np.linalg.norm(normal)

    # -------------------------------------------------
    # 5. Ensure distal side +Z
    # -------------------------------------------------
    if prox_cys_sg is not None:
        if np.dot(normal, prox_cys_sg.get_coord()) > 0:
            normal = -normal

    # -------------------------------------------------
    # 6. Rotate normal → +Z
    # -------------------------------------------------
    target = np.array([0, 0, 1])
    v = np.cross(normal, target)
    s = np.linalg.norm(v)
    c = np.dot(normal, target)

    if s < 1e-8:
        rot = np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]],
                       [v[2], 0, -v[0]],
                       [-v[1], v[0], 0]])
        rot = np.eye(3) + vx + (vx @ vx) * ((1 - c) / s**2)

    for atom in model.get_atoms():
        atom.set_coord(rot @ atom.get_coord())

    # -------------------------------------------------
    # 7. In-plane rotation (NA → +X)  <<< KUVEK STEP
    # -------------------------------------------------
    na_atom = n_atoms["NA"]
    na_vec = na_atom.get_coord()

    na_xy = np.array([na_vec[0], na_vec[1], 0.0])
    angle = np.arctan2(na_xy[1], na_xy[0])

    cos_a, sin_a = np.cos(-angle), np.sin(-angle)
    rot_z = np.array([[cos_a, -sin_a, 0],
                      [sin_a,  cos_a, 0],
                      [0, 0, 1]])

    for atom in model.get_atoms():
        atom.set_coord(rot_z @ atom.get_coord())

    # -------------------------------------------------
    # QC
    # -------------------------------------------------
    final_fe = fe_atom.get_coord()
    final_ns = np.array([a.get_coord() for a in n_atoms.values()])

    print("\n--- QC REPORT ---")
    print("FE position:", final_fe)
    print("Mean N Z:", np.mean(final_ns[:, 2]))
    print("NA direction:", na_atom.get_coord() / np.linalg.norm(na_atom.get_coord()))
    print("----------------\n")

    # Plot
    plt.figure(figsize=(5, 5))
    plt.scatter(final_ns[:, 0], final_ns[:, 1])
    plt.scatter(0, 0, marker="x", s=100)
    plt.axhline(0)
    plt.axvline(0)
    plt.gca().set_aspect("equal")
    plt.savefig("reference_qc_plot.png")

    # -------------------------------------------------
    # FINAL SAVE (must be LAST)
    # -------------------------------------------------
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python 03_prepare_reference.py <input.cif|pdb>")
    else:
        prepare_reference(sys.argv[1])
