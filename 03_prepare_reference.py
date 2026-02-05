#!/usr/bin/env python3

"""
USAGE
    python 03_prepare_reference.py data/4i3q.cif

DEPENDENCIES
    pip install biopython numpy

OUTPUT
    4I3Q_reference.pdb
"""

import sys
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO


# ---------------------------------------------------------
# helpers
# ---------------------------------------------------------

def rodrigues(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)

    if s < 1e-8:
        return np.eye(3)

    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s**2))


# ---------------------------------------------------------
# main
# ---------------------------------------------------------

def prepare_reference(infile, outfile="4I3Q_reference.pdb"):

    parser = MMCIFParser(QUIET=True) if infile.endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("REF", infile)
    model = structure[0]

    fe = None
    Ns = {}
    prox_sg = None
    best = 999.0

    for atom in model.get_atoms():
        res = atom.get_parent().get_resname().upper()
        name = atom.get_name().upper()

        if res == "HEM" and name == "FE":
            fe = atom

        if res == "HEM" and name in {"NA", "NB", "NC", "ND"}:
            Ns[name] = atom

    if fe is None or len(Ns) != 4:
        raise RuntimeError("Failed to find heme atoms")

    fe_pos = fe.get_coord()

    # find proximal cysteine
    for atom in model.get_atoms():
        if atom.get_parent().get_resname().upper() == "CYS" and atom.get_name().upper() == "SG":
            d = np.linalg.norm(atom.get_coord() - fe_pos)
            if d < 4.5 and d < best:
                best = d
                prox_sg = atom

    # ---------------------------------------------------------
    # 1. translate FE -> origin
    # ---------------------------------------------------------
    for a in model.get_atoms():
        a.set_coord(a.get_coord() - fe_pos)

    # ---------------------------------------------------------
    # 2. plane fit (CORRECT: centered SVD)
    # ---------------------------------------------------------
    n_coords = np.array([a.get_coord() for a in Ns.values()])
    centroid = n_coords.mean(axis=0)
    centered = n_coords - centroid

    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    normal /= np.linalg.norm(normal)

    # distal side should face +Z
    if prox_sg is not None and np.dot(normal, prox_sg.get_coord()) > 0:
        normal = -normal

    # ---------------------------------------------------------
    # 3. rotate normal -> +Z
    # ---------------------------------------------------------
    R = rodrigues(normal, np.array([0, 0, 1]))

    for a in model.get_atoms():
        a.set_coord(R @ a.get_coord())

    # ---------------------------------------------------------
    # 4. in-plane rotation (FE->NA -> +X)
    # ---------------------------------------------------------
    na = Ns["NA"].get_coord()
    na_xy = np.array([na[0], na[1], 0])

    theta = np.arctan2(na_xy[1], na_xy[0])

    cz = np.cos(-theta)
    sz = np.sin(-theta)

    Rz = np.array([
        [cz, -sz, 0],
        [sz,  cz, 0],
        [0,   0,  1]
    ])

    for a in model.get_atoms():
        a.set_coord(Rz @ a.get_coord())

    # ---- FINAL QC (correct metrics) ----
    Ns_after = np.array([a.get_coord() for a in Ns.values()])
    fe_after = fe.get_coord()

    # Recompute best-fit normal after all rotations (centered SVD)
    centroid2 = Ns_after.mean(axis=0)
    centered2 = Ns_after - centroid2
    _, _, vh2 = np.linalg.svd(centered2, full_matrices=False)
    normal2 = vh2[-1]
    normal2 /= np.linalg.norm(normal2)
    # force "up" for reporting
    if normal2[2] < 0:
        normal2 *= -1

    angle_to_z = np.degrees(np.arccos(np.clip(normal2[2], -1.0, 1.0)))

    # z spread (tilt/planarity proxy)
    z_mean = Ns_after[:, 2].mean()
    z_std  = Ns_after[:, 2].std()
    z_max_dev = np.max(np.abs(Ns_after[:, 2] - z_mean))

    # NA in-plane lock check
    na = Ns["NA"].get_coord()
    na_xy = np.array([na[0], na[1], 0.0])
    na_xy /= np.linalg.norm(na_xy)

    print("\n--- FINAL QC (Kuvek-frame) ---")
    print("FE:", fe_after, "(expect ~[0,0,0])")
    print("heme normal:", normal2, "angle_to_+Z(deg):", angle_to_z, "(expect ~0)")
    print("N z-mean:", z_mean, "(offset OK)")
    print("N z-std:", z_std, "N z-maxdev:", z_max_dev, "(expect small, e.g. <0.05–0.1 Å)")
    print("NA_xy:", na_xy, "(expect x>0, y~0)")
    print("-----------------------------\n")


    # ---------------------------------------------------------
    # save
    # ---------------------------------------------------------
    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile)


# ---------------------------------------------------------

if __name__ == "__main__":
    prepare_reference(sys.argv[1])
