#!/usr/bin/env python3
# Faithful PyMOL backbone alignment + Python-only QC
#
# Usage:
# python 04_batch_align.py \
#   4I3Q_reference.pdb \
#   data/structures \
#   data/04_out/aligned_struct \
#   data/04_out/qc

import os
import sys
import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser


REF = sys.argv[1]
IN_DIR = sys.argv[2]
OUT_DIR = sys.argv[3]
QC_DIR = sys.argv[4]

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)


# ---------------------------------------------------------
# EXACT PyMOL call (paper-faithful, reliable)
# ---------------------------------------------------------
PYMOL_SCRIPT = """
from pymol import cmd
import sys

ref, mob, outpdb, outtxt = sys.argv[1:5]

# The 20 conserved ranges from the author's shell script
ranges = "5-8+32-40+45-49+52-56+62-66+97-105+115-136+146-158+193-197+219-230+265-277+283-295+301-309+320-324+330-339+346-349+366-369+371-375+418-433+464-468"

cmd.load(ref, "ref")
cmd.load(mob, "mob")

# 1. Align using the author's specific conserved core ranges
sel_mob = f"mob and (resi {ranges}) and name N+CA+C+O"
sel_ref = f"ref and (resi {ranges}) and name N+CA+C+O"

# Using cealign as you did, but on the specific ranges
r = cmd.cealign(sel_mob, sel_ref)

# 2. THE FIX: Move the mobile HEME Iron exactly to (0,0,0)
# This prevents the "buried origin" 0.0 distance problem.
coords = cmd.get_coords("mob and resn HEM and name FE")
if coords is not None and len(coords) > 0:
    x, y, z = coords[0]
    cmd.translate([-x, -y, -z], "mob")
    print(f"   [PyMOL] Pocket Centered: FE moved from ({x:.2f}, {y:.2f}, {z:.2f}) to (0,0,0)")

cmd.save(outpdb, "mob")

with open(outtxt, "w") as f:
    if r and "RMSD" in r:
        f.write(str(r["RMSD"]))
    else:
        f.write("0.0")
"""


def run_align(ref, mob, outpdb, outtxt):
    tmp = "_pymol_align_tmp.py"

    with open(tmp, "w") as f:
        f.write(PYMOL_SCRIPT)

    subprocess.run(
        ["pymol", "-cq", tmp, "--", ref, mob, outpdb, outtxt],
        check=True
    )

    os.remove(tmp)


# ---------------------------------------------------------
# QC helpers (pure Python)
# ---------------------------------------------------------
parser = PDBParser(QUIET=True)


def backbone_centroid_radius(pdb_path):
    s = parser.get_structure("X", pdb_path)
    coords = [
        a.get_coord()
        for a in s.get_atoms()
        if a.get_name().strip().upper() in {"N", "CA", "C", "O"}
    ]
    coords = np.array(coords)
    return np.linalg.norm(coords.mean(axis=0))


def closest_cys_sg_z(pdb_path):
    s = parser.get_structure("X", pdb_path)
    best = None
    best_d = 1e9

    for a in s.get_atoms():
        if a.get_parent().get_resname().upper() == "CYS" and a.get_name().strip().upper() == "SG":
            d = np.linalg.norm(a.get_coord())
            if d < best_d:
                best = a
                best_d = d

    return best.get_coord()[2] if best else np.nan


# ---------------------------------------------------------
# batch alignment
# ---------------------------------------------------------
files = sorted(glob.glob(os.path.join(IN_DIR, "*.cif")))

if len(files) == 0:
    raise RuntimeError(f"No .cif files found in {IN_DIR}")

results = []

for i, cif in enumerate(files, 1):

    name = os.path.splitext(os.path.basename(cif))[0]
    print(f"[{i}/{len(files)}] aligning {name}", flush=True)

    outpdb = os.path.join(OUT_DIR, f"{name}.aligned.pdb")
    outtxt = os.path.join(QC_DIR, f"{name}.rmsd.txt")

    run_align(REF, cif, outpdb, outtxt)

    rmsd = float(open(outtxt).read().strip())
    centroid = backbone_centroid_radius(outpdb)
    cysz = closest_cys_sg_z(outpdb)

    results.append((name, rmsd, centroid, cysz))


# ---------------------------------------------------------
# save summary
# ---------------------------------------------------------
summary_path = os.path.join(QC_DIR, "alignment_summary.tsv")

with open(summary_path, "w") as f:
    f.write("id\trmsd\tcentroid_dist\tcys_sg_z\tpass_7A\n")
    for name, r, c, z in results:
        f.write(f"{name}\t{r:.4f}\t{c:.3f}\t{z:.3f}\t{int(r<=7)}\n")


# ---------------------------------------------------------
# plots
# ---------------------------------------------------------
rmsds = np.array([x[1] for x in results])
cysz = np.array([x[3] for x in results])

plt.hist(rmsds, bins=40)
plt.axvline(7)
plt.xlabel("Backbone RMSD (Å)")
plt.ylabel("Count")
plt.savefig(os.path.join(QC_DIR, "rmsd_hist.png"))
plt.close()

plt.hist(cysz[~np.isnan(cysz)], bins=40)
plt.axvline(0)
plt.xlabel("Closest Cys SG z (Å)")
plt.ylabel("Count")
plt.savefig(os.path.join(QC_DIR, "cys_z_hist.png"))
plt.close()


# ---------------------------------------------------------
# console report
# ---------------------------------------------------------
print("\nAlignment QC summary")
print("-------------------")
print("Structures:", len(results))
print("Mean RMSD:", rmsds.mean())
print("Max RMSD :", rmsds.max())
print("Pass ≤7Å:", np.sum(rmsds <= 7), "/", len(rmsds))
print("Outputs written to:", QC_DIR)
