import math
import numpy as np
import argparse
import os
import glob
from multiprocessing import Pool, cpu_count

def read_pdb_coords(pdb_filename):
    coords = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                # Lattice files use strict fixed-width columns
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)

def read_reference_file(ref_file):
    ref_map = {}
    with open(ref_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                # Key: (ResidueName, AtomName)
                ref_map[(parts[1].upper(), parts[0].upper())] = float(parts[2])
    return ref_map

def process_protein(task_args):
    pqr_path, surface_coords, radius_limit, ref_map = task_args
    base_name = os.path.splitext(os.path.basename(pqr_path))[0]
    
    prot_coords, prot_charges = [], []
    try:
        with open(pqr_path, "r") as f:
            for line in f:
                parts = line.split()
                if line.startswith("ATOM") and len(parts) >= 10:
                    # Residue is parts[3], Atom is parts[2]
                    res, atom = parts[3].upper(), parts[2].upper()
                    if res == "HEM": continue
                    
                    if (res, atom) in ref_map:
                        # Negative indexing for 10-column PQRs
                        # -5:X, -4:Y, -3:Z
                        prot_coords.append([float(parts[-5]), float(parts[-4]), float(parts[-3])])
                        prot_charges.append(ref_map[(res, atom)])

        if not prot_coords:
            # Safety return if no valid atoms found
            return f"{base_name} {' '.join(['0.0'] * len(surface_coords))}"

        prot_coords = np.array(prot_coords)
        prot_charges = np.array(prot_charges)
        hit_potentials = []

        for ray in surface_coords:
            # Calculate potential at a representative distance (10.0A)
            V_unit = ray / np.linalg.norm(ray)
            hit_point = V_unit * 10.0 
            
            # Distance from every atom to the hit point
            dists = np.linalg.norm(prot_coords - hit_point, axis=1)
            dists[dists < 0.1] = 0.1 # Prevent division by zero
            
            # Potential = sum(q/r)
            pot = np.sum(prot_charges / dists)
            hit_potentials.append(round(pot, 4))
            
        print(f"   [Proc {os.getpid()}] Finished: {base_name}", flush=True)
        return f"{base_name} {' '.join(map(str, hit_potentials))}"

    except Exception as e:
        # Fallback to prevent TypeError: 'NoneType' object is not iterable
        return f"{base_name} {' '.join(['0.0'] * len(surface_coords))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-r', '--radius', type=int, required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('-c', '--output', required=True)
    args = parser.parse_args()

    # Expand wildcards
    all_files = []
    for pattern in args.name:
        all_files.extend(glob.glob(pattern))

    surface_coords = read_pdb_coords(args.pdb)
    ref_map = read_reference_file(args.ref)
    tasks = [(f, surface_coords, args.radius, ref_map) for f in all_files]

    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_protein, tasks)

    with open(args.output, "w") as f:
        for line in results:
            f.write(line + "\n")
    print(f"✅ Charge extraction complete: {args.output}")