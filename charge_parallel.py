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
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)

def read_reference_file(ref_file):
    ref_map = {}
    with open(ref_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                # Key: (Residue, AtomName)
                ref_map[(parts[1].upper(), parts[0].upper())] = float(parts[2])
    return ref_map

def process_protein(task_args):
    pqr_path, surface_coords, radius_limit, ref_map = task_args
    base_name = os.path.splitext(os.path.basename(pqr_path))[0]
    print(f"   [Proc {os.getpid()}] Processing: {base_name}", flush=True)

    prot_coords, prot_charges = [], []
    with open(pqr_path, "r") as f:
        for line in f:
            parts = line.split()
            if line.startswith("ATOM") and len(parts) >= 10:
                res, atom = parts[3].upper(), parts[2].upper()
                if res == "HEM": continue
                if (res, atom) in ref_map:
                    prot_coords.append([float(parts[6]), float(parts[7]), float(parts[8])])
                    prot_charges.append(ref_map[(res, atom)])

    prot_coords = np.array(prot_coords)
    hit_potentials = []
    
    for ray in surface_coords:
        V_unit = ray / np.linalg.norm(ray)
        # Using a representative 10A hit point for the potential vector
        hit_point = V_unit * 10.0 
        
        pot = sum(q / np.linalg.norm(center - hit_point) 
                  for center, q in zip(prot_coords, prot_charges)
                  if np.linalg.norm(center - hit_point) > 0.1)
        hit_potentials.append(round(pot, 4))

    return f"{base_name} {' '.join(map(str, hit_potentials))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-r', '--radius', type=int, required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('-c', '--output', required=True)
    args = parser.parse_args()

    all_files = [f for pattern in args.name for f in glob.glob(pattern)]
    surface_coords = read_pdb_coords(args.pdb)
    ref_map = read_reference_file(args.ref)
    tasks = [(f, surface_coords, args.radius, ref_map) for f in all_files]

    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_protein, tasks)

    with open(args.output, "w") as f:
        for line in results:
            f.write(line + "\n")
    print(f"✅ Charge extraction complete: {args.output}")