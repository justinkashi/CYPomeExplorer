import math
import numpy as np
import argparse
import os
from multiprocessing import Pool, cpu_count

def read_pdb_coords(pdb_filename):
    coords = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                # Split-based reading for lattice coordinates
                parts = line.split()
                if len(parts) >= 6:
                    # Extracts X, Y, Z from the standard PDB/PQR columns
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)

def cavity_exact(name, surface_coords, radius_sphere):
    radius_limit = radius_sphere + 2
    prot_coords, prot_radii = [], []
    
    with open(name, "r") as f:
        for line in f:
            parts = line.split()
            # PQR lines must have at least 10 items (ATOM, Serial, Name, Res, Num, X, Y, Z, Q, R)
            if line.startswith("ATOM") and len(parts) >= 10:
                if parts[3] == "HEM": continue 
                
                try:
                    # Counting from the back: -1 is Radius, -2 is Charge, -3 is Z, -4 is Y, -5 is X
                    r = float(parts[-1])
                    z = float(parts[-3])
                    y = float(parts[-4])
                    x = float(parts[-5])
                    
                    # Filtering identical to your original surface.py
                    if z > -2:
                        d = math.sqrt(x*x + y*y + z*z)
                        if d < radius_limit:
                            prot_coords.append([x, y, z])
                            prot_radii.append(r)
                except (ValueError, IndexError):
                    continue
                    
    prot_coords = np.array(prot_coords)
    dist_results = []
def worker(task):
    filename, surface_coords, radius = task
    print(f"   [Proc {os.getpid()}] Processing: {os.path.basename(filename)}", flush=True)
    dist_results = cavity_exact(filename, surface_coords, radius)
    name = os.path.splitext(os.path.basename(filename))[0]
    return f"{name} {' '.join(map(str, dist_results))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-r','--radius', type=int, required=True)
    args = parser.parse_args()

    surface_coords = read_pdb_coords(args.pdb)
    tasks = [(f, surface_coords, args.radius) for f in args.name]
    
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(worker, tasks)

    with open(args.output, "w") as f:
        for line in results:
            f.write(line + "\n")
    print(f"✅ Surface extraction complete: {args.output}")