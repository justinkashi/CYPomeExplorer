import math
import numpy as np
import argparse
import os
from multiprocessing import Pool, cpu_count

def read_lattice_coords(pdb_file):
    """Reads XYZ coordinates from a standard PDB file."""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    # PDB fixed-width columns for coordinates
                    coords.append([
                        float(line[30:38]), 
                        float(line[38:46]), 
                        float(line[46:54])
                    ])
                except ValueError:
                    continue
    return np.array(coords)

def calculate_cavity_distances(pqr_path, lattice, radius_limit):
    """Calculates distances from origin to protein surface along lattice rays."""
    atoms = []
    radii = []
    
    with open(pqr_path, 'r') as f:
        for line in f:
            parts = line.split()
            # PQR lines usually have at least 10 columns (ATOM ... X Y Z Charge Radius)
            if line.startswith("ATOM") and len(parts) >= 10:
                # Exclude prosthetic groups or ligands if they are at the center
                if parts[3] == "HEM":
                    continue
                
                try:
                    # Negative indexing handles files with or without Chain IDs
                    # Radius is usually the last column, Coordinates precede Charge
                    r = float(parts[-1])
                    z = float(parts[-3])
                    y = float(parts[-4])
                    x = float(parts[-5])
                    
                    # Filter for atoms within a relevant vertical or radial range
                    if z > -2:
                        atoms.append([x, y, z])
                        radii.append(r)
                except ValueError:
                    continue
                    
    if not atoms:
        return [float(radius_limit)] * len(lattice)
        
    atoms = np.array(atoms)
    dist_results = []

    for ray in lattice:
        # Normalize to unit vector
        norm = np.linalg.norm(ray)
        if norm == 0:
            dist_results.append(float(radius_limit))
            continue
            
        V = ray / norm
        min_t = float(radius_limit)
        
        # Ray-sphere intersection check for each atom
        for center, r in zip(atoms, radii):
            t_proj = np.dot(center, V)
            if t_proj > 0:
                dist_sq = np.dot(center, center) - t_proj**2
                if dist_sq < r**2:
                    t_hit = t_proj - math.sqrt(max(0, r**2 - dist_sq))
                    if 0 <= t_hit < min_t:
                        min_t = t_hit
        dist_results.append(round(min_t, 3))
        
    return dist_results

def run_worker(task_args):
    """Multiprocessing worker to process individual files."""
    file_path, lattice, radius = task_args
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # A try-except block ensures the worker never returns None,
    # preventing 'NoneType' iteration errors in the main pool.
    try:
        results = calculate_cavity_distances(file_path, lattice, radius)
        if results is None:
            results = [float(radius)] * len(lattice)
    except Exception:
        results = [float(radius)] * len(lattice)
        
    return f"{base_name} {' '.join(map(str, results))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-r', '--radius', type=int, default=15)
    
    args = parser.parse_args()

    # Load shared lattice once
    lattice_coords = read_lattice_coords(args.pdb)
    
    # Map tasks to worker pool
    tasks = [(f, lattice_coords, args.radius) for f in args.name]
    
    with Pool(processes=cpu_count()) as pool:
        final_data = pool.map(run_worker, tasks)
        
    # Write aggregated results to output file
    with open(args.output, 'w') as out_file:
        for entry in final_data:
            out_file.write(entry + '\n')