import math
import numpy as np
import argparse
import os
import sys
from multiprocessing import Pool, cpu_count

def read_pdb_coords(pdb_filename):
    coords = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
    return np.array(coords)

def cavity_exact(name, surface_coords, radius_sphere):
    radius_limit = radius_sphere + 2
    protein_coords = []
    atom_radius = []
    
    # Filtering identical to surface.py
    with open(name, "r") as f:
        for line in f:
            if line[0:4] == "ATOM" and float(line[46:53]) > -2 and line[17:20] not in ["HEM"]:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:53])
                r = float(line[69:75])
                d = math.sqrt(x * x + y * y + z * z)
                if d < radius_limit:
                    protein_coords.append([x, y, z])
                    atom_radius.append(r)
                    
    protein_coords = np.array(protein_coords)
    atom_radius = np.array(atom_radius)

    def new_vector(V, t1):
        norm_V = np.linalg.norm(V)
        V_unit = V / norm_V
        return t1 * V_unit    

    distance_vectors = []
    for vector in surface_coords:
        # Replicating vector mutation
        V = np.array(vector)
        min_dist = float('inf')
        t1 = float('inf')
        temp_vector = V 
        
        for atom, r in zip(protein_coords, atom_radius):
            sphere_touching = True
            S = np.array(atom)
            d = np.round(np.linalg.norm(S), 3)
            
            # Replicating cumulative rounding drift
            cos_alpha = np.dot(S, V) / (d * np.round(np.linalg.norm(V), 3))
            
            if np.array_equal(S, V):
                t1 = d - r
                V = new_vector(V, t1)
            elif 0 <= cos_alpha <= 1:
                y = round(d * math.sqrt(1.0 - cos_alpha * cos_alpha), 3)
                if y < r:
                    x = math.sqrt(r * r - y * y)
                    t = np.round(abs(np.dot(S, V)) / np.linalg.norm(V), 3)
                    t1 = t - x
                    
                    hit_dist = t1 if t1 < radius_sphere else radius_sphere
                    temp_vector = new_vector(V, hit_dist)

                    # Replication of Hemisphere Gate
                    if temp_vector[2] >= 0:
                        V = temp_vector
                elif y == r:
                    t1 = np.dot(S, V) / np.linalg.norm(V)
                    hit_dist = t1 if t1 < radius_sphere else radius_sphere
                    temp_vector = new_vector(V, hit_dist)
                else:
                    sphere_touching = False
            else:
                sphere_touching = False
                
            if sphere_touching and t1 < min_dist:
                min_dist = t1
        
        if min_dist > radius_sphere:
            min_dist = radius_sphere
        distance_vectors.append(np.round(min_dist, 3))
    
    return distance_vectors

def worker(task_info):
    filename, surface_coords, radius = task_info
    if os.path.isdir(filename):
        print(f"   [!] Skipping directory: {filename}")
        return None
        
    # Print statement for start
    print(f"   [Proc {os.getpid()}] Starting: {os.path.basename(filename)}")
    
    dist_results = cavity_exact(filename, surface_coords, radius)
    name = os.path.splitext(os.path.basename(filename))[0]
    
    # Print statement for finish
    print(f"   [Proc {os.getpid()}] Finished: {os.path.basename(filename)}")
    return f"{name} {' '.join(map(str, dist_results))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-r','--radius', type=int, required=True)
    args = parser.parse_args()

    # Initial Progress Prints
    print(f"--- Starting Parallel Surface Extraction ---")
    surface_coords = read_pdb_coords(args.pdb)
    print(f"Lattice points: {len(surface_coords)}")
    print(f"Using {cpu_count()} CPU cores to process {len(args.name)} files...")

    tasks = [(f, surface_coords, args.radius) for f in args.name]
    
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(worker, tasks)
    
    # Filter results and write to output
    with open(args.output, 'w') as f:
        f.write('\n'.join(filter(None, results)) + '\n')
        
    print(f"✅ All tasks complete! Results saved to {args.output}")