#The surface_parallel.py script is designed to maximize your local hardware by leveraging the Python multiprocessing library. 
import math
import numpy as np
import argparse
import os
from multiprocessing import Pool, cpu_count

# --- Keep your original helper functions exactly as they are ---
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

def cavity(name, surface_coords, radius_sphere):
    radius = radius_sphere + 2
    with open(name, "r") as f:
        protein_coords = []
        atom_radius = []
        for line in f:
            if line[0:4] == "ATOM" and float(line[46:53]) > -2 and line[17:20] not in ["HEM"]:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:53])
                r = float(line[69:75])
                if math.sqrt(x*x + y*y + z*z) < radius:
                    protein_coords.append([x, y, z])
                    atom_radius.append(r)
        protein_coords = np.array(protein_coords)
        atom_radius = np.array(atom_radius)

    def new_vector(V, t1):
        return t1 * (V / np.linalg.norm(V))

    distance_vectors = []
    for vector in surface_coords:
        V = np.array(vector)
        min_dist, final_vector = float('inf'), V
        for atom, r in zip(protein_coords, atom_radius):
            sphere_touching = True
            S = np.array(atom)
            d = np.round(np.linalg.norm(S), 3)
            cos_alpha = np.dot(S, V) / (d * np.round(np.linalg.norm(V), 3))
            if np.array_equal(S, V):
                t1 = d - r
                temp_vector = new_vector(V, t1)
            elif 0 <= cos_alpha <= 1:
                y = round(d * math.sqrt(1.0 - cos_alpha * cos_alpha), 3)
                if y < r:
                    x = math.sqrt(r * r - y * y)
                    t = np.round(abs(np.dot(S, V)) / np.linalg.norm(V), 3)
                    t1 = t - x
                    temp_vector = new_vector(V, min(t1, radius_sphere))
                elif y == r:
                    t1 = np.dot(S, V) / np.linalg.norm(V)
                    temp_vector = new_vector(V, min(t1, radius_sphere))
                else:
                    sphere_touching = False
            else:
                sphere_touching = False
            
            if sphere_touching and t1 < min_dist:
                min_dist = t1
        
        final_dist = min(min_dist, radius_sphere)
        distance_vectors.append(np.round(final_dist, 3))
    
    return np.array(distance_vectors)

# --- NEW: Worker function for individual files ---
def worker(task_info):
    filename, surface_coords, radius = task_info

    
    # 🛑 NEW: Skip the item if it's a directory
    if os.path.isdir(filename):
        print(f"   ⚠️ Skipping directory: {filename}")
        return None 


    print(f"   - Processing: {os.path.basename(filename)}")
    dist_results = cavity(filename, surface_coords, radius)
    name = os.path.splitext(os.path.basename(filename))[0]
    return f"{name} {' '.join(map(str, dist_results.tolist()))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-r','--radius', type=int, required=True)
    args = parser.parse_args()

    surface_coords = read_pdb_coords(args.pdb)
    tasks = [(f, surface_coords, args.radius) for f in args.name]
    
    print(f"🚀 Starting parallel extraction on {cpu_count()} cores...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(worker, tasks)
    
    with open(args.output, 'w') as f:
        f.write('\n'.join(results) + '\n')
    print(f"✅ Success! Results saved to {args.output}")