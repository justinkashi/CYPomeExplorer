import math
import numpy as np
import argparse
import os
from multiprocessing import Pool, cpu_count

def read_pdb_coords(pdb_filename):
    coords = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)

def process_protein(task_info):
    pqr_path, surface_coords, radius_limit = task_info
    base_name = os.path.splitext(os.path.basename(pqr_path))[0]
    
    prot_coords, prot_radii, prot_charges = [], [], []
    with open(pqr_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and float(line[46:53]) > -2 and line[17:20] not in ["HEM"]:
                # Coordinates (Cols 31-54)
                prot_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                # Charge (Col 9 / bits 55-62)
                prot_charges.append(float(line[54:62].strip()))
                # Radius (Col 10 / bits 69-75)
                prot_radii.append(float(line[68:75].strip()))

    prot_coords = np.array(prot_coords)
    prot_radii = np.array(prot_radii)
    
    final_distances, final_charges = [], []

    for ray in surface_coords:
        V = ray / np.linalg.norm(ray) # Lattice ray direction
        min_t = float(radius_limit)
        
        # 1. Fixed Surface Calculation
        for center, r in zip(prot_coords, prot_radii):
            S = center
            d_origin = np.linalg.norm(S)
            cos_alpha = np.dot(S, V) / d_origin
            if 0 <= cos_alpha <= 1:
                y = d_origin * math.sqrt(1.0 - cos_alpha**2)
                if y < r:
                    t_hit = (d_origin * cos_alpha) - math.sqrt(r**2 - y**2)
                    if t_hit < min_t and (t_hit * V[2]) >= 0:
                        min_t = max(0, t_hit)

        # 2. Charge Calculation at the Surface Hit Point
        hit_point = V * min_t
        potential = 0.0
        for center, q in zip(prot_coords, prot_charges):
            dist_to_hit = np.linalg.norm(center - hit_point)
            if dist_to_hit > 0.1: 
                potential += q / dist_to_hit
        
        final_distances.append(round(min_t, 3))
        final_charges.append(round(potential, 4))

    return f"{base_name} {' '.join(map(str, final_distances))}", f"{base_name} {' '.join(map(str, final_charges))}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--names', nargs='+', required=True)
    parser.add_argument('-pdb', '--lattice', required=True)
    parser.add_argument('-r', '--radius', type=int, default=15)
    args = parser.parse_args()

    surface_coords = read_pdb_coords(args.lattice)
    tasks = [(f, surface_coords, args.radius) for f in args.names]

    with Pool(processes=max(1, cpu_count() - 1)) as pool:
        results = pool.map(process_protein, tasks)

    with open("surface_raw.txt", "w") as f_s, open("charge_raw.txt", "w") as f_c:
        for s_line, c_line in results:
            f_s.write(s_line + "\n")
            f_c.write(c_line + "\n")
    print("✅ Generated surface_raw.txt and charge_raw.txt using direct PQR data.")