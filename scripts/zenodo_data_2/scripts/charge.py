import math
import numpy as np
import argparse
import os
import glob # Added for wildcard handling

# Read coordinates from PDB
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

# Read new reference file format
def read_reference_file(ref_file):
    ref_map = {}
    with open(ref_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) != 4:
                continue
            atom_name, res_name = parts[0].strip(), parts[1].strip()
            charge, radius = float(parts[2]), float(parts[3])
            key = (res_name.upper(), atom_name.upper())
            ref_map[key] = (charge, radius)
    return ref_map

def new_vector(V, t1):
    V_unit = V / np.linalg.norm(V)
    return t1 * V_unit

# Main cavity calculation
def cavity(name_file, surface_coords, radius_sphere, ref_map):
    radius_limit = radius_sphere + 2
    protein_coords, atom_radii, charges, original_keys = [], [], [], []

    with open(name_file, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name, res_name = line[12:16].strip(), line[17:21].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:53])

            ref_key = (res_name.upper(), atom_name.upper())
            if ref_key not in ref_map:
                continue

            d = math.sqrt(x * x + y * y + z * z)
            if d < radius_limit:
                charge, radius = ref_map[ref_key]
                protein_coords.append([x, y, z])
                atom_radii.append(radius)
                charges.append(charge)
                original_keys.append(ref_key)

    if not protein_coords:
        return None

    protein_coords, atom_radii = np.array(protein_coords), np.array(atom_radii)
    distance_vectors, hit_keys = [], []

    for vector in surface_coords:
        V = np.array(vector)
        min_dist, final_index = float('inf'), -1

        for i, (S, r) in enumerate(zip(protein_coords, atom_radii)):
            d = np.linalg.norm(S)
            cos_alpha = np.dot(S, V) / (d * np.linalg.norm(V))

            if np.array_equal(S, V):
                t1 = d - r
            elif 0 <= cos_alpha <= 1:
                y = d * math.sqrt(1.0 - cos_alpha**2)
                if y < r:
                    x = math.sqrt(r * r - y * y)
                    t1 = (abs(np.dot(S, V)) / np.linalg.norm(V)) - x
                elif y == r:
                    t1 = np.dot(S, V) / np.linalg.norm(V)
                else: continue
            else: continue

            if t1 < min_dist:
                min_dist, final_index = t1, i

        distance_vectors.append(np.round(min(min_dist, radius_sphere), 3))
        hit_keys.append(original_keys[final_index] if final_index != -1 else None)

    return [ref_map[k][0] if k else 0.0 for k in hit_keys]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate partial charge vectors.")
    parser.add_argument('-n', '--name', nargs='+', type=str, required=True)
    parser.add_argument('-pdb', '--pdb', type=str, required=True)
    parser.add_argument('-r', '--radius', type=int, required=True)
    parser.add_argument('--ref', type=str, required=True)
    parser.add_argument('-c', '--charge_output', type=str, required=True)
    args = parser.parse_args()

    # Expand wildcards
    all_files = []
    for pattern in args.name:
        all_files.extend(glob.glob(pattern))
    
    if not all_files:
        print(f"❌ Error: No files found matching: {args.name}")
        exit(1)

    print(f"🚀 Found {len(all_files)} files. Loading lattice and reference map...")
    surface_coords = read_pdb_coords(args.pdb)
    ref_map = read_reference_file(args.ref)

    with open(args.charge_output, 'w') as charge_file:
        for idx, protein_file in enumerate(all_files):
            base_name = os.path.splitext(os.path.basename(protein_file))[0]
            
            # Progress print
            print(f"   [{idx+1}/{len(all_files)}] Processing: {base_name}...", end="\r")
            
            hit_charges = cavity(protein_file, surface_coords, args.radius, ref_map)
            
            if hit_charges is not None:
                charge_file.write(base_name + " " + ' '.join(map(str, np.round(hit_charges, 4))) + "\n")
            else:
                print(f"\n   ⚠️ Warning: No hits found for {base_name}")

    print(f"\n✅ Success! Charge vectors saved to: {args.charge_output}")