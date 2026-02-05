import math
import numpy as np
import argparse
import os

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
                continue  # skip malformed lines
            atom_name = parts[0].strip()
            res_name = parts[1].strip()
            charge = float(parts[2])
            radius = float(parts[3])
            key = (res_name.upper(), atom_name.upper())
            ref_map[key] = (charge, radius)
    return ref_map

# Normalize vector and scale by t1
def new_vector(V, t1):
    V_unit = V / np.linalg.norm(V)
    return t1 * V_unit

# Main cavity calculation
def cavity(name_file, surface_coords, radius_sphere, ref_map):
    radius_limit = radius_sphere + 2
    protein_coords = []
    atom_radii = []
    charges = []
    original_keys = []

    # Read and filter atoms based on reference file
    with open(name_file, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            atom_name = parts[2].strip()
            res_name = parts[3].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:53])

            # Only keep atoms that exist in the reference file
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
        print(f"⚠️ No hits found in: {name_file}")
        return [], [], [], [], []

    protein_coords = np.array(protein_coords)
    atom_radii = np.array(atom_radii)

    surface_vectors = []
    distance_vectors = []
    hit_keys = []

    for vector in surface_coords:
        V = np.array(vector)
        min_dist = float('inf')
        final_vector = V
        final_index = -1

        for i, (S, r) in enumerate(zip(protein_coords, atom_radii)):
            d = np.linalg.norm(S)
            cos_alpha = np.dot(S, V) / (d * np.linalg.norm(V))

            if np.array_equal(S, V):
                t1 = d - r
                temp_vector = new_vector(V, t1)
            elif 0 <= cos_alpha <= 1:
                y = d * math.sqrt(1.0 - cos_alpha**2)
                if y < r:
                    x = math.sqrt(r * r - y * y)
                    t = abs(np.dot(S, V)) / np.linalg.norm(V)
                    t1 = t - x
                    temp_vector = new_vector(V, t1 if t1 < radius_sphere else radius_sphere)
                elif y == r:
                    t1 = np.dot(S, V) / np.linalg.norm(V)
                    temp_vector = new_vector(V, t1 if t1 < radius_sphere else radius_sphere)
                else:
                    continue
            else:
                continue

            if t1 < min_dist:
                final_vector = temp_vector
                min_dist = t1
                final_index = i

        if min_dist > radius_sphere:
            min_dist = radius_sphere

        distance_vectors.append(np.round(min_dist, 3))
        surface_vectors.append(final_vector)
        hit_keys.append(original_keys[final_index] if final_index != -1 else None)

    hit_charges = [ref_map[k][0] if k else 0.0 for k in hit_keys]
    hit_atom_names = [k[1] if k else "UNK" for k in hit_keys]
    hit_residue_names = [k[0] if k else "UNK" for k in hit_keys]

    return np.array(distance_vectors), np.array(surface_vectors), hit_charges, hit_atom_names, hit_residue_names


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate vector distances and partial charges.")
    parser.add_argument('-n', '--name', nargs='+', type=str, required=True, help="Protein PDB file(s)")
    parser.add_argument('-pdb', '--pdb', type=str, required=True, help="PDB file with surface coordinates")
    parser.add_argument('-r', '--radius', type=int, required=True, help="Max sphere radius")
    parser.add_argument('--ref', type=str, required=True, help="Reference file with radii and charges")
    parser.add_argument('-c', '--charge_output', type=str, required=True, help="Output file for hit charges")

    args = parser.parse_args()

    surface_coords = read_pdb_coords(args.pdb)
    ref_map = read_reference_file(args.ref)

    with open(args.charge_output, 'w') as charge_file:
        for protein_file in args.name:
            base_name = os.path.splitext(os.path.basename(protein_file))[0]
            dist_vals, surface_vectors, hit_charges, hit_atom_names, hit_residue_names = cavity(
                protein_file, surface_coords, args.radius, ref_map
            )
            charge_file.write(base_name + " " + ' '.join(map(str, np.round(hit_charges, 4))) + "\n")

