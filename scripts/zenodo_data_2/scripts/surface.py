import math
import numpy as np
import argparse
import os

# Function to read coordinates from a PDB file
def read_pdb_coords(pdb_filename):
    coords = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract the x, y, z coordinates from the PDB line
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
    return np.array(coords)

# Function to process the protein structure and filter atoms
def cavity(name, surface_coords, radius_sphere):

    # Removing the atoms which are below the heme and outside the cutoff, radius is cutoff
    radius=radius_sphere+2
    with open(name, "r") as f:
        protein_coords = []
        atom_radius = []
        for line in f:
            atoms = []
            if line[0:4] == "ATOM" and float(line[46:53]) > -2 and line[17:20] not in ["HEM"]:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:53])
                r = float(line[69:75])
                d = math.sqrt(x * x + y * y + z * z)
                # if the sphere center is more away from iron than the cutoff we don't want these coordinates
                if d < radius:
                    atoms.append(x)
                    atoms.append(y)
                    atoms.append(z)
                    atom_radius.append(r)
                    protein_coords.append(atoms)
                    
        protein_coords = np.array(protein_coords)
        atom_radius = np.array(atom_radius)

    def new_vector(V, t1):
        # Normalize the vector V
        norm_V = np.linalg.norm(V)
        V_unit = V / norm_V
        
        # Scale the unit vector by distance t
        V_new = t1 * V_unit
        return V_new    

    # Vector generation for each connection point
    surface_vectors = []
    distance_vectors = []
    for vector in surface_coords:
        V = np.array(vector)
        min_dist = float('inf')
        t1 = float('inf')
        final_vector = V
        i = 0
        for atom, r in zip(protein_coords, atom_radius):
            i += 1
            sphere_touching = True
            S = np.array(atom)
            d = np.round(np.linalg.norm(S), 3)
            cos_alpha = np.dot(S, V) / (d * np.round(np.linalg.norm(V), 3))
            if np.array_equal(S, V):
                t1 = d - r
                V = new_vector(V, t1)
            elif cos_alpha >= 0 and cos_alpha <= 1:
                y = round(d * math.sqrt(1.0 - cos_alpha * cos_alpha), 3)
                if y < r:
                    x = math.sqrt(r * r - y * y)
                    t = np.round(abs(np.dot(S, V)) / np.linalg.norm(V), 3)
                    t1 = t - x
                    t2 = t + x
                    if (t1<radius_sphere):
                        temp_vector = new_vector(V, t1)
                    else:
                        temp_vector = new_vector(V, radius_sphere)

                    if temp_vector[2] >= 0:  # Ensure it stays in the upper hemisphere
                        V = temp_vector
                    if temp_vector[2] < 0:  # Ensure it stays in the upper hemisphere
                        # Debug: Print current vectors and parameters
                        print(f"Processing vector: {V}, Atom: {S}, Distance d: {d}")
                        print(f"Temp vector: {temp_vector}, Z: {temp_vector[2]}")
                        print(i, cos_alpha, t)
                elif y == r:
                    t1 = np.dot(S, V) / np.linalg.norm(V)
                    temp_vector = new_vector(V, t1)
                    if (t1<radius_sphere):
                        temp_vector = new_vector(V, t1)
                    else:
                        temp_vector = new_vector(V, radius_sphere)
                else:
                    sphere_touching = False
                
            if sphere_touching == True and t1 < min_dist:
                final_vector = temp_vector
                min_dist = t1
        if (min_dist > radius_sphere):
            min_dist=radius_sphere       
        distance_vectors.append(np.round(min_dist, 3)) #rounding the values
        surface_vectors.append(final_vector)
    
    distance_vectors = np.array(distance_vectors)
    for i in range(len(distance_vectors)):
        if distance_vectors[i] == float('inf'):
            distance_vectors[i] = radius_sphere
    surface_vectors = np.array(surface_vectors)
    '''with open("surface.pdb", "w") as file:
        atom_count = 1
        for line in surface_vectors:
            x, y, z = line  # Unpack the coordinates
            file.write(f"ATOM  {atom_count:5d}  H   XXX     1    {x:8.3f}{y:8.3f}{z:8.3f}\n")
            atom_count += 1
        file.write("END\n")'''

    return distance_vectors, surface_vectors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script")
    parser.add_argument('-n', '--name', nargs='+', type=str, required=True, help="Name of the pqr file")
    parser.add_argument('-pdb', '--pdb', type=str, required=True, help="PDB file with surface coordinates")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output name")
    parser.add_argument('-r','--radius', type=int, required=True, help="Max sphere radius")
    
    args = parser.parse_args()

    # Read the coordinates from the provided PDB file
    surface_coords = read_pdb_coords(args.pdb)
    
    with open(args.output, 'w') as output_file:
        for filename in args.name:
            distance_results, surface_vectors = cavity(filename, surface_coords,args.radius)
            file_name_without_ext = os.path.splitext(os.path.basename(filename))[0]
            result_with_filename = [file_name_without_ext] + distance_results.tolist()
            output_file.write(' '.join(map(str, result_with_filename)) + '\n')

