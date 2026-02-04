import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

# --- REPLICATION PARAMETERS ---
MAX_DIST = 20.0  # Å
NUM_SPHERE_POINTS = 492
VDW_RADII = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2} # Standard vdW
# Mock GROMOS 54a8 partial charges (In a real run, use a full mapping file)
GROMOS_CHARGES = {'N': -0.28, 'CA': 0.0, 'C': 0.38, 'O': -0.38} 

def generate_fibonacci_sphere(n=492):
    """Generates uniform points on a sphere and keeps the distal hemisphere."""
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
    for i in range(n):
        y = 1 - (i / float(n - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        # Kuvek Hemispheric Restriction: Keep only vectors pointing 'up' (+z)
        if z > 0:
            points.append(np.array([x, y, z]))
    return points

def get_intersection(vector, atoms):
    """Finds the first intersection of a ray with atom vdW spheres."""
    best_d = MAX_DIST
    hit_charge = 0.0
    
    for atom in atoms:
        # HEME EXCLUSION: Ignore heme atoms to see the pocket walls
        if atom.get_parent().get_resname() in ['HEM', 'HEC']: continue
        
        P = atom.get_coord()
        r = VDW_RADII.get(atom.element, 1.7)
        
        # Ray-Sphere Intersection math: |d*V - P|^2 = r^2
        v_dot_p = np.dot(vector, P)
        p_dot_p = np.dot(P, P)
        discriminant = (v_dot_p**2) - (p_dot_p - r**2)
        
        if discriminant >= 0:
            d = v_dot_p - np.sqrt(discriminant)
            if 0 < d < best_d:
                best_d = d
                # Assign GROMOS charge (Simplified lookup)
                hit_charge = GROMOS_CHARGES.get(atom.get_name(), 0.0)
                
    return best_d, hit_charge

def process_aligned_pdbs(input_dir):
    vectors = generate_fibonacci_sphere(NUM_SPHERE_POINTS)
    parser = PDBParser(QUIET=True)
    all_data = []

    for filename in os.listdir(input_dir):
        if not filename.endswith(".pdb"): continue
        
        struct_id = filename.split('_')[0]
        structure = parser.get_structure(struct_id, os.path.join(input_dir, filename))
        atoms = list(structure.get_atoms())
        
        lengths = []
        charges = []
        
        for v in vectors:
            d, c = get_intersection(v, atoms)
            lengths.append(d)
            charges.append(c)
        
        # Combine lengths and charges into one flat vector [L1...L260, C1...C260]
        fingerprint = lengths + charges
        all_data.append([struct_id] + fingerprint)
        print(f"Generated BSV for {struct_id}")

    # Save to Master BSV Matrix
    cols = ['UniProt'] + [f'L_{i}' for i in range(260)] + [f'C_{i}' for i in range(260)]
    df = pd.DataFrame(all_data, columns=cols)
    df.to_csv("tacca_bsv_matrix.csv", index=False)

if __name__ == "__main__":
    process_aligned_pdbs("aligned_pdbs")