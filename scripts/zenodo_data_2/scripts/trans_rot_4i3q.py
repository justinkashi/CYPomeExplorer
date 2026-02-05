import pymol
import numpy as np
import argparse
import sys

# Standardize PyMOL initialization for headless execution
pymol.finish_launching(['pymol', '-cq'])

def tranlate_to_iron(resn):
    print(f"[1/3] Translating {resn} Iron to origin...")
    pymol.cmd.select('fe_sel', f'resn {resn} and name FE')
    if pymol.cmd.count_atoms('fe_sel') == 0:
        print(f"ERROR: Could not find FE atom in residue {resn}")
        sys.exit(1)
        
    fe = pymol.cmd.get_coords('fe_sel')[0]
    # Explicitly cast to standard Python floats to avoid API deadlock
    translation_vector = [-float(fe[0]), -float(fe[1]), -float(fe[2])]
    
    print(f"      Targeting origin from current FE: {fe}")
    pymol.cmd.translate(translation_vector, 'all')
    print("      ✅ Translation finished.")

def unit_vector(vector):
    v = np.array(vector)
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else v

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def norm_cross_product(v1, v2):
    return unit_vector(np.cross(v1, v2))

def rad_to_degree(rad_angle):
    return float(np.rad2deg(rad_angle))

def rotate_around_iron(resn, axis):
    print(f"[2/3] Rotating heme into {axis.upper()} plane...")
    
    # Extract Nitrogen coordinates for plane definition
    pymol.cmd.select('na_sel', f'resn {resn} and name NA')
    pymol.cmd.select('nb_sel', f'resn {resn} and name NB')
    
    if pymol.cmd.count_atoms('na_sel') == 0 or pymol.cmd.count_atoms('nb_sel') == 0:
        print(f"ERROR: NA or NB atoms not found in {resn}")
        sys.exit(1)

    v1 = np.array(pymol.cmd.get_coords('na_sel')[0])
    v2 = np.array(pymol.cmd.get_coords('nb_sel')[0])
    
    # Calculate Heme normal
    current_norm = unit_vector(np.cross(v1, v2))
    target_norm = np.array([0., 0., 1.]) if axis == "xy" else \
                  np.array([0., 1., 0.]) if axis == "xz" else \
                  np.array([1., 0., 0.])
    
    angle_deg = rad_to_degree(angle_between(current_norm, target_norm))
    rot_axis = norm_cross_product(current_norm, target_norm)
    
    # CRITICAL: Convert everything to pure Python types before calling PyMOL
    rot_axis_list = [float(x) for x in rot_axis]
    angle_float = float(angle_deg)
    
    if angle_float > 0.001:
        print(f"      HEARTBEAT: Sending rotation order to PyMOL engine...")
        print(f"      Parameters: {angle_float:.4f} degrees around {rot_axis_list}")
        
        # Suspend updates to prevent engine hang during transformation
        pymol.cmd.set('suspend_updates', 'on')
        pymol.cmd.rotate(rot_axis_list, angle=angle_float, selection="all", origin=[0, 0, 0])
        pymol.cmd.set('suspend_updates', 'off')
        
        print("      ✅ Rotation finished.")
    else:
        print("      Plane is already aligned. Skipping rotation.")

def additional_rotation(resn, output_name):
    print(f"[3/3] Finalizing orientation (Fe-NA to X-axis)...")
    
    na_coords = np.array(pymol.cmd.get_coords(f'resn {resn} and name NA')[0])
    # Iron is at 0,0,0
    fe_na_vector = na_coords 
    reference_vector = np.array([1, 0, 0])
    
    angle_deg = rad_to_degree(angle_between(fe_na_vector, reference_vector))
    
    if angle_deg > 0.001:
        rot_axis = [float(x) for x in norm_cross_product(fe_na_vector, reference_vector)]
        pymol.cmd.rotate(rot_axis, angle=float(angle_deg), selection="all", origin=[0, 0, 0])
    
    print(f"      ✅ Orientation locked. Saving to: {output_name}")
    pymol.cmd.save(output_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pdb', required=True)
    parser.add_argument('-a', '--axis', default="xy")
    parser.add_argument('-n', '--name', default="HEM")
    parser.add_argument('-o', '--output', default="4i3q_std.pdb")
    args = parser.parse_args()
    
    print(f"\n--- REPRODUCTION PIPELINE: {args.pdb} ---")
    pymol.cmd.load(args.pdb, 'target')
    
    tranlate_to_iron(args.name)
    rotate_around_iron(args.name, args.axis)
    additional_rotation(args.name, args.output)
    
    # FINAL QC PRINT
    final_fe = pymol.cmd.get_coords(f'resn {args.name} and name FE')[0]
    print(f"\n--- FINAL COORDINATE CHECK ---")
    print(f"FE Position: {final_fe}")
    print(f"Targeting:   [0.0, 0.0, 0.0]")
    
    pymol.cmd.quit()