import pymol
import numpy as np
import argparse

def tranlate_to_iron(resn):
    pymol.cmd.select('select', f'resn {resn} and name FE')  
    fe = pymol.cmd.get_coords('select')
    translation_vector = [-fe[0][0], -fe[0][1], -fe[0][2]]
    pymol.cmd.translate(translation_vector, 'all')
    #pymol.cmd.save('fe.pdb')

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def norm_cross_product(v1, v2):
    return list(unit_vector(np.cross(v1, v2)))

def rad_to_degree(rad_angle):
    return np.round(np.rad2deg(rad_angle), 6)

def get_rotation_angle_and_axis(resn, axis):
    pymol.cmd.select('select', f'resn {resn} and name NA')
    v1 = tuple(pymol.cmd.get_coords('select')[0].tolist())
    pymol.cmd.select('select', f'resn {resn} and name NB')
    v2 = tuple(pymol.cmd.get_coords('select')[0].tolist())
      
    v = (0., 1., 0.) if axis == "xz" else (0., 0., 1.) if axis == "xy" else (1., 0., 0.)
    norm = norm_cross_product(v1, v2)
    angle_rad = angle_between(norm, v)
    angle_degree = rad_to_degree(angle_rad)
    vector_to_rotate_around = norm_cross_product(norm, v)
    return angle_degree, vector_to_rotate_around

def rotate_around_iron(resn, axis):
    angle_degree, vector_to_rotate_around = get_rotation_angle_and_axis(resn, axis)
    if angle_degree != "Done":
        pymol.cmd.rotate(vector_to_rotate_around, angle=angle_degree, selection="all", origin=[0, 0, 0])
        #pymol.cmd.save('fe_rotate.pdb')

def additional_rotation(resn, output_name):
    pymol.cmd.select('select', f'resn {resn} and name FE')
    fe = pymol.cmd.get_coords('select')[0]
    pymol.cmd.select('select', f'resn {resn} and name NA')
    na = pymol.cmd.get_coords('select')[0]
    
    fe_na_vector = np.array(na) - np.array(fe)
    reference_vector = np.array([1, 0, 0])
    
    angle_rad = angle_between(fe_na_vector, reference_vector)
    angle_degree = rad_to_degree(angle_rad)
    rotation_axis = norm_cross_product(fe_na_vector, reference_vector)
    
    pymol.cmd.rotate(rotation_axis, angle=angle_degree, selection="all", origin=[0, 0, 0])
    pymol.cmd.save(output_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align and rotate molecular structures.")
    parser.add_argument('-p', '--pdb', type=str, required=True, help="Path to pdb file")
    parser.add_argument('-a', '--axis', type=str, default="xy", required=False, help="Axis for angle calculation.")
    parser.add_argument('-n', '--name', type=str, default="HEMO", required=False, help="Residue name [upper case]")
    parser.add_argument('-o', '--output', type=str, default="fe_rotate.pdb", required=False, help="Output file name")
    args = parser.parse_args()
    
    pymol.cmd.load(args.pdb, 'ime')
    tranlate_to_iron(args.name)
    if get_rotation_angle_and_axis(args.name, args.axis) != "Done":
        rotate_around_iron(args.name, args.axis)
    additional_rotation(args.name, args.output)

