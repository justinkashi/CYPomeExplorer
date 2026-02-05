import numpy as np
import argparse
import os

def icosahedron_vertices():
    """Generate the vertices of a regular icosahedron centered at the origin."""
    phi = (1 + np.sqrt(5)) / 2  # The golden ratio

    vertices = np.array([
        [-1, phi, 0], [1, phi, 0], [-1, -phi, 0], [1, -phi, 0],
        [0, -1, phi], [0, 1, phi], [0, -1, -phi], [0, 1, -phi],
        [phi, 0, -1], [phi, 0, 1], [-phi, 0, -1], [-phi, 0, 1]
    ])

    # Normalize vertices to lie on the unit sphere
    vertices /= np.linalg.norm(vertices[0])
    return vertices

def icosahedron_faces():
    """Generate the faces of a regular icosahedron."""
    return np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
    ])

def subdivide_triangle(v1, v2, v3, n):
    """Subdivide a triangle into smaller triangles projected onto a sphere."""
    points = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            u = i / n
            v = j / n
            w = 1 - u - v
            point = u * v1 + v * v2 + w * v3
            points.append(point / np.linalg.norm(point))
    return np.array(points)

def triangular_lattice_on_sphere(radius, subdivisions, hemisphere=False):
    """
    Generate points on a full sphere or hemisphere using a triangular lattice.

    Parameters:
        radius (float): Radius of the sphere.
        subdivisions (int): Number of subdivisions per triangle.
        hemisphere (bool): If True, keep only z >= 0 (upper hemisphere).
    """
    vertices = icosahedron_vertices()
    faces = icosahedron_faces()

    points = []
    for face in faces:
        v1, v2, v3 = vertices[face]
        subdivided_points = subdivide_triangle(v1, v2, v3, subdivisions)
        points.append(subdivided_points)

    points = np.concatenate(points)
    rounded_points = np.round(points, decimals=8)
    unique_points = np.unique(rounded_points, axis=0)

    # Filter for hemisphere if requested
    if hemisphere:
        unique_points = unique_points[unique_points[:, 2] >= 0]

    unique_points *= radius
    return unique_points

def write_pdb(coords, filename):
    """Write 3D points to a PDB file."""
    with open(filename, "w") as pdb_file:
        pdb_file.write("REMARK   Generated PDB of triangular lattice points on a sphere\n")
        for i, (x, y, z) in enumerate(coords, start=1):
            pdb_file.write(f"ATOM  {i:5d}  C   SPH A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        pdb_file.write("TER\nEND\n")

def remove_duplicate_atoms(input_file, output_file):
    """Remove duplicate atom coordinates and renumber atoms sequentially."""
    unique_coordinates = set()
    unique_lines = []

    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("ATOM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                x, y, z = (0.0 if abs(coord) < 1e-6 else coord for coord in (x, y, z))
                if (x, y, z) not in unique_coordinates:
                    unique_coordinates.add((x, y, z))
                    unique_lines.append(line)

    with open(output_file, 'w') as outfile:
        for i, line in enumerate(unique_lines, start=1):
            new_line = f"{line[:6]}{i:5d}{line[11:]}"
            outfile.write(new_line)

# =============================
# Main execution block
# =============================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a triangular lattice (sphere or hemisphere) and remove duplicate points."
    )
    parser.add_argument('-r', '--radius', type=float, required=True, help="Radius of the sphere")
    parser.add_argument('-s', '--subdivisions', type=int, required=True, help="Number of subdivisions per triangle")
    parser.add_argument('-o', '--output_name', type=str, required=True, help="Name of the final PDB output file")
    parser.add_argument('--hemisphere', action='store_true', default=False,
                        help="Generate only the upper hemisphere (z >= 0). Default: full sphere")

    args = parser.parse_args()

    temp_file = "temp_lattice.pdb"

    # Generate lattice
    coords = triangular_lattice_on_sphere(args.radius, args.subdivisions, hemisphere=args.hemisphere)
    write_pdb(coords, temp_file)

    # Remove duplicates
    remove_duplicate_atoms(temp_file, args.output_name)

    # Clean up
    if os.path.exists(temp_file):
        os.remove(temp_file)

    print(f"\n✅ PDB file '{args.output_name}' created successfully.")
    print(f"   Radius: {args.radius}")
    print(f"   Subdivisions: {args.subdivisions}")
    print(f"   Hemisphere mode: {args.hemisphere}")
    print(f"   Total points generated: {len(coords)}")

