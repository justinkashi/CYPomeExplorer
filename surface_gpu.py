import cupy as cp
import numpy as np
import argparse
import os
import glob
import time

# --- PDB/PQR Loading Logic ---
def read_pdb_coords(pdb_filename):
    coords = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords, dtype=np.float64)

def load_pqr_atoms(pqr_filename, radius_sphere):
    radius_limit = radius_sphere + 2
    coords, radii = [], []
    with open(pqr_filename, "r") as f:
        for line in f:
            # Replicating exact surface.py filters
            if line[0:4] == "ATOM" and float(line[46:53]) > -2 and line[17:20] not in ["HEM"]:
                x, y, z, r = float(line[30:38]), float(line[38:46]), float(line[46:53]), float(line[69:75])
                if (x**2 + y**2 + z**2)**0.5 < radius_limit:
                    coords.append([x, y, z])
                    radii.append(r)
    return np.array(coords, dtype=np.float64), np.array(radii, dtype=np.float64)

# --- CUDA Kernel for Exact Parity ---
cuda_source = r'''
extern "C" __global__
void surface_kernel(const double* L, const double* P, const double* R, double* out, int nL, int nP, double limit) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nL) return;

    double Vx = L[i*3], Vy = L[i*3+1], Vz = L[i*3+2];
    double min_d = 1e18, t1 = 1e18;

    for (int j = 0; j < nP; j++) {
        double Sx = P[j*3], Sy = P[j*3+1], Sz = P[j*3+2], r = R[j];
        double d = round(sqrt(Sx*Sx + Sy*Sy + Sz*Sz) * 1000.0) / 1000.0;
        double norm_V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
        double r_norm_V = round(norm_V * 1000.0) / 1000.0;
        double cos_a = (Sx*Vx + Sy*Vy + Sz*Vz) / (d * r_norm_V);

        if (Sx == Vx && Sy == Vy && Sz == Vz) {
            t1 = d - r;
            Vx = (t1 * Vx) / norm_V; Vy = (t1 * Vy) / norm_V; Vz = (t1 * Vz) / norm_V;
        } else if (cos_a >= 0.0 && cos_a <= 1.0) {
            double y = round(d * sqrt(1.0 - cos_a * cos_a) * 1000.0) / 1000.0;
            if (y < r) {
                double t = round((abs(Sx*Vx + Sy*Vy + Sz*Vz) / norm_V) * 1000.0) / 1000.0;
                t1 = t - sqrt(r*r - y*y);
                double h = (t1 < limit) ? t1 : limit;
                if ((h * Vz / norm_V) >= 0.0) { 
                    Vx = (h * Vx) / norm_V; Vy = (h * Vy) / norm_V; Vz = (h * Vz) / norm_V;
                }
            } else if (y == r) {
                t1 = (Sx*Vx + Sy*Vy + Sz*Vz) / norm_V;
                double h = (t1 < limit) ? t1 : limit;
                Vx = (h * Vx) / norm_V; Vy = (h * Vy) / norm_V; Vz = (h * Vz) / norm_V;
            }
        }
        if (t1 < min_d) min_d = t1;
    }
    out[i] = round(fmin(min_d, limit) * 1000.0) / 1000.0;
}
'''
module = cp.RawModule(code=cuda_source)
kernel = module.get_function('surface_kernel')

def run_gpu(pqr_files, lattice_pdb, output_file, radius):
    L = cp.array(read_pdb_coords(lattice_pdb))
    with open(output_file, 'w') as f:
        for pqr in pqr_files:
            P_cpu, R_cpu = load_pqr_atoms(pqr, radius)
            if len(P_cpu) == 0: continue
            
            P, R, out = cp.array(P_cpu), cp.array(R_cpu), cp.zeros(len(L))
            grid = (int(np.ceil(len(L)/256)), 1, 1)
            kernel(grid, (256, 1, 1), (L, P, R, out, len(L), len(P), float(radius)))
            
            name = os.path.splitext(os.path.basename(pqr))[0]
            f.write(f"{name} {' '.join(map(str, out.get().tolist()))}\n")
            print(f"Finished: {name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-r', '--radius', type=int, required=True)
    args = parser.parse_args()
    run_gpu(args.name, args.pdb, args.output, args.radius)