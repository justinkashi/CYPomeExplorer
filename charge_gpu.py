import cupy as cp
import numpy as np
import argparse
import os
import glob
import time

# --- CUDA Kernel: Replicates charge.py static ray logic exactly ---
cuda_source = r'''
extern "C" __global__
void charge_kernel(
    const double* L,    // Lattice points [nL * 3]
    const double* P,    // Protein coords [nP * 3]
    const double* R,    // Atom radii [nP]
    const double* Q,    // Atom charges [nP]
    double* out,        // Result charges [nL]
    int nL, int nP, double radius_sphere) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nL) return;

    double Vx = L[i*3], Vy = L[i*3+1], Vz = L[i*3+2];
    double norm_V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
    
    double min_t1 = 1e18;
    double final_charge = 0.0;

    for (int j = 0; j < nP; j++) {
        double Sx = P[j*3], Sy = P[j*3+1], Sz = P[j*3+2];
        double r = R[j];
        double q = Q[j];

        double d = sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
        double dot_SV = Sx*Vx + Sy*Vy + Sz*Vz;
        double cos_a = dot_SV / (d * norm_V);

        double t1 = 1e18;

        if (Sx == Vx && Sy == Vy && Sz == Vz) {
            t1 = d - r;
        } else if (cos_a >= 0.0 && cos_a <= 1.0) {
            double y = d * sqrt(1.0 - cos_a * cos_a);
            if (y < r) {
                double x = sqrt(r * r - y * y);
                double t = abs(dot_SV) / norm_V;
                t1 = t - x;
            } else if (y == r) {
                t1 = dot_SV / norm_V;
            }
        }

        if (t1 < min_t1) {
            min_t1 = t1;
            final_charge = q;
        }
    }

    // Radius cutoff and assignment
    if (min_t1 > radius_sphere) {
        out[i] = 0.0;
    } else {
        out[i] = final_charge;
    }
}
'''

def read_reference(ref_file):
    ref_map = {}
    with open(ref_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                ref_map[(parts[1].upper(), parts[0].upper())] = (float(parts[2]), float(parts[3]))
    return ref_map

def load_pqr(pqr_file, ref_map, radius_sphere):
    limit = radius_sphere + 2
    coords, radii, charges = [], [], []
    with open(pqr_file, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            # Exact column parsing from charge.py
            atom_name, res_name = line[12:16].strip(), line[17:21].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:53])
            
            key = (res_name.upper(), atom_name.upper())
            if key in ref_map:
                if (x**2 + y**2 + z**2)**0.5 < limit:
                    q, r = ref_map[key]
                    coords.append([x, y, z]); radii.append(r); charges.append(q)
    return np.array(coords), np.array(radii), np.array(charges)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', nargs='+', required=True)
    parser.add_argument('-pdb', '--pdb', required=True)
    parser.add_argument('-r', '--radius', type=int, required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('-c', '--output', required=True)
    args = parser.parse_args()

    # Expand Wildcards
    all_files = [f for pattern in args.name for f in glob.glob(pattern)]
    ref_map = read_reference(args.ref)
    
    # Load and move Lattice to GPU
    def read_lattice(path):
        coords = []
        with open(path, 'r') as f:
            for l in f:
                if l.startswith("ATOM"): coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        return np.array(coords)

    L_gpu = cp.array(read_lattice(args.pdb), dtype=cp.float64)
    nL = len(L_gpu)

    # Compile Kernel
    mod = cp.RawModule(code=cuda_source)
    ker = mod.get_function('charge_kernel')

    print(f"🚀 Processing {len(all_files)} files on GPU...")
    with open(args.output, 'w') as out_f:
        for pqr in all_files:
            P_cpu, R_cpu, Q_cpu = load_pqr(pqr, ref_map, args.radius)
            if len(P_cpu) == 0: continue
            
            P_gpu, R_gpu, Q_gpu = cp.array(P_cpu), cp.array(R_cpu), cp.array(Q_cpu)
            res_gpu = cp.zeros(nL, dtype=cp.float64)
            
            grid = (int(np.ceil(nL/256)), 1, 1)
            ker(grid, (256, 1, 1), (L_gpu, P_gpu, R_gpu, Q_gpu, res_gpu, nL, len(P_gpu), float(args.radius)))
            
            name = os.path.splitext(os.path.basename(pqr))[0]
            out_f.write(f"{name} {' '.join(map(str, np.round(res_gpu.get(), 4).tolist()))}\n")
            print(f"   Finished: {name}")

if __name__ == "__main__":
    main()