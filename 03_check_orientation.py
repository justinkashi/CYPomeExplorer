# check_orientation.py
import numpy as np

def check_pdb(filename):
    fe_coords = None
    n_coords = []
    with open(filename, 'r') as f:
        for line in f:
            if "HEM" in line or "HEMO" in line:
                if " FE " in line:
                    fe_coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if " NA " in line or " NB " in line or " NC " in line or " ND " in line:
                    n_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    
    print(f"File: {filename}")
    print(f"  FE Position: {fe_coords}")
    if fe_coords is not None and np.allclose(fe_coords, [0,0,0], atol=0.1):
        print("  ✅ Origin Check: Iron is at (0,0,0).")
    else:
        print("  ❌ Origin Check: Iron is NOT at (0,0,0). You must use Trans-Rot.")

    if n_coords:
        z_vals = [c[2] for c in n_coords]
        if np.all(np.abs(z_vals) < 0.5):
            print("  ✅ Plane Check: Heme is in the XY-plane.")
        else:
            print(f"  ❌ Plane Check: Heme is tilted (Z-vals: {z_vals}). You must use Trans-Rot.")

if __name__ == "__main__":
    import sys
    check_pdb(sys.argv[1])
