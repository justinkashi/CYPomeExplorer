"""
USAGE: 
    python compare_batch.py --cpu surface_raw.txt --parallel surface_parallel.txt --gpu surface_gpu_results.txt
"""

import numpy as np
import argparse
import matplotlib.pyplot as plt

def load_data(filename):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            p = line.split()
            if p: data[p[0]] = np.array([float(x) for x in p[1:]])
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu", required=True)
    parser.add_argument("--parallel", required=True)
    parser.add_argument("--gpu", required=True)
    args = parser.parse_args()

    cpu, par, gpu = load_data(args.cpu), load_data(args.parallel), load_data(args.gpu)
    common = sorted(set(cpu.keys()) & set(par.keys()) & set(gpu.keys()))

    print(f"Comparing {len(common)} proteins...")

    # For the plot: Show average error across all proteins per lattice point
    all_diff_par = []
    all_diff_gpu = []

    for name in common:
        all_diff_par.append(np.abs(cpu[name] - par[name]))
        all_diff_gpu.append(np.abs(cpu[name] - gpu[name]))

    # Visualization
    plt.figure(figsize=(12, 6))
    plt.plot(np.mean(all_diff_par, axis=0), label="Avg Error: Parallel vs Ground Truth", color='blue', alpha=0.7)
    plt.plot(np.mean(all_diff_gpu, axis=0), label="Avg Error: GPU vs Ground Truth", color='red', alpha=0.7)
    
    plt.title("Batch Comparison: Mathematical Drift from Ground Truth")
    plt.xlabel("Lattice Point Index")
    plt.ylabel("Mean Absolute Difference (Å)")
    plt.legend()
    plt.yscale('log') # Log scale helps see tiny floating point differences
    plt.grid(True, which="both", ls="-", alpha=0.2)
    
    plt.savefig("batch_error_profile.png")
    print("Plot saved to batch_error_profile.png")
    plt.show()

if __name__ == "__main__":
    main()