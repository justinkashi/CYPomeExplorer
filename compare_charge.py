import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def load_results(filename):
    """Parses charge results into a dictionary {protein_name: np.array([charges])}"""
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts: continue
            data[parts[0]] = np.array([float(x) for x in parts[1:]])
    return data

def compare_results(ground_truth, parallel, gpu):
    common_keys = set(ground_truth.keys()) & set(parallel.keys()) & set(gpu.keys())
    print(f"✅ Comparing {len(common_keys)} common proteins...")

    # Lists to store mean absolute differences across the batch
    diff_parallel = []
    diff_gpu = []

    # Get the length of vectors from the first common key
    sample_key = list(common_keys)[0]
    num_points = len(ground_truth[sample_key])

    # Arrays to store point-by-point averages across the whole batch
    batch_parallel_err = np.zeros(num_points)
    batch_gpu_err = np.zeros(num_points)

    for key in common_keys:
        gt = ground_truth[key]
        pa = parallel[key]
        gp = gpu[key]

        # Calculate absolute difference for this protein
        err_pa = np.abs(gt - pa)
        err_gp = np.abs(gt - gp)

        batch_parallel_err += err_pa
        batch_gpu_err += err_gp

    # Average the error over the number of proteins
    batch_parallel_err /= len(common_keys)
    batch_gpu_err /= len(common_keys)

    return batch_parallel_err, batch_gpu_err

def plot_drift(err_pa, err_gp, output_image):
    plt.figure(figsize=(12, 6), dpi=150)
    
    indices = np.arange(len(err_pa))
    
    plt.plot(indices, err_pa, label='Avg Error: Parallel vs Ground Truth', color='blue', alpha=0.7)
    plt.plot(indices, err_gp, label='Avg Error: GPU vs Ground Truth', color='red', alpha=0.7)

    plt.yscale('log') # Use log scale because error should be tiny (1e-5 or smaller)
    plt.xlabel('Lattice Point Index')
    plt.ylabel('Mean Absolute Difference (Partial Charge)')
    plt.title('Batch Charge Comparison: Mathematical Parity Check')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_image)
    print(f"📊 Comparison plot saved to: {output_image}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--gt', required=True, help="Original charge_results.txt")
    parser.add_argument('--pa', required=True, help="Parallel charge_results.txt")
    parser.add_argument('--gp', required=True, help="GPU charge_results.txt")
    parser.add_argument('-o', '--output', default="charge_comparison.png")
    args = parser.parse_args()

    print("📖 Loading result files...")
    gt_data = load_results(args.gt)
    pa_data = load_results(args.pa)
    gp_data = load_results(args.gp)

    err_pa, err_gp = compare_results(gt_data, pa_data, gp_data)
    plot_drift(err_pa, err_gp, args.output)