import numpy as np

# === CONFIGURATION ===
input_file = "../charge.txt"    # <-- change this
output_file = "charge.txt"  # <-- change this

# === LOAD FILE & COMPUTE STD ===
all_values = []

with open(input_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        values = list(map(float, parts[1:]))  # skip snapshot ID
        all_values.extend(values)

all_values_np = np.array(all_values)
global_std = np.std(all_values_np)

print(f"File: {input_file}")
print(f"Global std: {global_std:.4f}")

# === NORMALIZE (values / global_std) & SAVE ===
with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    for line in fin:
        parts = line.strip().split()
        snapshot_id = parts[0]
        values = np.array(list(map(float, parts[1:])))
        norm_values = values / global_std
        norm_str = ' '.join(f"{v:.6f}" for v in norm_values)
        fout.write(f"{snapshot_id} {norm_str}\n")

print(f"\n✅ Normalization complete. Output saved to: {output_file}")

