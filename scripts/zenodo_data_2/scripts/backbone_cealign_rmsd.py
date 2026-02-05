#!/usr/bin/env python3
"""
parallel_pymol_rmsd_pdb.py

- Finds all PDB files in current directory
- Uses filename (without .pdb) as label
- Runs parallel PyMOL backbone RMSD calculations (align) on all pairs
- Outputs full RMSD matrix CSV

Usage:
  python parallel_pymol_rmsd_pdb.py --out_csv rmsd_matrix.csv --workers 32 --concurrency 8
"""

import os
import sys
import math
import argparse
import tempfile
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import cpu_count
from time import sleep

def gather_pdbs_in_folder(folder):
    pdb_files = [os.path.abspath(os.path.join(folder, f))
                 for f in os.listdir(folder) if f.lower().endswith('.pdb')]
    if not pdb_files:
        raise RuntimeError(f"No PDB files found in {folder}")
    labels = [os.path.splitext(os.path.basename(p))[0] for p in pdb_files]
    return pdb_files, labels

def build_pair_list(n):
    pairs = []
    for i in range(n):
        for j in range(i, n):
            pairs.append((i, j))
    return pairs

def chunk_pairs(pairs, workers):
    total = len(pairs)
    per = math.ceil(total / workers)
    chunks = [pairs[i:i+per] for i in range(0, total, per)]
    return chunks

def write_worker_script(tmpdir, idx, needed_files, labels_map, pairs, out_chunk_path, pdb_files):
    load_lines = []
    for pdb in needed_files:
        label = labels_map[pdb]
        pdb_escaped = pdb.replace("\\", "\\\\")
        load_lines.append(f'cmd.load(r"{pdb_escaped}", "{label}")\n')

    pairs_list_py = "[" + ", ".join([f'("{labels_map[pdb_files[i]]}", "{labels_map[pdb_files[j]]}")' 
                                     for (i, j) in pairs]) + "]"

    worker_code = f'''
from pymol import cmd, finish_launching
print("Worker {idx} PyMOL script started with offscreen mode", flush=True)
finish_launching(['pymol', '-cq'])

{"".join(load_lines)}

out_path = r"{out_chunk_path}"

with open(out_path, "w") as fout:
    for (lab_a, lab_b) in {pairs_list_py}:
        sel_a = f"{{lab_a}} and name N+CA+C+O"
        sel_b = f"{{lab_b}} and name N+CA+C+O"

        result = cmd.cealign(sel_a, sel_b)

        if result and "RMSD" in result:
            rms = result["RMSD"]
        else:
            rms = float("nan")

        fout.write(f"{{lab_a}},{{lab_b}},{{rms:.4f}}\\n")

cmd.quit()
'''
    script_path = os.path.join(tmpdir, f"worker_{idx:04d}.py")
    with open(script_path, "w") as fh:
        fh.write(worker_code)
    return script_path

def merge_chunks(chunk_files, labels, out_csv):
    n = len(labels)
    idx_map = {lab: i for i, lab in enumerate(labels)}
    mat = np.zeros((n, n), dtype=float)

    for cf in sorted(chunk_files):
        with open(cf, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                a, b, rms = line.split(",")
                i = idx_map[a]
                j = idx_map[b]
                val = float(rms)
                mat[i, j] = val
                mat[j, i] = val

    df = pd.DataFrame(mat, index=labels, columns=labels)
    df.to_csv(out_csv)
    print(f"Saved merged RMSD matrix to {out_csv}", flush=True)

def main():
    parser = argparse.ArgumentParser(description="Parallel PyMOL RMSD on all PDB files in current folder")
    parser.add_argument("--out_csv", default="rmsd_matrix_pymol.csv", help="Output CSV RMSD matrix file")
    parser.add_argument("--workers", type=int, default=32, help="Number of worker scripts (chunks)")
    parser.add_argument("--concurrency", type=int, default=min(8, cpu_count()), help="Number of pymol processes to run simultaneously")
    parser.add_argument("--folder", default=".", help="Folder containing PDB files")
    args = parser.parse_args()

    import shutil
    pymol_python = shutil.which("pymol")
    if pymol_python is None:
        raise RuntimeError("Could not find 'pymol' executable in PATH. Please ensure PyMOL is installed and accessible.")
    print(f"Using pymol executable at: {pymol_python}", flush=True)

    tmpdir = tempfile.mkdtemp(prefix="pymol_rmsd_")
    os.makedirs(tmpdir, exist_ok=True)
    print(f"Temporary working directory is: {tmpdir}", flush=True)

    pdb_files, labels = gather_pdbs_in_folder(args.folder)
    n = len(pdb_files)
    print(f"Found {n} PDB files: {labels}", flush=True)

    labels_map = {p: os.path.splitext(os.path.basename(p))[0] for p in pdb_files}

    pairs = build_pair_list(n)
    chunks = chunk_pairs(pairs, args.workers)
    print(f"Total pairs: {len(pairs)}, split into {len(chunks)} chunks", flush=True)

    worker_scripts = []
    chunk_out_files = []

    for idx, chunk in enumerate(chunks):
        needed_idx = set()
        for (i, j) in chunk:
            needed_idx.add(i)
            needed_idx.add(j)
        needed_files = [pdb_files[k] for k in needed_idx]
        out_chunk = os.path.join(tmpdir, f"chunk_{idx:04d}.csv")
        script_path = write_worker_script(tmpdir, idx, needed_files, labels_map, chunk, out_chunk, pdb_files)
        worker_scripts.append(script_path)
        chunk_out_files.append(out_chunk)

    running = []
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["PYMOL_NO_MAIN"] = "1"

    print(f"Launching up to {args.concurrency} PyMOL worker processes in offscreen mode.")
    for idx, script_path in enumerate(worker_scripts):
        while len(running) >= args.concurrency:
            for i, proc in enumerate(running):
                ret = proc.poll()
                if ret is not None:
                    stdout, stderr = proc.communicate(timeout=1)
                    if ret != 0:
                        print(f"[ERROR] Worker {i} exited with code {ret}")
                        print(f"----- stdout -----\n{stdout}\n----- stderr -----\n{stderr}\n------------------")
                    else:
                        print(f"[INFO] Worker {i} finished successfully")
                    running.pop(i)
                    break
            else:
                sleep(0.1)

        cmdline = [pymol_python, "-cq", script_path]
        print(f"Starting worker {idx} with command: {' '.join(cmdline)}", flush=True)
        proc = subprocess.Popen(cmdline,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                env=env,
                                universal_newlines=True,
                                bufsize=1)
        running.append(proc)

    while running:
        for i, proc in enumerate(running):
            ret = proc.poll()
            if ret is not None:
                stdout, stderr = proc.communicate(timeout=1)
                if ret != 0:
                    print(f"[ERROR] Worker {i} exited with code {ret}")
                    print(f"----- stdout -----\n{stdout}\n----- stderr -----\n{stderr}\n------------------")
                else:
                    print(f"[INFO] Worker {i} finished successfully")
                running.pop(i)
                break
        else:
            sleep(0.1)

    print("All workers completed. Merging chunk CSV files...")
    merge_chunks(chunk_out_files, labels, args.out_csv)

    print(f"Temporary directory used: {tmpdir}")
    print("Done.")

if __name__ == "__main__":
    main()

