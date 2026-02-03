#!/usr/bin/env python3
"""
01_fetch_seq_pdb.py

Download:
  • UniProt FASTA sequences
  • AlphaFold PDB structures

Input:
  plant_cyps_table_S2.csv  (columns: CYP,UniProt code)

Output:
  data/sequences/<CYP>.fasta
  data/structures/<CYP>.pdb

Usage:
  python 01_fetch_seq_pdb.py plant_cyps_table_S2.csv
"""

import os
import sys
import time
import requests
import pandas as pd

SEQ_DIR = "data/sequences"
PDB_DIR = "data/structures"

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{id}.fasta"
AF_PDB = "https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v4.pdb"


# ---------- utils ----------
def mkdir(p):
    os.makedirs(p, exist_ok=True)


def download(url, outpath):
    if os.path.exists(outpath):
        return "exists"

    r = requests.get(url, timeout=60)
    if r.status_code != 200:
        return "fail"

    with open(outpath, "wb") as f:
        f.write(r.content)
    return "ok"


# ---------- main ----------
def main(csv_path):

    df = pd.read_csv(csv_path)

    mkdir(SEQ_DIR)
    mkdir(PDB_DIR)

    total = len(df)
    print(f"\nFetching {total} plant CYPs\n")

    ok_seq = 0
    ok_pdb = 0

    for i, row in df.iterrows():
        cyp = row["CYP"]
        uid = row["UniProt code"]

        fasta_out = f"{SEQ_DIR}/{cyp}.fasta"
        pdb_out = f"{PDB_DIR}/{cyp}.pdb"

        s = download(UNIPROT_FASTA.format(id=uid), fasta_out)
        p = download(AF_PDB.format(id=uid), pdb_out)

        if s == "ok" or s == "exists":
            ok_seq += 1
        if p == "ok" or p == "exists":
            ok_pdb += 1

        print(f"[{i+1:03d}/{total}] {cyp:10s}  seq:{s:7s}  pdb:{p:7s}")
        time.sleep(0.1)  # polite rate limit

    print("\nQC")
    print("----")
    print("Sequences:", ok_seq)
    print("Structures:", ok_pdb)
    print("Expected:", total)

    if ok_seq != total or ok_pdb != total:
        print("\n⚠ Some downloads failed (likely deprecated UniProt/AF entries)")
    else:
        print("\n✓ All files fetched successfully")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 01_fetch_seq_pdb.py plant_cyps_table_S2.csv")
        sys.exit(1)

    main(sys.argv[1])
