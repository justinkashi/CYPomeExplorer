#!/usr/bin/env python3
"""
Download UniProt FASTAs + strict QC

Usage:
  python 01_fetch_sequences.py plant_cyps_table_S2.csv
"""

import os
import sys
import requests
import pandas as pd

OUTDIR = "data/sequences"
URL = "https://rest.uniprot.org/uniprotkb/{id}.fasta"

MIN_LEN = 350   # CYPs are ~480–520 aa; anything <350 is clearly broken


def fetch(uid, outpath):
    r = requests.get(URL.format(id=uid), timeout=60)
    if r.status_code != 200:
        return False

    text = r.text.strip()
    if not text.startswith(">"):
        return False

    seq = "".join(text.split("\n")[1:])

    if len(seq) < MIN_LEN:
        return False

    with open(outpath, "w") as f:
        f.write(text + "\n")

    return True


def main(csv):

    os.makedirs(OUTDIR, exist_ok=True)
    df = pd.read_csv(csv)

    ok = 0
    bad = []

    for _, r in df.iterrows():
        cyp = r["CYP"]
        uid = r["UniProt code"]

        out = f"{OUTDIR}/{cyp}.fasta"

        success = fetch(uid, out)

        if success:
            ok += 1
            print(f"{cyp:10s} ✓")
        else:
            bad.append((cyp, uid))
            print(f"{cyp:10s} ✗")

    print("\nQC SUMMARY")
    print("----------")
    print("valid sequences:", ok)
    print("expected:", len(df))

    if bad:
        print("\nFailures:")
        for x in bad:
            print(x)

    if ok != len(df):
        raise RuntimeError("Sequence QC failed")


if __name__ == "__main__":
    main(sys.argv[1])
