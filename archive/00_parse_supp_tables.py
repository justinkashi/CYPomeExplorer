#!/usr/bin/env python3
# pip install "camelot-py[cv]" 
# brew install ghostscript

# usage: python 00_parse_supp_tables.py kuvek_supp_2026.pdf .

# user specified var: EXPECTED = 343

import sys
from pathlib import Path
import camelot
import pandas as pd
import re

EXPECTED = 343


def parse_table(df):
    rows = []

    for _, r in df.iterrows():
        cells = [str(c).strip() for c in r if str(c).strip()]

        # each row contains repeating: CYP UniProt CYP UniProt CYP UniProt
        for i in range(0, len(cells) - 1, 2):
            cyp = cells[i]
            uni = cells[i + 1]

            if re.match(r"^\d+[A-Z]+\d+$", cyp) and re.match(r"^[A-Z0-9]{6,10}$", uni):
                rows.append((cyp, uni))

    df = pd.DataFrame(rows, columns=["CYP", "UniProt"])
    return df.drop_duplicates().sort_values("CYP").reset_index(drop=True)


def main():
    pdf = sys.argv[1]
    out = Path(sys.argv[2])
    out.mkdir(exist_ok=True)

    print("Reading tables with Camelot (lattice)...")

    tables = camelot.read_pdf(pdf, pages="all", flavor="lattice")

    plant = pd.DataFrame()

    for t in tables:
        df = parse_table(t.df)
        plant = pd.concat([plant, df])

    plant = plant.drop_duplicates().reset_index(drop=True)

    plant.to_csv(out / "plant_cyps_uniprot.csv", index=False)
    plant["UniProt"].to_csv(out / "plant_uniprot_ids.txt", index=False, header=False)

    n = len(plant)

    print("\nQC CHECK")
    print("rows:", n)
    print("expected:", EXPECTED)

    if n != EXPECTED:
        raise RuntimeError("QC FAILED — Camelot didn't capture full table")

    print("✓ QC PASSED")


if __name__ == "__main__":
    main()
