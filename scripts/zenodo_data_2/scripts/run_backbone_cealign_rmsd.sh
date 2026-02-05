#!/bin/bash
#SBATCH --job-name=rmsd_calc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32   # match --concurrency
#SBATCH --mem=32G            # adjust as needed
#SBATCH --time=72:00:00
#SBATCH --partition=MMS
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Initialize conda in this shell session (no need to run 'conda init' here)
source ~/miniconda3/etc/profile.d/conda.sh

# Activate your conda env
conda activate pymol_env

export QT_QPA_PLATFORM=offscreen

# Run your script with full paths and correct pymol_python executable
python3 /pool/teakuvek/new_vectors/pdb_humans/downloaded_pdb_cif/downloads_pdb_cifs_filtered_1/filtered_rmsd/fixed_pdbs/chainA_filtered/full_backbone_rmsd_cealign/script.py --out_csv /pool/teakuvek/new_vectors/pdb_humans/downloaded_pdb_cif/downloads_pdb_cifs_filtered_1/filtered_rmsd/fixed_pdbs/chainA_filtered/full_backbone_rmsd_cealign/rmsd_matrix.csv --workers 32 --concurrency 32 --folder /pool/teakuvek/new_vectors/pdb_humans/downloaded_pdb_cif/downloads_pdb_cifs_filtered_1/filtered_rmsd/fixed_pdbs/chainA_filtered/.
