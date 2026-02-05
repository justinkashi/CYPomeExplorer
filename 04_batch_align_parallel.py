import os
import sys
import glob
import subprocess
from multiprocessing import Pool, cpu_count

REF = sys.argv[1]
IN_DIR = sys.argv[2]
OUT_DIR = sys.argv[3]
QC_DIR = sys.argv[4]

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)

# The specific conserved ranges from the author's script
RANGES = "5-8+32-40+45-49+52-56+62-66+97-105+115-136+146-158+193-197+219-230+265-277+283-295+301-309+320-324+330-339+346-349+366-369+371-375+418-433+464-468"

def align_one(pdb_path):
    name = os.path.splitext(os.path.basename(pdb_path))[0]
    out_pdb = os.path.join(OUT_DIR, f"{name}.aligned.pdb")
    out_txt = os.path.join(QC_DIR, f"{name}_rmsd.txt")
    
    # Paper-faithful PyMOL command string
    pymol_cmd = f"""
from pymol import cmd
cmd.load('{REF}', 'ref')
cmd.load('{pdb_path}', 'mob')
sel_mob = 'mob and (resi {RANGES}) and name N+CA+C+O'
sel_ref = 'ref and (resi {RANGES}) and name N+CA+C+O'
r = cmd.cealign(sel_mob, sel_ref)
coords = cmd.get_coords('mob and resn HEM and name FE')
if coords is not None and len(coords) > 0:
    x, y, z = coords[0]
    cmd.translate([-x, -y, -z], 'mob')
cmd.save('{out_pdb}', 'mob')
with open('{out_txt}', 'w') as f:
    f.write(str(r.get('RMSD', 0.0)))
"""
    # Launching PyMOL in headless mode (-cq)
    subprocess.run(["pymol", "-cq", "-d", pymol_cmd], stdout=subprocess.DEVNULL)
    return name

if __name__ == "__main__":
    pdbs = glob.glob(os.path.join(IN_DIR, "*.pdb"))
    num_workers = max(1, cpu_count() - 1) # Leave one core free
    
    print(f"🚀 Aligning {len(pdbs)} proteins using {num_workers} cores...")
    
    with Pool(num_workers) as p:
        for i, name in enumerate(p.imap_unordered(align_one, pdbs), 1):
            print(f"   [{i}/{len(pdbs)}] Finished: {name}", end="\r")
            
    print("\n✅ Parallel alignment complete.")