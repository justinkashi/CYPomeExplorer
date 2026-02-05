**Feb5**
- 
- align = iterative least‑squares fit with outlier rejection. It can drop badly fitting regions (e.g., disordered tails). Good for robustness, but not exactly what the authors used for backbone‑based comparisons.
cealign = Combinatorial Extension alignment (CE). It aligns based on structural fragments and reports RMSD from that alignment. This is what the authors’ scripts use for backbone RMSD matrices, so it’s closer to paper‑faithful.

- mean N z is not a “plane is tilted” metric. It’s mostly measuring the plane offset along Z (i.e., the heme plane can be perfectly parallel to XY but sit at z≈0.40 after you center on FE). Kuvek’s PyMOL transform does not require the nitrogen plane to pass through z=0; it centers on FE, not the nitrogen centroid.What you actually need to verify (faithfully) is:

FE is at origin

the heme plane normal is aligned with +Z (angle ≈ 0°)

the heme nitrogens have small z-spread (stddev small) after rotation

in-plane lock: NA projected to XY points to +X (y≈0, x>0) 
- The 4 values in that list refer to the four pyrrole nitrogen atoms (labeled NA, NB, NC, and ND) that surround the central iron atom in the heme group.

In a Cytochrome P450, these four nitrogens form the "base" of the active site. To ensure your binding site vectors sample the pocket accurately, the entire heme macrocycle must be horizontal.
- Step3: Running transrot: ok i use conda venv and its running transrot but stuck during rotation The hang is occurring inside the PyMOL engine's rotate call. This often happens in headless environments on macOS when the PyMOL API receives NumPy float64 objects instead of standard Python float types, or when the geometry engine deadlocks while trying to update a large selection. The most reliable way to install PyMOL with its required binary libraries on macOS is through conda or mamba. These managers bundle the necessary .dylib files.. The error Library not loaded: @rpath/libGLEW.2.1.dylib indicates that the PyMOL installation in the current virtual environment is failing to locate its binary dependencies. This is a common issue with pip-installed PyMOL on macOS, as it often fails to correctly link system libraries like GLEW.The libGLEW.2.1.dylib error happened because your pip version of PyMOL was a "Wheel" compiled on a machine where the user was named Martin. By using conda-forge, the package manager will download a version specifically compiled for your macOS and link it to the libraries inside your kuvek_env folder instead of a random user's path.
- Step3: You translated FE → origin, but you did not (and should not, if you want to match Kuvek) force the nitrogen plane to pass through the origin. If FE is not perfectly coplanar with NA/NB/NC/ND in the deposited structure (common), then after you rotate the plane normal to +Z you will still see a non-zero mean N z (it’s basically the signed distance of the N-plane from FE along the normal). Kuvek’s approach (translate FE + rotate orientation) also does not guarantee “mean N z = 0”.

- aligned PDBs (UniProt)
   ↓
surface.py / charge.py
   ↓
combine.py
   ↓
bsv_combined.txt (UniProt labels)
   ↓
numeric remap
   ↓
first_clustering.py
   ↓
clusters
   ↓
map back to UniProt

- Feature	Your Zenodo Protocol	My Previous Protocol (CDPKit)
Logic Split	Two-pass: Separate scripts for Length (surface.py) and Charge (charge.py).	Unified: One script (gen_kuvek_bp_descr.py) does both in one pass.
Charge Source	Uses an external reference_charges.txt lookup.	Uses charges assigned directly to the atoms (via GROMACS/MOL2).
Radii Source	Expects PQR format (columns 70-75) for surface.py.	Uses the OpenBabel radii library internally.
Normalization	Requires manual path editing in normalization.py.	Automated via the cmp_kuvek_bp_descrs.py comparison tool.

**Feb 4/5**
- 01_topo/

Topology preparation

gmx/

*_gmx.top → main topology

mol_1.itp → molecule params

top_ff.itp → force field include

gromos/

.top files → GROMOS topologies

cyp_heme.top → heme parameters

linking scripts (*.arg)

Purpose:
Force-field parameterization

02_coord/

Starting coordinates

CYP72A188.pdb → starting structure

small helper scripts

Purpose:
Initial structure

03_min/

Energy minimization

em.imd

em.run

Purpose:
Relax structure

04_simbox/

Solvation box

sim_box_*

solvent configs

Purpose:
Add water box

05_ion/

Ion addition

ion.arg

Purpose:
Neutralize system

06_eq/

Equilibration

equilibration.imd

job scripts

Purpose:
Thermal/pressure equilibration

07_md/

Production MD

.gro

.mdp

index files

Purpose:
Trajectory generation

Critical conclusion

This entire tree is:

MD trajectory generation only

Not descriptor generation.
- These belong to:

CDPKit → compiled descriptor path

They are not part of the Zenodo pipeline.

They implement the descriptor inside CDPKit for speed.

They do the same job as:

surface.py + charge.py + normalization.py + combine.py

in one compiled call.
- WHEN IMPLEMENTING THE ZENOTO FOLDER SCRIPTS:  Looking at surface.py logic:

reads atom radius from fixed columns

expects PQR format or PDB with radii embedded

Consequence

If you feed plain PDB → all distances will be wrong (radii = 0 or garbage)

- Non-identical parameters: 54a7 vs 54a8 differ in residue/atom typing, partial charges, and/or LJ parameters for some functional groups; even “small” changes can affect electrostatics, H-bonding, salt-bridges, and stability on MD timescales.

Downstream comparability: your trajectories could drift differently, so any derived pocket descriptors/vectors can shift. If the paper’s conclusions depend on subtle binding-site electrostatics/shape distributions, this matters.

It won’t necessarily break the workflow, but it becomes “close reproduction” rather than “same conditions”.
- writing GROMOS in Rust instead of C++ ? 
- if kuvek method is looking at the binding site without ligands, an improvement would be to generate dynamic binding site vectors when a ligand is there also we could simulate the MD trajectory of a typical CYP catalytic reaction (PETase too)
- gromos force fields: 54a7 and 54a8 have identical protein partial charges, 54a8 mainly updated:lipids small molecules solvent tweaks
- Failed to kekulize aromatic bonds: That warning is OpenBabel telling you it couldn’t assign a specific single/double-bond “Kekulé” pattern for something it thinks is aromatic. In protein PDBs this usually happens because PDB files don’t contain explicit bond orders, and OpenBabel’s bond-order perception gets confused (common with histidine rings, aromatic residues, or odd/partial atoms in predicted structures).Kekulization is the chemical informatics process of converting aromatic, delocalized bond representations (e.g., in RDKit or CDK) into a specific Lewis structure featuring localized alternating single and double bonds. It assigns exact bond orders, replacing aromatic flags, which is often required for visualization, fingerprinting, or specific chemical property calculations.  It converts "aromatic" (e.g., c1ccccc1) SMILES into "kekule" (e.g., C1=CC=CC=C1) SMILES.
- is there a parallel way of doing obabel add hydrogen so its faster? 

**Feb 3**
- 1. The "Template-Based Affinity" Strategy (Similarity-to-Binder)Since you lack experimental affinity data for most Tacca CYPs, you can use the BSV similarity as a proxy.The Logic: If a well-characterized CYP (e.g., a known sterol 14$\alpha$-demethylase or a human CYP3A4 known to bind bulky steroids) has a specific BSV "fingerprint," and a Tacca CYP has a very low BSV-RMSD to that template, it is highly likely to accommodate the same steroid scaffold.Application: Identify a plant CYP with known high affinity for phytosterols (like CYP51). Generate its 260 vectors. Calculate the distance between your Tacca CYPs and this "Gold Standard" sterol binder.2. Utilizing "Ligand Vectors" (Occupancy Mapping)The authors briefly mention Ligand Vectors (as seen in Figure S11 of the supplement). You can apply this to your taccalonolide intermediates.The Logic: Instead of just measuring the pocket, you shoot vectors from the iron origin that terminate at the van der Waals surface of the ligand itself.The Comparison: You then calculate the "Vector Clearance"—the difference between the Binding Site Vector (length to protein) and the Ligand Vector (length to metabolite).The Insight: A small clearance across all 260 vectors indicates a "tight fit" or high shape complementarity, which is a strong prerequisite for high affinity. Large gaps suggest a loose fit or poor complementarity for bulky taccalonolide scaffolds.3. Integrating BSVs into a Docking WorkflowTo get actual affinity scores, you should use the BSVs to "gate" or inform molecular docking.Step A (BSV Filter): Use the BSV similarity trees to cluster your Tacca CYPs. Only pick one or two representative enzymes from clusters that "look" like steroid-binding pockets.Step B (Focused Docking): Since the BSV approach has already defined the $(0, 0, 0)$ origin and the 260-vector "envelope," you can use these coordinates to define a highly precise docking box for tools like AutoDock Vina.Step C (Validation): Compare the Vina binding energy (kcal/mol) with the BSV similarity. You are looking for a correlation where enzymes with "Steroid-like" BSV fingerprints also yield the lowest binding energies for phytosterol intermediates.Why this matters for Tacca CYPs:
Taccalonolides are complex, highly oxygenated steroids. Standard sequence-based phylogenetic trees often fail to predict which CYP handles which intermediate because a single loop change can drastically alter the pocket volume.

By using the BSV approach, you can identify if a Tacca CYP has the volume expansion required to transition from a simple phytosterol precursor to a bulky, acetylated taccalonolide intermediate. If the BSV lengths in the +z distal region are significantly longer in one Tacca CYP than another, that enzyme is your prime candidate for the later steps of the biosynthetic pathway.

Would you like a script to calculate the "Vector Clearance" between a metabolite and your aligned Tacca pockets? This would be the most direct way to estimate "fit" using the authors' vector logic.
- This part is still basically true: Bio.PDB.Superimposer expects you to provide paired atom lists (i.e., a correspondence). With 343 plant CYPs spanning many families, you usually don’t have a clean residue mapping “for free,” so you’d need an extra step (sequence alignment / mapping / anchor motifs) before you can use it. That’s not in the paper; their choice of PyMOL avoids you having to implement correspondence yourself.
- step 04: #,Method,Logic,Verdict for Replication
1,MSA-Guided,Align sequences first; map Cα​ to Cα​.,"The Winner. Matches the authors' ""All Cα​"" description."
2,Structural (Simple),Matches atoms by index (Residue 1 to 1).,Too risky; fails if plant CYPs have different lengths.
3,CEAlign,Fragments (helices) are matched geometrically.,"High quality, but often skips atoms the authors included."
4,TM-Align,"Weighting scores based on the ""folded core.""","Great for AlphaFold, but technically different from authors."
5,Local (Cys-Loop),Aligning only the 10 residues around the Heme.,"Too narrow; doesn't align the global ""envelope."""
6,Foldseek,"3D shape ""alphabet"" comparison.","Optimized for speed, not for precise origin anchoring."
7,USalign,Heuristic maximum overlapping substructures.,"The most rigorous, but likely ""overkill"" for this study."
- Here are the four most significant methods you’re missing:1. TM-align (Template Modeling)While CEAlign looks at fragments, TM-align is the industry standard for AlphaFold comparisons. It uses a TM-score $(0.0 \text{ to } 1.0)$ rather than RMSD.The Logic: RMSD is highly sensitive to local errors (like one floppy loop). TM-score weights smaller distances more heavily than large ones, meaning it prioritizes the "folded core" and ignores the "disordered tail" much better than a standard least-squares fit.Why use it? If you have plant CYPs with extremely long, unstructured regions that are causing "FAIL" statuses in your current scripts.2. Local Alignment to the Heme-Binding "Signature"The authors focused on the pocket. Instead of aligning the entire protein (Global), you can align specifically to the Cys-ligand loop (FXXGXRXCXG).The Logic: You identify the conserved Cysteine that holds the iron in 4I3Q and the corresponding Cysteine in your plant model. You align only those 10–15 residues.The Result: This "locks" the pocket origin $(0,0,0)$ with extreme precision, even if the rest of the protein is slightly tilted or distorted .3. Foldseek (Structural Alphabet)This is the newest "speed demon" in the field.The Logic: It converts 3D structure into a "structural alphabet" (strings of letters representing shapes). It then performs an alignment as if it were a sequence alignment but using 3D information.Why use it? If you eventually scale this from 343 proteins to 100,000, Foldseek is orders of magnitude faster than anything else.4. USalign (Universal Structural Alignment)This is the "successor" to TM-align. It can align monomers to monomers, or even monomers to complexes.The Logic: It uses a heuristic to find the maximum overlapping substructures between the plant CYP and 4I3Q.Kuvek et al. Relevance: The authors mentioned maintaining a "spatially consistent sampling." USalign is arguably the most mathematically rigorous way to ensure that the $(0,0,0)$ center isn't just "in the pocket," but is oriented the same way relative to the I-helix and the F-G loop.
-  Many structural biologists use a tool called cealign (Combinatorial Extension) or an iterative least-squares fit. This method:

Starts by aligning the whole thing.

Identifies which atoms are far apart (the flexible loops).

Throws them out and re-aligns using only the stable "Core."

Repeats until the RMSD is minimized for the conserved P450 fold.
-  Plant CYPs often have long, disordered N-terminal signal peptides (to anchor them to the ER membrane) that the human 4I3Q reference does not have. If the authors simply aligned "all atoms" by index, the signal peptide would pull the entire protein out of alignment.They likely used an algorithm that performs an initial sequence alignment to identify which residue in the plant CYP corresponds to which residue in the 4I3Q, and then aligned only those pairs. 
- To ensure your script is "faithful" to the research, it must handle the fact that AlphaFold models often have long, unstructured N-terminal "tails" (signal peptides) that 4I3Q doesn't have. If you align by raw index (Residue 1 of Plant vs. Residue 1 of 4I3Q), the alignment will fail.
- The Goal: Since AlphaFold models lack the Heme cofactor, the authors "inherited" the Heme's $(0,0,0)$ position from the 4I3Q template. By aligning the protein backbones, the 4I3Q Iron position effectively becomes the "ghost" origin for the plant CYP pocket.
- The "Ghost" Heme: Remember that your AlphaFold models are apoproteins. This alignment is the only way the computer knows where the heme "should" be. If the alignment is poor (RMSD > 7), the 320 rays will be projected from a nonsensical position in the protein.Coordinate Preservation: When you save the aligned files (the .pdb outputs), you must ensure your script doesn't "re-center" them upon saving. The coordinates in the files must remain relative to the 4I3Q origin so that Step 5 (Vector Generation) can simply read $(0,0,0)$ as the starting point for every ray.
- now at step 4
- now downloading the AF2 database on t7 ssd: 

- downloading the 3 uniprot databases on external t7 ssd:  

aria2c -x 8 -s 8 -c -d /Volumes/T7 https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

aria2c -x 8 -s 8 -c -d /Volumes/T7 https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

aria2c -x 8 -s 8 -c -d /Volumes/T7 https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_kb.fasta.gz
- 01_fetch_seq_pdb: can accelerate by downloading teh uniprot and AF2 database localy on my mac (or on my external SSD)
- would there be a way to accelerate tesseract/poppler with GPU, in the 00_parse.py section of going from table S2 to a csv file ? 
- both OCR and tesseract-CV way failed to parse table s2 into a csv......... because of how this specific pdf file is encoded — not because tools are weak. This is a document-format problem, not a technology problem.
- Why each method failed (and will always be flaky)
pdfplumber / text extraction Fails because: text order ≠ reading order
multi-column mixing

tokens interleaved

You get garbage streams.

Camelot lattice

Fails because:

needs borders

you have none

so it detects 0 tables

Working as designed.

OCR

Fails because:

dense 3-column layout

small font

column drift

characters lost

OCR is probabilistic → not deterministic.