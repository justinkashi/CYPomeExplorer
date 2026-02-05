*MD simulation files*

The directories human_cyps and plant_cyps include all the files needed to reproduce the MD simulation trajectories. Each CYP system contains seven numbered subdirectories, ordered from topology to md. They have a prefix number which indicates the usage order. Additional supporting files, such as ion topologies, .mtb, .ifp, and library files, can be found in the files/ directory

Useful info:
 -For each system, three GROMOS topologies are available:
 	- cyp.top — topology where the heme group is not bound to the protein.
 	- cyp_heme.top — topology where the heme group is bound to the cysteine residue (CYS).
 	- CYP90D1.top (example) — topology for the protein with bound heme and Na⁺/Cl⁻ ions.
	
Equilibration was run using GROMOS, and production MD simulations were run with GROMACS.

*Binding site vectors scripts*

	- triangular_latice_sphere.py: Generates a triangular lattice sphere used for binding site vectors generation. 
	Inputs: Sphere radius, number of triangle subdivisions, and output file name. 
	Includes an optional --hemisphere flag to output only the upper hemisphere (used for CYP work).

	- surface.py: Computes binding site vector lengths. 
	Inputs: PQR file, lattice sphere file, output text file name, and distance cutoff from the binding site center (which atoms to consider).
	Output: Text file containing vectors' length data.
	
	- charge.py: Computes binding site vector charges.
	Inputs: Same as surface.py, plus a reference charge file (reference_charges.txt, included in the scripts/ folder).
	Output: Text file containing vectors' charge data.
	The reference file contains partial charges for the GROMOS 54a8 force field that are parsed in the charge.py script.
	
	- normalization.py: Normalizes the outputs of surface.py or charge.py to enable direct comparison of binding site properties.
	
	- combine.py: Merges normalized vector length and charge data into a single file. It gives a full binding site description in one file. 
	
	- trans_rot_4i3q.py: Translates and rotates the reference 4I3Q PDB structure so that heme lies in the x–y plane, and the iron atom is at (0, 0, 0).
	
	- trans_first_frame.py: After alignment to the 4I3Q reference structure, translates the first frame of each trajectory so that the heme iron is positioned at (0, 0, 0).
	
	- align_selected_residues.py + run_align_selected_residues.sh: Aligns all trajectory snapshots to the first (aligned) snapshot using selected residues. 
	Executed via the accompanying Bash script.
	
	- first_clustering_vectors.py: Clusters trajectory snapshots of a single CYP system based on binding site vectors.
	Input: Text file retreived from combine.py (or any file containing vector lengths and/or charge data)
	Output: Text file with number of clusters, cluster centers, cluster population. 
	
	- post_first_clustering.py: Combines the exemplars of first clustering into a single text file, including their weights.
	
	- second_clustering_vectors.py: Clusters the exemplars from clustering round one. 
	Input: Text file retreived from post_first_clustering.py
	Output: Text file with number of clusters, cluster centers, cluster population.

	- MD_clustering.py: Clusters trajectory snapshots of a single CYP system based on backbone structural overlap.
	Inputs: Md trajectory in xtc format, corresponding topology in pdb format
	Output: Text file with number of clusters, cluster centers, cluster population.
	
	- backbone_cealign_rmsd.py + run_backbone_cealign_rmsd.sh: Generates pairwise RMSD matrices based on backbone similarity for all studied CYPs. 
	Executed via the provided Bash script.
	













