GOAL: Reproduce faithfully the results by Kuvek et al. for the plant cyp450 section (n=342), ignoring the human cyp section. 

1. The BSV Similarity Tree, (Figures S2 S3)
2. Functional clustering, (Case study 3) 
3. Comparing charge.txt  (FigureS10)
4. MD-Sim Insights (pocket fluctuations, representative conformers/exemplar from cluster centers in each sim)

PIPELINE
1. Fetching sequences 

2. Fetching structures via AF2 server 

3. 4I3Q Reference Preparation 
  - custom: 03_prepare_reference.py -> 4I3Q_reference.pdb
  - author: transrot_4i3q.py -> 4i3q_std.pdb 
  - test pymol align 2 results -> RMSD is 0.001 A

4. Batch align the 342 plant cyps via pymol cealign-> pass 
  - need to find correct residue ranges to align and not "bury" the vectors 

5. Generate hemisphere lattice
  - need -s 7 to generate 260 points 

6. Convert PDB -> PQR (openbabel, add radii)

7. Generate surface vectors (surface.py): pass  
  - bottleneck: uses one core -> tweaked it using python's multiprocessing (surface_parallel.py) -> pass -> issue: CPU 100% so its in bursts instead of smooth 
  - GPU/Cuda version: surface_gpu.py, pass 

8. Generate charge vectors (charge.py)
  - charge_parallel.py (100% CPU usage via multiprocessing)
  - charge_gpu.py (cuda) 

9. Normalize (normalization.py) 
  - discrepancy  

10. Combine to final BSV (combine.py)

11. compare to kuvek's vector_matrix.txt




For the MD-side: 
  - openbabel for vdW radii definitions
  - gromacs 54a8 force field
  - cdpkit library + CDPkit kuvek scripts 
  - MDAnalysis library (md trajectory handling)
  - sklearn AffinityPropagation 
