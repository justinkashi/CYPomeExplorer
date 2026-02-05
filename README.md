GOAL: Reproduce faithfully the results by Kuvek et al. for the plant cyp450 section (n=342), ignoring the human cyp section. 

1. The BSV Similarity Tree, (Figures S2 S3)
2. Functional clustering, (Case study 3) 
3. Comparing charge.txt  (FigureS10)
4. MD-Sim Insights (pocket fluctuations, representative conformers/exemplar from cluster centers in each sim)

PIPELINE
1. Fetching sequences 

2. Fetching structures via AF2 server 

3. 4I3Q Reference Preparation (trans_rot_4i3q.py)
  - I: 
  - O: 
  
4. 



For the MD-side: 
  - openbabel for vdW radii definitions
  - gromacs 54a8 force field
  - cdpkit library + CDPkit kuvek scripts 
  - MDAnalysis library (md trajectory handling)
  - sklearn AffinityPropagation 
