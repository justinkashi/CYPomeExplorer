GOAL: replicate faithfully the results by Kuvek et al. for the plant cyp450 section (n=343), ignoring the human cyp section.

1. Fetching 343 plant cyp seq across 45 distinct families 

2. Fetching 343 plant cyp .cif (via AF2 server, via uniprot IDs). Structures were validated by their high pLDDT scores (AlphaFold confidence)

3. PDB Refrence Preparation: Translate and rotate the 4I3Q reference so the heme iron sits at (0,0,0) and the macrocycle (large circular porphyrin made up of 4 pyrrole rings, holding iron atom in place, the chassis) lies in the xy-plane, orienting the proximal cysteine toward -z to ensure the distal pocket faces +z. This template serves as the master coordinate system for global C_alpha alignment of all 343 plant AlphaFold models, projecting the 4I3Q origin into their empty pockets; subsequently, mask all heme atoms during vector generation to ensure the 320 rays capture true pocket volume rather than the cofactor

  QC REPORT FOR: data/4I3Q.cif 
  Iron Position: [0. 0. 0.] (Expected: [0,0,0]) -> translation worked 
  Avg Nitrogen Z: 0.3998 (Expected: ~0.0) -> heme is leveled.
  Proximal Cys Z: -2.07 (Expected: Negative)
  And look at reference_qc_plot.png -> orientation works, distal pocket is +z 


4. Align: Via pymol-align, perform a global $C_{\alpha}$ alignment of 343 plant CYP models against the oriented 4I3Q reference to project the heme iron origin into the vacant plant active sites. Discard all models with a $C_{\alpha}$ RMSD $> 7.0$ Å. If the authors did global c_a alignment, in boinformatics this usually implies sequence-dependant alignment. Also filter for C_a RMSD < 7 A relative to 4I3Q template ( there will be 342 already kept cuz authors told us the n=342 as the ones who are below 7A) 

5. (to fill up)