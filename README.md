GOAL: replicate faithfully the results by Kuvek et al. for the plant cyp450 section (n=343), ignoring the human cyp section.

1. Fetching 343 plant cyp seq across 45 distinct families 

2. Fetching 343 plant cyp .cif (via AF2 server, via uniprot IDs). Structures were validated by their high pLDDT scores (AlphaFold confidence)

3. Prepare Reference PDB: and rotated reference template 4I3Q so that (1) the Heme Iron (Fe) atom is exactly at the coordinate origin (0, 0, 0) and (2) the Heme macrocycle sits perfectly in the xy-plane.

  * note: to optimize later on for novel plant CYP

4. Align: Global C_a alignment of 343 .cif against CYP3A4 (4I3Q.pdb) to physically move the plant CYP's predicted active site cavity to surround that same (0, 0, 0) point. 

5. Filter: models with C_a RMSD < 7 A relative to 4I3Q template were kept ()