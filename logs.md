Feb 3
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