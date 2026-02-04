Feb 3
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