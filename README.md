# Plant CYPome Explorer (CYPomeExplorer)

A comparative computational framework for **functional mapping** and **substrate triage** across **plant cytochrome P450s (CYPs)**, intended for academics in **plant CYP enzymology**, **synthetic biology**, **metabolic engineering**, and **protein engineering**. The core descriptor is **binding site vectors (BSVs)**: a high-resolution representation of active-site **shape** and **electrostatics** derived from fixed-index radial rays cast from a binding-site center to the first intersected atomic surface, recording **vector length** and **partial charge**. :contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1}

---

## Why this exists

CYP nomenclature is historically based on sequence identity thresholds (family/subfamily), but **sequence similarity often fails to predict substrate specificity and regioselectivity**, due in part to **CYP fold plasticity** and flexible access-channel/active-site loops. This project provides a structure/ensemble-aware alternative: **compare CYPs by their binding-site landscapes** rather than by sequence alone. :contentReference[oaicite:2]{index=2}

---

## What this repo does

### 1) **Map CYP functional landscapes** with BSV similarity
- Generate BSV descriptors for each CYP (static structures or MD exemplars).
- Compute pairwise distances (shape-only, charge-only, combined/weighted RMSD).
- Build similarity trees and clusters that reflect **active-site similarity**.

BSVs and their RMSD-style comparisons are described as a structured pocket similarity framework in the reference work. :contentReference[oaicite:3]{index=3} :contentReference[oaicite:4]{index=4}

### 2) **Ensemble-aware characterization** (optional but recommended for CYPs)
CYP pockets are dynamic. If you have trajectories:
- Cluster frames into representative **exemplars**.
- Weight distances by exemplar populations (paper example uses **N = 15,625** frames per simulation). :contentReference[oaicite:5]{index=5}
- Optionally use **Affinity Propagation** to discover shared cross-isoform “binding-site states.” :contentReference[oaicite:6]{index=6}

### 3) **Batch docking triage** for CYP × metabolite screening
This repo does **not** claim to compute binding free energies. Instead, it supports large-scale screening by combining:
- docking scores/poses (from your docking engine of choice), **plus**
- BSV-based pocket/pose compatibility signals (shape/electrostatics, ensemble support)

Use case: screening plant CYPomes against metabolite libraries (e.g., LC–MS-derived candidates) to prioritize likely enzyme–substrate pairs for validation.

---

## Method summary: Binding Site Vectors (BSVs)

### Vector construction
- Subdivide an **icosahedron**, project vertices to a sphere, and trace rays from the sphere center.
- The reference implementation uses **492** directions globally and restricts to a **260-vector hemisphere** for the catalytic pocket region. :contentReference[oaicite:7]{index=7}

### Descriptor per direction
For each ray:
- **Length (`l`)**: distance from center to first van der Waals surface intersection.
- **Charge (`q`)**: partial charge of the intersected atom (charge model must be defined consistently).

### Similarity metric (shape + charge)
Combined normalized RMSD (with optional weighting) is used to compare binding sites. :contentReference[oaicite:8]{index=8}  
Example plant normalization constants reported in the paper’s case study: **σ_l = 3.833 Å**, **σ_q = 0.213 e** (when plant CYPs analyzed as a set). :contentReference[oaicite:9]{index=9}

---

## Plant-specific context

Plant CYPs are frequently represented by **predicted structures** (e.g., AlphaFold-class models), often lacking cofactors/ligands and with sparse experimental structural coverage. This repo is designed to work in that regime by emphasizing:
- consistent alignment / reference framing,
- pocket geometry/electrostatics descriptors,
- ensemble sampling when available.

Support figures in the reference supplementary focus heavily on plant CYP similarity trees built from backbone RMSD vs BSV RMSD. :contentReference[oaicite:10]{index=10}

---

## Inputs

- **CYP structures**: PDB/mmCIF (single structures) and/or ensemble snapshots/exemplars.
- **CYP sequences**: FASTA (for phylogenetic comparison and annotation sanity checks).
- **Metabolites/ligands**:
  - SMILES supported for storage and enumeration, but docking/BSV pose triage requires **3D coordinates** (SDF/MOL2/PDB poses).
- **Docking poses** (optional): docked complexes for pose filtering/ranking.

---

## Outputs

- **BSV descriptors**: per CYP structure or per exemplar (vector lengths + charges).
- **Pairwise matrices**:
  - sequence similarity / phylogeny (optional),
  - backbone RMSD (optional),
  - BSV RMSD (shape, charge, combined/weighted).
- **Trees**: dendrograms comparing sequence/backbone/BSV landscapes. :contentReference[oaicite:11]{index=11}
- **Ensemble summaries**: exemplar weights, cross-isoform cluster memberships. :contentReference[oaicite:12]{index=12}
- **CYP × ligand ranking tables**: combined evidence from docking + BSV compatibility + ensemble support.

---

## Pipeline overview

1. **Standardize structures**
   - consistent residue naming, protonation conventions, and coordinate frames across the CYPome.

2. **Align to a reference frame**
   - ensure vectors are index-comparable across CYPs (alignment errors corrupt BSV comparisons).

3. **Generate BSV descriptors**
   - compute `(l, q)` over the 260-ray hemisphere per structure or exemplar. :contentReference[oaicite:13]{index=13}

4. **Compute distances and clusters**
   - build BSV RMSD matrices; cluster and visualize functional landscapes. :contentReference[oaicite:14]{index=14}

5. **(Optional) MD ensembles**
   - cluster trajectories into exemplars; weight distances by exemplar populations. :contentReference[oaicite:15]{index=15}

6. **Batch docking + pose triage**
   - dock candidate metabolites; prioritize CYP–ligand pairs using BSV-based compatibility and ensemble support.

---

## Scope and non-goals

- **Not** a binding free-energy engine; BSV similarity is a **comparative descriptor** for functional grouping and screening support, not a direct ΔG/Kd predictor. :contentReference[oaicite:16]{index=16}
- Intended to **prioritize** candidates for mechanistic follow-up and experimental validation.

---

## Citation

If you use BSV descriptors, similarity metrics, or replicate the comparative analyses, cite the reference work describing the BSV framework and CYP functional landscape mapping. :contentReference[oaicite:17]{index=17} :contentReference[oaicite:18]{index=18}
