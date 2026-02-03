Task 1 — Replicate the paper’s results (minimal reproducible plan)

Goal: reproduce the three “tree types” on static structures (Case Study 1), and the MD/ensemble clustering + weighted similarity trees + cross-isoform state overlap (Case Studies 2–3).

1A) Static-structure replication (Case Study 1)

Get the structure sets

Human CYP PDB IDs and plant UniProt IDs are listed in the SI tables. 

kuvek_supp_2026

Standardize and align all structures

Use a single reference: CYP3A4 (4I3Q), translated/rotated so heme plane = xy and Fe at (0,0,0). 

binding-site-vectors-enable-map…

Align all CYPs by backbone atoms to the reference. For plant AF2 models (heme absent), assume Fe at origin after alignment. 

binding-site-vectors-enable-map…

Generate BSV descriptors

Use 260 vectors (hemisphere above heme plane); assign vdW radii (OpenBabel) and partial charges (paper uses GROMOS 54a8). 

binding-site-vectors-enable-map…

Compute σl and σq per dataset

Compute separately for human and plant sets from all vector lengths and charges. 

binding-site-vectors-enable-map…

Build three similarity trees

Phylogeny: MSA (Clustal Omega) → sequence similarity matrix. 

binding-site-vectors-enable-map…

Backbone: PyMOL cealign pairwise RMSD. 

binding-site-vectors-enable-map…

BSV: normalized combined RMSD (eq 4) using σl, σq. 

binding-site-vectors-enable-map…

Dendrogram: SciPy hierarchical clustering (average linkage). 

binding-site-vectors-enable-map…

1B) Ensemble replication (Case Studies 2–3)

MD ensemble inputs

Paper uses 5×500 ns, 15,625 conformations/isoform. 

binding-site-vectors-enable-map…

For replication, you can start with smaller ensembles but keep the same analysis steps.

Per-isoform clustering

Compute RMSD matrices from BSVs (shape-only, charge-only, combined, weighted). 

binding-site-vectors-enable-map…

Cluster each isoform’s frames using Affinity Propagation, keep exemplars + weights. 

binding-site-vectors-enable-map…

Weighted similarity matrices between isoforms

Weight pairwise exemplar RMSDs by exemplar weights and total frames (eq 7) and sum across exemplar pairs (eq 8). 

binding-site-vectors-enable-map…

Cross-isoform conformational state overlap (Case Study 3)

Merge exemplars across isoforms; run a second Affinity Propagation clustering on BSV RMSD (eq 4), propagate weights to get per-isoform occupancy in each global state. 

binding-site-vectors-enable-map…

Reproduce SI figures

Trees: plant/human BSV dendrograms (S4/S6), weighted dendrograms (S8), seven-cluster representative pockets (S10), ligand-vectors proof-of-concept (S11). 

kuvek_supp_2026

 

kuvek_supp_2026

Task 2 — Apply your implementation to Tacca CYP450s + metabolites

Goal: produce (i) a functional landscape of the Tacca CYPome and (ii) a CYP×metabolite prioritization table, with optional ensemble support.

Inputs

CYP structures (heme-present preferred; if absent, you still can anchor by alignment and inferred Fe).

Metabolites: SMILES → 3D conformers/poses (SDF/MOL2) and/or docked poses.

Reference frame

Choose one heme-containing CYP template (can be 4I3Q or a representative plant CYP with heme modeled) and enforce the same origin and orientation across all Tacca structures.

BSV generation and clustering

Generate BSVs for all Tacca CYPs (static) → dendrogram/UMAP/cluster labels.

Docking + BSV “binary criterion” pose triage

For each docked pose, compute ligand vectors and apply: ligand_len(dir) ≤ pocket_len(dir) − margin (steric feasibility); optionally require local charge complementarity.

Summarize compatibility fraction across conformations (if you have ensembles).

Ranking

Per (CYP, metabolite): combine docking score + BSV compatibility (+ ensemble occupancy) into a single ranked output, keeping the components for interpretability.

Deliverables

tacca_bsv_tree.*, tacca_clusters.csv, cyp_metabolite_rankings.parquet/csv, optional state_overlap.csv.

Code for Task 1 and Task 2

Below is a self-contained Python implementation that:

generates 260-hemisphere BSVs (ray–sphere intersection against atom vdW spheres),

supports PQR inputs for per-atom charges (recommended; you can export PQR from your parametrization pipeline),

computes σl, σq,

computes BSV RMSD matrices, hierarchical trees, Affinity Propagation exemplars + weights,

computes weighted isoform–isoform RMSDs (eq 7–8),

computes ligand vectors and a binary docking compatibility check (as described conceptually in the paper). 

binding-site-vectors-enable-map…

Suggested repo layout
plant-cypome-explorer/
  cypomeexplorer/
    __init__.py
    vectors.py
    io.py
    bsv.py
    clustering.py
    docking.py
  scripts/
    01_bsv_static.py
    02_bsv_tree.py
    03_md_ap_exemplars.py
    04_weighted_isoform_matrix.py
    05_cross_isoform_states.py
    06_tacca_screen.py
  requirements.txt

requirements.txt
numpy
pandas
scipy
scikit-learn
biopython
rdkit-pypi

cypomeexplorer/io.py
# usage: from cypomeexplorer.io import read_pqr, read_pdb_atoms
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple, List
import re

@dataclass(frozen=True)
class Atom:
    x: float
    y: float
    z: float
    element: str
    charge: float
    vdw: float

_PDB_ATOM_RE = re.compile(r"^(ATOM  |HETATM)")

# Simple vdW radii (Å). If you want OpenBabel-perfect radii, precompute radii into PQR/B-factor or use an external mapping.
_VDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "FE": 1.80, "MG": 1.73, "ZN": 1.39, "CU": 1.40, "MN": 1.39,
}

def _guess_element(pdb_line: str) -> str:
    el = pdb_line[76:78].strip()
    if el:
        return el.upper()
    name = pdb_line[12:16].strip().upper()
    if len(name) >= 2 and name[0:2] in _VDW:
        return name[0:2]
    return name[0:1] if name else "C"

def vdw_radius(element: str) -> float:
    e = element.upper()
    return _VDW.get(e, _VDW.get(e[:1], 1.70))

def read_pqr(path: str | Path) -> List[Atom]:
    """
    Read PQR: columns include x,y,z, charge, radius.
    Recommended path to get charges consistent with your force field.
    """
    atoms: List[Atom] = []
    with open(path, "r") as f:
        for line in f:
            if not _PDB_ATOM_RE.match(line):
                continue
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            charge = float(line[54:62])
            radius = float(line[62:70])
            element = _guess_element(line)
            atoms.append(Atom(x, y, z, element, charge, radius))
    return atoms

def read_pdb_atoms(path: str | Path, default_charge: float = 0.0) -> List[Atom]:
    """
    Read PDB without charges; assigns default_charge and simple vdW radii.
    Use PQR if you want paper-like electrostatics.
    """
    atoms: List[Atom] = []
    with open(path, "r") as f:
        for line in f:
            if not _PDB_ATOM_RE.match(line):
                continue
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            element = _guess_element(line)
            atoms.append(Atom(x, y, z, element, default_charge, vdw_radius(element)))
    return atoms

cypomeexplorer/vectors.py
# usage: from cypomeexplorer.vectors import hemisphere_directions_icosa
from __future__ import annotations
import numpy as np
from functools import lru_cache
from typing import Tuple

def _icosahedron_vertices() -> np.ndarray:
    phi = (1.0 + 5.0 ** 0.5) / 2.0
    v = np.array([
        [-1,  phi, 0], [1,  phi, 0], [-1, -phi, 0], [1, -phi, 0],
        [0, -1,  phi], [0,  1,  phi], [0, -1, -phi], [0,  1, -phi],
        [ phi, 0, -1], [ phi, 0,  1], [-phi, 0, -1], [-phi, 0,  1],
    ], dtype=float)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v

def _icosahedron_faces() -> np.ndarray:
    # Standard icosahedron triangulation using vertex indices matching _icosahedron_vertices order.
    return np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ], dtype=int)

def _subdivide_triangle(a: np.ndarray, b: np.ndarray, c: np.ndarray, n: int) -> np.ndarray:
    # Barycentric subdivision of triangle into (n+1)(n+2)/2 points.
    pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            p = (i * a + j * b + k * c) / n
            pts.append(p)
    return np.array(pts, dtype=float)

@lru_cache(maxsize=8)
def sphere_directions_icosa(subdiv: int = 7) -> np.ndarray:
    """
    Returns unit vectors sampled by subdividing an icosahedron and projecting to the sphere.
    'subdiv' controls density. Increase if you want closer to paper's ~492 sphere directions.
    """
    v = _icosahedron_vertices()
    f = _icosahedron_faces()

    pts = []
    for tri in f:
        a, b, c = v[tri[0]], v[tri[1]], v[tri[2]]
        # subdiv >= 2 required
        sub_pts = _subdivide_triangle(a, b, c, subdiv)
        pts.append(sub_pts)

    pts = np.vstack(pts)
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)

    # Deduplicate (projection creates repeats across adjacent faces)
    key = np.round(pts, 6)
    _, idx = np.unique(key, axis=0, return_index=True)
    pts = pts[np.sort(idx)]
    return pts

def hemisphere_directions_icosa(subdiv: int = 7, zmin: float = 1e-9) -> np.ndarray:
    """
    Hemisphere above the heme plane: keep directions with z > 0 in the reference frame.
    Paper uses 260 rays; tune subdiv until len(...) ~= 260.
    """
    dirs = sphere_directions_icosa(subdiv=subdiv)
    hemi = dirs[dirs[:, 2] > zmin]
    return hemi

cypomeexplorer/bsv.py
# usage: from cypomeexplorer.bsv import bsv_from_atoms, bsv_rmsd, compute_sigmas
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, Dict, List
from .io import Atom

@dataclass(frozen=True)
class BSV:
    lengths: np.ndarray  # (Ndirs,)
    charges: np.ndarray  # (Ndirs,) charge at intersection atom

def _ray_sphere_first_intersection(u: np.ndarray, centers: np.ndarray, radii: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    For a set of spheres, compute the first intersection distance t along ray r(t)=t*u from origin.
    Returns (tmin, idx) where tmin is per-ray min distance, idx is argmin sphere index.
    """
    # centers: (Nat,3), radii: (Nat,)
    # For each atom: b = dot(c,u), d2 = |c|^2 - b^2, intersection t = b - sqrt(r^2 - d2)
    b = centers @ u  # (Nat,)
    mask_front = b > 0.0
    c2 = np.einsum("ij,ij->i", centers, centers)
    d2 = c2 - b * b
    r2 = radii * radii
    mask_hit = mask_front & (d2 <= r2)
    if not np.any(mask_hit):
        return np.inf, -1
    sqrt_term = np.sqrt(np.maximum(r2[mask_hit] - d2[mask_hit], 0.0))
    t = b[mask_hit] - sqrt_term
    # keep positive
    valid = t > 0.0
    if not np.any(valid):
        return np.inf, -1
    t_valid = t[valid]
    hit_indices = np.where(mask_hit)[0][valid]
    j = int(np.argmin(t_valid))
    return float(t_valid[j]), int(hit_indices[j])

def bsv_from_atoms(atoms: List[Atom], directions: np.ndarray) -> BSV:
    centers = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
    radii = np.array([a.vdw for a in atoms], dtype=float)
    charges_atoms = np.array([a.charge for a in atoms], dtype=float)

    N = directions.shape[0]
    lengths = np.full(N, np.inf, dtype=float)
    charges = np.zeros(N, dtype=float)

    for i in range(N):
        u = directions[i]
        t, idx = _ray_sphere_first_intersection(u, centers, radii)
        lengths[i] = t
        charges[i] = 0.0 if idx < 0 else charges_atoms[idx]

    # Replace inf (no hit) with max finite length (consistent cutoff) to keep RMSD defined.
    finite = np.isfinite(lengths)
    if np.any(finite):
        Lmax = float(np.max(lengths[finite]))
        lengths[~finite] = Lmax
    else:
        lengths[:] = 0.0
    return BSV(lengths=lengths, charges=charges)

def compute_sigmas(bsvs: List[BSV]) -> Tuple[float, float]:
    L = np.concatenate([b.lengths for b in bsvs]).astype(float)
    Q = np.concatenate([b.charges for b in bsvs]).astype(float)
    return float(np.std(L, ddof=0)), float(np.std(Q, ddof=0))

def bsv_rmsd(a: BSV, b: BSV, sigma_l: float, sigma_q: float, mode: str = "combined", d: float = 0.0) -> float:
    """
    mode: 'shape', 'charge', 'combined', 'weighted'
    weighted uses parameter d (paper-style shape vs charge emphasis).
    """
    dl = (a.lengths - b.lengths) / (sigma_l if sigma_l > 0 else 1.0)
    dq = (a.charges - b.charges) / (sigma_q if sigma_q > 0 else 1.0)

    if mode == "shape":
        x = dl
    elif mode == "charge":
        x = dq
    elif mode == "combined":
        x = np.concatenate([dl, dq])
    elif mode == "weighted":
        # paper-style weights: ws=2/(2+d), wc=2(1+d)/(2+d)
        ws = 2.0 / (2.0 + d)
        wc = 2.0 * (1.0 + d) / (2.0 + d)
        x = np.concatenate([np.sqrt(ws) * dl, np.sqrt(wc) * dq])
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return float(np.sqrt(np.mean(x * x)))

cypomeexplorer/clustering.py
# usage: from cypomeexplorer.clustering import ap_exemplars_from_rmsd, hierarchical_tree
from __future__ import annotations
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple
from sklearn.cluster import AffinityPropagation
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

def rmsd_matrix_from_pairwise(values: np.ndarray) -> np.ndarray:
    # expects square matrix already; helper kept for symmetry checks
    M = 0.5 * (values + values.T)
    np.fill_diagonal(M, 0.0)
    return M

def ap_exemplars_from_rmsd(rmsd: np.ndarray, preference: float | None = None, damping: float = 0.75) -> Tuple[np.ndarray, np.ndarray]:
    """
    AffinityPropagation expects similarities; we provide S = -rmsd.
    Returns (exemplar_indices, labels).
    """
    S = -rmsd
    if preference is None:
        preference = float(np.median(S))
    ap = AffinityPropagation(affinity="precomputed", preference=preference, damping=damping, random_state=0)
    ap.fit(S)
    exemplars = ap.cluster_centers_indices_
    labels = ap.labels_
    return exemplars, labels

def exemplar_weights(labels: np.ndarray) -> Dict[int, int]:
    uniq, cnt = np.unique(labels, return_counts=True)
    return {int(u): int(c) for u, c in zip(uniq, cnt)}

def hierarchical_tree(rmsd: np.ndarray, labels: List[str], method: str = "average") -> dict:
    """
    Returns dendrogram data; you can plot in your notebook/script.
    """
    condensed = squareform(rmsd, checks=False)
    Z = linkage(condensed, method=method)
    return {"Z": Z, "labels": labels}

cypomeexplorer/docking.py
# usage: from cypomeexplorer.docking import ligand_vectors, binary_compatibility
from __future__ import annotations
import numpy as np
from typing import List, Tuple, Optional
from .io import Atom
from .bsv import BSV, _ray_sphere_first_intersection

def ligand_vectors(ligand_atoms: List[Atom], directions: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Ligand vectors: for each direction, extend to ligand surface at the farthest intersection point.
    Returns (lens, charges, hitmask).
    """
    centers = np.array([[a.x, a.y, a.z] for a in ligand_atoms], dtype=float)
    radii = np.array([a.vdw for a in ligand_atoms], dtype=float)
    charges_atoms = np.array([a.charge for a in ligand_atoms], dtype=float)

    N = directions.shape[0]
    lens = np.full(N, np.nan, dtype=float)
    charges = np.zeros(N, dtype=float)
    hit = np.zeros(N, dtype=bool)

    for i in range(N):
        u = directions[i]
        # For ligand, we want farthest point along u among intersections (not nearest).
        b = centers @ u
        mask_front = b > 0.0
        c2 = np.einsum("ij,ij->i", centers, centers)
        d2 = c2 - b * b
        r2 = radii * radii
        mask_hit = mask_front & (d2 <= r2)
        if not np.any(mask_hit):
            continue
        sqrt_term = np.sqrt(np.maximum(r2[mask_hit] - d2[mask_hit], 0.0))
        # two intersections: b - sqrt, b + sqrt; farthest is b + sqrt
        t_far = b[mask_hit] + sqrt_term
        valid = t_far > 0.0
        if not np.any(valid):
            continue
        t_valid = t_far[valid]
        hit_indices = np.where(mask_hit)[0][valid]
        j = int(np.argmax(t_valid))
        lens[i] = float(t_valid[j])
        charges[i] = float(charges_atoms[int(hit_indices[j])])
        hit[i] = True

    return lens, charges, hit

def binary_compatibility(pocket: BSV, lig_lens: np.ndarray, lig_q: np.ndarray, hit: np.ndarray,
                         margin: float = 0.0, require_charge_complement: bool = False) -> bool:
    """
    Binary docking criterion (pose triage):
    - For each direction that intersects ligand: ligand length must be <= pocket length - margin.
    - Optional: require charge complementarity (sign-opposite) where ligand is near pocket.
    """
    ok_geom = np.all(lig_lens[hit] <= (pocket.lengths[hit] - margin))
    if not ok_geom:
        return False
    if not require_charge_complement:
        return True
    # Complementarity heuristic: product should be <= 0 (opposite sign or one near zero).
    return bool(np.all(lig_q[hit] * pocket.charges[hit] <= 0.0))

Scripts (Task 1 and Task 2)
scripts/01_bsv_static.py (Task 1A + Task 2 static)
# usage: python scripts/01_bsv_static.py --inputs data/plant_pqrs --out out/plant_bsv.parquet --subdiv 7
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

from cypomeexplorer.io import read_pqr, read_pdb_atoms
from cypomeexplorer.vectors import hemisphere_directions_icosa
from cypomeexplorer.bsv import bsv_from_atoms

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, help="Folder containing .pqr (preferred) or .pdb files already aligned to reference frame.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--subdiv", type=int, default=7, help="Icosa subdivision; tune until hemisphere ~260 directions.")
    args = ap.parse_args()

    inp = Path(args.inputs)
    dirs = hemisphere_directions_icosa(subdiv=args.subdiv)
    print(f"Using {dirs.shape[0]} hemisphere directions")

    rows = []
    for p in sorted(inp.glob("*.pqr")) + sorted(inp.glob("*.pdb")):
        atoms = read_pqr(p) if p.suffix.lower() == ".pqr" else read_pdb_atoms(p)
        b = bsv_from_atoms(atoms, dirs)
        rows.append({
            "id": p.stem,
            "lengths": b.lengths.astype(np.float32),
            "charges": b.charges.astype(np.float32),
        })

    df = pd.DataFrame(rows)
    df.to_parquet(args.out, index=False)
    print(f"Wrote {len(df)} BSVs -> {args.out}")

if __name__ == "__main__":
    main()

scripts/02_bsv_tree.py (Task 1A trees; Task 2 landscape)
# usage: python scripts/02_bsv_tree.py --bsv out/plant_bsv.parquet --out_prefix out/plant --mode combined
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd

from cypomeexplorer.bsv import BSV, compute_sigmas, bsv_rmsd
from cypomeexplorer.clustering import hierarchical_tree, rmsd_matrix_from_pairwise

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bsv", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--mode", choices=["shape","charge","combined","weighted"], default="combined")
    ap.add_argument("--d", type=float, default=0.0, help="Weighting parameter for mode=weighted")
    args = ap.parse_args()

    df = pd.read_parquet(args.bsv)
    bsvs = [BSV(np.array(r.lengths), np.array(r.charges)) for r in df.itertuples(index=False)]

    sigma_l, sigma_q = compute_sigmas(bsvs)
    print(f"sigma_l={sigma_l:.4f}, sigma_q={sigma_q:.4f}")

    n = len(bsvs)
    M = np.zeros((n,n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            M[i,j] = bsv_rmsd(bsvs[i], bsvs[j], sigma_l, sigma_q, mode=args.mode, d=args.d)
    M = rmsd_matrix_from_pairwise(M)

    tree = hierarchical_tree(M, labels=df["id"].tolist(), method="average")
    np.save(args.out_prefix + ".rmsd.npy", M)
    np.save(args.out_prefix + ".linkage.npy", tree["Z"])
    pd.Series({"sigma_l": sigma_l, "sigma_q": sigma_q}).to_csv(args.out_prefix + ".sigmas.csv")
    print(f"Wrote RMSD + linkage -> {args.out_prefix}.*")

if __name__ == "__main__":
    main()

scripts/03_md_ap_exemplars.py (Task 1B per-isoform AP exemplars)
# usage: python scripts/03_md_ap_exemplars.py --isoform_bsv out/md/CYP3A4_frames.parquet --out out/md/CYP3A4_ap.parquet
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd

from cypomeexplorer.bsv import BSV, compute_sigmas, bsv_rmsd
from cypomeexplorer.clustering import ap_exemplars_from_rmsd, exemplar_weights, rmsd_matrix_from_pairwise

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isoform_bsv", required=True, help="Parquet with one row per frame: lengths, charges")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["shape","charge","combined","weighted"], default="combined")
    ap.add_argument("--d", type=float, default=0.0)
    ap.add_argument("--max_frames", type=int, default=2000, help="Optional downsample for speed")
    args = ap.parse_args()

    df = pd.read_parquet(args.isoform_bsv)
    if len(df) > args.max_frames:
        df = df.sample(args.max_frames, random_state=0).reset_index(drop=True)

    bsvs = [BSV(np.array(r.lengths), np.array(r.charges)) for r in df.itertuples(index=False)]
    sigma_l, sigma_q = compute_sigmas(bsvs)

    n = len(bsvs)
    M = np.zeros((n,n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            M[i,j] = bsv_rmsd(bsvs[i], bsvs[j], sigma_l, sigma_q, mode=args.mode, d=args.d)
    M = rmsd_matrix_from_pairwise(M)

    exemplars, labels = ap_exemplars_from_rmsd(M)
    w = exemplar_weights(labels)

    out = pd.DataFrame({
        "frame_idx": np.arange(n),
        "cluster": labels,
    })
    out["cluster_weight"] = out["cluster"].map(w)
    out["is_exemplar"] = False
    out.loc[exemplars, "is_exemplar"] = True
    out.to_parquet(args.out, index=False)
    print(f"Exemplars={len(exemplars)} clusters={len(w)} -> {args.out}")

if __name__ == "__main__":
    main()

scripts/04_weighted_isoform_matrix.py (Task 1B eq 7–8)
# usage: python scripts/04_weighted_isoform_matrix.py --isoforms out/md_isoforms_list.txt --bsv_dir out/md_bsv --ap_dir out/md_ap --out out/isoform_weighted.npy
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from cypomeexplorer.bsv import BSV, bsv_rmsd
from cypomeexplorer.clustering import rmsd_matrix_from_pairwise

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isoforms", required=True, help="Text file: one isoform name per line")
    ap.add_argument("--bsv_dir", required=True, help="Folder with {isoform}.parquet (frames) or {isoform}_exemplars.parquet")
    ap.add_argument("--ap_dir", required=True, help="Folder with {isoform}_ap.parquet (frame->cluster, exemplar flags, cluster_weight)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sigma_l", type=float, required=True)
    ap.add_argument("--sigma_q", type=float, required=True)
    ap.add_argument("--mode", choices=["shape","charge","combined","weighted"], default="combined")
    ap.add_argument("--d", type=float, default=0.0)
    ap.add_argument("--N_total", type=int, default=15625, help="Total frames per isoform used for weighting (paper uses 15625).")
    args = ap.parse_args()

    isoforms = [l.strip() for l in open(args.isoforms) if l.strip()]
    bsv_dir = Path(args.bsv_dir)
    ap_dir = Path(args.ap_dir)

    # Load exemplar BSVs + weights per isoform
    exemplars = {}
    weights = {}
    for iso in isoforms:
        frames = pd.read_parquet(bsv_dir / f"{iso}.parquet")
        apmap = pd.read_parquet(ap_dir / f"{iso}_ap.parquet")
        ex_idx = apmap.index[apmap["is_exemplar"].values].to_numpy()
        # exemplar BSVs
        ex_bsv = [BSV(np.array(frames.loc[i, "lengths"]), np.array(frames.loc[i, "charges"])) for i in ex_idx]
        exemplars[iso] = ex_bsv
        # cluster weights correspond to exemplar’s cluster
        ex_clusters = apmap.loc[ex_idx, "cluster"].to_numpy()
        cluster_weight = apmap.drop_duplicates("cluster").set_index("cluster")["cluster_weight"].to_dict()
        ex_w = np.array([cluster_weight[int(c)] for c in ex_clusters], dtype=float)
        weights[iso] = ex_w

    n = len(isoforms)
    M = np.zeros((n,n), dtype=float)

    for i, isoA in enumerate(isoforms):
        for j, isoB in enumerate(isoforms[i+1:], start=i+1):
            rmsd_sum = 0.0
            for a_idx, a in enumerate(exemplars[isoA]):
                for b_idx, b in enumerate(exemplars[isoB]):
                    rij = bsv_rmsd(a, b, args.sigma_l, args.sigma_q, mode=args.mode, d=args.d)
                    wij = (weights[isoA][a_idx] / args.N_total) * (weights[isoB][b_idx] / args.N_total) * rij
                    rmsd_sum += wij
            M[i, j] = rmsd_sum

    M = rmsd_matrix_from_pairwise(M)
    np.save(args.out, M)
    print(f"Wrote weighted isoform RMSD matrix -> {args.out}")

if __name__ == "__main__":
    main()

scripts/06_tacca_screen.py (Task 2 docking triage + ranking)
# usage: python scripts/06_tacca_screen.py --pocket_bsv out/tacca_bsv.parquet --poses data/poses --out out/tacca_screen.csv
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from cypomeexplorer.io import read_pqr, read_pdb_atoms
from cypomeexplorer.vectors import hemisphere_directions_icosa
from cypomeexplorer.bsv import BSV
from cypomeexplorer.docking import ligand_vectors, binary_compatibility

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pocket_bsv", required=True, help="Parquet: id,lengths,charges for Tacca CYPs (static or exemplar).")
    ap.add_argument("--poses", required=True, help="Folder with docked ligand poses as PQR/PDB, named: {cyp_id}__{lig_id}.pqr")
    ap.add_argument("--out", required=True)
    ap.add_argument("--subdiv", type=int, default=7)
    ap.add_argument("--margin", type=float, default=0.25)
    ap.add_argument("--require_charge_complement", action="store_true")
    args = ap.parse_args()

    dirs = hemisphere_directions_icosa(subdiv=args.subdiv)
    pockets_df = pd.read_parquet(args.pocket_bsv)
    pockets = {r.id: BSV(np.array(r.lengths), np.array(r.charges)) for r in pockets_df.itertuples(index=False)}

    pose_dir = Path(args.poses)
    rows = []
    for pose in sorted(pose_dir.glob("*.pqr")) + sorted(pose_dir.glob("*.pdb")):
        stem = pose.stem
        if "__" not in stem:
            continue
        cyp_id, lig_id = stem.split("__", 1)
        if cyp_id not in pockets:
            continue

        lig_atoms = read_pqr(pose) if pose.suffix.lower() == ".pqr" else read_pdb_atoms(pose)
        lig_lens, lig_q, hit = ligand_vectors(lig_atoms, dirs)
        ok = binary_compatibility(
            pockets[cyp_id], lig_lens, lig_q, hit,
            margin=args.margin,
            require_charge_complement=args.require_charge_complement
        )
        hit_frac = float(hit.mean())
        rows.append({
            "cyp_id": cyp_id,
            "lig_id": lig_id,
            "compatible": bool(ok),
            "ligand_ray_hit_fraction": hit_frac,
        })

    out = pd.DataFrame(rows)
    out.sort_values(["compatible", "ligand_ray_hit_fraction"], ascending=[False, False], inplace=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote screening table -> {args.out}")

if __name__ == "__main__":
    main()

How you run this for Task 1 vs Task 2
Task 1 (paper replication)

Prepare aligned PQRs for the SI structure list (human and plant).

Run:

01_bsv_static.py to generate BSVs

02_bsv_tree.py to produce RMSD matrices + linkages (trees)

For MD:

export per-frame PQRs (or per-frame atom tables) → build per-isoform frames.parquet

run 03_md_ap_exemplars.py per isoform

compute weighted isoform matrix with 04_weighted_isoform_matrix.py using σl/σq appropriate to that ensemble (paper uses σl=4.183 Å, σq=0.214 e for that MD exemplar analysis context). 

binding-site-vectors-enable-map…

Task 2 (Tacca)

Align Tacca CYPs to your chosen reference frame.

Generate tacca_bsv.parquet with 01_bsv_static.py.

Dock ligands externally; export poses as {cyp_id}__{lig_id}.pqr (preferred).

Run 06_tacca_screen.py to get a first-pass compatibility table.

(Optional) upgrade to ensemble compatibility by generating BSVs per exemplar and repeating screening per exemplar, then aggregating compatibility fractions.

Notes you should not skip

The paper’s electrostatics assume GROMOS 54a8 charges. If you don’t feed PQR charges from an equivalent parametrization, your q channel won’t match their trees quantitatively (shape will still be meaningful). 

binding-site-vectors-enable-map…

Their alignment protocol includes a stable-secondary-structure alignment for MD frames. If you skip that and align all residues, you can blur binding-site comparisons. 

binding-site-vectors-enable-map…

If you paste your current folder conventions (where your Tacca PDB/PQRs and ligand poses live), I can adapt the scripts to your exact paths and naming so you can run them without edits.