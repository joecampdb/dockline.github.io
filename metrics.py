"""
Compute interaction metrics for DiffDock docking poses.

Per-pose metrics:
  - H-bonds (N/O atoms within 3.5 A)
  - Hydrophobic contacts (C-C pairs within 4.0 A)
  - Salt bridges (ligand N to ASP OD / GLU OE within 4.0 A)
  - Close contacts (unique residues within 4.0 A of any ligand atom)
  - Buried surface area estimate (ligand atoms with protein neighbor < 4.0 A * 15 A^2)
  - RMSD to crystal (when available, via rdMolAlign.GetBestRMS)

Ensemble metrics (all 10 poses per complex):
  - Confidence scores extracted from filenames
  - Pairwise RMSD matrix (mean / min / max)
"""

from __future__ import annotations

import os
import re
import glob
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from Bio.PDB import PDBParser

from crystal import get_crystal_mol

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
CLEAN_PDB_DIR = os.path.join(BASE_DIR, "pdb_clean")

# Known complex -> clean PDB filename
KNOWN_PROTEINS = {
    "aripiprazole_5HT1A": "5HT1A_7E2Z_clean.pdb",
    "aripiprazole_5HT2A": "5HT2A_7VOE_clean.pdb",
    "aripiprazole_DRD2": "DRD2_6CM4_clean.pdb",
}

_protein_cache: dict = {}

HBOND_ELEMENTS = {"N", "O"}
HYDROPHOBIC_ELEMENT = "C"
SALT_BRIDGE_RESIDUES = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
}
DISTANCE_HBOND = 3.5
DISTANCE_CONTACT = 4.0
BSA_PER_ATOM = 15.0  # approximate buried SA per buried ligand atom (A^2)


# ---------------------------------------------------------------------------
# Protein atom parsing
# ---------------------------------------------------------------------------

def _parse_protein_atoms(complex_name: str) -> list[dict]:
    """
    Parse protein heavy atoms from clean PDB into a list of dicts.
    Returns [{name, element, resname, resseq, chain, coord}, ...].
    Cached per complex name.
    """
    if complex_name in _protein_cache:
        return _protein_cache[complex_name]

    pdb_file = KNOWN_PROTEINS.get(complex_name)
    if pdb_file:
        pdb_path = os.path.join(CLEAN_PDB_DIR, pdb_file)
    else:
        pdb_path = os.path.join(RESULTS_DIR, complex_name, "protein.pdb")

    if not os.path.exists(pdb_path):
        _protein_cache[complex_name] = []
        return []

    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure("prot", pdb_path)

    atoms = []
    for atom in structure.get_atoms():
        elem = atom.element.strip().upper()
        if elem == "H":
            continue
        res = atom.get_parent()
        atoms.append({
            "name": atom.get_name(),
            "element": elem,
            "resname": res.get_resname(),
            "resseq": res.get_id()[1],
            "chain": res.get_parent().get_id(),
            "coord": atom.get_vector().get_array(),
        })
    _protein_cache[complex_name] = atoms
    return atoms


def _get_ligand_atoms(mol) -> list[dict]:
    """Extract heavy atom info from an RDKit Mol with 3D coordinates."""
    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        elem = atom.GetSymbol().upper()
        if elem == "H":
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append({
            "idx": atom.GetIdx(),
            "element": elem,
            "coord": np.array([pos.x, pos.y, pos.z]),
        })
    return atoms


# ---------------------------------------------------------------------------
# Distance matrix helpers
# ---------------------------------------------------------------------------

def _distance_matrix(lig_atoms: list[dict], prot_atoms: list[dict]) -> np.ndarray:
    """Return (n_lig, n_prot) distance matrix."""
    lig_coords = np.array([a["coord"] for a in lig_atoms])
    prot_coords = np.array([a["coord"] for a in prot_atoms])
    from scipy.spatial.distance import cdist
    return cdist(lig_coords, prot_coords)


# ---------------------------------------------------------------------------
# Per-pose metrics
# ---------------------------------------------------------------------------

def compute_pose_metrics(complex_name: str, rank: int) -> dict:
    """
    Compute all interaction metrics for a specific pose.

    Returns a dict with:
      confidence, rmsd_to_crystal, hbonds, salt_bridges,
      hydrophobic_contacts, close_contacts, buried_surface_area
    """
    # Find the SDF file
    sdf_path = _find_sdf(complex_name, rank)
    if sdf_path is None:
        return {"error": f"Pose rank {rank} not found for {complex_name}"}

    # Parse confidence from filename
    confidence = _parse_confidence(sdf_path)

    # Load ligand mol
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next(iter(supplier), None)
    if mol is None:
        return {"error": f"Could not read molecule from {sdf_path}"}

    lig_atoms = _get_ligand_atoms(mol)
    prot_atoms = _parse_protein_atoms(complex_name)

    if not prot_atoms:
        return {
            "confidence": confidence,
            "rmsd_to_crystal": None,
            "hbonds": [],
            "salt_bridges": [],
            "hydrophobic_contacts": 0,
            "close_contacts": [],
            "buried_surface_area": 0.0,
            "error": "Protein atoms not available",
        }

    dist = _distance_matrix(lig_atoms, prot_atoms)

    # H-bonds: ligand N/O to protein N/O within 3.5 A
    hbonds = _find_hbonds(lig_atoms, prot_atoms, dist)

    # Salt bridges: ligand N to ASP(OD1/OD2) or GLU(OE1/OE2) within 4.0 A
    salt_bridges = _find_salt_bridges(lig_atoms, prot_atoms, dist)

    # Hydrophobic contacts: C-C within 4.0 A (count unique ligand carbons)
    hydrophobic = _count_hydrophobic(lig_atoms, prot_atoms, dist)

    # Close contacts: unique residues with any atom within 4.0 A
    close = _find_close_contacts(lig_atoms, prot_atoms, dist)

    # Buried surface area estimate
    bsa = _estimate_buried_sa(lig_atoms, prot_atoms, dist)

    # RMSD to crystal
    rmsd = _compute_rmsd_to_crystal(complex_name, mol)

    return {
        "confidence": confidence,
        "rmsd_to_crystal": rmsd,
        "hbonds": hbonds,
        "salt_bridges": salt_bridges,
        "hydrophobic_contacts": hydrophobic,
        "close_contacts": close,
        "buried_surface_area": round(bsa, 1),
    }


def _find_hbonds(lig_atoms, prot_atoms, dist) -> list[dict]:
    """Find H-bond donor/acceptor pairs (N/O on both sides within 3.5 A)."""
    results = []
    seen = set()
    for i, la in enumerate(lig_atoms):
        if la["element"] not in HBOND_ELEMENTS:
            continue
        for j, pa in enumerate(prot_atoms):
            if pa["element"] not in HBOND_ELEMENTS:
                continue
            if dist[i, j] <= DISTANCE_HBOND:
                key = (pa["resname"], pa["resseq"], pa["name"])
                if key not in seen:
                    seen.add(key)
                    results.append({
                        "protein_residue": f"{pa['resname']}{pa['resseq']}",
                        "protein_atom": pa["name"],
                        "distance": round(float(dist[i, j]), 2),
                    })
    return results


def _find_salt_bridges(lig_atoms, prot_atoms, dist) -> list[dict]:
    """Find salt bridges: ligand N atoms to ASP/GLU carboxylate oxygens."""
    results = []
    seen = set()
    for i, la in enumerate(lig_atoms):
        if la["element"] != "N":
            continue
        for j, pa in enumerate(prot_atoms):
            resname = pa["resname"]
            if resname not in SALT_BRIDGE_RESIDUES:
                continue
            valid_atoms = SALT_BRIDGE_RESIDUES[resname]
            if pa["name"] not in valid_atoms:
                continue
            if dist[i, j] <= DISTANCE_CONTACT:
                key = (resname, pa["resseq"])
                if key not in seen:
                    seen.add(key)
                    results.append({
                        "residue": f"{resname}{pa['resseq']}",
                        "distance": round(float(dist[i, j]), 2),
                    })
    return results


def _count_hydrophobic(lig_atoms, prot_atoms, dist) -> int:
    """Count unique ligand carbons with a protein carbon within 4.0 A."""
    lig_c_indices = set()
    for i, la in enumerate(lig_atoms):
        if la["element"] != HYDROPHOBIC_ELEMENT:
            continue
        for j, pa in enumerate(prot_atoms):
            if pa["element"] != HYDROPHOBIC_ELEMENT:
                continue
            if dist[i, j] <= DISTANCE_CONTACT:
                lig_c_indices.add(i)
                break
    return len(lig_c_indices)


def _find_close_contacts(lig_atoms, prot_atoms, dist) -> list[str]:
    """Find unique residues with any heavy atom within 4.0 A of any ligand atom."""
    residues = set()
    for i in range(len(lig_atoms)):
        for j, pa in enumerate(prot_atoms):
            if dist[i, j] <= DISTANCE_CONTACT:
                residues.add(f"{pa['resname']}{pa['resseq']}")
    return sorted(residues)


def _estimate_buried_sa(lig_atoms, prot_atoms, dist) -> float:
    """Estimate buried surface area: count ligand atoms near protein * BSA_PER_ATOM."""
    buried = 0
    for i in range(len(lig_atoms)):
        if np.any(dist[i, :] < DISTANCE_CONTACT):
            buried += 1
    return buried * BSA_PER_ATOM


def _compute_rmsd_to_crystal(complex_name: str, pose_mol) -> float | None:
    """Compute RMSD between pose and crystal ligand. Returns None if unavailable."""
    crystal = get_crystal_mol(complex_name)
    if crystal is None:
        return None
    try:
        # Remove Hs for consistent comparison
        pose_noH = Chem.RemoveHs(pose_mol)
        crystal_noH = Chem.RemoveHs(crystal)
        rmsd = rdMolAlign.GetBestRMS(crystal_noH, pose_noH)
        return round(rmsd, 2)
    except Exception as exc:
        print(f"[metrics] RMSD calculation failed for {complex_name}: {exc}")
        return None


# ---------------------------------------------------------------------------
# Ensemble metrics
# ---------------------------------------------------------------------------

def compute_ensemble_metrics(complex_name: str) -> dict:
    """
    Compute ensemble-level metrics for all poses of a complex.

    Returns:
      confidence_scores: [{rank, confidence, filename}]
      pairwise_rmsd: {mean, min, max} (None if < 2 poses)
    """
    complex_dir = os.path.join(RESULTS_DIR, complex_name)
    if not os.path.isdir(complex_dir):
        return {"error": f"Complex directory not found: {complex_name}"}

    # Gather all ranked poses with confidence scores
    pattern = os.path.join(complex_dir, "rank*_confidence*.sdf")
    files = sorted(glob.glob(pattern), key=_rank_sort_key)

    scores = []
    mols = []
    for f in files:
        rank = _parse_rank(f)
        conf = _parse_confidence(f)
        scores.append({
            "rank": rank,
            "confidence": conf,
            "filename": os.path.basename(f),
        })
        supplier = Chem.SDMolSupplier(f, removeHs=False)
        mol = next(iter(supplier), None)
        if mol is not None:
            mols.append(mol)

    # Pairwise RMSD
    pairwise = None
    if len(mols) >= 2:
        rmsds = []
        for i in range(len(mols)):
            for j in range(i + 1, len(mols)):
                try:
                    m1 = Chem.RemoveHs(mols[i])
                    m2 = Chem.RemoveHs(mols[j])
                    r = rdMolAlign.GetBestRMS(m1, m2)
                    rmsds.append(r)
                except Exception:
                    continue
        if rmsds:
            pairwise = {
                "mean": round(float(np.mean(rmsds)), 2),
                "min": round(float(np.min(rmsds)), 2),
                "max": round(float(np.max(rmsds)), 2),
            }

    return {
        "confidence_scores": scores,
        "pairwise_rmsd": pairwise,
    }


# ---------------------------------------------------------------------------
# SDF file helpers
# ---------------------------------------------------------------------------

def _find_sdf(complex_name: str, rank: int) -> str | None:
    """Find the SDF file for a given complex and rank."""
    complex_dir = os.path.join(RESULTS_DIR, complex_name)
    if not os.path.isdir(complex_dir):
        return None

    # Try ranked file with confidence first
    pattern = os.path.join(complex_dir, f"rank{rank}_confidence*.sdf")
    matches = glob.glob(pattern)
    if matches:
        return matches[0]

    # Fall back to bare rank file (rank1.sdf)
    bare = os.path.join(complex_dir, f"rank{rank}.sdf")
    if os.path.exists(bare):
        return bare

    return None


def _parse_confidence(sdf_path: str) -> float | None:
    """Extract confidence score from filename like rank1_confidence-1.22.sdf."""
    fname = os.path.basename(sdf_path)
    match = re.search(r"confidence(-?\d+\.\d+)", fname)
    if match:
        return float(match.group(1))
    return None


def _parse_rank(sdf_path: str) -> int:
    """Extract rank number from filename."""
    fname = os.path.basename(sdf_path)
    match = re.search(r"rank(\d+)", fname)
    return int(match.group(1)) if match else 0


def _rank_sort_key(path: str) -> int:
    return _parse_rank(path)
