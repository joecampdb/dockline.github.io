"""
Extract crystal-structure ligands from raw PDB files for RMSD comparison.

Maps known DiffDock complex names to the PDB entry, chain, and residue ID
of the co-crystallised ligand. The extracted HETATM block is converted to
an RDKit Mol with correct bond orders using the ligand SMILES as a template.

Supported ligands:
- Caffeine (CFF) in 5MZP (A2A receptor)
- Imatinib (STI) in 2HYY (c-Abl) and 1T46 (c-KIT)
- Quercetin (QUE) in 1E8W (PI3K) and 2O3P (Pim1)
"""

from __future__ import annotations

import os
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_PDB_DIR = os.path.join(BASE_DIR, "pdb_raw")

# Ligand SMILES for bond order assignment
LIGAND_SMILES = {
    "caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(n4)c5cccnc5",
    "quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
}

# Complex name -> (pdb_filename, chain_id, residue_name, ligand_key)
CRYSTAL_MAP = {
    # Caffeine - A2A has caffeine co-crystallized
    "caffeine_A2A": ("5MZP.pdb", "A", "CFF", "caffeine"),
    # Caffeine - A1 does not have caffeine in 5N2S (has PSB36)
    # "caffeine_A1": None,

    # Imatinib
    "imatinib_Abl": ("2HYY.pdb", "A", "STI", "imatinib"),
    "imatinib_cKIT": ("1T46.pdb", "A", "STI", "imatinib"),

    # Quercetin
    "quercetin_PI3K": ("1E8W.pdb", "A", "QUE", "quercetin"),
    "quercetin_Pim1": ("2O3P.pdb", "A", "QUE", "quercetin"),
}

_cache: dict = {}


def _extract_hetatm_block(pdb_path: str, chain_id: str, resname: str) -> str:
    """Return the HETATM lines for a given chain and residue from a PDB file."""
    lines = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            file_chain = line[21]
            file_resname = line[17:20].strip()
            if file_chain == chain_id and file_resname == resname:
                lines.append(line)
    if not lines:
        raise ValueError(
            f"No HETATM records for chain={chain_id} resname={resname} in {pdb_path}"
        )
    lines.append("END\n")
    return "".join(lines)


def get_crystal_mol(complex_name: str):
    """
    Return an RDKit Mol for the crystal ligand, or None if unavailable.

    Results are cached for the lifetime of the process.
    """
    if complex_name in _cache:
        return _cache[complex_name]

    entry = CRYSTAL_MAP.get(complex_name)
    if entry is None:
        _cache[complex_name] = None
        return None

    pdb_file, chain_id, resname, ligand_key = entry
    pdb_path = os.path.join(RAW_PDB_DIR, pdb_file)
    if not os.path.exists(pdb_path):
        _cache[complex_name] = None
        return None

    smiles = LIGAND_SMILES.get(ligand_key)
    if smiles is None:
        _cache[complex_name] = None
        return None

    try:
        block = _extract_hetatm_block(pdb_path, chain_id, resname)
        raw_mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
        if raw_mol is None:
            raise ValueError("MolFromPDBBlock returned None")

        template = Chem.MolFromSmiles(smiles)
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
        _cache[complex_name] = mol
        return mol
    except Exception as exc:
        print(f"[crystal] Could not extract crystal ligand for {complex_name}: {exc}")
        _cache[complex_name] = None
        return None
