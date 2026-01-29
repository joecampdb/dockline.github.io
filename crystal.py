"""
Extract crystal-structure ligands from raw PDB files for RMSD comparison.

Maps known DiffDock complex names to the PDB entry, chain, and residue ID
of the co-crystallised aripiprazole (residue name 9SC).  The extracted
HETATM block is converted to an RDKit Mol with correct bond orders using
the aripiprazole SMILES as a template.

Only 5-HT1A (7E2Z) and 5-HT2A (7VOE) have a co-crystallised aripiprazole;
DRD2 (6CM4) does not, so RMSD-to-crystal is unavailable for that target.
"""

import os
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_PDB_DIR = os.path.join(BASE_DIR, "pdb_raw")

ARIPIPRAZOLE_SMILES = "O=C1CCc2ccc(OCCCCN3CCN(c4cccc(Cl)c4Cl)CC3)cc2N1"

# Complex name -> (pdb_filename, chain_id, residue_name)
CRYSTAL_MAP = {
    "aripiprazole_5HT1A": ("7E2Z.pdb", "R", "9SC"),
    "aripiprazole_5HT2A": ("7VOE.pdb", "A", "9SC"),
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

    pdb_file, chain_id, resname = entry
    pdb_path = os.path.join(RAW_PDB_DIR, pdb_file)
    if not os.path.exists(pdb_path):
        _cache[complex_name] = None
        return None

    try:
        block = _extract_hetatm_block(pdb_path, chain_id, resname)
        raw_mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
        if raw_mol is None:
            raise ValueError("MolFromPDBBlock returned None")

        template = Chem.MolFromSmiles(ARIPIPRAZOLE_SMILES)
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
        _cache[complex_name] = mol
        return mol
    except Exception as exc:
        print(f"[crystal] Could not extract crystal ligand for {complex_name}: {exc}")
        _cache[complex_name] = None
        return None
