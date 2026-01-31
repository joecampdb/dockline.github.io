#!/usr/bin/env python
"""
Download and prepare receptor PDB structures for the new ligands.

Caffeine targets (Blue):
- 5MZP: Adenosine A2A receptor with caffeine (chain A)
- 5N2S: Adenosine A1 receptor (chain A)

Imatinib targets (Orange):
- 2HYY: c-Abl kinase with imatinib (chain A)
- 1T46: c-KIT kinase with imatinib (chain A)

Quercetin targets (Green):
- 1E8W: PI3Kγ with quercetin (chain A)
- 2O3P: Pim1 kinase with quercetin (chain A)
"""

from __future__ import annotations

import os
import urllib.request
from Bio.PDB import PDBParser, PDBIO, Select

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(BASE_DIR, "pdb_raw")
CLEAN_DIR = os.path.join(BASE_DIR, "pdb_clean")

# Receptor definitions: (PDB ID, chain, ligand name, crystal ligand residue)
RECEPTORS = {
    # Caffeine targets
    "caffeine_A2A": {
        "pdb_id": "5MZP",
        "chain": "A",
        "label": "Adenosine A2A (5MZP)",
        "ligand": "caffeine",
        "crystal_ligand": "CFF",  # Caffeine residue name in PDB
        "color": "#3b82f6",
    },
    "caffeine_A1": {
        "pdb_id": "5N2S",
        "chain": "A",
        "label": "Adenosine A1 (5N2S)",
        "ligand": "caffeine",
        "crystal_ligand": None,  # No caffeine in this structure
        "color": "#3b82f6",
    },
    # Imatinib targets
    "imatinib_Abl": {
        "pdb_id": "2HYY",
        "chain": "A",
        "label": "c-Abl Kinase (2HYY)",
        "ligand": "imatinib",
        "crystal_ligand": "STI",  # Imatinib/STI-571 residue name
        "color": "#f97316",
    },
    "imatinib_cKIT": {
        "pdb_id": "1T46",
        "chain": "A",
        "label": "c-KIT Kinase (1T46)",
        "ligand": "imatinib",
        "crystal_ligand": "STI",
        "color": "#f97316",
    },
    # Quercetin targets
    "quercetin_PI3K": {
        "pdb_id": "1E8W",
        "chain": "A",
        "label": "PI3K gamma (1E8W)",
        "ligand": "quercetin",
        "crystal_ligand": "QUE",  # Quercetin residue name
        "color": "#22c55e",
    },
    "quercetin_Pim1": {
        "pdb_id": "2O3P",
        "chain": "A",
        "label": "Pim1 Kinase (2O3P)",
        "ligand": "quercetin",
        "crystal_ligand": "QUE",
        "color": "#22c55e",
    },
}


class ChainSelect(Select):
    """Select only specified chain and exclude water/heteroatoms."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # Keep only standard amino acids (exclude water and heteroatoms)
        hetfield = residue.id[0]
        return hetfield == " "  # Standard residue


def download_pdb(pdb_id: str, output_path: str) -> bool:
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except Exception as e:
        print(f"  Error downloading {pdb_id}: {e}")
        return False


def clean_pdb(input_path: str, output_path: str, chain_id: str) -> bool:
    """Extract single chain and remove water/heteroatoms."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", input_path)
    except Exception as e:
        print(f"  Error parsing PDB: {e}")
        return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, ChainSelect(chain_id))
    return True


def main():
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(CLEAN_DIR, exist_ok=True)

    print("=" * 60)
    print("Downloading and preparing receptor structures")
    print("=" * 60)

    for name, info in RECEPTORS.items():
        pdb_id = info["pdb_id"]
        chain = info["chain"]

        print(f"\n{name}")
        print("-" * 40)
        print(f"  PDB: {pdb_id}, Chain: {chain}")
        print(f"  Label: {info['label']}")

        # Download
        raw_path = os.path.join(RAW_DIR, f"{pdb_id}.pdb")
        if not os.path.exists(raw_path):
            print(f"  Downloading {pdb_id}...")
            if not download_pdb(pdb_id, raw_path):
                continue
        else:
            print(f"  Already downloaded: {pdb_id}.pdb")

        # Clean
        clean_path = os.path.join(CLEAN_DIR, f"{name}_clean.pdb")
        print(f"  Cleaning (chain {chain})...")
        if clean_pdb(raw_path, clean_path, chain):
            # Get file size
            size_kb = os.path.getsize(clean_path) / 1024
            print(f"  Saved: {clean_path} ({size_kb:.1f} KB)")
        else:
            print(f"  FAILED to clean {name}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
