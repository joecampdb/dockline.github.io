#!/usr/bin/env python
"""
Generate 3D structures for caffeine, imatinib, and quercetin.
Save as SDF files for DiffDock docking.
"""

from __future__ import annotations

import os
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Ligand definitions with SMILES
LIGANDS = {
    "caffeine": {
        "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        "color": "#3b82f6",  # Blue
    },
    "imatinib": {
        "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(n4)c5cccnc5",
        "color": "#f97316",  # Orange
    },
    "quercetin": {
        "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "color": "#22c55e",  # Green
    },
}


def generate_3d_sdf(smiles: str, name: str, output_path: str) -> bool:
    """Generate 3D conformer and save as SDF."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Error: Invalid SMILES for {name}")
        return False

    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        print(f"  Warning: Could not embed {name}, trying with random coords")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    # Optimize with MMFF
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        print(f"  Warning: MMFF optimization failed for {name}")

    # Set molecule name
    mol.SetProp("_Name", name)

    # Write SDF
    writer = Chem.SDWriter(output_path)
    writer.write(mol)
    writer.close()

    return True


def main():
    print("=" * 60)
    print("Generating ligand 3D structures")
    print("=" * 60)

    for name, info in LIGANDS.items():
        print(f"\n{name.upper()}")
        print("-" * 40)

        smiles = info["smiles"]
        print(f"  SMILES: {smiles}")

        # Generate 3D structure
        output_path = os.path.join(BASE_DIR, f"{name}.sdf")
        if generate_3d_sdf(smiles, name, output_path):
            print(f"  Saved: {output_path}")
        else:
            print(f"  FAILED to generate {name}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
