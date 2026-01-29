#!/usr/bin/env python
"""
Export pre-computed docking data for GitHub Pages static demo.

Creates docs/data/ with:
- complexes.json (list of complexes)
- {complex}/protein.pdb
- {complex}/pose_{rank}.sdf
- {complex}/metrics_{rank}.json
- {complex}/ensemble.json
- comparison.json
"""

from __future__ import annotations

import json
import os
import shutil
import sys

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import metrics as metrics_mod

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
CLEAN_PDB_DIR = os.path.join(BASE_DIR, "pdb_clean")
DOCS_DIR = os.path.join(BASE_DIR, "docs")
DATA_DIR = os.path.join(DOCS_DIR, "data")

# Known complexes for the demo
DEMO_COMPLEXES = {
    "aripiprazole_5HT1A": {
        "label": "5-HT1A Receptor (7E2Z)",
        "pdb_id": "7E2Z",
        "pdb_file": "5HT1A_7E2Z_clean.pdb",
    },
    "aripiprazole_5HT2A": {
        "label": "5-HT2A Receptor (7VOE)",
        "pdb_id": "7VOE",
        "pdb_file": "5HT2A_7VOE_clean.pdb",
    },
    "aripiprazole_DRD2": {
        "label": "D2 Receptor DRD2 (6CM4)",
        "pdb_id": "6CM4",
        "pdb_file": "DRD2_6CM4_clean.pdb",
    },
}


def find_sdf(complex_name: str, rank: int) -> str | None:
    """Find SDF file for a given complex and rank."""
    import glob
    cdir = os.path.join(RESULTS_DIR, complex_name)
    pattern = os.path.join(cdir, f"rank{rank}_confidence*.sdf")
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    bare = os.path.join(cdir, f"rank{rank}.sdf")
    if os.path.exists(bare):
        return bare
    return None


def export_complex(name: str, info: dict):
    """Export all data for a single complex."""
    print(f"  Exporting {name}...")

    complex_dir = os.path.join(DATA_DIR, name)
    os.makedirs(complex_dir, exist_ok=True)

    # Copy protein PDB
    pdb_src = os.path.join(CLEAN_PDB_DIR, info["pdb_file"])
    pdb_dst = os.path.join(complex_dir, "protein.pdb")
    if os.path.exists(pdb_src):
        shutil.copy2(pdb_src, pdb_dst)
        print(f"    [OK] protein.pdb")
    else:
        print(f"    [SKIP] protein.pdb not found")

    # Copy poses and compute metrics
    num_poses = 0
    for rank in range(1, 11):
        sdf_src = find_sdf(name, rank)
        if sdf_src:
            # Copy SDF
            sdf_dst = os.path.join(complex_dir, f"pose_{rank}.sdf")
            shutil.copy2(sdf_src, sdf_dst)
            num_poses += 1

            # Compute and save metrics
            try:
                m = metrics_mod.compute_pose_metrics(name, rank)
                metrics_path = os.path.join(complex_dir, f"metrics_{rank}.json")
                with open(metrics_path, "w") as f:
                    json.dump(m, f, indent=2)
            except Exception as e:
                print(f"    [WARN] metrics for rank {rank}: {e}")

    print(f"    [OK] {num_poses} poses")

    # Compute ensemble metrics
    try:
        ensemble = metrics_mod.compute_ensemble_metrics(name)
        ensemble_path = os.path.join(complex_dir, "ensemble.json")
        with open(ensemble_path, "w") as f:
            json.dump(ensemble, f, indent=2)
        print(f"    [OK] ensemble.json")
    except Exception as e:
        print(f"    [WARN] ensemble: {e}")

    return num_poses


def main():
    print("=" * 60)
    print("Exporting static demo data for GitHub Pages")
    print("=" * 60)

    # Ensure data directory exists
    os.makedirs(DATA_DIR, exist_ok=True)

    # Build complexes list
    complexes = []
    comparison = []

    for name, info in DEMO_COMPLEXES.items():
        num_poses = export_complex(name, info)

        complexes.append({
            "name": name,
            "label": info["label"],
            "pdb_id": info["pdb_id"],
            "num_poses": num_poses,
        })

        # Add to comparison (rank-1 metrics)
        try:
            m = metrics_mod.compute_pose_metrics(name, 1)
            m["complex_name"] = name
            m["label"] = info["label"]
            comparison.append(m)
        except Exception as e:
            print(f"    [WARN] comparison for {name}: {e}")

    # Write complexes.json
    complexes_path = os.path.join(DATA_DIR, "complexes.json")
    with open(complexes_path, "w") as f:
        json.dump(complexes, f, indent=2)
    print(f"\n[OK] complexes.json")

    # Write comparison.json
    comparison_path = os.path.join(DATA_DIR, "comparison.json")
    with open(comparison_path, "w") as f:
        json.dump(comparison, f, indent=2)
    print(f"[OK] comparison.json")

    print("\n" + "=" * 60)
    print(f"Export complete! Data written to: {DATA_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
