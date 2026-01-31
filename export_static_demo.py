#!/usr/bin/env python
"""
Export pre-computed docking data for GitHub Pages static demo.

Creates docs/data/ with:
- complexes.json (list of complexes with color info)
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

# Demo complexes with color coding by ligand
DEMO_COMPLEXES = {
    # Caffeine targets (Blue)
    "caffeine_A2A": {
        "label": "Adenosine A2A (5MZP)",
        "pdb_id": "5MZP",
        "pdb_file": "caffeine_A2A_clean.pdb",
        "ligand": "caffeine",
        "color": "#3b82f6",
    },
    "caffeine_A1": {
        "label": "Adenosine A1 (5N2S)",
        "pdb_id": "5N2S",
        "pdb_file": "caffeine_A1_clean.pdb",
        "ligand": "caffeine",
        "color": "#3b82f6",
    },
    # Imatinib targets (Orange)
    "imatinib_Abl": {
        "label": "c-Abl Kinase (2HYY)",
        "pdb_id": "2HYY",
        "pdb_file": "imatinib_Abl_clean.pdb",
        "ligand": "imatinib",
        "color": "#f97316",
    },
    "imatinib_cKIT": {
        "label": "c-KIT Kinase (1T46)",
        "pdb_id": "1T46",
        "pdb_file": "imatinib_cKIT_clean.pdb",
        "ligand": "imatinib",
        "color": "#f97316",
    },
    # Quercetin targets (Green)
    "quercetin_PI3K": {
        "label": "PI3K gamma (1E8W)",
        "pdb_id": "1E8W",
        "pdb_file": "quercetin_PI3K_clean.pdb",
        "ligand": "quercetin",
        "color": "#22c55e",
    },
    "quercetin_Pim1": {
        "label": "Pim1 Kinase (2O3P)",
        "pdb_id": "2O3P",
        "pdb_file": "quercetin_Pim1_clean.pdb",
        "ligand": "quercetin",
        "color": "#22c55e",
    },
}


def find_sdf(complex_name: str, rank: int) -> str | None:
    """Find SDF file for a given complex and rank."""
    import glob
    cdir = os.path.join(RESULTS_DIR, complex_name)
    if not os.path.isdir(cdir):
        return None
    pattern = os.path.join(cdir, f"rank{rank}_confidence*.sdf")
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    bare = os.path.join(cdir, f"rank{rank}.sdf")
    if os.path.exists(bare):
        return bare
    return None


def export_complex(name: str, info: dict) -> int:
    """Export all data for a single complex. Returns number of poses."""
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
        return 0

    # Check if results exist
    results_dir = os.path.join(RESULTS_DIR, name)
    if not os.path.isdir(results_dir):
        print(f"    [SKIP] No docking results for {name}")
        return 0

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
    if num_poses > 0:
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

        if num_poses > 0:
            complexes.append({
                "name": name,
                "label": info["label"],
                "pdb_id": info["pdb_id"],
                "ligand": info["ligand"],
                "color": info["color"],
                "num_poses": num_poses,
            })

            # Add to comparison (rank-1 metrics)
            try:
                m = metrics_mod.compute_pose_metrics(name, 1)
                m["complex_name"] = name
                m["label"] = info["label"]
                m["ligand"] = info["ligand"]
                m["color"] = info["color"]
                comparison.append(m)
            except Exception as e:
                print(f"    [WARN] comparison for {name}: {e}")

    # Write complexes.json
    complexes_path = os.path.join(DATA_DIR, "complexes.json")
    with open(complexes_path, "w") as f:
        json.dump(complexes, f, indent=2)
    print(f"\n[OK] complexes.json ({len(complexes)} complexes)")

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
