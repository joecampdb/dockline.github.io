"""
Dockline – Flask web application for molecular docking visualisation and job submission.

Routes:
  GET  /                              -> index.html
  GET  /api/complexes                 -> JSON list of discovered complexes
  GET  /api/complex/<name>/protein    -> PDB text for 3Dmol.js
  GET  /api/complex/<name>/pose/<rank>-> SDF text for a specific rank
  GET  /api/complex/<name>/metrics/<rank> -> per-pose interaction metrics
  GET  /api/complex/<name>/ensemble   -> ensemble-level metrics
  GET  /api/comparison                -> rank-1 metrics for all complexes
  POST /api/dock                      -> submit a new docking job
  GET  /api/dock/<job_id>/status      -> poll job progress
"""

import glob
import os
import re

from flask import Flask, jsonify, render_template, request, Response

import metrics as metrics_mod
import docking_runner

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
CLEAN_PDB_DIR = os.path.join(BASE_DIR, "pdb_clean")

# Known complexes: complex_name -> {label, pdb_id, pdb_file, ligand, color}
KNOWN_COMPLEXES = {
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

app = Flask(__name__)


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    return render_template("index.html")


# ---------------------------------------------------------------------------
# Complex discovery
# ---------------------------------------------------------------------------

@app.route("/api/complexes")
def api_complexes():
    """Return a JSON list of complexes found in results/."""
    complexes = []
    if not os.path.isdir(RESULTS_DIR):
        return jsonify(complexes)

    for name in sorted(os.listdir(RESULTS_DIR)):
        cdir = os.path.join(RESULTS_DIR, name)
        if not os.path.isdir(cdir):
            continue
        # Must contain at least rank1.sdf or rank1_confidence*.sdf
        has_rank1 = (
            os.path.exists(os.path.join(cdir, "rank1.sdf"))
            or glob.glob(os.path.join(cdir, "rank1_confidence*.sdf"))
        )
        if not has_rank1:
            continue

        known = KNOWN_COMPLEXES.get(name)
        label = known["label"] if known else name
        pdb_id = known["pdb_id"] if known else None

        # Count poses
        poses = glob.glob(os.path.join(cdir, "rank*_confidence*.sdf"))
        num_poses = len(poses) if poses else 1

        complexes.append({
            "name": name,
            "label": label,
            "pdb_id": pdb_id,
            "num_poses": num_poses,
        })

    return jsonify(complexes)


# ---------------------------------------------------------------------------
# Protein & pose serving
# ---------------------------------------------------------------------------

@app.route("/api/complex/<name>/protein")
def api_protein(name):
    """Serve the protein PDB text for 3Dmol.js."""
    pdb_path = _resolve_protein_path(name)
    if pdb_path is None or not os.path.exists(pdb_path):
        return Response("Protein PDB not found", status=404)
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        return Response(f.read(), mimetype="text/plain")


@app.route("/api/complex/<name>/pose/<int:rank>")
def api_pose(name, rank):
    """Serve the SDF text for a specific pose rank."""
    sdf_path = _find_sdf(name, rank)
    if sdf_path is None:
        return Response(f"Pose rank {rank} not found", status=404)
    with open(sdf_path, "r", encoding="utf-8", errors="replace") as f:
        return Response(f.read(), mimetype="text/plain")


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

@app.route("/api/complex/<name>/metrics/<int:rank>")
def api_metrics(name, rank):
    """Per-pose interaction metrics."""
    result = metrics_mod.compute_pose_metrics(name, rank)
    return jsonify(result)


@app.route("/api/complex/<name>/ensemble")
def api_ensemble(name):
    """Ensemble-level metrics (confidence scores + pairwise RMSD)."""
    result = metrics_mod.compute_ensemble_metrics(name)
    return jsonify(result)


@app.route("/api/comparison")
def api_comparison():
    """Rank-1 metrics for every discovered complex."""
    rows = []
    if not os.path.isdir(RESULTS_DIR):
        return jsonify(rows)

    for name in sorted(os.listdir(RESULTS_DIR)):
        cdir = os.path.join(RESULTS_DIR, name)
        if not os.path.isdir(cdir):
            continue
        has_rank1 = (
            os.path.exists(os.path.join(cdir, "rank1.sdf"))
            or glob.glob(os.path.join(cdir, "rank1_confidence*.sdf"))
        )
        if not has_rank1:
            continue

        m = metrics_mod.compute_pose_metrics(name, 1)
        known = KNOWN_COMPLEXES.get(name)
        m["complex_name"] = name
        m["label"] = known["label"] if known else name
        rows.append(m)

    return jsonify(rows)


# ---------------------------------------------------------------------------
# Docking job submission
# ---------------------------------------------------------------------------

@app.route("/api/dock", methods=["POST"])
def api_dock():
    """Submit a new docking job."""
    complex_name = request.form.get("complex_name", "").strip()
    if not complex_name:
        return jsonify({"error": "complex_name is required"}), 400

    # Protein: file upload or PDB ID
    protein_path = None
    pdb_id = request.form.get("protein_pdb_id", "").strip() or None

    protein_file = request.files.get("protein_file")
    if protein_file and protein_file.filename:
        dest = os.path.join(docking_runner.UPLOADS_DIR, protein_file.filename)
        protein_file.save(dest)
        protein_path = dest

    if not protein_path and not pdb_id:
        return jsonify({"error": "Provide protein_file or protein_pdb_id"}), 400

    # Ligand: SMILES or file upload
    ligand_smiles = request.form.get("ligand_smiles", "").strip() or None
    ligand_file_path = None
    ligand_file = request.files.get("ligand_file")
    if ligand_file and ligand_file.filename:
        dest = os.path.join(docking_runner.UPLOADS_DIR, ligand_file.filename)
        ligand_file.save(dest)
        ligand_file_path = dest

    if not ligand_smiles and not ligand_file_path:
        return jsonify({"error": "Provide ligand_smiles or ligand_file"}), 400

    samples = int(request.form.get("samples", 10))

    job = docking_runner.submit_job(
        complex_name=complex_name,
        ligand_smiles=ligand_smiles,
        ligand_file_path=ligand_file_path,
        protein_path=protein_path,
        pdb_id=pdb_id,
        samples=samples,
    )

    return jsonify({"job_id": job.job_id, "status": job.status})


@app.route("/api/dock/<job_id>/status")
def api_dock_status(job_id):
    """Poll docking job progress."""
    status = docking_runner.get_job_status(job_id)
    if status is None:
        return jsonify({"error": "Job not found"}), 404
    return jsonify(status)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_protein_path(complex_name: str):
    """Find the protein PDB file for a given complex."""
    known = KNOWN_COMPLEXES.get(complex_name)
    if known:
        return os.path.join(CLEAN_PDB_DIR, known["pdb_file"])
    # User-submitted job: protein.pdb in results dir
    candidate = os.path.join(RESULTS_DIR, complex_name, "protein.pdb")
    if os.path.exists(candidate):
        return candidate
    return None


def _find_sdf(complex_name: str, rank: int):
    """Find the SDF file for a given complex and rank."""
    cdir = os.path.join(RESULTS_DIR, complex_name)
    if not os.path.isdir(cdir):
        return None
    # Prefer file with confidence
    pattern = os.path.join(cdir, f"rank{rank}_confidence*.sdf")
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    bare = os.path.join(cdir, f"rank{rank}.sdf")
    if os.path.exists(bare):
        return bare
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app.run(debug=True, threaded=True, port=5000)
