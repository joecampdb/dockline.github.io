"""
Asynchronous DiffDock job submission and progress tracking.

Manages a background thread per docking job.  The thread launches
``python -m inference`` as a subprocess inside the DiffDock/ directory
and monitors stdout for progress keywords.

Jobs are stored in an in-memory dict keyed by UUID.
"""

import os
import re
import shutil
import subprocess
import threading
import uuid
from dataclasses import dataclass, field
from typing import Optional

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DIFFDOCK_DIR = os.path.join(BASE_DIR, "DiffDock")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
UPLOADS_DIR = os.path.join(BASE_DIR, "uploads")


def _find_python_executable() -> str:
    """
    Find the Python executable for running DiffDock.

    Priority:
    1. DIFFDOCK_PYTHON environment variable (explicit override)
    2. Current Python interpreter (sys.executable) - works if running in correct env
    3. Common Conda paths as fallback
    """
    import sys

    # Check environment variable first
    env_python = os.environ.get("DIFFDOCK_PYTHON")
    if env_python and os.path.isfile(env_python):
        return env_python

    # Use current Python if we're already in the right environment
    # This works when user activates conda env before running
    current_python = sys.executable
    if current_python and os.path.isfile(current_python):
        return current_python

    # Fallback: try common Conda installation paths
    home = os.path.expanduser("~")
    candidates = [
        os.path.join(home, "anaconda3", "envs", "diffdock", "python.exe"),
        os.path.join(home, "miniconda3", "envs", "diffdock", "python.exe"),
        os.path.join(home, "mambaforge", "envs", "diffdock", "python.exe"),
        os.path.join(home, "Anaconda3", "envs", "diffdock", "python.exe"),
        os.path.join(home, "Miniconda3", "envs", "diffdock", "python.exe"),
        # Linux/Mac paths
        os.path.join(home, "anaconda3", "envs", "diffdock", "bin", "python"),
        os.path.join(home, "miniconda3", "envs", "diffdock", "bin", "python"),
    ]

    for path in candidates:
        if os.path.isfile(path):
            return path

    # Last resort: hope 'python' is on PATH and correct
    return "python"


DIFFDOCK_PYTHON = _find_python_executable()

os.makedirs(UPLOADS_DIR, exist_ok=True)

# Progress keywords spotted in DiffDock stdout/stderr
_PROGRESS_MAP = [
    (re.compile(r"ESM", re.IGNORECASE), "embeddings"),
    (re.compile(r"loading", re.IGNORECASE), "loading_model"),
    (re.compile(r"sampl", re.IGNORECASE), "diffusion"),
    (re.compile(r"confidence", re.IGNORECASE), "scoring"),
    (re.compile(r"writ|sav", re.IGNORECASE), "saving"),
]


@dataclass
class DockingJob:
    job_id: str
    complex_name: str
    status: str = "queued"          # queued | preparing | running | complete | failed
    progress_step: str = "queued"
    error: Optional[str] = None
    _thread: Optional[threading.Thread] = field(default=None, repr=False)


_jobs: dict[str, DockingJob] = {}
_lock = threading.Lock()


def submit_job(
    complex_name: str,
    ligand_smiles: Optional[str] = None,
    ligand_file_path: Optional[str] = None,
    protein_path: Optional[str] = None,
    pdb_id: Optional[str] = None,
    samples: int = 10,
) -> DockingJob:
    """
    Submit a new docking job.

    Provide either ``protein_path`` (local PDB) or ``pdb_id`` (download from
    RCSB).  Provide either ``ligand_smiles`` or ``ligand_file_path``.
    """
    job_id = uuid.uuid4().hex[:12]
    job = DockingJob(job_id=job_id, complex_name=complex_name)

    with _lock:
        _jobs[job_id] = job

    t = threading.Thread(
        target=_run_job,
        args=(job, complex_name, ligand_smiles, ligand_file_path,
              protein_path, pdb_id, samples),
        daemon=True,
    )
    job._thread = t
    t.start()
    return job


def get_job_status(job_id: str) -> Optional[dict]:
    with _lock:
        job = _jobs.get(job_id)
    if job is None:
        return None
    return {
        "job_id": job.job_id,
        "complex_name": job.complex_name,
        "status": job.status,
        "progress_step": job.progress_step,
        "error": job.error,
    }


# ---------------------------------------------------------------------------
# Background job logic
# ---------------------------------------------------------------------------

def _run_job(
    job: DockingJob,
    complex_name: str,
    ligand_smiles: Optional[str],
    ligand_file_path: Optional[str],
    protein_path: Optional[str],
    pdb_id: Optional[str],
    samples: int,
):
    try:
        job.status = "preparing"
        job.progress_step = "preparing"

        # Resolve protein
        if protein_path:
            abs_protein = os.path.abspath(protein_path)
        elif pdb_id:
            abs_protein = _download_pdb(pdb_id)
        else:
            raise ValueError("Either protein_path or pdb_id must be provided")

        if not os.path.exists(abs_protein):
            raise FileNotFoundError(f"Protein file not found: {abs_protein}")

        # Resolve ligand
        if ligand_file_path:
            ligand_desc = os.path.abspath(ligand_file_path)
        elif ligand_smiles:
            ligand_desc = ligand_smiles
        else:
            raise ValueError("Either ligand_smiles or ligand_file_path must be provided")

        # DiffDock creates a <complex_name>/ subdir inside out_dir,
        # so point out_dir at results/ (not results/<name>/).
        out_dir = os.path.abspath(RESULTS_DIR)
        os.makedirs(out_dir, exist_ok=True)

        # Build command using the diffdock env's Python directly
        cmd = [
            DIFFDOCK_PYTHON, "-m", "inference",
            "--protein_path", abs_protein,
            "--ligand_description", ligand_desc,
            "--complex_name", complex_name,
            "--out_dir", out_dir,
            "--samples_per_complex", str(samples),
        ]

        job.status = "running"
        job.progress_step = "starting"

        env = os.environ.copy()
        env["KMP_DUPLICATE_LIB_OK"] = "TRUE"

        proc = subprocess.Popen(
            cmd,
            cwd=DIFFDOCK_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )

        # Track last lines for error reporting
        last_lines = []
        for line in proc.stdout:
            _update_progress(job, line)
            last_lines.append(line.rstrip())
            if len(last_lines) > 20:
                last_lines.pop(0)

        proc.wait()

        if proc.returncode != 0:
            job.status = "failed"
            tail = "\n".join(last_lines[-5:])
            job.error = f"DiffDock exited with code {proc.returncode}\n{tail}"
        else:
            # Copy protein into results/<complex_name>/ so the app can serve it
            result_complex_dir = os.path.join(out_dir, complex_name)
            if os.path.isdir(result_complex_dir):
                protein_dest = os.path.join(result_complex_dir, "protein.pdb")
                if not os.path.exists(protein_dest):
                    shutil.copy2(abs_protein, protein_dest)

            job.status = "complete"
            job.progress_step = "complete"

    except Exception as exc:
        job.status = "failed"
        job.error = str(exc)


def _update_progress(job: DockingJob, line: str):
    """Update job progress based on stdout line content."""
    for pattern, step in _PROGRESS_MAP:
        if pattern.search(line):
            job.progress_step = step
            break


def _download_pdb(pdb_id: str) -> str:
    """Download a PDB file from RCSB into uploads/."""
    import urllib.request

    pdb_id = pdb_id.strip().upper()
    dest = os.path.join(UPLOADS_DIR, f"{pdb_id}.pdb")
    if not os.path.exists(dest):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, dest)
    return dest
