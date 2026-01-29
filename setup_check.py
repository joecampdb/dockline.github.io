#!/usr/bin/env python
"""
Dockline Setup Checker

Run this script to verify your installation is complete:
    python setup_check.py

This checks:
1. Required Python packages
2. CUDA/GPU availability
3. DiffDock directory and models
4. Directory structure
"""

import sys
import os
import importlib.util

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def check_mark(ok: bool) -> str:
    return "[OK]" if ok else "[MISSING]"

def check_package(name: str, import_name: str = None) -> bool:
    """Check if a Python package is installed."""
    import_name = import_name or name
    try:
        spec = importlib.util.find_spec(import_name)
        return spec is not None
    except ModuleNotFoundError:
        return False

def main():
    print("=" * 60)
    print("Dockline Setup Checker")
    print("=" * 60)
    print()

    all_ok = True

    # 1. Check Python version
    print("1. Python Version")
    py_version = sys.version_info
    py_ok = py_version >= (3, 9)
    print(f"   {check_mark(py_ok)} Python {py_version.major}.{py_version.minor}.{py_version.micro}")
    if not py_ok:
        print("      Requires Python 3.9+")
        all_ok = False
    print()

    # 2. Check required packages
    print("2. Required Packages")
    packages = [
        ("flask", "flask"),
        ("rdkit", "rdkit"),
        ("biopython", "Bio"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
    ]
    for name, import_name in packages:
        ok = check_package(name, import_name)
        print(f"   {check_mark(ok)} {name}")
        if not ok:
            all_ok = False
    print()

    # 3. Check DiffDock packages
    print("3. DiffDock Dependencies")
    diffdock_packages = [
        ("torch", "torch"),
        ("torch_geometric", "torch_geometric"),
        ("esm", "esm"),
        ("e3nn", "e3nn"),
        ("pytorch_lightning", "pytorch_lightning"),
    ]
    for name, import_name in diffdock_packages:
        ok = check_package(name, import_name)
        print(f"   {check_mark(ok)} {name}")
        if not ok:
            all_ok = False
    print()

    # 4. Check CUDA
    print("4. GPU/CUDA")
    try:
        import torch
        cuda_available = torch.cuda.is_available()
        print(f"   {check_mark(cuda_available)} CUDA available")
        if cuda_available:
            print(f"   {check_mark(True)} GPU: {torch.cuda.get_device_name(0)}")
            print(f"   {check_mark(True)} CUDA version: {torch.version.cuda}")
        else:
            print("      Warning: No GPU detected. DiffDock will be slow on CPU.")
    except ImportError:
        print(f"   {check_mark(False)} PyTorch not installed")
        all_ok = False
    print()

    # 5. Check DiffDock directory
    print("5. DiffDock Directory")
    diffdock_dir = os.path.join(BASE_DIR, "DiffDock")
    diffdock_exists = os.path.isdir(diffdock_dir)
    print(f"   {check_mark(diffdock_exists)} DiffDock/")

    if diffdock_exists:
        # Check for model files
        score_model = os.path.join(diffdock_dir, "workdir", "v1.1", "score_model", "best_ema_inference_epoch_model.pt")
        conf_model = os.path.join(diffdock_dir, "workdir", "v1.1", "confidence_model", "best_model_epoch75.pt")

        score_ok = os.path.isfile(score_model)
        conf_ok = os.path.isfile(conf_model)

        print(f"   {check_mark(score_ok)} Score model")
        print(f"   {check_mark(conf_ok)} Confidence model")

        if not score_ok or not conf_ok:
            print("      Models will be downloaded on first run (~2GB)")
    else:
        print("      Clone DiffDock: git clone https://github.com/gcorso/DiffDock.git")
        all_ok = False
    print()

    # 6. Check directory structure
    print("6. Directory Structure")
    dirs = ["results", "uploads", "pdb_clean", "pdb_raw", "static", "templates"]
    for d in dirs:
        path = os.path.join(BASE_DIR, d)
        ok = os.path.isdir(path)
        print(f"   {check_mark(ok)} {d}/")
        if not ok:
            os.makedirs(path, exist_ok=True)
            print(f"      Created {d}/")
    print()

    # 7. Check web app files
    print("7. Web Application Files")
    files = [
        "app.py",
        "metrics.py",
        "crystal.py",
        "docking_runner.py",
        "templates/index.html",
        "static/js/app.js",
        "static/css/style.css",
    ]
    for f in files:
        path = os.path.join(BASE_DIR, f)
        ok = os.path.isfile(path)
        print(f"   {check_mark(ok)} {f}")
        if not ok:
            all_ok = False
    print()

    # Summary
    print("=" * 60)
    if all_ok:
        print("Setup complete! Run 'python app.py' to start Dockline.")
    else:
        print("Some components are missing. See above for details.")
    print("=" * 60)

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
