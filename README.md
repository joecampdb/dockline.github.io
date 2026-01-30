# Dockline

A web-based GUI for [DiffDock](https://github.com/gcorso/DiffDock) molecular docking with interactive 3D visualization and quantitative binding metrics.

![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![Flask](https://img.shields.io/badge/flask-2.3+-green.svg)
![License](https://img.shields.io/badge/license-MIT-blue.svg)

## Live Demo

**[View the interactive demo](https://joecampdb.github.io/dockline.github.io)** - Browse pre-computed docking results for aripiprazole binding to serotonin and dopamine receptors.

## Features

- **Interactive 3D Viewer**: Visualize protein-ligand complexes using 3Dmol.js
- **Pose Navigation**: Browse through all ranked docking poses with a slider
- **Binding Metrics**: H-bonds, salt bridges, hydrophobic contacts, buried surface area
- **RMSD Validation**: Compare poses to crystal structures (when available)
- **Confidence Analysis**: Chart.js visualization of DiffDock confidence scores
- **Cross-Receptor Comparison**: Side-by-side metrics table for multiple targets
- **Job Submission**: Submit new docking jobs via web interface with progress tracking

## Screenshots

*Coming soon*

## Requirements

- **GPU**: NVIDIA GPU with CUDA 12.1+ support (RTX 30/40 series recommended)
- **RAM**: 16GB+ recommended
- **Storage**: ~5GB for DiffDock models and dependencies

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/dockline.git
cd dockline
```

### 2. Clone DiffDock

```bash
git clone https://github.com/gcorso/DiffDock.git
```

### 3. Create Conda Environment

```bash
# Create environment from file
conda env create -f environment.yml

# Activate
conda activate diffdock
```

**Alternative (manual setup):**

```bash
# Create base environment
conda create -n diffdock python=3.10
conda activate diffdock

# Install PyTorch with CUDA 12.1
conda install pytorch=2.1.2 pytorch-cuda=12.1 -c pytorch -c nvidia

# Install PyTorch Geometric
pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-2.1.2+cu121.html
pip install torch-geometric

# Install other dependencies
pip install -r requirements.txt
pip install fair-esm e3nn pytorch-lightning biopandas
```

### 4. Download DiffDock Models

The first time you run DiffDock, it will automatically download the required model weights (~2GB) to the `DiffDock/workdir/` directory.

Alternatively, download manually:
```bash
cd DiffDock
python -m inference --help  # Triggers model download
```

## Usage

### Start the Web Server

```bash
conda activate diffdock
python app.py
```

Open http://127.0.0.1:5000 in your browser.

### Submit a Docking Job

1. Enter a **Job Name** (e.g., `aspirin_COX2`)
2. Provide a **Protein**:
   - Enter a PDB ID (e.g., `7E2Z`) - will be downloaded from RCSB
   - Or upload a `.pdb` file
3. Provide a **Ligand**:
   - Enter a SMILES string (e.g., `CC(=O)Oc1ccccc1C(=O)O`)
   - Or upload a `.sdf`, `.mol`, or `.mol2` file
4. Adjust **Samples** (number of poses to generate, default: 10)
5. Click **Run Docking**

### Browse Results

- Select a complex from the **Browse Results** dropdown
- Use the **Pose Rank** slider to navigate through poses
- View metrics in the sidebar panels
- Click **Comparison** to see all complexes side-by-side

## Project Structure

```
dockline/
├── app.py                 # Flask web application
├── metrics.py             # Binding interaction calculations
├── crystal.py             # Crystal ligand extraction for RMSD
├── docking_runner.py      # Async job submission backend
├── templates/
│   └── index.html         # Single-page web interface
├── static/
│   ├── css/style.css      # Custom styles
│   └── js/app.js          # Frontend logic
├── DiffDock/              # DiffDock repository (clone separately)
├── results/               # Docking output (generated)
├── uploads/               # User uploads (generated)
├── pdb_clean/             # Cleaned protein structures
├── pdb_raw/               # Raw crystal structures
├── requirements.txt       # Python dependencies
└── environment.yml        # Conda environment
```

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | Main web interface |
| GET | `/api/complexes` | List discovered complexes |
| GET | `/api/complex/<name>/protein` | Protein PDB data |
| GET | `/api/complex/<name>/pose/<rank>` | Pose SDF data |
| GET | `/api/complex/<name>/metrics/<rank>` | Pose interaction metrics |
| GET | `/api/complex/<name>/ensemble` | Ensemble metrics |
| GET | `/api/comparison` | Cross-complex comparison |
| POST | `/api/dock` | Submit new docking job |
| GET | `/api/dock/<job_id>/status` | Poll job progress |

## Configuration

### Environment Variables

- `DIFFDOCK_PYTHON`: Path to Python executable for DiffDock (auto-detected by default)

### Customizing Known Complexes

Edit `KNOWN_COMPLEXES` in `app.py` to add pre-configured receptor mappings:

```python
KNOWN_COMPLEXES = {
    "my_complex": {
        "label": "My Receptor (PDB ID)",
        "pdb_id": "XXXX",
        "pdb_file": "my_receptor_clean.pdb",
    },
}
```

## Metrics Explained

- **H-bonds**: Hydrogen bonds between N/O atoms (< 3.5 Å)
- **Salt Bridges**: Ionic interactions with ASP/GLU residues (< 4.0 Å)
- **Hydrophobic Contacts**: Carbon-carbon contacts (< 4.0 Å)
- **Close Contacts**: Residues within 4.0 Å of ligand
- **Buried Surface Area**: Estimated buried surface (~15 Å² per buried atom)
- **RMSD to Crystal**: Root-mean-square deviation from co-crystallized ligand

## Troubleshooting

### "DiffDock exited with code 1"

- Ensure CUDA is properly installed: `nvidia-smi`
- Check PyTorch CUDA: `python -c "import torch; print(torch.cuda.is_available())"`
- Verify DiffDock models exist in `DiffDock/workdir/v1.1/`

### "Protein PDB not found"

- For known complexes, ensure `pdb_clean/` contains the expected files
- For user jobs, check that the protein was uploaded correctly

### NumPy version errors

```bash
pip install "numpy<2"
```

### CUDA out of memory

- Reduce the number of samples (poses) per job
- Close other GPU-intensive applications

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [DiffDock](https://github.com/gcorso/DiffDock) - Diffusion-based molecular docking
- [3Dmol.js](https://3dmol.csb.pitt.edu/) - WebGL molecular visualization
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
- [BioPython](https://biopython.org/) - Biological computation library

## Citation

If you use Dockline in your research, please cite DiffDock:

```bibtex
@article{corso2023diffdock,
  title={DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking},
  author={Corso, Gabriele and St{\"a}rk, Hannes and Jing, Bowen and Barzilay, Regina and Jaakkola, Tommi},
  journal={International Conference on Learning Representations (ICLR)},
  year={2023}
}
```
