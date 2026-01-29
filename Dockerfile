# Dockline Dockerfile
#
# Build:
#   docker build -t dockline .
#
# Run (with GPU):
#   docker run --gpus all -p 5000:5000 dockline
#
# Run (CPU only - slower):
#   docker run -p 5000:5000 dockline

FROM nvidia/cuda:12.1.0-runtime-ubuntu22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.10 \
    python3-pip \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3.10 as default
RUN ln -sf /usr/bin/python3.10 /usr/bin/python

# Set working directory
WORKDIR /app

# Copy requirements first (for caching)
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir flask rdkit biopython scipy "numpy<2"

# Install PyTorch with CUDA
RUN pip install --no-cache-dir torch==2.1.2 --index-url https://download.pytorch.org/whl/cu121

# Install PyTorch Geometric
RUN pip install --no-cache-dir \
    torch-scatter torch-sparse torch-cluster \
    -f https://data.pyg.org/whl/torch-2.1.2+cu121.html
RUN pip install --no-cache-dir torch-geometric

# Install DiffDock dependencies
RUN pip install --no-cache-dir fair-esm e3nn pytorch-lightning biopandas

# Clone DiffDock
RUN git clone https://github.com/gcorso/DiffDock.git

# Copy application code
COPY app.py metrics.py crystal.py docking_runner.py ./
COPY templates/ templates/
COPY static/ static/

# Create runtime directories
RUN mkdir -p results uploads pdb_clean pdb_raw

# Copy example data if available
COPY aripiprazole.sdf ./

# Expose port
EXPOSE 5000

# Set environment variables
ENV FLASK_APP=app.py
ENV FLASK_ENV=production
ENV KMP_DUPLICATE_LIB_OK=TRUE

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s \
    CMD wget --no-verbose --tries=1 --spider http://localhost:5000/ || exit 1

# Run the application
CMD ["python", "app.py"]
