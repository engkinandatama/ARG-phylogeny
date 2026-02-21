#!/bin/bash
# =============================================================================
# Google Cloud Shell / Codespace / WSL Setup Script
# =============================================================================
# Sets up the entire pipeline environment in one command.
#
# Usage:
#   bash setup_env.sh
#
# Estimated:
#   - Download: ~2 GB (Conda + tools)
#   - Disk after install: ~3.5 GB
#   - Time: 10-20 minutes
# =============================================================================

set -euo pipefail

echo "=============================================="
echo " ARG Phylogenomics Toolkit - Environment Setup"
echo "=============================================="

# ─── 1. Check if conda/mamba is available ─────────────────────────────
if command -v mamba &> /dev/null; then
    PKG_MGR="mamba"
elif command -v conda &> /dev/null; then
    PKG_MGR="conda"
else
    echo "[setup] Installing Miniforge (lightweight conda)..."
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
    rm Miniforge3-Linux-x86_64.sh
    export PATH="$HOME/miniforge3/bin:$PATH"
    eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
    PKG_MGR="mamba"
fi

echo "[setup] Using: $PKG_MGR"

# ─── 2. Install Nextflow ─────────────────────────────────────────────
if ! command -v nextflow &> /dev/null; then
    echo "[setup] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    mkdir -p "$HOME/bin"
    mv nextflow "$HOME/bin/"
    export PATH="$HOME/bin:$PATH"
    echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
fi
echo "[setup] Nextflow: $(nextflow -version 2>&1 | head -3 | tail -1)"

# ─── 3. Create Conda environment ─────────────────────────────────────
ENV_NAME="arg-phylo"
if ! conda env list | grep -q "$ENV_NAME"; then
    echo "[setup] Creating conda environment '$ENV_NAME'..."
    echo "[setup] This will install: MAFFT, IQ-TREE2, HyPhy, CD-HIT, TrimAl, etc."
    $PKG_MGR create -n "$ENV_NAME" -y \
        -c bioconda -c conda-forge \
        mafft iqtree cd-hit trimal hmmer blast \
        hyphy=2.5 \
        python=3.10 biopython pandas numpy scipy \
        statsmodels matplotlib plotly seaborn networkx \
        ete3 pyyaml jinja2 requests tqdm kaleido
else
    echo "[setup] Environment '$ENV_NAME' already exists"
fi

# ─── 4. Activate and verify ──────────────────────────────────────────
echo ""
echo "[setup] Activating environment..."
conda activate "$ENV_NAME" 2>/dev/null || source activate "$ENV_NAME" 2>/dev/null || true

echo ""
echo "=============================================="
echo " Verification"
echo "=============================================="
echo -n "Python:    "; python --version 2>&1
echo -n "MAFFT:     "; mafft --version 2>&1 | head -1 || echo "not found"
echo -n "IQ-TREE2:  "; iqtree2 --version 2>&1 | head -1 || echo "not found"
echo -n "HyPhy:     "; hyphy --version 2>&1 | head -1 || echo "not found"
echo -n "CD-HIT:    "; cd-hit-est 2>&1 | head -1 || echo "not found"
echo -n "TrimAl:    "; trimal --version 2>&1 | head -1 || echo "not found"
echo -n "Nextflow:  "; nextflow -version 2>&1 | head -3 | tail -1

echo ""
echo "=============================================="
echo " Setup Complete!"
echo "=============================================="
echo ""
echo "To activate:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To run mini test:"
echo "  nextflow run main.nf -params-file conf/test_mini.yaml -stub"
echo ""
echo "To run with real data:"
echo "  nextflow run main.nf -params-file conf/test_mini.yaml -profile local"
echo ""
