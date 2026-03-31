#!/bin/bash
# Quantum Explorer — One-Shot Setup
# Run this ONCE inside WSL2 to set up the entire environment
#
# Prerequisites:
#   - WSL2 with Ubuntu installed
#   - NVIDIA drivers installed on Windows (NOT inside WSL)
#   - Claude Code installed and authenticated
#
# Usage: bash setup.sh

set -e

echo "=== Quantum Explorer Setup ==="
echo ""

# 1. Create directory structure
echo ">>> Creating project structure..."
mkdir -p results sprints logs
git init 2>/dev/null || true

# 2. Check for conda
if ! command -v conda &> /dev/null; then
    echo ">>> Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b
    eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
    conda init
    echo "!!! Close and reopen your terminal, then run this script again."
    exit 0
fi

# 3. Create conda environment
echo ">>> Creating conda environment 'quantum'..."
conda create -n quantum python=3.12 -y 2>/dev/null || echo "Environment 'quantum' already exists"
eval "$(conda shell.bash hook)"
conda activate quantum

# 4. Install Qiskit stack
echo ">>> Installing Qiskit and dependencies..."
pip install -q qiskit qiskit-ibm-runtime qiskit-aer
pip install -q pennylane pennylane-qiskit
pip install -q matplotlib scipy pandas numpy
pip install -q pyphi 2>/dev/null || echo "Note: pyphi install may need manual attention"

# 5. Try GPU-accelerated Aer (may fail without CUDA — that's ok)
echo ">>> Attempting GPU-accelerated simulator install..."
pip install -q qiskit-aer-gpu 2>/dev/null || echo "Note: GPU Aer not available — using CPU (still works fine)"

# 6. Check IBM Quantum credentials
echo ">>> Checking IBM Quantum credentials..."
python -c "
from qiskit_ibm_runtime import QiskitRuntimeService
try:
    service = QiskitRuntimeService()
    backends = service.backends(operational=True, simulator=False)
    print(f'Connected! Available backends: {[b.name for b in backends[:5]]}...')
except Exception as e:
    print(f'Not configured yet. Run this in Python:')
    print(f'  from qiskit_ibm_runtime import QiskitRuntimeService')
    print(f'  QiskitRuntimeService.save_account(token=\"YOUR-IBM-TOKEN\", set_as_default=True)')
    print(f'Get your token from: https://quantum.cloud.ibm.com/')
"

# 7. Make loop script executable
chmod +x loop.sh

# 8. Initial git commit
git add -A
git commit -m "Initial setup — Quantum Explorer" 2>/dev/null || true

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Next steps:"
echo "  1. If IBM credentials aren't set, follow the instructions above"
echo "  2. Test it interactively first:"
echo "     cd $(pwd)"
echo "     conda activate quantum"
echo "     claude"
echo "     > Follow the instructions in TASK.md."
echo ""
echo "  3. When ready for autonomous overnight runs:"
echo "     ./loop.sh 10"
echo ""
echo "  4. Check results in the morning:"
echo "     cat CHANGELOG.md"
echo "     ls sprints/"
echo ""
