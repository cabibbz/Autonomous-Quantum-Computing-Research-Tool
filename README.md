# Quantum Explorer

Autonomous AI-driven quantum computing research using Claude Code + IBM Quantum.

Based on [chanind's autonomous SAE research](https://github.com/chanind/claude-auto-research-synthsaebench) and [Karpathy's autoresearch](https://github.com/karpathy/autoresearch) patterns.

## How it works

Claude Code runs in a loop. Each iteration is a "research sprint" — Claude reads its notes from last time, generates an idea, implements it, runs the experiment, analyzes results, and writes up a report. Then the next sprint picks up where it left off.

No predefined objectives. The agent decides what to explore.

## Current status

108 sprints completed. Research arc: entanglement archetypes → quantum error correction → phase transitions → Potts model CFT → walking regime diagnostics. See `STATE.md` for current position and `KNOWLEDGE.md` for accumulated findings.

## Files

- `TASK.md` — The agent's instructions (it can edit these)
- `STATE.md` — Current position (rewritten each sprint)
- `KNOWLEDGE.md` — Accumulated findings, organized by topic
- `CHANGELOG.md` — Sprint log (recent detailed, older compressed)
- `loop.sh` — The loop (runs sprints autonomously)
- `gpu_utils.py` — Drop-in GPU eigensolver (CuPy, auto-selects GPU/CPU)
- `db_utils.py` — SQLite measurements database
- `results.db` — Queryable database of key measurements
- `exp_NNN*.py` — Standalone experiment scripts (one per experiment)
- `sprints/` — Individual sprint reports
- `results/` — Raw experiment data (JSON)

## Quick start

```bash
# 1. Set up environment (run once)
bash setup.sh

# 2. Set IBM Quantum credentials (run once)
source ~/quantum-env/bin/activate
python -c "
from qiskit_ibm_runtime import QiskitRuntimeService
QiskitRuntimeService.save_account(token='YOUR-TOKEN', set_as_default=True)
"

# 3. Test interactively
claude
> Follow the instructions in TASK.md.

# 4. Run overnight
./loop.sh 10

# 5. Check results
cat STATE.md
ls sprints/
```

## Requirements

- Windows + WSL2 (Ubuntu)
- Claude Code (Max plan recommended)
- IBM Quantum account (free)
- Python 3.12+, Qiskit 2.x
- NVIDIA GPU + CuPy (optional, for larger exact diag)
- physics-tenpy (for DMRG)
