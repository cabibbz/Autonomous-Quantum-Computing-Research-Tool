# Quantum Explorer

Autonomous AI-driven quantum computing research using Claude Code + IBM Quantum.

Based on [chanind's autonomous SAE research](https://github.com/chanind/claude-auto-research-synthsaebench) and [Karpathy's autoresearch](https://github.com/karpathy/autoresearch) patterns.

## How it works

Claude Code runs in a loop. Each iteration is a "research sprint" — Claude reads its notes from last time, generates an idea, implements it, runs the experiment, analyzes results, and writes up a report. Then the next sprint picks up where it left off.

No predefined objectives. The agent decides what to explore.

## Files

- `TASK.md` — The agent's instructions (it can edit these)
- `CHANGELOG.md` — The agent's long-term memory (append-only)
- `explorer.py` — The experiment script (modified each sprint)
- `loop.sh` — The Ralph Wiggum loop (runs sprints autonomously)
- `setup.sh` — One-shot environment setup
- `sprints/` — Individual sprint reports
- `results/` — Raw experiment data
- `logs/` — Loop execution logs

## Quick start

```bash
# 1. Set up environment (run once)
bash setup.sh

# 2. Set IBM Quantum credentials (run once)
conda activate quantum
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
cat CHANGELOG.md
ls sprints/
```

## Requirements

- Windows + WSL2 (Ubuntu)
- Claude Code (Max plan recommended)
- IBM Quantum account (free)
- Python 3.12+, Qiskit 2.3+
