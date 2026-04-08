#!/bin/bash
# Ralph Wiggum Loop — Autonomous Quantum Research
# Based on chanind's SAE research setup
#
# Usage: ./loop.sh [number_of_sprints]
# Default: 10 sprints
#
# Each sprint: Claude reads TASK.md, does a research sprint,
# logs results, and exits. Next iteration picks up from the notes.

SPRINTS=${1:-10}
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
TIMESTAMP=$(date +%Y%m%d-%H%M%S)

cd "$PROJECT_DIR"

echo "============================================"
echo "  Quantum Explorer — Ralph Wiggum Loop"
echo "  Starting $SPRINTS sprints at $(date)"
echo "  Project: $PROJECT_DIR"
echo "============================================"

for i in $(seq 1 $SPRINTS); do
    echo ""
    echo ">>> Sprint $i/$SPRINTS — $(date) <<<"
    echo ""
    
    claude -p "Follow the instructions in TASK.md." \
        --dangerously-skip-permissions \
        --verbose \
        --max-turns 30 \
        2>&1 | tee -a "logs/loop-$TIMESTAMP.log"
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -ne 0 ]; then
        echo "!!! Sprint $i exited with code $EXIT_CODE — pausing 60s before retry"
        sleep 60
    fi
    
    # Brief pause between sprints to avoid rate limits
    echo "--- Cooling down 30s before next sprint ---"
    sleep 10
done

echo ""
echo "============================================"
echo "  Loop complete: $SPRINTS sprints finished"
echo "  $(date)"
echo "  Check CHANGELOG.md for results"
echo "============================================"
