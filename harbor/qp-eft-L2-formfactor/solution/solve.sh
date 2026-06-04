#!/bin/bash
# Oracle: install the reference solution as the agent's /app/run_qp.py (+ its
# Julia helpers). Used by Harbor to confirm the task is solvable.
set -e
cp /solution/run_qp.py /solution/gold_runner.jl /solution/atomic_hf.jl /app/
echo "oracle reference solution installed at /app/run_qp.py"
