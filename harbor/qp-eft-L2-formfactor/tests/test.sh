#!/bin/bash
# Harbor verifier: run the agent's /app/run_qp.py on the concealed metals,
# score against held-out ARPES (flooding guard + KS-baseline gate), and write
# a binary reward to /logs/verifier/reward.txt (1 = overall PASS).
set -uo pipefail
mkdir -p /logs/verifier /tmp/preds

fail=0
for el in K Mg; do
  python3 /app/run_qp.py \
    --element-config /tests/hidden/$el/element_config.json \
    --grid /tests/hidden/$el/grid.csv \
    --out /tmp/preds/${el}_out.csv \
    || { echo "[verifier] runner failed for $el"; fail=1; }
done

python3 /tests/score.py \
  --pred-dir /tmp/preds --hidden /tests/hidden --gold /tests/gold \
  --json /logs/verifier/result.json
pass=$?

if [ "$fail" -eq 0 ] && [ "$pass" -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
echo "[verifier] reward=$(cat /logs/verifier/reward.txt)"
