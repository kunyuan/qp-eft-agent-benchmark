#!/bin/bash
# Harbor verifier: run the agent's /app/run_qp.py on the concealed metals, score
# against held-out ARPES (flooding guard + KS-baseline gate), and write a binary
# reward to /logs/verifier/reward.txt (1 = overall PASS).
set -uo pipefail
mkdir -p /logs/verifier /tmp/preds /tmp/runner_inputs /tmp/julia_depot
chmod -R a+rwX /tmp/preds /tmp/runner_inputs /tmp/julia_depot
chmod -R go-rwx /tests/hidden /tests/gold
chmod go-rwx /tests/test.sh /tests/score.py
chmod -R a+rX /app
export JULIA_DEPOT_PATH=/tmp/julia_depot:/opt/julia_depot
cd /app

sanitize_case() {
  local el="$1"
  local case_name="$2"
  local dst="/tmp/runner_inputs/${case_name}"
  rm -rf "$dst"
  mkdir -p "$dst"
  python3 - "$el" "$dst" <<'PY'
import json
import shutil
import sys
from pathlib import Path

el = sys.argv[1]
dst = Path(sys.argv[2])
src = Path("/tests/hidden") / el
for path in src.iterdir():
    if path.name == "arpes_reference.csv":
        continue
    target = dst / path.name
    if path.is_dir():
        shutil.copytree(path, target)
    elif path.suffix == ".json":
        data = json.loads(path.read_text())
        if isinstance(data, dict):
            data.pop("element", None)
        target.write_text(json.dumps(data, indent=2) + "\n")
    else:
        shutil.copy2(path, target)
PY
  chmod -R a+rX "$dst"
}

run_as_nobody() {
  python3 - "$@" <<'PY'
import os
import pwd
import subprocess
import sys

pw = pwd.getpwnam("nobody")
os.setgid(pw.pw_gid)
os.setuid(pw.pw_uid)
raise SystemExit(subprocess.run(sys.argv[1:]).returncode)
PY
}

fail=0
idx=0
for el in K Mg; do
  idx=$((idx + 1))
  case_name=$(printf "case_%03d" "$idx")
  sanitize_case "$el" "$case_name"
  case_out="/tmp/preds/${case_name}_out.csv"
  run_as_nobody python3 /app/run_qp.py \
    --element-config /tmp/runner_inputs/$case_name/element_config.json \
    --grid /tmp/runner_inputs/$case_name/grid.csv \
    --out "$case_out" \
    || { echo "[verifier] runner failed for $el"; fail=1; }
  if [ -f "$case_out" ]; then
    cp "$case_out" /tmp/preds/${el}_out.csv
  fi
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
exit $((1 - $(cat /logs/verifier/reward.txt)))
