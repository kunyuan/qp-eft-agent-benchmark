#!/bin/bash
# Harbor verifier. Defense-in-depth against hidden-answer leaks:
#  - hand the runner a SANITIZED input dir with NO arpes_reference.csv;
#  - lock /tests/hidden + /tests/gold to root-only and run the agent's code as
#    the unprivileged `nobody`, so it cannot read the answer key even if it tries;
#  - score (flooding guard, exact band-count, KS-baseline gate) and write a
#    binary reward to /logs/verifier/reward.txt (1 = overall PASS).
set -uo pipefail
mkdir -p /logs/verifier /tmp/preds /tmp/runner_inputs /tmp/julia_depot
chmod -R a+rwX /tmp/preds /tmp/runner_inputs /tmp/julia_depot
chmod -R go-rwx /tests/hidden /tests/gold
chmod go-rwx /tests/test.sh /tests/score.py
chmod -R a+rX /app
export JULIA_DEPOT_PATH=/tmp/julia_depot:/opt/julia_depot

# copy the hidden inputs MINUS the ARPES answer key (element identity is kept;
# the runner needs it to load the pseudopotential, and Z_nuclear reveals it anyway)
sanitize_case() {
  local el="$1" dst="/tmp/runner_inputs/$1"
  rm -rf "$dst"; mkdir -p "$dst"
  for f in /tests/hidden/"$el"/*; do
    [ "$(basename "$f")" = "arpes_reference.csv" ] && continue
    cp -r "$f" "$dst"/
  done
  chmod -R a+rX "$dst"
}

run_as_nobody() {
  python3 - "$@" <<'PY'
import os, pwd, subprocess, sys
pw = pwd.getpwnam("nobody")
os.setgid(pw.pw_gid); os.setuid(pw.pw_uid)
raise SystemExit(subprocess.run(sys.argv[1:]).returncode)
PY
}

fail=0
for el in K Mg; do
  sanitize_case "$el"
  run_as_nobody python3 /app/run_qp.py \
    --element-config /tmp/runner_inputs/$el/element_config.json \
    --grid /tmp/runner_inputs/$el/grid.csv \
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
exit $((1 - $(cat /logs/verifier/reward.txt)))
