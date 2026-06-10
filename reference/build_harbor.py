#!/usr/bin/env python3
"""Assemble Harbor tasks (one per difficulty level) from the built packet +
evaluator + reference oracle. Level-independent files (Dockerfile, oracle,
verifier) are identical across levels; only instruction.md, task.toml metadata,
and the packet/hidden data differ.

Usage:
  python build_harbor.py            # builds L1, L2, L3, L4
  python build_harbor.py 2          # just L2
"""
from __future__ import annotations
import shutil, stat, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

LEVELS = {
    1: dict(slug="L1-apply", doc="THEORY.md", tag="level-1",
            difficulty=("Apply the given frozen-core quasiparticle correction. The "
                        "formula AND the core form factors are provided; wire up DFTK "
                        "with the pinned setup, extract the Bloch coefficients, assemble "
                        "the state-dependent renormalization, and beat bare KS against "
                        "held-out ARPES on concealed metals. Tests correct DFT plumbing "
                        "and generalization.")),
    2: dict(slug="L2-formfactor", doc="THEORY.md", tag="level-2",
            difficulty=("Reconstruct and implement the frozen-core quasiparticle "
                        "correction. The closed-form correction is given but the core "
                        "form factor is NOT; compute it from the provided atomic core "
                        "data, assemble the renormalization on DFTK Kohn-Sham bands, and "
                        "beat bare KS against held-out ARPES on concealed metals.")),
    3: dict(slug="L3-derive", doc="SETUP.md", tag="level-3",
            difficulty=("Derive and implement the leading frozen-core quasiparticle "
                        "correction from the physical setup alone (no formula given) plus "
                        "the provided atomic core data, then beat bare KS against held-out "
                        "ARPES on concealed metals. Frontier rung: amounts to "
                        "reconstructing the paper's derivation.")),
    4: dict(slug="L4-frontier", doc="SETUP.md", tag="level-4",
            difficulty=("Open frontier: derive the frozen-core quasiparticle correction "
                        "from first principles with NO formula, NO structural ansatz, and "
                        "NO atomic data — compute your own atomic inputs, declare and "
                        "justify every approximation, and beat bare KS against held-out "
                        "ARPES on concealed metals. The published leading-order treatment "
                        "is the baseline to match or beat.")),
}

PREAMBLE = """\
# Working environment

You are in a container with Julia + DFTK 0.7.25 + PseudoPotentialData and Python
(numpy/scipy) preinstalled; the pinned Julia project is at `$JULIA_PROJECT`.
Public development data is in `./packet/Na/` and `./packet/Al/`. Write your
solution as **`run_qp.py` in the working directory (`/app`)**. The verifier will
invoke it as

```
python run_qp.py --element-config <config.json> --grid <grid.csv> --out <out.csv>
```

on concealed held-out metals (not Na/Al). Each hidden element's data files
({HIDDEN_FILES})
are placed in one directory; your code is given the config and grid paths and
should read the rest relative to the config.

---

"""

DOCKERFILE = """\
# Pinned environment for the QP-EFT benchmark: Julia + DFTK 0.7.25 + Python.
# No PackageCompiler sysimage — for a single-agent exploration harness the
# one-time precompile/JIT cost is acceptable. DFTK is version-pinned via
# Manifest.toml (gold/thresholds are calibrated to 0.7.25) and the GTH pseudo
# artifact is pre-fetched so the runtime needs no network.
FROM julia:1.12.1-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \\
        python3 python3-numpy python3-scipy ca-certificates && \\
    rm -rf /var/lib/apt/lists/*

ENV JULIA_PROJECT=/opt/jenv \\
    JULIA_DEPOT_PATH=/opt/julia_depot

WORKDIR /opt/jenv
COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' && \\
    julia -e 'using PseudoPotentialData; fam = PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"); foreach(el -> fam[el], [:Na, :Al, :K, :Mg, :Ca, :Li, :Si])' && \\
    chmod -R a+rX /opt/jenv /opt/julia_depot   # so the unprivileged `nobody` verifier user can read the baked env

WORKDIR /app
COPY packet/ /app/packet/
"""

ORACLE_RUNQP = '''\
#!/usr/bin/env python3
"""Oracle reference solution (installed by solve.sh as /app/run_qp.py).
Reads the element from the config and drives the pinned coherent QP pipeline in
gold_runner.jl (next to this file). JULIA_PROJECT is set in the image."""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--element-config", required=True)
    p.add_argument("--grid", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    element = json.loads(Path(a.element_config).read_text())["element"]
    r = subprocess.run(["julia", str(HERE / "gold_runner.jl"), element, a.grid, a.out],
                       text=True, capture_output=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout + "\\n" + r.stderr)
        return r.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''

SOLVE_SH = """\
#!/bin/bash
# Oracle: install the reference solution as the agent's /app/run_qp.py (+ Julia
# helpers). Used by Harbor to confirm the task is solvable.
set -e
cp /solution/run_qp.py /solution/gold_runner.jl /solution/atomic_hf.jl /app/
echo "oracle reference solution installed at /app/run_qp.py"
"""

TEST_SH = """\
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
  run_as_nobody python3 /app/run_qp.py \\
    --element-config /tmp/runner_inputs/$el/element_config.json \\
    --grid /tmp/runner_inputs/$el/grid.csv \\
    --out /tmp/preds/${el}_out.csv \\
    || { echo "[verifier] runner failed for $el"; fail=1; }
done

python3 /tests/score.py \\
  --pred-dir /tmp/preds --hidden /tests/hidden --gold /tests/gold \\
  --json /logs/verifier/result.json
pass=$?

if [ "$fail" -eq 0 ] && [ "$pass" -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
echo "[verifier] reward=$(cat /logs/verifier/reward.txt)"
exit $((1 - $(cat /logs/verifier/reward.txt)))
"""

# self-contained scorer (mirrors evaluator/validate_submission.py); captured up
# front because building a level deletes its dir before we copy this back in.
SCORE_PY_TEXT = (ROOT / "harbor" / "qp-eft-L2-formfactor" / "tests" / "score.py").read_text()


def task_toml(meta) -> str:
    return f'''version = "1.0"

[metadata]
author_name = "QP-EFT Benchmark"
author_email = "noreply@example.com"
difficulty_explanation = "{meta['difficulty']}"
category = "scientific-computing"
tags = ["physics", "dft", "dftk", "quasiparticle", "frozen-core", "julia", "{meta['tag']}"]

# DFTK SCF + bands on the two hidden metals (K, Mg) take a few minutes each.
[verifier]
timeout_sec = 2400.0

# Generous: exploration involves repeated DFTK runs (each fresh Julia process
# pays ~tens of seconds of load/JIT on top of the SCF).
[agent]
timeout_sec = 7200.0

# Building the image instantiates + precompiles DFTK and pre-fetches artifacts.
[environment]
build_timeout_sec = 3600.0
'''


def _exe(path: Path):
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def build_level(level: int):
    meta = LEVELS[level]
    dst = ROOT / "harbor" / f"qp-eft-{meta['slug']}"
    if dst.exists():
        shutil.rmtree(dst)
    (dst / "environment" / "packet").mkdir(parents=True)
    (dst / "solution").mkdir()
    (dst / "tests" / "hidden").mkdir(parents=True)
    (dst / "tests" / "gold").mkdir()

    lvl = ROOT / "agent_packet" / "levels" / f"L{level}"

    # instruction.md = preamble + level README + theory/setup
    hidden_files = ("`element_config.json` and `grid.csv` only -- no atomic data"
                    if level == 4 else
                    "`element_config.json`, `grid.csv`, `core_model.json`, and this level's data")
    instruction = (PREAMBLE.replace("{HIDDEN_FILES}", hidden_files)
                   + (lvl / "README.md").read_text()
                   + "\n\n---\n\n" + (lvl / meta["doc"]).read_text())
    (dst / "instruction.md").write_text(instruction)
    (dst / "task.toml").write_text(task_toml(meta))

    # environment
    (dst / "environment" / "Dockerfile").write_text(DOCKERFILE)
    for f in ("Project.toml", "Manifest.toml"):
        shutil.copy(ROOT / "environment" / f, dst / "environment" / f)
    for el in ("Na", "Al"):
        shutil.copytree(lvl / el, dst / "environment" / "packet" / el)

    # solution (oracle)
    (dst / "solution" / "run_qp.py").write_text(ORACLE_RUNQP)
    solve = dst / "solution" / "solve.sh"
    solve.write_text(SOLVE_SH); _exe(solve)
    gr = (ROOT / "reference" / "gold_runner.jl").read_text().replace(
        'include(joinpath(@__DIR__, "gold", "sodium", "atomic_hf.jl"))',
        'include(joinpath(@__DIR__, "atomic_hf.jl"))')
    (dst / "solution" / "gold_runner.jl").write_text(gr)
    shutil.copy(ROOT / "reference" / "gold" / "sodium" / "atomic_hf.jl",
                dst / "solution" / "atomic_hf.jl")

    # tests (verifier)
    test = dst / "tests" / "test.sh"
    test.write_text(TEST_SH); _exe(test)
    (dst / "tests" / "score.py").write_text(SCORE_PY_TEXT)
    for el in ("K", "Mg"):
        shutil.copytree(ROOT / "evaluator" / "hidden" / f"L{level}" / el,
                        dst / "tests" / "hidden" / el)
        shutil.copy(ROOT / "evaluator" / "gold" / f"{el}_gold.csv",
                    dst / "tests" / "gold" / f"{el}_gold.csv")
    print(f"built {dst.relative_to(ROOT)}")


if __name__ == "__main__":
    levels = [int(sys.argv[1])] if len(sys.argv) > 1 else [1, 2, 3, 4]
    for lv in levels:
        build_level(lv)
