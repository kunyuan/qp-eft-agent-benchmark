#!/usr/bin/env python3
"""Reference submission (gold). Thin Python entry over reference/gold_runner.jl.

Reads the element identity from --element-config and drives the pinned coherent
QP pipeline. A real agent submission would implement the physics itself; this
exists to (a) calibrate, (b) end-to-end test the evaluator.

  python run_qp.py --element-config <cfg.json> --grid <grid.csv> --out <out.csv>
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENV = HERE.parent / "environment"
GOLD = HERE / "gold_runner.jl"
ELEMENT_BY_Z = {3: "Li", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 19: "K", 20: "Ca"}


def element_from_config(path: str) -> str:
    cfg = json.loads(Path(path).read_text())
    if cfg.get("element"):
        return cfg["element"]
    return ELEMENT_BY_Z[int(cfg["Z_nuclear"])]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--element-config", required=True)
    p.add_argument("--grid", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    element = element_from_config(a.element_config)
    cmd = ["julia", f"--project={ENV}", str(GOLD), element, a.grid, a.out]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + "\n" + proc.stderr)
        return proc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
