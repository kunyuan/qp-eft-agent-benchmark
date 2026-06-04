#!/usr/bin/env python3
"""Oracle reference solution (installed by solve.sh as /app/run_qp.py).
Reads the element from the config and drives the pinned coherent QP pipeline in
gold_runner.jl (next to this file). JULIA_PROJECT is set in the image."""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ELEMENT_BY_Z = {3: "Li", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 19: "K", 20: "Ca"}


def element_from_config(path):
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
    r = subprocess.run(["julia", str(HERE / "gold_runner.jl"), element, a.grid, a.out],
                       text=True, capture_output=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout + "\n" + r.stderr)
        return r.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
