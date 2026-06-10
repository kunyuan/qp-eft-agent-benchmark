#!/usr/bin/env python3
"""
run_qp.py — quasiparticle band predictions for simple metals from the
frozen-core dynamical correction (Level 3).

Usage:
    python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv

Pipeline:
  1. Run the pinned Kohn-Sham DFT setup (Julia + DFTK 0.7.25, LDA, GTH
     cp2k.nc.sr.lda.v0_1.largecore.gth, Ecut/kgrid/smearing from the config)
     and obtain, at every grid k-point k = t * endpoint: KS eigenvalues,
     the Fermi level, and the plane-wave coefficients of the Bloch states.
  2. For every band, compute the frozen-core dynamical renormalization
         z_core(nu,k) = 1 / (1 + w),
         w = (Z^2/z_val) sum_{c,atoms} |<phi_c| 1/r |psi~_{nu,k}>|^2 / DeltaE_c^2
     (core s-channels c from core_model.json; see method.md).
  3. Emit, per grid point, every occupied band plus the single lowest
     unoccupied band:  E_pred = z_core * (eps_KS - eps_F)  in eV.

No network access, no fitted parameters, single code path for all elements.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

SOLUTION_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SOLUTION_DIR)

from qp_correction import (HARTREE_EV, load_core_channels, w_dynamical)
from ks_io import read_dump

# Pinned Julia project (offline, pre-installed).  Can be overridden via env.
DEFAULT_JULIA_PROJECT = "/Users/kunchen/project/qp_eft_L3_fable5_run/environment"


def run_ks_dump(config_path, grid_path, prefix):
    julia_project = os.environ.get("QP_JULIA_PROJECT", DEFAULT_JULIA_PROJECT)
    script = os.path.join(SOLUTION_DIR, "ks_dump.jl")
    cmd = ["julia", f"--project={julia_project}", script,
           os.path.abspath(config_path), os.path.abspath(grid_path), prefix]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(prefix + "_meta.json"):
        sys.stderr.write(res.stdout[-4000:] + "\n" + res.stderr[-4000:] + "\n")
        raise RuntimeError("KS (DFTK) step failed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.element_config) as f:
        cfg = json.load(f)
    element = cfg["element"]
    Z = float(cfg["Z_nuclear"])
    z_val = float(cfg["dft"]["z_valence"])
    cfg_dir = os.path.dirname(os.path.abspath(args.element_config))
    core_model_path = os.path.join(
        cfg_dir, cfg.get("frozen_core", {}).get("core_model_file", "core_model.json"))
    channels = load_core_channels(core_model_path)

    workdir = tempfile.mkdtemp(prefix=f"qp_{element}_")
    prefix = os.path.join(workdir, "ks")
    run_ks_dump(args.element_config, args.grid, prefix)

    meta, kdata = read_dump(prefix)
    eF = meta["eF_Ha"]
    B = np.array(meta["recip_lattice_cols_invbohr"]).T  # columns b1,b2,b3
    vol = float(meta["unit_cell_volume_bohr3"])
    tau = np.array(meta["atom_positions_frac"], dtype=float).reshape(-1, 3)
    point_ids = meta["point_ids"]
    ts = meta["ts"]
    evals = [np.array(e) for e in meta["eigenvalues_Ha"]]

    rows = []
    for ik in range(len(point_ids)):
        kfrac = np.array(meta["kcoords_frac"][ik])
        G = kdata[ik]["G"]
        psi = kdata[ik]["psi"]                      # (nG, nb)
        q = (G + kfrac) @ B.T                       # cartesian k+G (Bohr^-1)
        eks_eV = (evals[ik] - eF) * HARTREE_EV      # ascending (DFTK ordering)
        order = np.argsort(eks_eV)
        eks_eV = eks_eV[order]
        psi = psi[:, order]
        n_occ = int(np.sum(eks_eV < 0.0))
        n_emit = min(n_occ + 1, len(eks_eV))
        w = w_dynamical(channels, Z, z_val, q, G, psi[:, :n_emit], tau, vol)
        for ib in range(n_emit):
            z_core = 1.0 / (1.0 + w[ib])
            rows.append((element, point_ids[ik], ts[ik], z_core * eks_eV[ib]))

    with open(args.out, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["element", "point_id", "t", "E_pred_eV"])
        for r in rows:
            wr.writerow([r[0], r[1], f"{r[2]:.6g}", f"{r[3]:.6f}"])
    print(f"wrote {args.out}: {len(rows)} rows, {len(point_ids)} grid points")


if __name__ == "__main__":
    main()
