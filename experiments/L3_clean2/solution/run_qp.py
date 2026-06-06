#!/usr/bin/env python3
"""
run_qp.py -- Frozen-core dynamical quasiparticle correction to LDA Kohn-Sham bands.

CLI:
    python run_qp.py --element-config <config.json> --grid <grid.csv> --out <out.csv>

Physics (derived in method.md, summary here)
--------------------------------------------
A frozen-core pseudopotential integrates the core out STATICALLY: the LDA/GTH
Kohn-Sham band is that static reduction.  The leading piece it discards is the
core's *dynamical* response -- a valence electron virtually excites a core
particle-hole pair (energy DeltaE_c) and the core relaxes back.  This is the
second-order, core-induced self-energy.  Integrating the high-energy core mode
out by closure at the single scale DeltaE_c, the consequence for the valence
quasiparticle is a reduction of its spectral weight,

      z = 1 - lambda,        lambda = sum_c  <u_c^2 V_H_c> / DeltaE_c^2 ,

where, per core s-channel c,
   <u_c^2 V_H_c> = INT_0^inf u_c(r)^2 V_H_c(r) dr   (the core orbital's Hartree
                   self-energy, in Hartree; u_c normalized to INT u_c^2 dr = 1),
   DeltaE_c      = the core excitation energy (Hartree).
The weight renormalizes the band energy measured FROM the Fermi level:

      E_QP - E_F = z * (E_KS - E_F).

z is a single number per element (a local, k-independent dynamical self-energy,
the eDMFT picture): the core fluctuation is on-site, so it acts equally on every
Bloch state -- it narrows the occupied bandwidth without moving E_F.  The
correction is parameter-free: every quantity comes from the atomic core data
(u_c, V_H_c, DeltaE_c) and the converged KS run; there is no fitted constant and
no per-element branch.

Inputs read relative to the config path:
   core_model.json     -> per-channel DeltaE_c and the atomic_core_<c>.csv name
   atomic_core_<c>.csv  -> columns r_bohr, u_c, V_H_c
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HARTREE_EV = 27.211386245988


def trapz(y, x):
    return float(np.trapezoid(y, x)) if hasattr(np, "trapezoid") else float(np.trapz(y, x))


def compute_lambda(config_path):
    """lambda = sum_c <u_c^2 V_H_c> / DeltaE_c^2  (dimensionless, parameter-free)."""
    cfgdir = os.path.dirname(os.path.abspath(config_path))
    cfg = json.load(open(config_path))
    cm_name = cfg["frozen_core"]["core_model_file"]
    cm = json.load(open(os.path.join(cfgdir, cm_name)))

    lam = 0.0
    per_channel = {}
    for ch in cm["core_s_channels"]:
        dE = float(ch["DeltaE_c_Ha"])
        data = np.loadtxt(os.path.join(cfgdir, ch["atomic_core_file"]),
                          delimiter=",", skiprows=1)
        r = data[:, 0]
        u = data[:, 1]       # u_c(r) = r R_c(r), INT u^2 dr = 1
        vh = data[:, 2]      # single-core-orbital Hartree potential (Hartree)
        # core orbital Hartree self-energy <u_c^2 V_H_c> (Hartree)
        S = trapz(u * u * vh, r)
        lam_c = S / (dE * dE)
        per_channel[ch["channel"]] = (S, dE, lam_c)
        lam += lam_c
    return lam, per_channel


def compute_ks_bands(config_path, grid_path):
    """Run one DFTK SCF + bands at the grid k-points; return list of (point_id, t, bands_eV)."""
    here = os.path.dirname(os.path.abspath(__file__))
    helper = os.path.join(here, "_ks_bands.jl")

    # Locate the pinned Julia project: $JULIA_PROJECT, else a local jenv/.
    jproj = os.environ.get("JULIA_PROJECT", "")
    if not jproj:
        for cand in (os.path.join(here, "jenv"), here):
            if os.path.exists(os.path.join(cand, "Project.toml")):
                jproj = cand
                break

    env = dict(os.environ)
    # Be a good citizen on a shared machine.
    env.setdefault("OPENBLAS_NUM_THREADS", "4")
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("JULIA_NUM_THREADS", "1")

    with tempfile.NamedTemporaryFile(suffix="_ks.csv", delete=False) as tf:
        ks_out = tf.name
    try:
        cmd = ["julia"]
        if jproj:
            cmd.append(f"--project={jproj}")
        cmd += [helper, config_path, grid_path, ks_out]
        subprocess.run(cmd, env=env, check=True)

        rows = []
        with open(ks_out) as fh:
            for row in csv.DictReader(fh):
                bands = [float(x) for x in row["bands_eV"].split(";") if x != ""]
                rows.append((int(row["point_id"]), float(row["t"]), bands))
        return rows
    finally:
        try:
            os.remove(ks_out)
        except OSError:
            pass


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = json.load(open(args.element_config))
    element = cfg["element"]

    lam, per_channel = compute_lambda(args.element_config)
    z = 1.0 - lam

    sys.stderr.write(f"[run_qp] element={element}  lambda={lam:.5f}  z={z:.5f}\n")
    for chn, (S, dE, lc) in per_channel.items():
        sys.stderr.write(f"[run_qp]   {chn}: <u^2 V_H>={S:.4f} Ha  DeltaE={dE:.4f} Ha"
                         f"  lambda_c={lc:.5f}\n")

    ks_rows = compute_ks_bands(args.element_config, args.grid)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "point_id", "t", "E_pred_eV"])
        for point_id, t, bands in ks_rows:
            for e_ks in bands:                       # already E_KS - E_F (eV), ascending
                e_qp = z * e_ks                      # E_QP - E_F = z (E_KS - E_F)
                w.writerow([element, point_id, f"{t:.6f}", f"{e_qp:.6f}"])

    sys.stderr.write(f"[run_qp] wrote {args.out}\n")


if __name__ == "__main__":
    main()
