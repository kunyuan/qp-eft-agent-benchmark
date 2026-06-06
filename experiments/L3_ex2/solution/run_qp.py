#!/usr/bin/env python3
"""
run_qp.py -- frozen-core dynamical quasiparticle bands (Level 3).

CLI:
    python3 run_qp.py --element-config <config.json> --grid <grid.csv> --out <out.csv>

Pipeline:
  1. Run ONE DFTK SCF + bands at the requested explicit k-points (dft_bands.jl),
     using the pinned setup from the config (LDA, GTH large-core, Ecut, kgrid,
     Fermi-Dirac smearing).  Per k-point we read eps_KS, eps_F, |k+G| and the
     plane-wave coefficients c_nk(G).
  2. Apply the derived frozen-core dynamical renormalization (qp_physics.py):
         z_core_{nu,k} = 1 / (1 + sum_c |g_{c,nk}|^2 / DeltaE_c^2)
         eps_QP        = eps_F + z_core (eps_KS - eps_F)
     g_{c,nk} = 2 sqrt(4pi) int u_c V_H^core psi0 r dr  (see method.md).
  3. Emit, per grid point, every occupied band (E_KS < E_F) plus the single
     lowest unoccupied band, with E_pred_eV = eps_QP - eps_F.

Parameter-free; no per-element branches; no ARPES used.  No network.
All per-element data (core_model.json, atomic_core_*.csv) read relative to the
config path.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

import qp_physics as qp

HA_EV = qp.HA_EV
HERE = os.path.dirname(os.path.abspath(__file__))


def load_core_channels(config_dir, core_model_path):
    """Read core_model.json (relative to config dir) and the atomic_core_*.csv."""
    with open(core_model_path) as f:
        cm = json.load(f)
    channels = []
    for ch in cm["core_s_channels"]:
        csv_path = os.path.join(config_dir, ch["atomic_core_file"])
        d = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        r, u = d[:, 0], d[:, 1]
        channels.append(
            dict(
                name=ch["channel"],
                r=r,
                u=u,
                dE=float(ch["DeltaE_c_Ha"]),
                # closed s-shell: 2 electrons (the monopole fluctuation amplitude).
                occ=2,
            )
        )
    return qp.build_core(channels)


def run_dft(config_path, grid_path, n_bands):
    """Shell out to dft_bands.jl; return the parsed JSON result."""
    julia = os.environ.get("JULIA_BIN", "julia")
    out_json = tempfile.NamedTemporaryFile(
        suffix=".json", delete=False, dir=tempfile.gettempdir()
    )
    out_json.close()
    script = os.path.join(HERE, "dft_bands.jl")
    env = dict(os.environ)
    # keep the shared machine calm
    env.setdefault("OPENBLAS_NUM_THREADS", "4")
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("JULIA_NUM_THREADS", "1")
    # Ensure DFTK is found.  Honour an existing JULIA_PROJECT (the pinned env in
    # the grading container); otherwise fall back to a local jenv next to this
    # script if one exists.
    if not env.get("JULIA_PROJECT"):
        local_jenv = os.path.join(HERE, "jenv")
        if os.path.isdir(local_jenv):
            env["JULIA_PROJECT"] = local_jenv
    cmd = [julia, script, config_path, grid_path, out_json.name, str(n_bands)]
    # IMPORTANT: DFTK prints a long SCF log; capturing it through a fixed-size
    # PIPE would deadlock once the buffer fills.  Stream it to a log file instead.
    log = tempfile.NamedTemporaryFile(
        suffix=".log", delete=False, dir=tempfile.gettempdir()
    )
    with open(log.name, "w") as logf:
        proc = subprocess.run(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        with open(log.name) as f:
            sys.stderr.write(f.read())
        raise RuntimeError("DFTK run failed")
    with open(out_json.name) as f:
        res = json.load(f)
    os.unlink(out_json.name)
    os.unlink(log.name)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    config_path = os.path.abspath(args.element_config)
    grid_path = os.path.abspath(args.grid)
    config_dir = os.path.dirname(config_path)

    with open(config_path) as f:
        cfg = json.load(f)
    element = cfg["element"]
    lv = np.array(cfg["structure"]["lattice_vectors_bohr"], dtype=float)
    # lattice vectors are columns; cell volume = |det|
    Omega = abs(np.linalg.det(lv.T))
    n_atom = len(cfg["structure"]["atom_positions_frac"])
    Omega_at = Omega / n_atom

    z_val = int(cfg["dft"]["z_valence"])
    # need occupied set + first unoccupied; occupied count ~ z_val * n_atom / 2 (spin)
    # ascending KS order. Compute a few extra bands above the (spin-paired) count.
    n_occ_max = max(1, int(np.ceil(z_val * n_atom / 2.0)))
    n_bands = n_occ_max + 5

    core_model_path = os.path.join(
        config_dir, cfg["frozen_core"]["core_model_file"]
    )
    channels = load_core_channels(config_dir, core_model_path)

    res = run_dft(config_path, grid_path, n_bands)
    epsF = float(res["epsF"])

    rows = []
    for kp in res["kpoints"]:
        pid = kp["point_id"]
        t = kp["t"]
        kpg = np.array(kp["kpg"], dtype=float)
        # bands in ascending energy order
        bands = kp["bands"]
        eps = np.array([b["eps"] for b in bands], dtype=float)
        order = np.argsort(eps)
        # occupied = E_KS < E_F; plus the single lowest unoccupied
        occ_idx = [i for i in order if eps[i] < epsF]
        unocc_idx = [i for i in order if eps[i] >= epsF]
        emit_idx = occ_idx + (unocc_idx[:1] if unocc_idx else [])
        if not emit_idx:  # safety: emit lowest band
            emit_idx = [order[0]]
        for i in emit_idx:
            b = bands[i]
            c_re = np.array(b["c_re"], dtype=float)
            c_im = np.array(b["c_im"], dtype=float)
            z = qp.z_core_state(channels, c_re, c_im, kpg, Omega_at)
            eps_qp = epsF + z * (eps[i] - epsF)
            E_pred_eV = (eps_qp - epsF) * HA_EV
            rows.append((element, pid, t, E_pred_eV))

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["element", "point_id", "t", "E_pred_eV"])
        for el, pid, t, e in rows:
            w.writerow([el, pid, t, f"{e:.6f}"])

    sys.stderr.write(f"wrote {len(rows)} rows to {args.out}\n")


if __name__ == "__main__":
    main()
