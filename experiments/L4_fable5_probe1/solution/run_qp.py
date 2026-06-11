#!/usr/bin/env python3
"""Quasiparticle band energies for simple metals.

    python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv

Pipeline (all computed at run time, element blind):
 1. Pinned KS step (Julia/DFTK, exact setup from element_config.json) at the
    explicit path k-points; dumps eigenvalues, Fermi level and the Bloch
    plane-wave coefficients, plus the SCF-grid eigenvalues for the
    Fermi-level refill.
 2. All-electron spherical KS-LDA atomic solver (same functional set) gives
    the core orbitals and the KS excitation spectrum; the core occupation
    follows from Z_nuclear and dft.z_valence.
 3. The derived dynamical core-polarization (retardation-only) self-energy is
    evaluated on-shell for every needed Bloch state (see corepol.py and
    method.md for the derivation).
 4. The Fermi level is recomputed from the corrected spectrum; energies are
    emitted relative to it following the band-emission rules.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corepol import CoreChannels, SelfEnergy, HARTREE_EV  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def find_julia_project():
    for cand in (os.environ.get("QP_JULIA_PROJECT"),
                 os.environ.get("JULIA_PROJECT"),
                 os.path.abspath(os.path.join(HERE, "..", "environment"))):
        if cand and os.path.isfile(os.path.join(cand, "Project.toml")):
            return cand
    return None


def run_ks(config_path, grid_path, workdir):
    dump = os.path.join(workdir, "ks_dump.json")
    cmd = ["julia"]
    proj = find_julia_project()
    if proj:
        cmd.append(f"--project={proj}")
    cmd += [os.path.join(HERE, "ks_dump.jl"),
            os.path.abspath(config_path), os.path.abspath(grid_path), dump]
    print("[run_qp] running pinned KS step:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    with open(dump) as fh:
        return json.load(fh)


def fermi_refill(scf_grid, n_electrons, T, shift_of_eps):
    """Recompute the Fermi level from corrected eigenvalues on the SCF grid
    (spin-degenerate occupation 2, FD smearing T)."""
    eps, wgt = [], []
    for kp in scf_grid:
        for e in kp["eps_Ha"]:
            eps.append(e + shift_of_eps(e))
            wgt.append(kp["kweight"])
    eps, wgt = np.array(eps), np.array(wgt)

    def count(mu):
        x = np.clip((eps - mu) / max(T, 1e-12), -60, 60)
        return np.sum(wgt * 2.0 / (1.0 + np.exp(x)))

    lo, hi = eps.min() - 1.0, eps.max() + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if count(mid) < n_electrons:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-cross", action="store_true",
                    help="skip G!=G' interference terms (diagnostic)")
    args = ap.parse_args()

    with open(args.element_config) as fh:
        cfg = json.load(fh)
    elem = cfg["element"]
    Z = int(cfg["Z_nuclear"])
    zval = int(round(float(cfg["dft"]["z_valence"])))

    reuse = os.environ.get("QP_KS_DUMP")   # development shortcut only
    if reuse and os.path.isfile(reuse):
        with open(reuse) as fh:
            ks = json.load(fh)
    else:
        with tempfile.TemporaryDirectory() as workdir:
            ks = run_ks(args.element_config, args.grid, workdir)

    EF = ks["EF_Ha"]
    T = ks["temperature_Ha"]

    print(f"[run_qp] atomic solver: Z={Z}, z_valence={zval} "
          f"(core = {Z - zval} electrons)", flush=True)
    channels = CoreChannels(Z, zval)
    print(f"[run_qp] {channels.nchannels} core fluctuation channels kept; "
          f"alpha0(dipole,KS) = {channels.alpha0_dipole():.4f} a.u.", flush=True)

    se = SelfEnergy(channels, EF=EF, volume=ks["volume_bohr3"],
                    n_electrons=ks["n_electrons"], n_atoms=ks["n_atoms"],
                    positions_frac=ks["positions_frac"],
                    recip_cols=ks["recip_cols"], T=T)

    # --- corrections for every needed path state ---
    # need: every band with eps < EF plus two above (emission needs +1;
    # +2 gives the dE(eps) interpolant coverage above EF)
    records = []   # (point_id, t, band_idx, eps, dE)
    cache = {}
    for pt in ks["path"]:
        epslist = [b["eps_Ha"] for b in pt["bands"]]
        order = np.argsort(epslist)
        nocc = sum(1 for e in epslist if e < EF)
        need = order[:min(len(order), nocc + 2)]
        for n in need:
            b = pt["bands"][n]
            key = (round(pt["t"], 9), n)
            dE = se.correction(b["eps_Ha"], pt["k_frac"], b["G"],
                               b["c_re"], b["c_im"],
                               do_cross=not args.no_cross)
            cache[key] = dE
            records.append((pt["point_id"], pt["t"], n, b["eps_Ha"], dE))
        print(f"[run_qp] point {pt['point_id']:3d} t={pt['t']:.4f} done "
              f"({len(need)} states)", flush=True)

    # --- smooth dE(eps) for the Fermi refill ---
    eps_arr = np.array([r[3] for r in records])
    dE_arr = np.array([r[4] for r in records])
    srt = np.argsort(eps_arr)
    xs, ys = eps_arr[srt], dE_arr[srt]
    # window-average duplicates / near-duplicates for a monotone abscissa
    xu, yu = [], []
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] - xs[i] < 2e-3:
            j += 1
        xu.append(np.mean(xs[i:j + 1]))
        yu.append(np.mean(ys[i:j + 1]))
        i = j + 1
    xu, yu = np.array(xu), np.array(yu)

    def shift_of_eps(e):
        return np.interp(e, xu, yu)

    EF_corr = fermi_refill(ks["scf_grid"], ks["n_electrons"], T, shift_of_eps)
    print(f"[run_qp] EF(KS) = {EF:.6f} Ha ; EF(corrected) = {EF_corr:.6f} Ha "
          f"(shift {(EF_corr - EF) * HARTREE_EV:+.4f} eV)", flush=True)

    # --- emission ---
    lines = ["element,point_id,t,E_pred_eV"]
    for pt in ks["path"]:
        epslist = np.array([b["eps_Ha"] for b in pt["bands"]])
        order = np.argsort(epslist)
        nocc = sum(1 for e in epslist if e < EF)
        emit = order[:nocc + 1]   # every occupied + lowest unoccupied
        for n in emit:
            key = (round(pt["t"], 9), n)
            dE = cache.get(key)
            if dE is None:
                dE = shift_of_eps(epslist[n])
            E = (epslist[n] + dE - EF_corr) * HARTREE_EV
            lines.append(f"{elem},{pt['point_id']},{pt['t']},{E:.6f}")
    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[run_qp] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
