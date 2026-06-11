#!/usr/bin/env python3
"""Quasiparticle bands of simple metals from a pinned KS-LDA calculation plus a
derived dynamical frozen-core correction.

    python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv

Method (see method.md for the derivation and the consistency ledger):
  1. Pinned Kohn-Sham step (LDA / GTH largecore / settings from the config),
     run via DFTK from ks_band.jl, at the explicit grid k-points.
  2. All-electron radial LDA atom (own solver, element-blind aufbau from
     Z_nuclear and z_valence) -> core orbitals chi_a, eigenvalues eps_a, and
     the atomic valence reference eps_v.
  3. Dynamical core correction: the static pseudopotential removes the core
     levels; the pseudo Bloch state |psi~_kn> is not orthogonal to them.
     Integrating the core fields out leaves the hybridization-type self-energy
         Sigma_kn(omega) = sum_{R,a,m} |Xi <chi_a^{R,m}|psi~_kn>|^2 / (omega - eps_a),
     whose frequency dependence the static PSP omits.  Linearized at the band
     energy and anchored at E_F:
         E_QP - E_F = (eps_KS - E_F) / (1 + lambda_kn),
         lambda_kn  = Xi^2 sum_{R,a,m} |<chi_a|psi~_kn>|^2 / Delta_a^2,
         Delta_a    = eps_v(atom) - eps_a(atom)  (core -> valence-level excitation),
     with the universal coupling Xi^2 = 9*pi Ha^2 fixed by closing Stage 1 on
     Li (see method.md; the Li closed-form anchor gives Xi^2 = 28.25, 9*pi =
     28.27).  Same code path for every element; no per-element constants.
  4. Output, per grid point: every occupied band (E_KS < E_F) plus the single
     lowest unoccupied band, as E_pred_eV = E_QP - E_F.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from scipy.special import spherical_jn

HA_EV = 27.211386245988
XI2 = 9.0 * np.pi  # Ha^2, universal core-hybridization coupling (method.md)

# ----------------------------------------------------------------------------
# all-electron radial LDA (PZ81) atom, element-blind
# ----------------------------------------------------------------------------
AUFBAU = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1),
          (5, 0), (4, 2), (5, 1), (6, 0), (4, 3), (5, 2), (6, 1), (7, 0)]


def occupations(Z):
    occ, n_e = [], int(round(Z))
    for (n, l) in AUFBAU:
        if n_e <= 0:
            break
        f = min(2 * (2 * l + 1), n_e)
        occ.append((n, l, float(f)))
        n_e -= f
    return occ


def core_and_valence(Z, zval):
    """Core shells (n,l,f) and valence levels, removing z_valence electrons
    from the top of the aufbau filling."""
    occ = occupations(Z)
    rem, shells, vlev = zval, [], []
    for (n, l, f) in occ[::-1]:
        take = min(f, rem)
        if take > 0:
            vlev.append((n, l))
        if f - take > 0:
            shells.append((n, l, f - take))
        rem -= take
    return sorted(shells), vlev


def vxc_pz81(n):
    n = np.maximum(n, 1e-30)
    rs = (3.0 / (4.0 * np.pi * n)) ** (1.0 / 3.0)
    vx = -(3.0 / np.pi * n) ** (1.0 / 3.0)
    vc = np.empty_like(n)
    hi = rs >= 1.0
    lo = ~hi
    g, b1, b2 = -0.1423, 1.0529, 0.3334
    sq = np.sqrt(rs[hi])
    ec = g / (1.0 + b1 * sq + b2 * rs[hi])
    vc[hi] = ec * (1.0 + 7.0 / 6.0 * b1 * sq + 4.0 / 3.0 * b2 * rs[hi]) \
        / (1.0 + b1 * sq + b2 * rs[hi])
    A, Bc, C, D = 0.0311, -0.048, 0.0020, -0.0116
    lnr = np.log(rs[lo])
    vc[lo] = A * lnr + (Bc - A / 3.0) + 2.0 / 3.0 * C * rs[lo] * lnr \
        + (2.0 * D - C) / 3.0 * rs[lo]
    return vx + vc


def hartree(r, rho4, h):
    q_in = np.concatenate(([0.0], np.cumsum((rho4[:-1] + rho4[1:]) * 0.5 * h)))
    f = rho4 / r
    tot = np.concatenate(([0.0], np.cumsum((f[:-1] + f[1:]) * 0.5 * h)))
    return q_in / r + (tot[-1] - tot)


def solve_atom(Z, h=0.001, rmax=35.0, mixing=0.35, tol=1e-7, maxiter=300):
    Z = float(Z)
    N = int(rmax / h) - 1
    r = h * np.arange(1, N + 1)
    occ = occupations(Z)
    lmax = max(l for (_, l, _) in occ)
    rho4 = 4.0 * Z ** 3 * r ** 2 * np.exp(-2.0 * Z * r)
    rho4 *= Z / np.trapezoid(rho4, r)
    Vold, eps_out, u_out, drho = None, {}, {}, 1.0
    for it in range(maxiter):
        VH = hartree(r, rho4, h)
        V = -Z / r + VH + vxc_pz81(rho4 / (4.0 * np.pi * r ** 2))
        if Vold is not None:
            V = mixing * V + (1.0 - mixing) * Vold
        rho_new = np.zeros_like(rho4)
        eps_new, u_new = {}, {}
        for l in range(lmax + 1):
            states_l = [(n, f) for (n, ll, f) in occ if ll == l]
            if not states_l:
                continue
            nmaxl = max(n for (n, _) in states_l) - l
            d = 1.0 / h ** 2 + V + l * (l + 1) / (2.0 * r ** 2)
            e = -0.5 / h ** 2 * np.ones(N - 1)
            w, v = eigh_tridiagonal(d, e, select='i', select_range=(0, nmaxl - 1))
            for (n, f) in states_l:
                u = v[:, n - l - 1] / np.sqrt(h)
                rho_new += f * u ** 2
                eps_new[(n, l)] = w[n - l - 1]
                u_new[(n, l)] = u
        drho = np.trapezoid(np.abs(rho_new - rho4), r)
        rho4, Vold, eps_out, u_out = rho_new, V, eps_new, u_new
        if drho < tol:
            break
    if drho >= tol:
        print(f"WARNING: atomic SCF not fully converged (drho={drho:.2e})",
              file=sys.stderr)
    return dict(Z=Z, r=r, h=h, eps=eps_out, u=u_out)


# ----------------------------------------------------------------------------
# KS step via DFTK
# ----------------------------------------------------------------------------
def julia_project():
    for var in ("QP_JULIA_PROJECT", "JULIA_PROJECT"):
        if os.environ.get(var):
            return os.environ[var]
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "environment"),
                 os.path.join(os.path.dirname(here), "environment")):
        if os.path.isfile(os.path.join(cand, "Project.toml")):
            return cand
    return None


def run_ks(config_path, grid_path, out_json):
    here = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(here, "ks_band.jl")
    cmd = ["julia"]
    proj = julia_project()
    if proj:
        cmd.append(f"--project={proj}")
    cmd += [script, config_path, grid_path, out_json]
    print("[run_qp] running KS step:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


# ----------------------------------------------------------------------------
# core-overlap lambda
# ----------------------------------------------------------------------------
def legendre(l, x):
    if l == 0:
        return np.ones_like(x)
    if l == 1:
        return x
    if l == 2:
        return 1.5 * x ** 2 - 0.5
    if l == 3:
        return 2.5 * x ** 3 - 1.5 * x
    pm2, pm1 = np.ones_like(x), x
    for ll in range(2, l + 1):
        pm2, pm1 = pm1, ((2 * ll - 1) * x * pm1 - (ll - 1) * pm2) / ll
    return pm1


def compute(config_path, grid_path, out_path):
    cfg = json.load(open(config_path))
    element = cfg["element"]
    Z = int(cfg["Z_nuclear"])
    zval = int(cfg["dft"]["z_valence"])

    with tempfile.TemporaryDirectory() as td:
        ksj = os.path.join(td, "ks.json")
        run_ks(config_path, grid_path, ksj)
        d = json.load(open(ksj))

    EF = d["EF_Ha"]
    B = np.array(d["recip_lattice"]).T          # columns b1,b2,b3
    A = np.array(d["lattice"]).T                # columns a1,a2,a3
    Om = d["unit_cell_volume"]
    pos_frac = np.array(d["atom_positions_frac"], dtype=float).reshape(-1, 3)
    Rsites = pos_frac @ A.T                     # cartesian site positions

    atom = solve_atom(Z)
    shells, vlev = core_and_valence(Z, zval)
    eps_v = max(atom["eps"][(n, l)] for (n, l) in vlev)
    r = atom["r"]
    print(f"[run_qp] {element}: core shells {shells}; "
          f"eps_a = { {(n,l): round(atom['eps'][(n,l)],4) for (n,l,_) in shells} }; "
          f"eps_v = {eps_v:.4f} Ha", flush=True)

    qmax = 0.0
    for kp in d["kpoints"]:
        G = np.array(kp["G"], dtype=float)
        kf = np.array(kp["kfrac"])
        qmax = max(qmax, np.linalg.norm((G + kf) @ B.T, axis=1).max())
    qs = np.linspace(1e-4, qmax * 1.02, 600)
    fts = {}
    for (n, l, f) in shells:
        u = atom["u"][(n, l)]
        fts[(n, l)] = CubicSpline(
            qs, [np.trapezoid(u * spherical_jn(l, q * r) * r, r) for q in qs])

    rows = []
    for kp in d["kpoints"]:
        t, pid = kp["t"], kp["point_id"]
        kf = np.array(kp["kfrac"])
        G = np.array(kp["G"], dtype=float)
        kpG = (G + kf) @ B.T
        q = np.linalg.norm(kpG, axis=1)
        qh = np.where(q[:, None] > 1e-12, kpG / np.maximum(q, 1e-12)[:, None], 0.0)
        # Legendre matrices only for l-channels actually coupled (s-wave
        # selection rule => none today; kept general for completeness)
        cosm = {}
        lset = sorted(set(l for (_, l, _) in shells if l == 0))
        if lset and max(lset) > 0:
            cth = qh @ qh.T
            for l in lset:
                if l > 0:
                    cosm[l] = legendre(l, cth)
        phases = [np.exp(1j * (kpG @ Rs)) for Rs in Rsites]

        bands = sorted(kp["bands"], key=lambda b: b["eps_Ha"])
        evs = [b["eps_Ha"] for b in bands]
        # occupied + first unoccupied
        nocc = sum(1 for e in evs if e < EF)
        nsel = min(nocc + 1, len(evs))
        if nocc + 1 > len(evs):
            print(f"WARNING: not enough bands at point {pid}", file=sys.stderr)
        for b in bands[:nsel]:
            eps = b["eps_Ha"]
            c = np.array(b["c_re"]) + 1j * np.array(b["c_im"])
            lam = 0.0
            for (n, l, f) in shells:
                if l != 0:
                    # zero-range (point-core) coupling scatters s-waves only:
                    # higher-l core shells enter at the next order in
                    # (core radius x valence momentum) and are dropped
                    # (method.md, ledger L1.2b).
                    continue
                Rt = fts[(n, l)](q)
                D = eps_v - atom["eps"][(n, l)]
                for ph in phases:
                    w = c * ph * Rt
                    if l == 0:
                        s2 = (4.0 * np.pi / Om) * abs(np.sum(w)) ** 2
                    else:
                        s2 = (4.0 * np.pi * (2 * l + 1) / Om) * \
                            float(np.real(np.conj(w) @ (cosm[l] @ w)))
                    lam += XI2 * max(s2, 0.0) / D ** 2
            E_pred = (eps - EF) / (1.0 + lam) * HA_EV
            rows.append((pid, t, E_pred))

    with open(out_path, "w") as fo:
        fo.write("element,point_id,t,E_pred_eV\n")
        for (pid, t, E) in rows:
            fo.write(f"{element},{pid},{t},{E:.6f}\n")
    print(f"[run_qp] wrote {out_path} ({len(rows)} rows)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    compute(args.element_config, args.grid, args.out)


if __name__ == "__main__":
    main()
