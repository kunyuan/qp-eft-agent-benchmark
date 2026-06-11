#!/usr/bin/env python3
"""Quasiparticle band prediction from the pinned KS-LDA pipeline plus a
parameter-free dynamical frozen-core correction (see method.md).

    python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv

Element-blind: all atomic inputs (core shells, orbitals, eigenvalues, form
factors) are computed at run time by the embedded radial LDA solver from
Z_nuclear and dft.z_valence in the config. No network access.

Physics (derived in method.md):
  E_QP(nk) - E_F = (eps_nk - E_F) / (1 + lambda_nk)
  lambda_nk = sum_{core shells c} E0^2 * Q_{c,nk} / (eps_nk - eps_c)^2
  Q_{c,nk}  = n_atoms * (2l+1)/(4*pi*Omega) * sum_G |c_nk(G)|^2 * cf_c(|k+G|)^2
  cf_c(kappa) = 4*pi int j_l(kappa r) u_c(r) r dr      (core form factor)
  E0 = pi * E_h  (universal closure scale of the core-projection channel)
"""
import argparse, json, math, os, subprocess, sys, tempfile, shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from atom_solver import RadialAtom, occupations, split_core_valence  # noqa: E402

HA_EV = 27.211386245988
E0SQ = math.pi**2          # [Ha^2] universal closure scale, derived in method.md

PINNED_ENV = "/Users/kunchen/project/qp_eft_L4_fable5_run3/environment"


def julia_project():
    """Pinned Julia project (DFTK 0.7.25). Overridable via env; falls back to
    the default Julia environment if the pinned path is absent."""
    env = os.environ.get("QP_JULIA_PROJECT") or os.environ.get("JULIA_PROJECT")
    if env and os.path.exists(os.path.join(env, "Project.toml")):
        return env
    if os.path.exists(os.path.join(PINNED_ENV, "Project.toml")):
        return PINNED_ENV
    return None


def run_dftk(cfg_path, ts, nbands, outdir):
    tpath = os.path.join(outdir, "tlist.txt")
    with open(tpath, "w") as f:
        f.write("\n".join(f"{t:.10f}" for t in ts))
    proj = julia_project()
    cmd = ["julia"] + ([f"--project={proj}"] if proj else []) + [
           os.path.join(HERE, "dftk_bands.jl"), cfg_path, tpath, outdir, str(nbands)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if "DFTK_BANDS_OK" not in res.stdout:
        sys.stderr.write(res.stdout[-4000:] + "\n" + res.stderr[-4000:] + "\n")
        raise RuntimeError("DFTK band run failed")


def core_form_factors(Z, zv, kmax=14.0, nk=700):
    """Radial-LDA core shells: (l, eps_c, kappa-grid, cf(kappa))."""
    from scipy.special import spherical_jn
    at = RadialAtom(int(round(Z)))
    occs = occupations(int(round(Z)))
    res = at.scf(occs)
    core = split_core_valence(occs, zv)
    ks = np.linspace(1e-4, kmax, nk)
    shells = []
    for (n, l, occ) in core:
        u = res["u"][(n, l)]
        cf = np.array([4*np.pi*np.trapezoid(spherical_jn(l, k*at.r)*u*at.r, at.r)
                       for k in ks])
        shells.append(dict(n=n, l=l, eps=res["eps"][(n, l)], ks=ks, cf=cf))
    return shells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = json.load(open(args.element_config))
    element = cfg["element"]
    Z = cfg["Z_nuclear"]
    zv = cfg["dft"]["z_valence"]
    natoms = len(cfg["structure"]["atom_positions_frac"])

    # grid
    pts = []
    import csv as _csv
    with open(args.grid) as f:
        for row in _csv.DictReader(f):
            pts.append((row["point_id"], float(row["t"])))

    nbands = int(math.ceil(zv*natoms/2.0)) + 4

    workdir = tempfile.mkdtemp(prefix="qp_")
    try:
        # 1. pinned KS step
        cfg_abs = os.path.abspath(args.element_config)
        run_dftk(cfg_abs, [t for (_, t) in pts], nbands, workdir)
        eF = float(open(os.path.join(workdir, "ef.txt")).read().strip())
        cell = open(os.path.join(workdir, "cell.txt")).read().split()
        Omega = float(cell[0])

        # 2. self-computed atomic core inputs
        shells = core_form_factors(Z, zv)

        # 3. assemble lambda_nk and predictions
        import pandas as pd
        bands = pd.read_csv(os.path.join(workdir, "bands.csv"))
        psi = pd.read_csv(os.path.join(workdir, "psi.csv"))

        lam_map = {}
        for (ik, n), grp in psi.groupby(["ik", "n"]):
            kg = grp["kg"].values
            w = grp["w"].values
            eps = bands[(bands.ik == ik) & (bands.n == n)].eps_Ha.values[0]
            lam = 0.0
            for sh in shells:
                cfk = np.interp(kg, sh["ks"], sh["cf"])
                lam += (E0SQ * natoms * (2*sh["l"]+1)/(4*np.pi)
                        * np.sum(w*cfk**2)/Omega/(eps - sh["eps"])**2)
            lam_map[(ik, n)] = lam

        rows = []
        for ik, (pid, t) in enumerate(pts, start=1):
            sub = bands[bands.ik == ik].sort_values("eps_Ha")
            eps = sub.eps_Ha.values
            ns = sub.n.values
            nocc = int(np.sum(eps < eF))
            nemit = min(nocc + 1, len(eps))
            for j in range(nemit):
                lam = lam_map[(ik, ns[j])]
                Epred = (eps[j] - eF)/(1.0 + lam)*HA_EV
                rows.append((element, pid, t, Epred))

        with open(args.out, "w") as f:
            f.write("element,point_id,t,E_pred_eV\n")
            for (el, pid, t, E) in rows:
                f.write(f"{el},{pid},{t:.6f},{E:.6f}\n")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
