#!/usr/bin/env python3
"""
run_qp.py -- Level-3 frozen-core dynamical quasiparticle correction.

Usage:
    python run_qp.py --element-config <config.json> --grid <grid.csv> --out <out.csv>

Derivation (see method.md): the large-core pseudopotential already contains the
STATIC core-valence interaction (the core Hartree potential V_H_c is inside the
PSP). What it omits is the DYNAMICAL response of the core: a valence electron
virtually excites a core s-shell (c -> c', cost DeltaE_c) and the core relaxes
back. This is a second-order self-energy in the valence-core Coulomb coupling.

Closure (Unsold) over the excited core manifold at the single scale DeltaE_c
turns the c'-sum into a variance of the Coulomb field of one core electron seen
at the valence position r:

    F_c(r) = <c| 1/|r-r'|^2 |c> - ( <c| 1/|r-r'| |c> )^2          (>= 0)
           = W2_c(r) - V_H_c(r)^2,

    W2_c(r) = int u_c(r')^2 (1/(2 r r')) ln|(r+r')/(r-r')| dr'.

The dynamical self-energy couples the valence electron to the core fluctuation
F_c through the valence resolvent (intermediate plane-wave states). Because the
electron gas is Galilean invariant over the Fermi ball, the energy- and
momentum-derivatives of Sigma are locked (Ward identity), so the only effect on
the band measured from E_F is a POLE RESCALING (mass renormalization):

    E_QP(nk) - E_F = Z_nk * (eps_nk - E_F),   Z_nk = 1 / (1 + lambda_nk) <= 1,

with the dimensionless frozen-core coupling (see method.md for the coefficient)

    lambda_nk = (1/pi) sum_c N_c * W_c(nk) / DeltaE_c^2 ,
    W_c(nk)   = (Omega/n_atom) * <psi_nk| F_c |psi_nk>  (per-atom on-site
                fluctuation strength; Omega = cell volume, n_atom = atoms/cell).

N_c = 2 (closed s-shell). The 1/DeltaE_c^2 is the scale-separation (retardation)
suppression; 1/pi is the closure phase-space coefficient of the second-order
self-energy. W_c is a near-universal atomic quantity (~0.8 Ha), so the element
trend (strong narrowing in Na, weak in Al) is carried by 1/DeltaE_c^2.
Parameter-free: Omega, n_atom, F_c, DeltaE_c are all fixed by the inputs.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HARTREE_TO_EV = 27.211386245988
HERE = os.path.dirname(os.path.abspath(__file__))


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# F_c(r): variance of the Coulomb field of one core electron at valence point r.
# ---------------------------------------------------------------------------
def compute_Fc_on_grid(r_core, u_core, V_H_core, r_out):
    """
    F_c(r_out) = W2_c(r_out) - V_H_c(r_out)^2.
    W2_c(r) = int_0^inf u_c(r')^2 * (1/(2 r r')) ln|(r+r')/(r-r')| dr'.
    """
    rp = r_core
    u2 = u_core ** 2
    F = np.empty_like(r_out)
    W2 = np.empty_like(r_out)
    for i, r in enumerate(r_out):
        # angular-averaged 1/|r-r'|^2 kernel
        denom = 2.0 * r * rp
        ratio = np.abs((r + rp) / (r - rp + 1e-300))
        kern = np.log(ratio) / denom
        # remove the integrable log singularity at r'=r by limiting: as r'->r,
        # ln|(r+r')/(r-r')| -> diverges but integrable; trapezoid on fine grid ok.
        kern = np.where(np.isfinite(kern), kern, 0.0)
        W2[i] = np.trapezoid(u2 * kern, rp)
    # V_H_c(r_out) interpolated from the provided V_H_c
    VH = np.interp(r_out, r_core, V_H_core, left=V_H_core[0], right=None)
    # beyond tabulated range V_H ~ 1/r (single electron); enforce tail
    tail = r_out > r_core[-1]
    VH = np.where(tail, 1.0 / r_out, VH)
    F = W2 - VH ** 2
    F = np.maximum(F, 0.0)  # variance is non-negative; clip numerical noise
    return F


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-config", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg_path = os.path.abspath(args.element_config)
    cfg_dir = os.path.dirname(cfg_path)
    with open(cfg_path) as f:
        cfg = json.load(f)

    element = cfg["element"]
    struct = cfg["structure"]
    dft = cfg["dft"]
    endpoint = [float(x) for x in struct["path"]["endpoint_frac"]]
    z_val = int(dft["z_valence"])
    n_atoms = len(struct["atom_positions_frac"])

    # grid (t values)
    grid = np.genfromtxt(args.grid, delimiter=",", names=True)
    ts = np.atleast_1d(grid["t"]).astype(float)
    point_ids = np.atleast_1d(grid["point_id"]).astype(int)

    # core model
    core_model_file = cfg["frozen_core"]["core_model_file"]
    with open(os.path.join(cfg_dir, core_model_file)) as f:
        core_model = json.load(f)

    # Radial BIN EDGES for the nearest-atom distance histogram of |psi_nk|^2.
    # F_c decays as 1/r^4 so a few Bohr around each atom suffices; bin finely near
    # the core where F_c is largest.
    rbins = np.concatenate([
        np.linspace(0.0, 2.0, 80),
        np.linspace(2.0, 6.0, 40)[1:],
    ])
    r_centers = 0.5 * (rbins[:-1] + rbins[1:])

    # Per-channel F_c(r_centers) and DeltaE_c
    channels = []
    for ch in core_model["core_s_channels"]:
        dat = np.loadtxt(os.path.join(cfg_dir, ch["atomic_core_file"]),
                         delimiter=",", skiprows=1)
        r_c, u_c, vh_c = dat[:, 0], dat[:, 1], dat[:, 2]
        Fc = compute_Fc_on_grid(r_c, u_c, vh_c, r_centers)
        channels.append({
            "channel": ch["channel"],
            "DeltaE_c": float(ch["DeltaE_c_Ha"]),
            "Fc": Fc,
            "Nc": 2.0,  # closed s-shell: 2 electrons
        })

    # number of bands to compute: occupied + a few buffer
    n_occ = z_val * n_atoms
    n_bands = n_occ + 6

    # Build spec for Julia and run ONE SCF + bands.
    spec = {
        "element": element,
        "lattice_vectors_bohr": struct["lattice_vectors_bohr"],
        "atom_positions_frac": struct["atom_positions_frac"],
        "ecut_Ha": float(dft["ecut_Ha"]),
        "kgrid": [int(x) for x in dft["kgrid"]],
        "smearing_Ha": float(dft["smearing_Ha"]),
        "pseudopotential_family": dft["pseudopotential_family"],
        "endpoint_frac": endpoint,
        "ts": ts.tolist(),
        "n_bands": int(n_bands),
        "rbins": rbins.tolist(),
    }

    tmpdir = tempfile.mkdtemp(prefix="qp_")
    spec_file = os.path.join(tmpdir, "spec.json")
    out_json = os.path.join(tmpdir, "bands.json")
    spec["out_json"] = out_json
    with open(spec_file, "w") as f:
        json.dump(spec, f)

    env = dict(os.environ)
    env.setdefault("OPENBLAS_NUM_THREADS", "2")
    env.setdefault("OMP_NUM_THREADS", "2")
    env.setdefault("JULIA_NUM_THREADS", "1")
    env.setdefault("JULIA_PROJECT", os.path.join(HERE, "jenv"))

    jl = os.path.join(HERE, "compute_bands.jl")
    log("Running DFTK SCF + bands ...")
    subprocess.run(["julia", jl, spec_file], check=True, env=env)

    with open(out_json) as f:
        bands = json.load(f)

    eF = float(bands["eF_Ha"])
    Omega = float(bands["Omega"])           # cell volume (Bohr^3)
    Omega_atom = Omega / n_atoms            # volume per atom (per-site density norm)
    inv_pi = 1.0 / np.pi                    # closure phase-space coefficient

    rows = []
    for kp in bands["kpoints"]:
        t = float(kp["t"])
        eigs = np.array(kp["eigenvalues_Ha"])
        prof = np.array(kp["prof"])  # (nb, nbin): prob. mass of |psi|^2 per shell bin
        # which bands to emit: every occupied (eps<eF) + the single lowest unocc
        order = np.argsort(eigs)
        eigs_s = eigs[order]
        prof_s = prof[order]
        occ_mask = eigs_s < eF
        n_occ_here = int(np.sum(occ_mask))
        emit_idx = list(range(n_occ_here + 1))  # occupied + first unoccupied
        emit_idx = [i for i in emit_idx if i < len(eigs_s)]

        # match the point_id by t
        pid = point_ids[np.argmin(np.abs(ts - t))]

        for bi in emit_idx:
            eps = eigs_s[bi]
            p = prof_s[bi]
            # <psi|F_c|psi>_nk = sum_bins P(r_bin) * F_c(r_bin)  (P = int_{bin}|psi|^2 dV,
            # summed over ALL atomic sites in the cell).  Per-atom on-site
            # fluctuation W_c = Omega_atom * (per-site overlap) = Omega_atom * <F>/n_atom.
            # lambda = (1/pi) sum_c N_c W_c / DeltaE_c^2.
            lam = 0.0
            for ch in channels:
                F_cell = float(np.sum(p * ch["Fc"]))       # summed over sites
                W_c = Omega_atom * (F_cell / n_atoms)      # per-atom on-site strength
                lam += inv_pi * ch["Nc"] * W_c / ch["DeltaE_c"] ** 2
            Z = 1.0 / (1.0 + lam)
            E_qp = eF + Z * (eps - eF)
            E_pred_eV = (E_qp - eF) * HARTREE_TO_EV
            rows.append((element, int(pid), t, E_pred_eV))

    with open(args.out, "w") as f:
        f.write("element,point_id,t,E_pred_eV\n")
        for el, pid, t, e in rows:
            f.write(f"{el},{pid},{t:.6f},{e:.6f}\n")
    log(f"Wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
