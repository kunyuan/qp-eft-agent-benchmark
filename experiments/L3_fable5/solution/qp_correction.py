"""
qp_correction.py — frozen-core dynamical quasiparticle correction.

Derived correction (see method.md for the step-by-step derivation):

Integrating the core electrons out of the all-electron action leaves, besides
the static pseudopotential V_PSP already inside H_KS, a dynamical
(frequency-dependent) induced potential delta_V_pp(omega).  Its physical
content: the core's screening of the bare nuclear attraction is retarded --
the core electrons in shell c can only follow the valence motion on the time
scale 1/DeltaE_c.  At frequencies far below DeltaE_c the valence sees the
fully screened (pseudized) ion; the leading omitted dynamics is carried by the
core s-shell (monopole) channels c, each contributing a causal one-pole term

    Sigma_c(i omega) = |Lambda_c(nu,k)|^2 / (i omega - eps_c),   eps_c = -DeltaE_c,

with the channel coupling (bare-nucleus vertex, normalized per valence
electron fed into the channel)

    |Lambda_c(nu,k)|^2 = (Z^2 / z_val) * sum_atoms |<phi_c | 1/r | psi~_{nu,k}>|^2 .

The static part of Sigma_c (level repulsion |Lambda|^2/DeltaE_c) is part of
V_PSP and must not be re-added; the leading dynamical remainder is the
linear-in-omega term, which renormalizes the i*omega coefficient:

    z_core(nu,k) = 1 / (1 + w_{nu,k}),
    w_{nu,k} = sum_c |Lambda_c(nu,k)|^2 / DeltaE_c^2
             = (Z^2/z_val) sum_{c,atoms} |<phi_c|1/r|psi~>|^2 / DeltaE_c^2 ,

    eps_QP - eps_F = z_core * (eps_KS - eps_F).

All ingredients come from the packet data (u_c(r), DeltaE_c, Z_nuclear,
z_valence) and the DFTK Bloch coefficients.  No fitted parameters.
"""

import json
import os
import numpy as np

HARTREE_EV = 27.211386245988


def load_core_channels(core_model_path):
    """Return list of dicts: name, r (Bohr), u (radial u_c=r*R_c), dE (Ha)."""
    base = os.path.dirname(os.path.abspath(core_model_path))
    with open(core_model_path) as f:
        cm = json.load(f)
    channels = []
    for ch in cm["core_s_channels"]:
        dat = np.loadtxt(os.path.join(base, ch["atomic_core_file"]),
                         delimiter=",", skiprows=1)
        r, u = dat[:, 0], dat[:, 1]
        # trapezoid weights on the native grid
        wts = np.zeros_like(r)
        wts[1:-1] = 0.5 * (r[2:] - r[:-2])
        wts[0] = 0.5 * (r[1] - r[0]) + r[0]  # include [0, r0] sliver
        wts[-1] = 0.5 * (r[-1] - r[-2])
        channels.append({"name": ch["channel"], "r": r, "u": u,
                         "w": wts, "dE": float(ch["DeltaE_c_Ha"])})
    return channels


def coupling_form_factor(channel, qnorms):
    """I_c(q) = int dr u_c(r) j0(q r)  for an array of |k+G| values.

    This equals <phi_c| 1/r |e^{iqr}> / sqrt(4 pi) with phi_c = (u_c/r) Y00:
        <phi_c|1/r|PW> = sqrt(4 pi) * int dr u_c(r) j0(qr).
    """
    r, u, w = channel["r"], channel["u"], channel["w"]
    out = np.empty(len(qnorms))
    # chunk to bound memory: (nq x nr)
    chunk = max(1, int(4e7 / len(r)))
    uw = u * w
    for i in range(0, len(qnorms), chunk):
        qr = np.outer(qnorms[i:i + chunk], r)
        out[i:i + chunk] = np.sinc(qr / np.pi) @ uw
    return out


def w_dynamical(channels, Z, z_val, qvecs_cart, Gint, coeffs, tau_frac, volume):
    """w_{nu,k} for every band at one k-point.

    qvecs_cart : (nG,3) cartesian k+G in Bohr^-1
    Gint       : (nG,3) integer G vectors
    coeffs     : (nG, nbands) complex plane-wave coefficients (DFTK normalized)
    tau_frac   : (natoms,3) fractional atomic positions
    volume     : unit-cell volume in Bohr^3
    Returns (nbands,) array of w.
    """
    qn = np.linalg.norm(qvecs_cart, axis=1)
    nb = coeffs.shape[1]
    w = np.zeros(nb)
    sq4pi = np.sqrt(4.0 * np.pi)
    for ch in channels:
        Ic = coupling_form_factor(ch, qn)
        for ta in tau_frac:
            phase = np.exp(2j * np.pi * (Gint @ np.asarray(ta)))
            # <phi_c|1/r|psi~_b> = sq4pi/sqrt(Om) * sum_G c_Gb phase_G I_c(q_G)
            amp = (coeffs * (phase * Ic)[:, None]).sum(axis=0)
            M2 = (sq4pi ** 2 / volume) * np.abs(amp) ** 2
            w += (Z ** 2 / z_val) * M2 / ch["dE"] ** 2
    return w
