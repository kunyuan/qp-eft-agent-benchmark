"""
qp_physics.py -- frozen-core dynamical quasiparticle renormalization z_core.

Closed form derived in method.md (parameter-free, prefactor 1):

    z_core_{nu,k} = 1 / ( 1 + sum_c  g_{c,nk}^2 / DeltaE_c^2 )

with the valence<->core dynamical coupling (monopole / s-shell channels)

    g_{c,nk} = 2 * < phi_c | V_H^core | psi_{nk} >                       (Eq. A)

  V_H^core(r) : the radial monopole Hartree potential built from the FULL core
                density n_core(r) = sum_{c'} occ_{c'} |phi_{c'}(r)|^2 (occ=2 per
                closed s-shell) -- i.e. the Coulomb field the core leaves behind.
  phi_c(r)    = u_c(r) / (r sqrt(4pi))      (atomic core s-orbital, monopole)
  psi_{nk}    : Bloch state; only its s-wave (monopole) component at the atom
                couples to the s-core, psi0(r) = (1/sqrt(Omega_at)) sum_G c_G j0(q_G r),
                q_G = |k+G| (Bohr^-1), Omega_at = cell volume / n_atom.
  factor 2    : the closed s-shell carries two electrons (the antisymmetrized
                direct+exchange monopole fluctuation; see method.md).

In radial form, with phi_c r^2 4pi -> u_c r sqrt(4pi):

    g_{c,nk} = 2 * sqrt(4pi) * int_0^inf u_c(r) V_H^core(r) psi0(r) r dr.

z_core is STATE-DEPENDENT (it varies with n,k through psi0, i.e. through the
valence-core overlap), < 1, and reduces the occupied bandwidth.  The QP energy is

    eps_QP = eps_F + z_core (eps_KS - eps_F).
"""

import numpy as np
from scipy.special import spherical_jn

HA_EV = 27.211386245988


def radial_hartree(r, dens):
    """Monopole Hartree potential V(r) = int dens(r')/max(r,r') dr', where
    dens(r') is a radial CHARGE density already carrying the 4 pi r'^2 factor
    (so that total charge = int dens dr').  V(r) = (1/r) int_0^r dens dr'
    + int_r^inf dens/r' dr'.  Uniform or non-uniform grid (trapezoid)."""
    seg = 0.5 * (dens[1:] + dens[:-1]) * np.diff(r)
    cum_in = np.concatenate([[0.0], np.cumsum(seg)])           # int_0^r dens
    dover = dens / r
    seg2 = 0.5 * (dover[1:] + dover[:-1]) * np.diff(r)
    cum_out = np.concatenate([[0.0], np.cumsum(seg2)])         # int_0^r dens/r'
    return cum_in / r + (cum_out[-1] - cum_out)


def build_core(channels):
    """channels: list of dicts {name, r, u, dE, occ}.  Returns the per-channel
    data plus a callable V_H^core interpolated onto each channel's r grid."""
    # full core density on the longest grid
    rfull = max((c["r"] for c in channels), key=len)
    ncore = np.zeros_like(rfull)
    for c in channels:
        ui = np.interp(rfull, c["r"], c["u"], right=0.0)
        ncore += c["occ"] * ui ** 2                            # 4pi r^2 n = sum occ u^2
    VHfull = radial_hartree(rfull, ncore)
    for c in channels:
        c["VH"] = np.interp(c["r"], rfull, VHfull)
    return channels


def coupling_g(channel, c_re, c_im, kpg, Omega_at):
    """g_{c,nk} = 2 sqrt(4pi) int u_c V_H^core psi0 r dr  (Eq. A, radial)."""
    r, u, VH = channel["r"], channel["u"], channel["VH"]
    c = c_re + 1j * c_im
    j0 = spherical_jn(0, np.outer(kpg, r))                     # (nG, nr)
    psi0 = (c @ j0) / np.sqrt(Omega_at)                        # monopole Bloch amp
    g = 2.0 * np.sqrt(4.0 * np.pi) * np.trapz(u * VH * psi0 * r, r)
    return g                                                   # complex; |g| used


def z_core_state(channels, c_re, c_im, kpg, Omega_at):
    """z_core for one Bloch state given all core s-channels."""
    s = 0.0
    for ch in channels:
        g = coupling_g(ch, c_re, c_im, kpg, Omega_at)
        s += abs(g) ** 2 / ch["dE"] ** 2
    return 1.0 / (1.0 + s)
