#!/usr/bin/env python3
"""All-electron radial atomic LDA (Dirac-exchange) solver.

Produces the all-electron core orbitals u_c(r)=r R_c(r) needed for the
frozen-core form factor f_c(K) (paper Eq. 7). Finite-difference radial
eigensolver on a uniform grid; self-consistent Hartree + Dirac (Slater)
exchange, optional PW92 correlation.

Stage 1 (this file): radial eigensolver, validated on hydrogen.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import eigh_tridiagonal


def radial_eigenstates(r, h, Veff, l, nstates):
    """Lowest `nstates` eigenpairs of -1/2 u'' + [V + l(l+1)/2r^2] u = e u
    on uniform grid r (u(0)=u(rmax)=0 Dirichlet). Returns (energies, U[:,k]).
    U columns are normalized so that sum(u^2) h = 1 (i.e. ∫|R|^2 r^2 dr = 1)."""
    N = len(r)
    centrifugal = l * (l + 1) / (2.0 * r**2)
    diag = 1.0 / h**2 + Veff + centrifugal
    off = -0.5 / h**2 * np.ones(N - 1)
    evals, evecs = eigh_tridiagonal(diag, off, select="i",
                                    select_range=(0, nstates - 1))
    # normalize: ∫ u^2 dr = 1  (trapz on uniform grid ≈ sum*h)
    for k in range(evecs.shape[1]):
        nrm = np.sqrt(np.sum(evecs[:, k] ** 2) * h)
        evecs[:, k] /= nrm
        if evecs[1, k] < 0:  # fix sign: u>0 near origin
            evecs[:, k] *= -1
    return evals, evecs


# neutral-atom electron configurations: list of (n, l, occupation)
CONFIGS = {
    "Li": [(1, 0, 2), (2, 0, 1)],
    "Na": [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 1)],
    "Mg": [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2)],
    "Al": [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 1)],
    "Si": [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 2)],
    "K":  [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (4, 0, 1)],
    "Ca": [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (4, 0, 2)],
}
Z_OF = {"Li": 3, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "K": 19, "Ca": 20}


def hartree_potential(rho_u, r, h):
    """V_H(r) from rho_u(r) = 4π n r^2 = sum_i occ_i u_i^2.
    V_H = (1/r)∫_0^r rho_u dr' + ∫_r^∞ rho_u/r' dr'."""
    inner = np.cumsum(rho_u) * h            # ∫_0^r rho_u
    tail = np.cumsum((rho_u / r)[::-1])[::-1] * h  # ∫_r^∞ rho_u/r'
    return inner / r + tail


def dirac_exchange_potential(rho_u, r):
    """Slater/Dirac exchange V_x = -(3/π n)^{1/3}, n = rho_u/(4π r^2)."""
    n = np.maximum(rho_u / (4.0 * np.pi * r**2), 1e-30)
    return -(3.0 / np.pi * n) ** (1.0 / 3.0)


def solve_atom(symbol, h=0.002, rmax=30.0, alpha=0.3, tol=1e-6, maxit=300,
               verbose=False):
    """Self-consistent radial LDA (exchange-only / Dirac) atom.
    Returns dict: r, h, orbitals {(n,l): (energy, u_array)}, V_total, V_H."""
    Z = Z_OF[symbol]
    cfg = CONFIGS[symbol]
    r = np.arange(h, rmax + h, h)
    lmax = max(l for _, l, _ in cfg)

    V = -Z / r  # initial: bare nucleus
    orbitals = {}
    eprev = None
    for it in range(maxit):
        rho_u = np.zeros_like(r)
        orbitals = {}
        for l in range(lmax + 1):
            ns_l = [n for (n, ll, _) in cfg if ll == l]
            if not ns_l:
                continue
            kmax = max(n - l - 1 for n in ns_l)
            ev, U = radial_eigenstates(r, h, V, l, kmax + 1)
            for (n, ll, occ) in cfg:
                if ll != l:
                    continue
                k = n - l - 1
                u = U[:, k]
                orbitals[(n, l)] = (ev[k], u)
                rho_u += occ * u**2
        VH = hartree_potential(rho_u, r, h)
        Vx = dirac_exchange_potential(rho_u, r)
        Vnew = -Z / r + VH + Vx
        V = (1 - alpha) * V + alpha * Vnew
        es = np.array([orbitals[(n, l)][0] for (n, l, _) in cfg])
        if eprev is not None and np.max(np.abs(es - eprev)) < tol:
            if verbose:
                print(f"  {symbol}: converged in {it+1} iters")
            break
        eprev = es
    return {"r": r, "h": h, "Z": Z, "orbitals": orbitals,
            "V_total": V, "V_H": VH, "rho_u": rho_u, "config": cfg}


def core_form_factor(symbol, n, l, Kgrid, h=0.002, rmax=30.0):
    """f_c(K) for core orbital (n,l), paper Eq. 7 (s-wave core):
        f_c^K = (sqrt(4π)/K) ∫_0^∞ u_c(r) [V_{H,c}(r) - J_c] sin(Kr) dr
    where V_{H,c} is the single-orbital Hartree potential of u_c (∫u_c^2 dr=1)
    and J_c = ∫ u_c^2 V_{H,c} dr the self-Coulomb integral. Returns f_c on Kgrid."""
    res = solve_atom(symbol, h=h, rmax=rmax)
    r = res["r"]
    _, u = res["orbitals"][(n, l)]
    u2 = u**2
    VHc = hartree_potential(u2, r, h)          # single-orbital Hartree
    Jc = np.sum(u2 * VHc) * h                   # self-Coulomb average
    w = u * (VHc - Jc)                          # radial weight
    f = np.empty_like(Kgrid, dtype=float)
    pref = np.sqrt(4.0 * np.pi)
    for i, K in enumerate(Kgrid):
        if K < 1e-8:
            f[i] = pref * np.sum(w * r) * h     # sin(Kr)/K -> r as K->0
        else:
            f[i] = pref / K * np.sum(w * np.sin(K * r)) * h
    return f, Jc


def _test_form_factor():
    print("\nCore form factor f_c(K) (Na 2s), paper Fig.2 plots f_c/DeltaE_c:")
    dE = 2.70  # Na 2s, from element config
    K = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0])
    f, Jc = core_form_factor("Na", 2, 0, K)
    print(f"  Na 2s: J_c={Jc:.4f} Ha")
    print("   K(1/Bohr) :", "  ".join(f"{k:5.1f}" for k in K))
    print("   f_c (Ha)  :", "  ".join(f"{v:5.2f}" for v in f))
    print("   f_c/DeltaE:", "  ".join(f"{v/dE:5.2f}" for v in f))
    print(f"  f_c(0)={f[0]:.3f} Ha;  (f_c(0)/DeltaE)^2 = {(f[0]/dE)**2:.3f}"
          f"  -> z_core(0-th approx, G=0 only) = {1/(1+(f[0]/dE)**2):.3f}  [target 0.80]")


def _test_atoms():
    print("\nSelf-consistent atomic LDA (Dirac exchange) orbital energies (Ha):")
    for sym in ["Li", "Na", "Mg", "Al", "K"]:
        res = solve_atom(sym, verbose=True)
        line = "  ".join(f"{n}{'spdf'[l]}={res['orbitals'][(n,l)][0]:.3f}"
                         for (n, l, _) in res["config"])
        print(f"  {sym:3s}: {line}")


def _test_hydrogen():
    h = 0.002
    rmax = 60.0
    r = np.arange(h, rmax + h, h)
    V = -1.0 / r  # bare Coulomb, Z=1
    print("Hydrogen radial eigensolver test (Ha):")
    print(f"  grid: h={h}, rmax={rmax}, N={len(r)}")
    for l, label, exact in [(0, "1s/2s/3s", [-0.5, -0.125, -0.05556]),
                            (1, "2p/3p/4p", [-0.125, -0.05556, -0.03125])]:
        ev, _ = radial_eigenstates(r, h, V, l, 3)
        print(f"  l={l} {label}: computed {np.round(ev, 5)}  exact {exact}")


if __name__ == "__main__":
    _test_hydrogen()
    _test_atoms()
