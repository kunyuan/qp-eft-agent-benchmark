#!/usr/bin/env python3
"""Self-contained spherical all-electron LDA (PZ81) radial atomic solver.

Element-blind: input Z (nuclear charge) and z_valence; core = aufbau filling of
Z electrons with the z_valence outermost removed. Outputs core orbitals u_nl(r),
LDA eigenvalues eps_nl, and spherical-Bessel form factors chi~_nl(kappa).

Log-radial grid; generalized tridiagonal eigenproblem per l channel:
  x = ln r, u = r R, v = u / sqrt(r):
  [-d2/dx2 + 2 r^2 V(r) + (l+1/2)^2] v = 2 r^2 E v
Solved with scipy sparse shift-invert (fallback: dense eigh).
Atomic units throughout.
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh

AUFBAU = [(1,0),(2,0),(2,1),(3,0),(3,1),(4,0),(3,2),(4,1),(5,0),(4,2),(5,1),(6,0),(4,3),(5,2),(6,1),(7,0)]

def occupations(Z):
    """Aufbau occupation list [(n,l,occ)] for neutral atom with Z electrons."""
    occs, nel = [], 0
    for (n, l) in AUFBAU:
        if nel >= Z:
            break
        cap = 2*(2*l+1)
        o = min(cap, Z - nel)
        occs.append([n, l, float(o)])
        nel += o
    return occs

def split_core_valence(occs, z_valence):
    """Remove z_valence electrons from the latest-filled shells; rest = core."""
    occs = [list(o) for o in occs]
    rem = z_valence
    for o in reversed(occs):
        take = min(o[2], rem)
        o[2] -= take
        rem -= take
        if rem <= 0:
            break
    core = [o for o in occs if o[2] > 0]
    return core

def pz81_vxc(n):
    """PZ81 LDA xc potential and energy density (unpolarized). n: density array."""
    n = np.maximum(n, 1e-30)
    rs = (3.0/(4.0*np.pi*n))**(1.0/3.0)
    # exchange
    ex = -0.75*(3.0/np.pi)**(1.0/3.0)*n**(1.0/3.0)
    vx = -(3.0/np.pi)**(1.0/3.0)*n**(1.0/3.0)
    # correlation PZ81
    A, B, C, D = 0.0311, -0.048, 0.0020, -0.0116
    g, b1, b2 = -0.1423, 1.0529, 0.3334
    ec = np.empty_like(rs); vc = np.empty_like(rs)
    lo = rs < 1.0
    lrs = np.log(rs[lo])
    ec[lo] = A*lrs + B + C*rs[lo]*lrs + D*rs[lo]
    vc[lo] = A*lrs + (B - A/3.0) + (2.0/3.0)*C*rs[lo]*lrs + (2.0*D - C)/3.0*rs[lo]
    hi = ~lo
    sq = np.sqrt(rs[hi])
    den = 1.0 + b1*sq + b2*rs[hi]
    ec[hi] = g/den
    vc[hi] = ec[hi]*(1.0 + 7.0/6.0*b1*sq + 4.0/3.0*b2*rs[hi])/den
    return vx + vc, ex + ec

class RadialAtom:
    def __init__(self, Z, N=3200, rmin=None, rmax=60.0):
        self.Z = float(Z)
        rmin = rmin or 1e-7/Z
        self.x = np.linspace(np.log(rmin), np.log(rmax), N)
        self.h = self.x[1] - self.x[0]
        self.r = np.exp(self.x)
        self.N = N

    def solve_channel(self, V, l, nstates):
        """Lowest nstates eigenpairs of the radial problem in channel l.
        Returns (eps[], u[][]) with u = r*R normalized: int u^2 dr = 1."""
        r, h, N = self.r, self.h, self.N
        w = 2.0*r*r*V + (l + 0.5)**2
        main = 2.0/h**2 + w
        off = -1.0/h**2*np.ones(N-1)
        A = diags([off, main, off], [-1, 0, 1], format='csc')
        Bd = 2.0*r*r
        B = diags([Bd], [0], format='csc')
        try:
            sigma = -0.6*self.Z**2 - 1.0
            vals, vecs = eigsh(A, k=nstates, M=B, sigma=sigma, which='LM')
            order = np.argsort(vals)
            vals, vecs = vals[order], vecs[:, order]
        except Exception:
            Ad = A.toarray(); Bdense = np.diag(Bd)
            vals, vecs = eigh(Ad, Bdense)
            vals, vecs = vals[:nstates], vecs[:, :nstates]
        us = []
        for j in range(vecs.shape[1]):
            u = vecs[:, j]*np.sqrt(r)          # u = v*sqrt(r)
            norm = np.sqrt(np.trapezoid(u*u, r))
            u = u/norm
            if u[np.argmax(np.abs(u))] < 0:
                u = -u
            us.append(u)
        return vals, us

    def scf(self, occs, mix=0.35, tol=1e-8, maxiter=200, verbose=False):
        """occs: [[n,l,occ]]; returns dict with eps, u per shell and potentials."""
        r = self.r
        nelec = sum(o[2] for o in occs)
        # initial guess: Thomas-Fermi-ish screened density via hydrogenic 1s
        zeta0 = max(self.Z - 0.3, 1.0)
        n_e = nelec*(zeta0**3/np.pi)*np.exp(-2*zeta0*r)
        eps_prev = None
        for it in range(maxiter):
            # Hartree potential of density n_e (electrons, positive density)
            rho_r2 = 4*np.pi*n_e*r*r
            q_in = np.concatenate([[0.0], np.cumsum(0.5*(rho_r2[1:]+rho_r2[:-1])*np.diff(r))])
            # outer integral int_r^inf rho*r dr
            rho_r = 4*np.pi*n_e*r
            o_int = np.concatenate([[0.0], np.cumsum(0.5*(rho_r[1:]+rho_r[:-1])*np.diff(r))])
            o_out = o_int[-1] - o_int
            VH = q_in/r + o_out
            vxc, _ = pz81_vxc(n_e)
            V = -self.Z/r + VH + vxc
            # solve channels
            lmax = max(o[1] for o in occs)
            new_n = np.zeros_like(r)
            eps_all, u_all = {}, {}
            for l in range(lmax+1):
                # number of states needed in this channel
                ns = [o for o in occs if o[1] == l]
                if not ns:
                    continue
                nstates = max(o[0]-l for o in ns)
                vals, us = self.solve_channel(V, l, nstates)
                for o in ns:
                    idx = o[0]-l-1                     # node count = n-l-1
                    eps_all[(o[0], l)] = vals[idx]
                    u_all[(o[0], l)] = us[idx]
                    new_n += o[2]*us[idx]**2/(4*np.pi*r*r)
            # convergence on eigenvalues
            if eps_prev is not None:
                d = max(abs(eps_all[k]-eps_prev[k]) for k in eps_all)
                if verbose:
                    print(f"  it={it} d={d:.2e}")
                if d < tol:
                    n_e = (1-mix)*n_e + mix*new_n
                    break
            eps_prev = eps_all
            n_e = (1-mix)*n_e + mix*new_n
        return {"eps": eps_all, "u": u_all, "V": V, "n": n_e, "r": r, "occs": occs}

def total_energy(atom, res):
    """LDA total energy from converged result (for Delta-SCF uses)."""
    r = res["r"]; n_e = res["n"]
    rho_r2 = 4*np.pi*n_e*r*r
    q_in = np.concatenate([[0.0], np.cumsum(0.5*(rho_r2[1:]+rho_r2[:-1])*np.diff(r))])
    rho_r = 4*np.pi*n_e*r
    o_int = np.concatenate([[0.0], np.cumsum(0.5*(rho_r[1:]+rho_r[:-1])*np.diff(r))])
    o_out = o_int[-1] - o_int
    VH = q_in/r + o_out
    vxc, exc = pz81_vxc(n_e)
    Eband = sum(o[2]*res["eps"][(o[0], o[1])] for o in res["occs"])
    EH = 0.5*np.trapezoid(4*np.pi*r*r*n_e*VH, r)
    Exc = np.trapezoid(4*np.pi*r*r*n_e*exc, r)
    Evxc = np.trapezoid(4*np.pi*r*r*n_e*vxc, r)
    return Eband - EH + Exc - Evxc

def form_factor(r, u, l, kappas):
    """chi~_nl(kappa) = 4pi int j_l(kr) R_nl(r) r^2 dr,  R = u/r."""
    from scipy.special import spherical_jn
    out = np.empty(len(kappas))
    for i, k in enumerate(kappas):
        jl = spherical_jn(l, k*r)
        out[i] = 4*np.pi*np.trapezoid(jl*u*r, r)
    return out

if __name__ == "__main__":
    import sys
    Z = float(sys.argv[1]) if len(sys.argv) > 1 else 3
    zv = float(sys.argv[2]) if len(sys.argv) > 2 else 1
    at = RadialAtom(Z)
    occs = occupations(int(Z))
    res = at.scf(occs, verbose=False)
    print(f"Z={Z} occs={occs}")
    for (n, l), e in sorted(res["eps"].items()):
        print(f"  n={n} l={l} eps={e:.5f} Ha")
    core = split_core_valence(occs, zv)
    print("core:", core)
    print("Etot:", total_energy(at, res))
