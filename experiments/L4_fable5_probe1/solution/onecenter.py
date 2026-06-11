"""One-center (all-electron) evaluation of the dynamical core-polarization
self-energy.

Same derived diagram as corepol.py (Fan-Migdal self-energy from the core
fluctuation propagator, static closure limit subtracted), but the Bloch-state
contraction is done in the all-electron local representation:

  dE_nk = sum_sites sum_l W_l^{site}(nk) * delta_eps_l(eps_nk)

  W_l(nk): weight of the local angular-momentum-l component of the (pseudo)
  Bloch state inside the site's atomic (WS) sphere.  By norm conservation the
  pseudo state carries the correct *norm* per channel; its core-region SHAPE
  is wrong, so the radial wave used in the matrix elements is the regular
  all-electron solution u_l(r; eps) of the atomic KS potential at the band
  energy, unit-normalized in the sphere (this restores the core-orthogonality
  lobes that dominate the coupling integrals; see method.md audit).

  delta_eps_l(eps) = sum_b sum_{l',j'} P(b,l',l) Rad^2 * K(eps, eps'_j', D_b)

  Rad = int_0^... u_{l'j'}(r) w_L^{(b)}(r) u_l(r; eps) dr
  w_L^{(b)}(r) = int dr' u_x(r') u_c(r') r_<^L / r_>^{L+1}   (radial Coulomb
                 potential of the channel transition density)
  P(b,l',l) = W_bL * (2l'+1) threej(l,L,l';000)^2 / ((2L+1)^2 (2l+1))
  W_bL = 2 (2l_c+1)(2l_x+1)(2L+1) threej(l_c,L,l_x;000)^2

  K(e,e',D) = (1-f')/(e-e'-D) + f'/(e-e'+D) + 1/D   (retardation-only kernel)

  Intermediates |l'j'> are the eigenstates of the same atomic KS potential
  (bound virtuals + box-discretized continuum), aligned to the crystal via
  eps' = eps_{l'j'} + U0,  U0 = EF - kF^2/2  (free-dispersion band offset),
  occupied core states excluded, Fermi sea handled by f'(eps').

Element blind: inputs are Z_nuclear, z_valence, and crystal data only.
"""
import numpy as np
from scipy.special import spherical_jn

from atom import Atom
from corepol import threej000_sq, HARTREE_EV


# ---------------------------------------------------------------- channels
class AtomChannels:
    """Atomic KS reference: core shells, excited spectra, channel groups."""

    def __init__(self, Z, z_valence, emax=60.0, Lmax=2, rmax=50.0):
        self.Z, self.zval, self.Lmax = Z, z_valence, Lmax
        self.at = Atom(Z, rmax=rmax).scf()
        self.core, self.val = self.at.split_core_valence(z_valence)
        self.occupied_nl = {c["nl"] for c in self.core}
        self.mesh = self.at.mesh
        self.emax = emax
        lmax_x = max(c["nl"][1] for c in self.core) + Lmax
        lmax_inter = 4 + Lmax  # local waves up to l=4 scatter to l' <= l+L
        self.spec = {l: self.at.spectrum(l, emax)
                     for l in range(0, max(lmax_x, lmax_inter) + 1)}

        # channel groups: (core shell c, lx, L) with all x vectorized
        r, dr = self.mesh.r, self.mesh.dr
        self.groups = []
        for c in self.core:
            n_c, l_c = c["nl"]
            for lx in range(0, l_c + Lmax + 1):
                eps_x, U_x = self.spec[lx]
                for L in range(Lmax + 1):
                    if (l_c + lx + L) % 2 or not (abs(l_c - lx) <= L <= l_c + lx):
                        continue
                    tj2 = threej000_sq(l_c, L, lx)
                    if tj2 == 0.0:
                        continue
                    keep = [j for j in range(U_x.shape[1])
                            if (j + lx + 1, lx) not in self.occupied_nl
                            and eps_x[j] - c["eps"] > 0]
                    if not keep:
                        continue
                    W = 2.0 * (2 * l_c + 1) * (2 * lx + 1) * (2 * L + 1) * tj2
                    Delta = eps_x[keep] - c["eps"]
                    # radial Coulomb potentials w_L(r) for all x in the group
                    f = U_x[:, keep] * c["u"][:, None]          # (nr, nx)
                    rl = r ** L
                    rml = r ** (-(L + 1.0))
                    a = f * rl[:, None]
                    g1 = np.zeros_like(f)
                    g1[1:] = np.cumsum(0.5 * (a[1:] + a[:-1]) * dr, axis=0)
                    b = f * rml[:, None]
                    g2tot = np.trapezoid(b, dx=dr, axis=0)
                    g2 = g2tot[None, :] - np.concatenate(
                        [np.zeros((1, f.shape[1])),
                         np.cumsum(0.5 * (b[1:] + b[:-1]) * dr, axis=0)])
                    wL = g1 * rml[:, None] + g2 * rl[:, None]    # (nr, nx)
                    self.groups.append({
                        "shell": (n_c, l_c), "lx": lx, "L": L, "W": W,
                        "Delta": Delta, "wL": wL})

    def alpha0_dipole(self):
        """Static dipole polarizability 2 sum |<z>|^2/Delta (m,spin summed)."""
        r, dr = self.mesh.r, self.mesh.dr
        alpha = 0.0
        core = {c["nl"]: c for c in self.core}
        for g in self.groups:
            if g["L"] != 1:
                continue
            c = core[g["shell"]]
            eps_x, U_x = self.spec[g["lx"]]
            keep = [j for j in range(U_x.shape[1])
                    if (j + g["lx"] + 1, g["lx"]) not in self.occupied_nl
                    and eps_x[j] - c["eps"] > 0]
            d = (U_x[:, keep] * (r * c["u"])[:, None]).sum(axis=0) * dr
            alpha += np.sum(2.0 * g["W"] * d ** 2 / (9.0 * g["Delta"]))
        return float(alpha)


# ---------------------------------------------------------------- kernel
def kernel(eps, epsp, Delta, EF, T):
    x = np.clip((epsp - EF) / max(T, 1e-12), -40.0, 40.0)
    f = 1.0 / (1.0 + np.exp(x))
    de = eps - epsp
    return (1.0 - f) / (de - Delta) + f / (de + Delta) + 1.0 / Delta


# ---------------------------------------------------------------- evaluator
class OneCenterSigma:
    def __init__(self, channels, EF, volume, n_electrons, n_atoms,
                 positions_frac, recip_cols, T=0.001, lmax_loc=4):
        self.chn = channels
        self.EF, self.Omega, self.T = EF, volume, T
        self.n_atoms = n_atoms
        self.tau_frac = np.array(positions_frac)
        self.B = np.array(recip_cols).T   # columns b1,b2,b3
        self.lmax_loc = lmax_loc
        n_dens = n_electrons / volume
        self.kF = (3.0 * np.pi ** 2 * n_dens) ** (1.0 / 3.0)
        self.U0 = EF - 0.5 * self.kF ** 2
        self.R_at = (3.0 * volume / (4.0 * np.pi * n_atoms)) ** (1.0 / 3.0)
        mesh = channels.mesh
        self.iR = min(int(self.R_at / mesh.dr), mesh.N - 1)
        # exclude occupied (core) intermediates per l'
        self.inter = {}
        for l, (eps_p, U_p) in channels.spec.items():
            mask = np.ones(len(eps_p), bool)
            for (nn, ll) in channels.occupied_nl:
                if ll == l and nn - ll - 1 < len(eps_p):
                    mask[nn - ll - 1] = False
            self.inter[l] = (eps_p[mask], U_p[:, mask])
        self._dnodes = None

    # ---- local radial wave at energy eps (unit norm in sphere) ----
    def local_wave(self, l, eps, mode="ae"):
        """mode='ae': regular all-electron solution of the atomic KS
        potential at energy eps (Numerov).  mode='ps': smooth pseudo-shape
        reference, the free wave r*j_l(kappa r) at the pseudo crystal
        momentum kappa(eps) = sqrt(2 max(eps-U0, ~0)) -- the local shape the
        plane-wave evaluation effectively uses."""
        mesh = self.chn.mesh
        r, dr = mesh.r, mesh.dr
        if mode == "ps":
            kappa = np.sqrt(2.0 * max(eps - self.U0, 1e-8))
            u = r * spherical_jn(l, kappa * r)
            u = u.copy()
            u[self.iR + 10:] = 0.0
        else:
            V = self.chn.at.V
            n = self.iR + 10
            u = np.zeros(mesh.N)
            u[0] = r[0] ** (l + 1)
            u[1] = r[1] ** (l + 1)
            f = 2.0 * (V[:n + 1] + l * (l + 1) / (2 * r[:n + 1] ** 2) - eps)
            h2 = dr * dr
            w = 1.0 - h2 / 12.0 * f
            g = np.zeros(n + 1)
            g[0], g[1] = u[0] * w[0], u[1] * w[1]
            for i in range(1, n):
                g[i + 1] = 2 * g[i] - g[i - 1] + h2 * f[i] * u[i]
                u[i + 1] = g[i + 1] / w[i + 1]
            u[n + 1:] = 0.0
        nrm = np.sqrt(np.sum(u[:self.iR] ** 2) * dr)
        return u / max(nrm, 1e-300)

    # ---- delta_eps_l(eps): one-center correction for unit local wave ----
    def delta_eps(self, l_i, eps, u_i):
        dr = self.chn.mesh.dr
        tot = 0.0
        for g in self.chn.groups:
            L = g["L"]
            for lp in range(abs(l_i - L), l_i + L + 1):
                if (l_i + L + lp) % 2:
                    continue
                tjp2 = threej000_sq(l_i, L, lp)
                if tjp2 == 0.0:
                    continue
                # Wigner-Eckart reduction (audited against the exact closure
                # identity and the plane-wave/free limit of corepol.py):
                # sum over m_c,m_x,M,m' of |<lp m'|V_b|l_i m_i>|^2
                #   = Rad^2 * W * (2lp+1) * 3j(l_i,L,lp;000)^2 / (2L+1)^2,
                # independent of m_i.
                pref = g["W"] * (2 * lp + 1) * tjp2 / (2 * L + 1) ** 2
                eps_p, U_p = self.inter[lp]
                Rad = (g["wL"] * u_i[:, None]).T @ U_p * dr   # (nx, nj')
                K = kernel(eps, (eps_p + self.U0)[None, :],
                           g["Delta"][:, None], self.EF, self.T)
                tot += pref * np.sum(Rad ** 2 * K)
        return tot

    def tabulate(self, eps_lo, eps_hi, n_nodes=6):
        """Chebyshev-node tabulation of delta_eps_l^AE and delta_eps_l^PS
        on [eps_lo, eps_hi]."""
        x = np.cos(np.pi * (np.arange(n_nodes) + 0.5) / n_nodes)[::-1]
        nodes = 0.5 * (eps_hi - eps_lo) * x + 0.5 * (eps_hi + eps_lo)
        tab = {}
        for mode in ("ae", "ps"):
            for l in range(self.lmax_loc + 1):
                vals = []
                for e in nodes:
                    u_i = self.local_wave(l, e, mode=mode)
                    vals.append(self.delta_eps(l, e, u_i))
                tab[(mode, l)] = np.array(vals)
        self._dnodes = (nodes, tab)
        return nodes, tab

    def delta_interp(self, l, eps, mode="ae"):
        nodes, tab = self._dnodes
        e = np.clip(eps, nodes[0], nodes[-1])
        return np.interp(e, nodes, tab[(mode, l)])

    # ---- local weights of a Bloch state in each site sphere ----
    def local_weights(self, kfrac, Gints, cre, cim, wmin=1e-5):
        c = np.asarray(cre) + 1j * np.asarray(cim)
        G = np.asarray(Gints, dtype=float)
        keep = np.abs(c) ** 2 > wmin
        c, G = c[keep], G[keep]
        kcart = self.B @ np.asarray(kfrac, dtype=float)
        pvec = kcart[None, :] + (self.B @ G.T).T          # (nG,3)
        pn = np.linalg.norm(pvec, axis=1)
        phat = pvec / np.maximum(pn, 1e-30)[:, None]
        cosGG = np.clip(phat @ phat.T, -1.0, 1.0)
        r = self.chn.mesh.r[:self.iR]
        dr = self.chn.mesh.dr
        # radial overlap integrals S_l[G,G'] = int j_l(p r) j_l(p' r) r^2 dr
        W = np.zeros((self.n_atoms, self.lmax_loc + 1))
        for ia in range(self.n_atoms):
            phase = np.exp(2j * np.pi * (G @ self.tau_frac[ia]))
            cc = c * phase
            for l in range(self.lmax_loc + 1):
                jl = spherical_jn(l, np.outer(pn, r))      # (nG, nr)
                S = (jl * r[None, :] ** 2) @ jl.T * dr     # (nG,nG)
                Pl = np.polynomial.legendre.legval(
                    cosGG, [0] * l + [1])
                M = np.real(np.outer(cc, np.conj(cc)) * S * Pl)
                W[ia, l] = (4 * np.pi * (2 * l + 1) / self.Omega) * np.sum(M)
        return W

    # ---- public: the short-range AE-minus-pseudo one-center correction ----
    def correction(self, eps, kfrac, Gints, cre, cim):
        """sum_sites sum_l W_l [delta_l^AE(eps) - delta_l^PS(eps)] --
        to be ADDED to the plane-wave (pseudo) evaluation of the same
        diagram (corepol.SelfEnergy)."""
        W = self.local_weights(kfrac, Gints, cre, cim)
        dE = 0.0
        for ia in range(self.n_atoms):
            for l in range(self.lmax_loc + 1):
                dE += W[ia, l] * (self.delta_interp(l, eps, "ae")
                                  - self.delta_interp(l, eps, "ps"))
        return dE, W
