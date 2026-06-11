"""Dynamical core-polarization quasiparticle correction (retardation only).

Derivation summary (full derivation in method.md):

Starting from the all-electron action, the core electrons at each lattice
site form a gapped, localized subsystem.  Integrating them out of the path
integral yields, beyond the static (mean-field + adiabatic-polarization)
pieces that the pinned frozen-core pseudopotential already contains, a
retarded valence-valence kernel mediated by the core's connected
density-density fluctuation propagator

    chi_c(r,r'; i nu) = - sum_b rho_b(r) rho_b*(r') * 2 Delta_b / (nu^2 + Delta_b^2)

with rho_b = phi_x* phi_c the transition density of core excitation
b = (c -> x) and Delta_b = eps_x - eps_c its excitation energy (KS
eigenvalue difference of the same converged neutral-atom KS potential).
The valence electron couples to channel b through the fluctuation potential

    V_b(r) = int d^3r' rho_b(r') / |r - r'|        (the DERIVED coupling)

The resulting second-order (Fan-Migdal) on-shell self-energy, with the full
static (adiabatic closure) limit subtracted because every static piece is
already counted inside the pinned pseudopotential, gives the quasiparticle
correction

    dE_nk = sum_{b,R} sum_{m'k'} |<m'k'| V_b(r-R) |nk>|^2 *
            [ (1-f')/(eps_nk - eps' - Delta_b)
              + f'/(eps_nk - eps' + Delta_b)
              + 1/Delta_b ]                                  (*)

The +1/Delta_b term removes the adiabatic core-polarization potential
(closure sum over ALL intermediates, no Fermi blocking - exactly the object
an atomic pseudopotential construction can absorb); what remains vanishes
identically when all valence energy differences are negligible against
Delta_b, i.e. it is pure retardation + Fermi-sea blocking, which no static
pseudopotential can contain.

Evaluation: intermediate valence states are plane waves with the
free-electron dispersion eps'(kappa) = EF + (kappa^2 - kF^2)/2 (NFE metals),
initial states are the pinned DFTK Bloch states |nk> = sum_G c_G e^{i(k+G)r}.
Using <x m_x| e^{-ip.r} |c m_c> summed over m's and spin,

    sum_{m,spin} V~_b(p) V~_b*(p') = (16 pi^2/(p^2 p'^2)) *
        sum_L W_bL F_bL(p) F_bL(p') P_L(p^.p'^)
    W_bL = 2 (2l_c+1)(2l_x+1)(2L+1) threej(l_c,L,l_x;0,0,0)^2
    F_bL(p) = int dr u_x(r) j_L(p r) u_c(r)

so that with p_G = q - G (all in atomic units, energies Ha)

    dE_nk = (1/Omega) int d^3q/(2pi)^3 sum_{G,G'} c_G c_G'* S_GG'
            (16 pi^2/(p_G^2 p_G'^2)) sum_{b,L} W_bL F_bL(p_G) F_bL(p_G')
            P_L(p^_G . p^_G') K(eps_nk, eps'(|k+q|), Delta_b)

    K(e, e', D) = (1-f')/(e-e'-D) + f'/(e-e'+D) + 1/D
    S_GG' = sum_atoms e^{i(G-G').tau_j}    (S_GG = n_atoms)

Everything below is element blind: inputs are Z_nuclear, z_valence and the
crystal data from element_config.json only.
"""
import numpy as np
from scipy.special import spherical_jn
from atom import Atom

HARTREE_EV = 27.211386245988


# ---------------- Wigner 3j (m=0) squared ----------------
def threej000_sq(l1, l2, l3):
    """(l1 l2 l3; 0 0 0) 3j symbol squared."""
    J = l1 + l2 + l3
    if J % 2 == 1:
        return 0.0
    if (abs(l1 - l2) > l3) or (l3 > l1 + l2):
        return 0.0
    g = J // 2
    from math import lgamma
    def lf(n):  # log n!
        return lgamma(n + 1)
    logv = (lf(J - 2 * l1) + lf(J - 2 * l2) + lf(J - 2 * l3) - lf(J + 1)) \
        + 2 * (lf(g) - lf(g - l1) - lf(g - l2) - lf(g - l3))
    return np.exp(logv)


# ---------------- channel construction ----------------
class CoreChannels:
    """All core fluctuation channels b=(c->x,L) of the neutral atom,
    with multipole form factors F_bL(p) tabulated on a p-grid."""

    def __init__(self, Z, z_valence, emax=60.0, Lmax=2, rmax=50.0,
                 pmax=18.0, npgrid=451):
        self.Z, self.zval = Z, z_valence
        self.at = Atom(Z, rmax=rmax).scf()
        core, val = self.at.split_core_valence(z_valence)
        self.core_shells = core
        occupied_core_nl = {c["nl"] for c in core}
        mesh = self.at.mesh
        r, dr = mesh.r, mesh.dr
        self.pgrid = np.linspace(0.0, pmax, npgrid)
        self.dp = self.pgrid[1] - self.pgrid[0]

        lx_needed = sorted({lx for c in core
                            for lx in range(0, c["nl"][1] + Lmax + 1)
                            if any((abs(c["nl"][1] - lx) <= L <= c["nl"][1] + lx)
                                   and ((c["nl"][1] + lx + L) % 2 == 0)
                                   for L in range(Lmax + 1))})
        spectra = {lx: self.at.spectrum(lx, emax) for lx in lx_needed}

        # channels grouped by multipole L: lists of (Delta, W, F_row_index)
        chans = {L: {"Delta": [], "W": [], "F": []} for L in range(Lmax + 1)}
        diag = []
        for L in range(Lmax + 1):
            # J_L(p,r) once per L
            JL = spherical_jn(L, np.outer(self.pgrid, r))
            for c in core:
                n_c, l_c = c["nl"]
                u_c = c["u"]
                for lx in range(abs(0), l_c + Lmax + 1):
                    if (l_c + lx + L) % 2 == 1:
                        continue
                    if not (abs(l_c - lx) <= L <= l_c + lx):
                        continue
                    tj2 = threej000_sq(l_c, L, lx)
                    if tj2 == 0.0:
                        continue
                    W = 2.0 * (2 * l_c + 1) * (2 * lx + 1) * (2 * L + 1) * tj2
                    eps_x, U_x = spectra[lx]
                    # F rows for all x at once: (npgrid, nx)
                    F = (JL * u_c[None, :]) @ U_x * dr
                    for j in range(len(eps_x)):
                        nl_x = (j + lx + 1, lx)
                        if nl_x in occupied_core_nl:
                            continue
                        Delta = eps_x[j] - c["eps"]
                        if Delta <= 0:
                            continue
                        Frow = F[:, j].copy()
                        if L == 0:
                            # exact: monopole transition density carries zero
                            # net charge (orthogonal orbitals) => F_0(0)=0;
                            # remove the numerical orthogonality defect.
                            Frow[0] = 0.0
                        chans[L]["Delta"].append(Delta)
                        chans[L]["W"].append(W)
                        chans[L]["F"].append(Frow)
                        diag.append(((n_c, l_c), lx, j, L, Delta, W))
        self.channel_meta = diag
        self.ch = {}
        for L, d in chans.items():
            if not d["Delta"]:
                continue
            self.ch[L] = {"Delta": np.array(d["Delta"]),
                          "W": np.array(d["W"]),
                          "F": np.array(d["F"])}  # (nch, npgrid)
        # prune negligible channels: importance ~ W * int F^2/p^2 dp / Delta^2
        self._prune(keep=0.9995)

    def _prune(self, keep):
        scores, tags = [], []
        p2 = np.maximum(self.pgrid, self.dp) ** 2
        for L, d in self.ch.items():
            s = d["W"] * np.sum(d["F"] ** 2 / p2[None, :], axis=1) * self.dp \
                / d["Delta"] ** 2
            scores.append(s)
            tags.append(np.full(len(s), L))
        allsc = np.concatenate(scores)
        order = np.argsort(allsc)[::-1]
        cum = np.cumsum(allsc[order])
        ncut = int(np.searchsorted(cum, keep * cum[-1])) + 1
        chosen = set(order[:ncut])
        # map back
        idx0 = 0
        for L in list(self.ch.keys()):
            n = len(self.ch[L]["Delta"])
            sel = [i for i in range(n) if (idx0 + i) in chosen]
            if not sel:
                del self.ch[L]
            else:
                for key in ("Delta", "W"):
                    self.ch[L][key] = self.ch[L][key][sel]
                self.ch[L]["F"] = self.ch[L]["F"][sel]
            idx0 += n
        self.nchannels = sum(len(d["Delta"]) for d in self.ch.values())

    # ---- diagnostics ----
    def alpha0_dipole(self):
        """Static dipole core polarizability from the channel set.
        alpha_zz = 2 sum_b sum_{m,spin} |<x|z|c>|^2 / Delta_b, and the
        m,spin-summed |<z>|^2 equals lim_{p->0} W_b1 F_1(p)^2/p^2 = W_b1 d^2/9
        (F_1 -> p d/3). So alpha = sum_b 2 W_b1 d_b^2/(9 Delta_b)."""
        if 1 not in self.ch:
            return 0.0
        d = self.ch[1]
        p1 = self.pgrid[1]
        dip = 3.0 * d["F"][:, 1] / p1
        return float(np.sum(2.0 * d["W"] * dip ** 2 / (9.0 * d["Delta"])))

    def interp_F(self, L, p):
        """Linear interpolation of F_bL onto arbitrary |p| array.
        Returns (nch, *p.shape)."""
        d = self.ch[L]
        x = np.clip(p / self.dp, 0, len(self.pgrid) - 1 - 1e-9)
        i0 = x.astype(int)
        frac = x - i0
        return d["F"][:, i0] * (1 - frac) + d["F"][:, i0 + 1] * frac


# ---------------- kernel ----------------
def kernel(eps, epsp, Delta, EF, T):
    """K(e,e',D) with FD occupation at temperature T (Ha). Shapes:
    eps scalar, epsp (...), Delta (nch,1...) broadcastable."""
    x = (epsp - EF) / max(T, 1e-12)
    f = np.where(x > 40, 0.0, np.where(x < -40, 1.0, 1.0 / (1.0 + np.exp(np.clip(x, -40, 40)))))
    de = eps - epsp
    return (1.0 - f) / (de - Delta) + f / (de + Delta) + 1.0 / Delta


# ---------------- integration grids ----------------
def radial_grid(pmax):
    """Composite Gauss-Legendre radial grid on [0,pmax]."""
    edges = [0.0, 0.4, 1.0, 2.0, 3.5, 6.0, 10.0, pmax]
    nseg = [10, 12, 14, 14, 12, 10, 8]
    xs, ws = [], []
    for (a, b, n) in zip(edges[:-1], edges[1:], nseg):
        x, w = np.polynomial.legendre.leggauss(n)
        xs.append(0.5 * (b - a) * x + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(xs), np.concatenate(ws)


class SelfEnergy:
    def __init__(self, channels, EF, volume, n_electrons, n_atoms,
                 positions_frac, recip_cols, T=0.001, nmu=48, chunk=400):
        self.chn = channels
        self.EF, self.Omega, self.T = EF, volume, T
        self.n_atoms = n_atoms
        self.tau_frac = np.array(positions_frac)      # (nat,3)
        self.B = np.array(recip_cols).T               # columns b1,b2,b3
        n_dens = n_electrons / volume
        self.kF = (3.0 * np.pi ** 2 * n_dens) ** (1.0 / 3.0)
        self.pr, self.pw = radial_grid(channels.pgrid[-1])
        self.mu, self.muw = np.polynomial.legendre.leggauss(nmu)
        self.chunk = chunk

    def epsp(self, kappa):
        return self.EF + 0.5 * (kappa ** 2 - self.kF ** 2)

    # ----- diagonal G = G' term (2D integral, azimuthal symmetry) -----
    def _diag_one(self, eps, a):
        """(1/Omega)*n_atoms * 4 * int dp dmu (1/p^2) sum_bL W F^2 K  for
        |c_G|^2 = 1; a = |k+G|."""
        p = self.pr[:, None]
        mu = self.mu[None, :]
        kappa = np.sqrt(np.maximum(a * a + p * p + 2 * a * p * mu, 0.0))
        ep = self.epsp(kappa)                              # (np, nmu)
        total = 0.0
        for L, d in self.chn.ch.items():
            F = self.chn.interp_F(L, self.pr)              # (nch, np)
            WF2 = d["W"][:, None] * F ** 2                 # (nch, np)
            nch = len(d["Delta"])
            for i0 in range(0, nch, self.chunk):
                sl = slice(i0, min(i0 + self.chunk, nch))
                K = kernel(eps, ep[None, :, :], d["Delta"][sl][:, None, None],
                           self.EF, self.T)                # (nc, np, nmu)
                total += np.einsum('cp,cpm,p,m->', WF2[sl], K,
                                   self.pw / self.pr ** 2, self.muw)
        return 4.0 * self.n_atoms / self.Omega * total

    # ----- cross terms G != G' (3D integral) -----
    def _cross_one(self, eps, avec, dG, nrad=40, nth=20, nph=20):
        """Integral for pair (G,G'): grid centered at G (p = q-G);
        avec = k+G (cart), dG = G-G' (cart, so p' = p + dG).
        Returns I = (1/Omega) int d3q/(2pi)^3 (16pi^2/(p^2 p'^2))
                    sum W F(p)F(p') P_L K   (without c_G c_G'* S factor)."""
        x, w = np.polynomial.legendre.leggauss(nth)
        th_mu, th_w = x, w
        ph = (np.arange(nph) + 0.5) * (2 * np.pi / nph)
        phw = 2 * np.pi / nph
        # radial: reuse composite grid but coarser
        redges = [0.0, 0.5, 1.5, 3.0, 6.0, 12.0]
        rn = [8, 10, 10, 8, 6]
        rs, rw = [], []
        for (a, b, n) in zip(redges[:-1], redges[1:], rn):
            xx, ww = np.polynomial.legendre.leggauss(n)
            rs.append(0.5 * (b - a) * xx + 0.5 * (a + b))
            rw.append(0.5 * (b - a) * ww)
        rs, rw = np.concatenate(rs), np.concatenate(rw)
        # axis: z along dG, need avec too -> build full vectors
        ez = dG / np.linalg.norm(dG) if np.linalg.norm(dG) > 1e-12 else np.array([0, 0, 1.0])
        tmp = np.array([1.0, 0, 0]) if abs(ez[0]) < 0.9 else np.array([0, 1.0, 0])
        ex = np.cross(tmp, ez); ex /= np.linalg.norm(ex)
        ey = np.cross(ez, ex)
        st = np.sqrt(1 - th_mu ** 2)
        # p vectors: (nr, nth, nph, 3)
        dirs = (st[:, None, None] * np.cos(ph)[None, :, None] * ex[None, None, :]
                + st[:, None, None] * np.sin(ph)[None, :, None] * ey[None, None, :]
                + th_mu[:, None, None] * ez[None, None, :])      # (nth, nph, 3)
        pvec = rs[:, None, None, None] * dirs[None, :, :, :]
        pnorm = rs[:, None, None]
        ppvec = pvec + dG[None, None, None, :]
        ppn = np.linalg.norm(ppvec, axis=-1)
        cosang = np.einsum('rtpx,rtpx->rtp', pvec, ppvec) / \
            np.maximum(pnorm * ppn, 1e-30)
        kq = pvec + avec[None, None, None, :]
        kappa = np.linalg.norm(kq, axis=-1)
        ep = self.epsp(kappa)
        wgt = (rw[:, None, None] * th_w[None, :, None] * phw)    # measure /  (r^2 cancels)
        # full weight: p^2 dp dmu dph / p^2 / p'^2 -> dp dmu dph / p'^2
        wgt = wgt / np.maximum(ppn, 1e-30) ** 2
        total = 0.0
        for L, d in self.chn.ch.items():
            Fp = self.chn.interp_F(L, rs)                        # (nch, nr)
            Fpp = self.chn.interp_F(L, ppn)                      # (nch, nr,nth,nph)
            # Legendre P_L(cosang)
            if L == 0:
                PL = np.ones_like(cosang)
            elif L == 1:
                PL = cosang
            else:
                PL = 0.5 * (3 * cosang ** 2 - 1)
            nch = len(d["Delta"])
            for i0 in range(0, nch, 200):
                sl = slice(i0, min(i0 + 200, nch))
                K = kernel(eps, ep[None], d["Delta"][sl][:, None, None, None],
                           self.EF, self.T)
                contrib = np.einsum('cr,crtp,crtp,rtp->',
                                    d["W"][sl][:, None] * Fp[sl], Fpp[sl],
                                    K, wgt * PL)
                total += contrib
        # prefactor: (1/Omega)(1/(2pi)^3) * 16 pi^2 = 2/(pi Omega)
        return 2.0 / (np.pi * self.Omega) * total

    # ----- public: correction for one Bloch state -----
    def correction(self, eps, kfrac, Gints, cre, cim, wmin_diag=1e-4,
                   wmin_cross=2e-3, max_cross=36, do_cross=True):
        c = np.asarray(cre) + 1j * np.asarray(cim)
        G = np.asarray(Gints, dtype=float)                # (nG,3) frac
        w = np.abs(c) ** 2
        kcart = self.B @ np.asarray(kfrac, dtype=float)
        dE = 0.0
        keep = np.where(w > wmin_diag)[0]
        for i in keep:
            a = np.linalg.norm(kcart + self.B @ G[i])
            dE += w[i] * self._diag_one(eps, a)
        if do_cross:
            # pairs by |c_i c_j|
            idx = np.argsort(w)[::-1][:12]
            pairs = [(i, j) for ii, i in enumerate(idx) for j in idx[ii + 1:]
                     if abs(c[i] * c[j]) > wmin_cross]
            pairs = sorted(pairs, key=lambda ij: -abs(c[ij[0]] * c[ij[1]]))[:max_cross]
            for (i, j) in pairs:
                Sgg = np.sum(np.exp(2j * np.pi * (self.tau_frac @ (G[i] - G[j]))))
                pref = c[i] * np.conj(c[j]) * Sgg
                if abs(pref) < 1e-12:
                    continue
                # symmetrize the two centerings
                I1 = self._cross_one(eps, kcart + self.B @ G[i],
                                     self.B @ (G[i] - G[j]))
                I2 = self._cross_one(eps, kcart + self.B @ G[j],
                                     self.B @ (G[j] - G[i]))
                dE += 2.0 * np.real(pref) * 0.5 * (I1 + I2)
        return dE
