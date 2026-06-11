"""Full-box one-center evaluation of the core-polarization Fan diagram for
actual Bloch states, with optional OPW-style all-electron augmentation of the
external state inside r_m.

dE_nk = sum_sites sum_l (4pi(2l+1)/Omega) sum_GG' Re(ct_G ct*_G') P_l(cos_GG')
        sum_{g,lp} pref sum_{x,j} R_G[x,j] R_G'[x,j] K(eps, eps'_j+U0, Delta_x)

R_G[x,j] = int_0^box  phi_l(p_G; r) w_Lx(r) U_j(r) dr
phi_l(p,r) = j_l(p r) r                      (mode 'pw')
           = j_l(p r) r            r > r_m   (mode 'ae')
             beta(p) uAE_l(epsf;r) r < r_m,  beta = j_l(p r_m) r_m / uAE(r_m)
epsf = eps_band - U0 (atomic frame), uAE = regular Numerov solution of the
atomic KS potential.  Intermediates U_j: free box states (mode 'pw' audit) or
atomic KS spectrum (mode 'ae').  Closure-exact handling of the +1/Delta
counterterm; dynamic kernel K - 1/Delta summed over the discretized box states.
"""
import sys, json
sys.path.insert(0, '/Users/kunchen/project/qp_eft_L4_fable5_run/solution')
import numpy as np
from scipy.special import spherical_jn
from onecenter import AtomChannels, threej000_sq
from atom import radial_states

HARTREE_EV = 27.211386245988


def kernel_dyn(eps, epsp, Delta, EF, T):
    """K - 1/Delta (decays for high eps')."""
    x = np.clip((epsp - EF) / max(T, 1e-12), -40.0, 40.0)
    f = 1.0 / (1.0 + np.exp(x))
    de = eps - epsp
    return (1.0 - f) / (de - Delta) + f / (de + Delta)


class FullBox:
    def __init__(self, Z, zval, ks, lmax=4, r_match=2.6, emax_int=None):
        self.ch = AtomChannels(Z, zval)
        self.lmax = lmax
        self.r_m = r_match
        self.EF, self.T = ks['EF_Ha'], ks['temperature_Ha']
        self.Omega = ks['volume_bohr3']
        self.n_atoms = ks['n_atoms']
        self.tau_frac = np.array(ks['positions_frac'])
        self.B = np.array(ks['recip_cols']).T
        kF = (3 * np.pi ** 2 * ks['n_electrons'] / self.Omega) ** (1. / 3.)
        self.kF = kF
        self.U0 = self.EF - 0.5 * kF ** 2
        mesh = self.ch.mesh
        self.im = int(self.r_m / mesh.dr)
        self.r_m = mesh.r[self.im]   # snap to mesh
        # downsampled radial grid (indices + weights) for Rad integrals
        dr = mesh.dr
        i_fine = np.arange(0, int(8.0 / dr))
        stride = 8
        i_coarse = np.arange(int(8.0 / dr), mesh.N - stride, stride)
        self.idx = np.concatenate([i_fine, i_coarse])
        self.wts = np.concatenate([np.full(len(i_fine), dr),
                                   np.full(len(i_coarse), dr * stride)])
        self.rs = mesh.r[self.idx]
        self.in_mask = self.idx < self.im
        # intermediates per lp: (eps, U) atomic and free
        emax = emax_int or self.ch.emax
        lmax_lp = lmax + self.ch.Lmax
        zeroV = np.zeros(mesh.N)
        self.int_at, self.int_free = {}, {}
        for lp in range(lmax_lp + 1):
            if lp in self.ch.spec:
                e, U = self.ch.spec[lp]
            else:
                e, U = radial_states(mesh, self.ch.at.V, lp, emax=emax)
            self.int_at[lp] = (e, U[self.idx, :])
            e2, U2 = radial_states(mesh, zeroV, lp, emax=emax)
            self.int_free[lp] = (e2, U2[self.idx, :])
        # channel radial couplings on the downsampled grid
        for g in self.ch.groups:
            g['wLd'] = g['wL'][self.idx, :]
        # core orbitals for closure exclusion (atomic intermediates)
        self.core_by_l = {}
        self.coreT_by_l = {}
        r = mesh.r
        for c in self.ch.core:
            self.core_by_l.setdefault(c['nl'][1], []).append(c['u'][self.idx])
            # kinetic centroid <T> = eps_c - <V> (free-energy centroid of the
            # core orbital's weight in the plane-wave continuum)
            Tc = c['eps'] - np.trapezoid(c['u'] ** 2 * self.ch.at.V, dx=dr)
            self.coreT_by_l.setdefault(c['nl'][1], []).append(
                (c['u'][self.idx], Tc))

    # ---- log-derivative of the AE regular solution at r_m ----
    def _ae_logderiv(self, l, epsf):
        u = self._uae_raw(l, epsf)
        # centered derivative at im (raw fine-mesh array)
        dr = self.ch.mesh.dr
        return (u[self.im + 1] - u[self.im - 1]) / (2 * dr * u[self.im]), u

    def _uae_raw(self, l, epsf):
        mesh = self.ch.mesh
        r, dr = mesh.r, mesh.dr
        n = self.im + 2
        V = self.ch.at.V
        u = np.zeros(n + 1)
        u[0], u[1] = r[0] ** (l + 1), r[1] ** (l + 1)
        f = 2.0 * (V[:n + 1] + l * (l + 1) / (2 * r[:n + 1] ** 2) - epsf)
        h2 = dr * dr
        w = 1.0 - h2 / 12.0 * f
        g = np.zeros(n + 1)
        g[0], g[1] = u[0] * w[0], u[1] * w[1]
        for i in range(1, n):
            g[i + 1] = 2 * g[i] - g[i - 1] + h2 * f[i] * u[i]
            u[i + 1] = g[i + 1] / w[i + 1]
        return u

    def match_epsf(self, l, nu_target, nnodes_req, elo=-3.0, ehi=2.5, nscan=56):
        """Find epsf with the AE regular solution's log-derivative at r_m equal
        to nu_target, on the branch with the required interior node count."""
        es = np.linspace(elo, ehi, nscan)
        best = None
        prev = None
        for e in es:
            nu, u = self._ae_logderiv(l, e)
            nn = int(np.sum(u[5:self.im - 1] * u[6:self.im] < 0))
            if prev is not None:
                e0, nu0, nn0 = prev
                # logderiv decreases with E between poles; root if target crossed
                if nn0 == nnodes_req and nn == nnodes_req and \
                        (nu0 - nu_target) * (nu - nu_target) <= 0 and nu0 > nu:
                    a, bb = e0, e
                    for _ in range(48):
                        m = 0.5 * (a + bb)
                        num, _ = self._ae_logderiv(l, m)
                        if num > nu_target:
                            a = m
                        else:
                            bb = m
                    best = 0.5 * (a + bb)
                    break
            prev = (e, nu, nn)
        return best

    # ---- AE inside wave (regular solution at atomic-frame energy) ----
    def uae(self, l, epsf):
        mesh = self.ch.mesh
        r, dr = mesh.r, mesh.dr
        n = self.im + 2
        V = self.ch.at.V
        u = np.zeros(n + 1)
        u[0], u[1] = r[0] ** (l + 1), r[1] ** (l + 1)
        f = 2.0 * (V[:n + 1] + l * (l + 1) / (2 * r[:n + 1] ** 2) - epsf)
        h2 = dr * dr
        w = 1.0 - h2 / 12.0 * f
        g = np.zeros(n + 1)
        g[0], g[1] = u[0] * w[0], u[1] * w[1]
        for i in range(1, n):
            g[i + 1] = 2 * g[i] - g[i - 1] + h2 * f[i] * u[i]
            u[i + 1] = g[i + 1] / w[i + 1]
        # rescale to avoid overflow issues: normalize value at r_m to 1
        val = u[self.im]
        full = np.zeros(mesh.N)
        full[:n + 1] = u / val
        return full[self.idx]          # downsampled, value 1 at r_m

    # ---- per-state evaluation ----
    def correction(self, eps, kfrac, Gints, cre, cim, mode='ae', nG=14,
                   verbose=False, with_ct=True, pauli=False):
        """with_ct: include the +1/Delta closure counterterm (adiabatic-CPP
        removal).  pauli: subtract the core-subspace weight from the
        intermediate k' sum (Pauli exclusion of core states), evaluating the
        removed weight at the core orbital's kinetic-centroid energy."""
        c = np.asarray(cre) + 1j * np.asarray(cim)
        G = np.asarray(Gints, dtype=float)
        order = np.argsort(np.abs(c))[::-1][:nG]
        c, G = c[order], G[order]
        kcart = self.B @ np.asarray(kfrac, dtype=float)
        pvec = kcart[None, :] + (self.B @ G.T).T
        pn = np.linalg.norm(pvec, axis=1)
        phat = pvec / np.maximum(pn, 1e-30)[:, None]
        cosGG = np.clip(phat @ phat.T, -1.0, 1.0)
        # modes: 'pw' (free ext + free int), 'ae' (AE ext + atomic int),
        #        'aefree' (AE ext + free int)
        inter = self.int_at if mode == 'ae' else self.int_free
        ae_ext = mode in ('ae', 'aefree')
        rs, wts = self.rs, self.wts
        dE = 0.0
        ncore_l = {}
        for cs in self.ch.core:
            ncore_l[cs['nl'][1]] = ncore_l.get(cs['nl'][1], 0) + 1
        for ia in range(self.n_atoms):
            ct = c * np.exp(2j * np.pi * (G @ self.tau_frac[ia]))
            cc = np.real(np.outer(ct, np.conj(ct)))
            for l in range(self.lmax + 1):
                # external profiles phi_l(p_G; r) (nG, nr)
                phi = spherical_jn(l, np.outer(pn, rs)) * rs[None, :]
                if ae_ext:
                    # target log-derivative of the (l,m)-aggregated pseudo
                    # profile r*sqrt(rho_l) at r_m: self-aligning matching
                    Pl_ = np.polynomial.legendre.legval(cosGG, [0] * l + [1])
                    h = 1e-3
                    rho = lambda rr: float(np.einsum(
                        'ab,a,b->', cc * Pl_, spherical_jn(l, pn * rr),
                        spherical_jn(l, pn * rr)))
                    r0 = self.r_m
                    rho0 = rho(r0)
                    nu_t = 1.0 / r0 + (rho(r0 + h) - rho(r0 - h)) / (4 * h * max(rho0, 1e-300))
                    nnod = ncore_l.get(l, 0)
                    epsf = self.match_epsf(l, nu_t, nnod)
                    if epsf is None:
                        epsf = eps - self.U0   # fallback
                    if verbose:
                        print(f"    l={l}: nu_target={nu_t:+.4f} -> epsf={epsf:+.4f}")
                    uin = self.uae(l, epsf)        # value 1 at r_m
                    jm = spherical_jn(l, pn * self.r_m) * self.r_m
                    phi = np.where(self.in_mask[None, :],
                                   jm[:, None] * uin[None, :], phi)
                Pl = np.polynomial.legendre.legval(cosGG, [0] * l + [1])
                wl = (4 * np.pi * (2 * l + 1) / self.Omega) * cc * Pl  # (nG,nG)
                for g in self.ch.groups:
                    L = g['L']
                    for lp in range(abs(l - L), l + L + 1):
                        if (l + L + lp) % 2:
                            continue
                        tjp2 = threej000_sq(l, L, lp)
                        if tjp2 == 0.0:
                            continue
                        pref = g['W'] * (2 * lp + 1) * tjp2 / (2 * L + 1) ** 2
                        eps_p, U_p = inter[lp]
                        # R_G[x,j]: (nG, nx, nj)
                        A = phi * wts[None, :]                  # (nG, nr)
                        nx = g['wLd'].shape[1]
                        Kd = kernel_dyn(eps, (eps_p + self.U0)[None, :],
                                        g['Delta'][:, None], self.EF, self.T)
                        # dynamic part
                        S = 0.0
                        Rg = np.empty((len(pn), nx, U_p.shape[1]))
                        for ix in range(nx):
                            M = (g['wLd'][:, ix][:, None] * U_p)  # (nr, nj)
                            Rg[:, ix, :] = A @ M
                        S += np.einsum('ab,axj,bxj,xj->', wl, Rg, Rg, Kd)
                        if pauli and mode != 'ae':
                            # remove core-subspace weight from the k' sum
                            for (uc, ek) in self.coreT_by_l.get(lp, []):
                                v = (phi * wts[None, :]) @ (uc[:, None] * g['wLd'])  # (nG,nx)
                                Kc = kernel_dyn(eps, ek + self.U0,
                                                g['Delta'], self.EF, self.T)  # (nx,)
                                S -= np.einsum('ab,ax,bx,x->', wl, v, v, Kc)
                        if with_ct:
                            # closure-exact +1/Delta piece
                            WW = (g['wLd'] ** 2 / g['Delta'][None, :])  # (nr, nx)
                            cw = WW.sum(axis=1)                          # (nr,)
                            Cab = (phi * cw[None, :] * wts[None, :]) @ phi.T
                            if mode == 'ae' or pauli:
                                for uc_ in self.core_by_l.get(lp, []):
                                    v = (phi * wts[None, :]) @ (uc_[:, None] * g['wLd'])
                                    Cab -= np.einsum('ax,bx->ab',
                                                     v / g['Delta'][None, :], v)
                            S += np.sum(wl * Cab)
                        dE += pref * S
        return dE


if __name__ == '__main__':
    import time
    ROOT = '/Users/kunchen/project/qp_eft_L4_fable5_run'
    elem = sys.argv[1] if len(sys.argv) > 1 else 'na'
    Z, zv = (11, 1) if elem == 'na' else (13, 3)
    ks = json.load(open(f'{ROOT}/scratch/{elem}_dump.json'))
    EF = ks['EF_Ha']
    fb = FullBox(Z, zv, ks)
    print(f"r_m={fb.r_m:.3f} U0={fb.U0:.4f} kF={fb.kF:.4f}")
    pts = {p['point_id']: p for p in ks['path']}
    sel = sorted(pts, key=lambda i: pts[i]['t'])
    tests = [sel[0], sel[len(sel) // 2], sel[-1]]
    for pid in tests:
        pt = pts[pid]
        b = min(pt['bands'], key=lambda bb: bb['eps_Ha'])
        eps = b['eps_Ha']
        for mode in ('pw', 'ae'):
            t0 = time.time()
            dE = fb.correction(eps, pt['k_frac'], b['G'], b['c_re'], b['c_im'],
                               mode=mode)
            print(f"pid={pid} t={pt['t']:.3f} eps-EF={(eps-EF)*HARTREE_EV:+7.3f}  "
                  f"{mode}: dE={dE*HARTREE_EV:+8.4f} eV ({time.time()-t0:.0f}s)")
