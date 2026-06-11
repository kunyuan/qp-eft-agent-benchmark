"""All-electron spherical LDA atom solver on a uniform radial mesh.

- KS-LDA, spin-restricted, Slater exchange + PW92 correlation: the same
  functional set as the pinned DFTK solid calculation [:lda_x, :lda_c_pw].
- Eigenstates per angular momentum l from a symmetric tridiagonal
  finite-difference Hamiltonian on u(r) = r R(r), u(0) = u(R_box) = 0.
  One diagonalization yields bound AND box-discretized continuum states on
  the same footing, mutually orthogonal, with consistent KS excitation
  energies  Delta_{c->x} = eps_x - eps_c.
- Element blind: neutral-atom configuration from Z_nuclear by aufbau;
  core = the innermost (Z - z_valence) electrons in energy order.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

# ---------------- LDA functional (Slater x + PW92 c) -----------------
def lda_x_potential(n):
    return -(3.0 / np.pi) ** (1.0 / 3.0) * n ** (1.0 / 3.0)

def lda_x_energy_density(n):
    return -0.75 * (3.0 / np.pi) ** (1.0 / 3.0) * n ** (1.0 / 3.0)

def _pw92_G(rs):
    A, a1, b1, b2, b3, b4 = 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294
    srs = np.sqrt(rs)
    Q0 = -2 * A * (1 + a1 * rs)
    Q1 = 2 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs)
    G = Q0 * np.log(1 + 1.0 / Q1)
    dQ0 = -2 * A * a1
    dQ1 = A * (b1 / srs + 2 * b2 + 3 * b3 * srs + 4 * b4 * rs)
    dG = dQ0 * np.log(1 + 1.0 / Q1) - Q0 * dQ1 / (Q1 * Q1 + Q1)
    return G, dG

def lda_c_pw92(n):
    n = np.maximum(n, 1e-30)
    rs = (3.0 / (4.0 * np.pi * n)) ** (1.0 / 3.0)
    ec, dec = _pw92_G(rs)
    vc = ec - rs * dec / 3.0
    return ec, vc

def vxc_of_n(n):
    ec, vc = lda_c_pw92(n)
    return lda_x_potential(n) + vc

# ---------------- mesh -----------------
class Mesh:
    def __init__(self, Z, rmax=50.0, dr=None):
        if dr is None:
            dr = min(0.004, 0.05 / max(Z, 1))
        self.dr = dr
        self.r = np.arange(1, int(round(rmax / dr)) + 1) * dr
        self.N = len(self.r)

    def integrate(self, f):
        return np.trapezoid(f, dx=self.dr)

def radial_states(mesh, V, l, n_states=None, emax=None):
    """-1/2 u'' + (V + l(l+1)/2r^2) u = e u ; returns eps, U (columns u_j,
    normalized int u^2 dr = 1)."""
    r, dr = mesh.r, mesh.dr
    main = 1.0 / dr ** 2 + V + l * (l + 1) / (2.0 * r * r)
    off = np.full(mesh.N - 1, -0.5 / dr ** 2)
    if emax is not None:
        lower = (V + l * (l + 1) / (2.0 * r * r)).min() - 1.0
        eps, W = eigh_tridiagonal(main, off, select="v",
                                  select_range=(lower, emax))
    else:
        eps, W = eigh_tridiagonal(main, off, select="i",
                                  select_range=(0, n_states - 1))
    U = W / np.sqrt(dr)  # mesh-orthonormal -> int u^2 dr = 1
    return eps, U

def hartree(mesh, n):
    """V_H for spherical density n(r)."""
    r, dr = mesh.r, mesh.dr
    f = 4 * np.pi * n * r * r
    q = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dr)])
    g = 4 * np.pi * n * r
    I2tot = np.trapezoid(g, dx=dr)
    I2 = I2tot - np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * dr)])
    return q / r + I2

AUFBAU = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1),
          (5, 0), (4, 2), (5, 1), (6, 0), (4, 3), (5, 2), (6, 1), (7, 0)]

def configuration(Z):
    cfg, ne = [], Z
    for (n, l) in AUFBAU:
        if ne <= 0:
            break
        occ = min(2 * (2 * l + 1), ne)
        cfg.append((n, l, occ))
        ne -= occ
    return cfg

class Atom:
    def __init__(self, Z, rmax=50.0, dr=None):
        self.Z = Z
        self.mesh = Mesh(Z, rmax=rmax, dr=dr)
        self.conf = configuration(Z)

    def scf(self, tol=1e-8, maxiter=300, mix=0.3, verbose=False):
        m, r = self.mesh, self.mesh.r
        n = self.Z ** 4 / np.pi * np.exp(-2.0 * self.Z * r) / self.Z  # guess
        lmax = max(l for (_, l, _) in self.conf)
        Vold = None
        for it in range(maxiter):
            V = -self.Z / r + hartree(m, n) + vxc_of_n(n)
            if Vold is not None:
                V = mix * V + (1 - mix) * Vold
            Vold = V
            nnew = np.zeros_like(n)
            self.levels = {}
            for l in range(lmax + 1):
                nsh = sum(1 for (_, ll, _) in self.conf if ll == l)
                eps, U = radial_states(m, V, l, n_states=nsh + 1)
                for (nn, ll, occ) in self.conf:
                    if ll != l:
                        continue
                    idx = nn - ll - 1
                    self.levels[(nn, ll)] = (eps[idx], U[:, idx], occ)
                    nnew += occ * (U[:, idx] / r) ** 2 / (4 * np.pi)
            dn = m.integrate(4 * np.pi * r * r * np.abs(nnew - n))
            n = nnew
            if verbose and it % 10 == 0:
                print("  scf iter", it, "dn =", dn)
            if dn < tol:
                break
        self.n = n
        self.V = -self.Z / r + hartree(m, n) + vxc_of_n(n)
        return self

    def split_core_valence(self, z_valence):
        ncore = self.Z - z_valence
        levels = sorted(self.levels.items(), key=lambda kv: kv[1][0])
        core, val, left = [], [], ncore
        for (nl, (e, u, occ)) in levels:
            if left >= occ:
                core.append({"nl": nl, "eps": e, "u": u, "occ": occ})
                left -= occ
            elif left == 0:
                val.append({"nl": nl, "eps": e, "u": u, "occ": occ})
            else:
                raise RuntimeError("core does not close at a shell boundary")
        return core, val

    def spectrum(self, l, emax):
        """All eigenstates of the converged KS potential with eps <= emax."""
        return radial_states(self.mesh, self.V, l, emax=emax)


if __name__ == "__main__":
    import sys, time
    Z = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    t0 = time.time()
    at = Atom(Z).scf()
    print("Z =", Z, " (%.1f s)" % (time.time() - t0))
    for (nl, (e, u, occ)) in sorted(at.levels.items(), key=lambda kv: kv[1][0]):
        rmean = at.mesh.integrate(u * u * at.mesh.r)
        print(f"  n={nl[0]} l={nl[1]} occ={occ}  eps = {e:.4f} Ha  <r> = {rmean:.3f} a0")
