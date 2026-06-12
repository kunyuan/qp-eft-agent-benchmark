#!/usr/bin/env python3
"""GUIDED SESSION: derived hybridization vertex from the two-body-induced
transfer-amplitude class (derivation_notes.md sec. 5).

Per core shell c:
  rho_spec,c(r) = (N_c-1) u_c(r)^2/(4 pi r^2)        (same-shell spectators)
  V_spec,c     = Poisson[rho_spec,c]                  (radial)
  U_c          = int V_spec,c u_c^2 dr                (orbital average / diagonal)
  G_c(kappa)   = 4 pi int j_l(kappa r) V_spec,c u_c r dr
  cf_c(kappa)  = 4 pi int j_l(kappa r) u_c r dr
  M_c(kappa)   = G_c - U_c cf_c                       (zero-diagonal vertex)

lambda_nk (production, incoherent-G, drop-in for E0*cf -> M_c):
  lambda_nk = sum_c n_at (2l+1)/(4 pi Omega) sum_G |c_G|^2 M_c(|k+G|)^2 / (eps-eps_c)^2

Outputs: Li closed-form cross-check, per-shell U_c and E0_eff at gate states,
Li Gamma gate, Na/Al ARPES RMSE with the derived vertex. NO calibration.
"""
import json, sys
import numpy as np
import pandas as pd
from scipy.special import spherical_jn

sys.path.insert(0, '/Users/kunchen/project/qp_eft_L4_fable5_run3/scratch')
from atom import RadialAtom, occupations, split_core_valence

HA_EV = 27.211386245988
PACKET = '/Users/kunchen/project/qp_eft_L4_fable5_run3/packet'
SCR = '/Users/kunchen/project/qp_eft_L4_fable5_run3/scratch'
DEMANDED = {('Li', 1, 0): 3.35, ('Na', 2, 0): 3.18, ('Al', 2, 0): 2.65}


def radial_hartree(r, rho):
    """V_H(r) for spherical density rho(r) [e/bohr^3]: q_in/r + int_r^inf rho r' dr'."""
    rho_r2 = 4*np.pi*rho*r*r
    q_in = np.concatenate([[0.0], np.cumsum(0.5*(rho_r2[1:]+rho_r2[:-1])*np.diff(r))])
    rho_r = 4*np.pi*rho*r
    o_int = np.concatenate([[0.0], np.cumsum(0.5*(rho_r[1:]+rho_r[:-1])*np.diff(r))])
    return q_in/r + (o_int[-1] - o_int)


def shell_vertex(r, u, l, occ, kmax=14.0, nk=900):
    """Derived per-shell vertex tables: ks, cf, G, M, U."""
    nspec = occ - 1.0                      # same-shell spectators
    rho_spec = nspec*u*u/(4*np.pi*r*r)
    Vs = radial_hartree(r, rho_spec)
    U = np.trapezoid(Vs*u*u, r)
    ks = np.linspace(1e-4, kmax, nk)
    cf = np.empty(nk); G = np.empty(nk)
    for i, k in enumerate(ks):
        jl = spherical_jn(l, k*r)
        cf[i] = 4*np.pi*np.trapezoid(jl*u*r, r)
        G[i] = 4*np.pi*np.trapezoid(jl*Vs*u*r, r)
    M = G - U*cf
    return dict(ks=ks, cf=cf, G=G, M=M, U=U)


def li_closed_form_check():
    zeta = 3.0 - 5.0/16.0
    r = np.geomspace(1e-7, 40.0, 6000)
    u = 2.0*zeta**1.5*r*np.exp(-zeta*r)            # u = r R_1s, normalized
    sv = shell_vertex(r, u, 0, 2.0)
    ks = sv['ks']
    G_cf = 4*np.pi*np.sqrt(zeta**3/np.pi)*(1.0/(ks**2+zeta**2) - 1.0/(ks**2+9*zeta**2)
                                           - 6*zeta**2/(ks**2+9*zeta**2)**2)
    chit = 8*np.sqrt(np.pi*zeta**5)/(ks**2+zeta**2)**2
    err_G = np.max(np.abs(sv['G']-G_cf))
    err_cf = np.max(np.abs(sv['cf']-chit))
    print(f"[Li closed-form check, hydrogenic zeta={zeta}]")
    print(f"  U = {sv['U']:.6f} Ha (closed form 5/8 zeta = {0.625*zeta:.6f})")
    print(f"  max|G_num - G_closed| = {err_G:.2e};  max|cf_num - chit| = {err_cf:.2e}")
    print(f"  E_eff(0) = {sv['M'][0]/sv['cf'][0]:+.4f} Ha "
          f"(closed form (33/81-5/8)zeta = {(33/81-0.625)*zeta:+.4f})")
    print(f"  E_eff(kmax) = {sv['M'][-1]/sv['cf'][-1]:+.4f} Ha (asymptote 3/8 zeta = {0.375*zeta:+.4f})")


def element_tables(el):
    cfg = json.load(open(f'{PACKET}/{el}/element_config.json'))
    Z = int(cfg['Z_nuclear']); zv = cfg['dft']['z_valence']
    at = RadialAtom(Z)
    occs = occupations(Z)
    res = at.scf(occs)
    core = split_core_valence(occs, zv)
    shells = []
    for (n, l, occ) in core:
        sv = shell_vertex(at.r, res['u'][(n, l)], l, occ)
        sv.update(n=n, l=l, occ=occ, eps=res['eps'][(n, l)])
        shells.append(sv)
    return cfg, shells


def lam_states(el, shells, vertex_key='M'):
    """Per-(ik,n) lambda with the chosen vertex table (production incoherent-G)."""
    d = f'{SCR}/dump_{el}'
    eF = float(open(f'{d}/ef.txt').read().strip())
    cell = open(f'{d}/cell.txt').read().split('\n')
    Omega = float(cell[0])
    n_at = len([ln for ln in cell[1:] if ln.strip()])
    bands = pd.read_csv(f'{d}/bands.csv')
    psi = pd.read_csv(f'{d}/psi.csv')
    out = []
    for (ik, n), grp in psi.groupby(['ik', 'n']):
        kg = grp[['kgx', 'kgy', 'kgz']].values
        w = grp['re'].values**2 + grp['im'].values**2
        kmag = np.linalg.norm(kg, axis=1)
        row = bands[(bands.ik == ik) & (bands.n == n)]
        eps = row.eps_Ha.values[0]; t = row.t.values[0]
        lam = 0.0; rec = dict(ik=ik, n=n, t=t, eps=eps)
        for sh in shells:
            Mv = np.interp(kmag, sh['ks'], sh[vertex_key])
            cfv = np.interp(kmag, sh['ks'], sh['cf'])
            den2 = (eps - sh['eps'])**2
            num_M = np.sum(w*Mv**2)
            num_cf = np.sum(w*cfv**2)
            lc = n_at*(2*sh['l']+1)/(4*np.pi)*num_M/Omega/den2
            lam += lc
            rec[f"lam_{sh['n']}{'spdf'[sh['l']]}"] = lc
            rec[f"E0eff_{sh['n']}{'spdf'[sh['l']]}"] = np.sqrt(num_M/num_cf) if num_cf > 0 else np.nan
        rec['lam'] = lam
        out.append(rec)
    return pd.DataFrame(out), eF


def gates(el, df, eF):
    rows = []
    for ik, grp in df.groupby('ik'):
        grp = grp.sort_values('eps')
        occ = grp[grp.eps < eF]
        sel = grp.iloc[:len(occ)+1]
        Epred = (sel.eps.values - eF)/(1.0 + sel.lam.values)*HA_EV
        Eks = (sel.eps.values - eF)*HA_EV
        rows.append(dict(ik=ik, t=grp.t.values[0], Epred=Epred, Eks=Eks))
    pred = pd.DataFrame(rows)
    try:
        arp = pd.read_csv(f'{PACKET}/{el}/arpes_reference.csv')
    except FileNotFoundError:
        return None, None
    errs, errs_ks = [], []
    for _, row in arp.iterrows():
        cand = pred.iloc[(pred.t - row.t).abs().argsort().iloc[0]]
        errs.append(np.min(np.abs(cand.Epred - row.E_expt_eV)))
        errs_ks.append(np.min(np.abs(cand.Eks - row.E_expt_eV)))
    return float(np.sqrt(np.mean(np.array(errs)**2))), float(np.sqrt(np.mean(np.array(errs_ks)**2)))


if __name__ == '__main__':
    li_closed_form_check()
    print()
    for el in (sys.argv[1:] or ['Li', 'Na', 'Al']):
        cfg, shells = element_tables(el)
        print(f"== {el} (Z={cfg['Z_nuclear']}, zv={cfg['dft']['z_valence']}) core shells:")
        for sh in shells:
            tag = f"{sh['n']}{'spdf'[sh['l']]}"
            dem = DEMANDED.get((el, sh['n'], sh['l']))
            print(f"   {tag}: occ={sh['occ']:.0f} eps={sh['eps']:+.4f} Ha U_c={sh['U']:.4f} Ha"
                  + (f"   [gate-demanded vertex: {dem} Ha]" if dem else ""))
        df, eF = lam_states(el, shells, 'M')
        df.to_csv(f'{SCR}/lamderived_{el}.csv', index=False)
        g1 = df[df.ik == 1].sort_values('eps').iloc[0]
        depth_ks = (g1.eps - eF)*HA_EV
        depth = depth_ks/(1.0 + g1.lam)
        print(f"   Gamma band1: eps-eF={depth_ks:+.3f} eV lam={g1.lam:.5f} -> depth {depth:+.3f} eV"
              f" z={1/(1+g1.lam):.4f}")
        for sh in shells:
            tag = f"{sh['n']}{'spdf'[sh['l']]}"
            print(f"     E0_eff[{tag}] at Gamma = {g1[f'E0eff_{tag}']:.4f} Ha"
                  f" (lam_{tag}={g1[f'lam_{tag}']:.2e})")
        rmse, rmse_ks = gates(el, df, eF)
        if rmse is not None:
            print(f"   ARPES RMSE derived-vertex: {rmse:.3f} eV (bare KS {rmse_ks:.3f})")
        # occupied-state-averaged effective vertex per shell (weighted comparison)
        occd = df[df.eps < eF]
        for sh in shells:
            tag = f"{sh['n']}{'spdf'[sh['l']]}"
            print(f"     <E0_eff[{tag}]> over occupied states = {occd[f'E0eff_{tag}'].mean():.4f} Ha")
        print()
