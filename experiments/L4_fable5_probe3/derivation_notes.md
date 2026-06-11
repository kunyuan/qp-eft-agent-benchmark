# Derivation notes (live) — L4 frozen-core dynamical correction

Atomic units (Hartree, bohr) throughout unless eV stated.

## 0. Strategy

Start from the all-electron action (packet SETUP). Choose a mean-field reference
determinant: core shells {chi_a} (localized, per site, fully occupied:
N_core = Z_nuclear - z_valence electrons) + valence Bloch sea. Split the single
antisymmetrized field psi = psi_c + psi_v by the projector P_c = sum_a |chi_a><chi_a|
(per site), psi_v = (1-P_c) psi. Wick-contract the ONE electron field; classify
the 2^4 = 16 leg assignments of the quartic Coulomb vertex
(1/2) int psibar(r) psibar(r') v(r-r') psi(r') psi(r) by Delta N_c = (#c) - (#cbar):

| class | legs | Delta N_c | fate |
|---|---|---|---|
| K0 | cccc (cbar cbar c c) | 0 | core-core Coulomb: in the atomic core reference energy / PSP construction. Static. Dropped (counted once there). |
| K1 | vvvv | 0 | valence-valence: handled by the pinned KS step (Hartree + LDA xc). Dropped (that is the KS baseline, not the core correction). |
| K2 | cbar c vbar v, direct (densities at separate vertices) | 0 | static core Hartree V_H[n_c] on the valence: contained in the local part of V_PSP. The FLUCTUATION part requires an occupied->occupied core transition, which is Pauli-forbidden in a full frozen core => no dynamical content. Dropped (counted in V_PSP). |
| K3 | cbar v vbar c, exchange assignment | 0 | static nonlocal core-valence exchange K_c: contained in the PSP at the atomic-LDA level (GTH is fitted to all-electron atomic valence levels). Static. Dropped (counted, with ledger error = LDA-vs-Fock core exchange difference). |
| K4 | vbar vbar v c + h.c. | -1/+1 | ONE-CORE-EXCITATION vertex W1: a valence density scatters a core electron into the non-core space. -> 2nd-order channels A and B below. KEPT. |
| K5 | vbar cbar c c + h.c. | -1/+1 | three-core vertex W3: core pair interaction ejecting one core electron. At 2nd order it dresses the CORE propagator only; enters the valence self-energy first at 3rd order ~ (W1 W3 W1). Dropped with error estimate O(W3^2/Delta_core) relative shift of I_c, absorbed into the Delta-SCF excitation energies (which already include core relaxation). |
| K6 | vbar vbar c c + h.c. | -2/+2 | CORE-PAIR vertex W2 (does NOT conserve core particle number): both electrons of a core pair scatter into the valence space. -> 2nd-order channel C. KEPT. This is the assignment an enumeration that conserves core number separately would miss. |

Second-order valence self-energy = all MP2-like terms with at least one core
label. Complete list (i,j occupied labels; a,b unoccupied; external p = valence k):

- 2p1h, i = core:        channel A (core polarization). From <W1 W1^dag>.
- 2h1p, (i,j) = (core, valence-occ): channel B. From <W1^dag W1>.
- 2h1p, (i,j) = (core, core):        channel C (pair channel). From <W2^dag W2>.
- 2p1h/2h1p with no core label: valence correlation = LDA's job (K1). Excluded.
- Cross terms <W1 W2> etc.: intermediate states have different core-hole number
  => orthogonal => vanish at 2nd order. (Noted in ledger.)
- First-order terms (K2,K3 diagonals): static, in V_PSP; this is the ledger-item-3
  pointer: the diagonal <chi|e^{iq.r}|chi> = F_c(q) Coulomb-dressed IS the core
  Hartree form factor in V_PSP's local part, and the orthogonalized vertex below
  subtracts exactly that diagonal.

## 1. Self-energy expressions (per atom; external = plane wave k, spin 1/2)

Notation: v_q = 4pi/q^2; n_at = atoms/volume; chi-tilde(kappa) = FT of core
orbital; F(q) = <chi|e^{iq.r}|chi>. Intermediate excited electrons = plane waves
ORTHOGONALIZED to the core (OPW) — orthogonalization is what removes the static
(diagonal) core component from the vertex and makes the q->0 behavior dipole-like:

  T(p,q) = <OPW_p| e^{iq.r} |chi> = chi-tilde(|p-q|) - chi-tilde(p) F(q),  T -> O(q) as q->0.

Channel A (2p1h, core hole; spin-summed, direct + exchange):
  Sigma_A(k,w) = n_at Int d3q/(2pi)^3 d3p/(2pi)^3  Theta(|k-q|-k_F) Theta(p-k_F)
                 [ v_q T(p,q) ( 2 v_q T(p,q) - v_{|k-p|} T(k-q, k-p) ) ]
                 / ( w - I_c - p^2/2 - eps_{k-q} + eps-floor-consistent )
with denominator D_A = w_kin - ( -I_c ... ) written concretely below; on shell
(w = k^2/2, all kinetic energies from the same free-electron zero):
  D_A = k^2/2 - (k-q)^2/2 - p^2/2 - I_c   (always < 0; |D_A| >= I_c).
I_c = Delta-SCF core ionization energy (explicit expression, sec. 2).

Channel C (2h1p, core singlet pair removed; one spatial core orbital chi):
  Sigma_C(k,w) = n_at Int d3p/(2pi)^3 Theta(p-k_F) |A_orth(k,p)|^2 / D_C,
  D_C = k^2/2 + p^2/2 + I_pair  ( > 0 ),
  A(k,p) = <k p| v |chi chi> = (1/(2pi)^3) Int d3u v_u chi-tilde(|k-u|) chi-tilde(|p+u|),
  A_orth = A - chit(k) G(p) - chit(p) G(k) + chit(k) chit(p) U,
  G(p) = <chi p|v|chi chi> = Int e^{-ip.r} chi(r) V_H[chi^2](r) d3r,  U = <chi chi|v|chi chi>.
Spin algebra for the doubly occupied spatial orbital gives weight 1*|A_orth|^2
(the same-spin combination Pauli-cancels). I_pair = Delta-SCF double ionization.

Channel B (2h1p, one core + one valence hole): suppressed by the occupied-valence
phase space (k_F^3 Omega_cell /(2pi)^3 ~ z_val); evaluated explicitly on Li to
bound it (sec. 4); expression analogous with M = <k p|v|chi j>, |j|<k_F.

## 2. Li closed forms (Stage 1)

Core = He-like 1s^2, variational hydrogenic orbital:
  chi(r) = (zeta^3/pi)^(1/2) e^(-zeta r),  zeta = Z - 5/16 = 2.6875 (Z=3).
  chi-tilde(kappa) = 8 sqrt(pi zeta^5) / (kappa^2 + zeta^2)^2.
  F(q) = [1 + (q/(2 zeta))^2]^(-2).
  U = <chi chi|v|chi chi> = (5/8) zeta = 1.6797 Ha.
  E(1s^2) = zeta^2 - 2 Z zeta + (5/8) zeta = -zeta^2 = -7.22266 Ha.
  E(1s^1) = -Z^2/2 = -4.5 Ha (exact hydrogenic).
  Koopmans orbital energy: eps_1s = zeta^2/2 - Z zeta + (5/8)zeta = -2.77149 Ha.
  Delta-SCF: I_1 = E(1s^1) - E(1s^2) = zeta^2 - Z^2/2 = 2.72266 Ha  (USED in D_A).
             I_pair = -E(1s^2) = zeta^2 = 7.22266 Ha                (USED in D_C).
  V_H[chi^2](r) = (1/r)[1 - e^{-2 zeta r}(1 + zeta r)]  (closed form),
  G(p) = 4pi sqrt(zeta^3/pi) * [ 1/(zeta^2+p^2)... ] (radial integral, closed form,
         evaluated by 1D quadrature in code; pure 1s objects).

Excitation-energy explicit expressions (ledger item 2):
  Delta E_A(p; k,q) = I_1 + p^2/2 + ((k-q)^2 - k^2)/2   (core hole + OPW electron p
      + scattered valence; free intermediate electron, core-hole attraction in the
      intermediate state neglected -> D overestimated -> correction UNDERestimated;
      bound: exciton binding <~ Z_eff^2/8 relative to I_1, declared in ledger).
  Delta E_C(p; k) = I_pair + p^2/2 + k^2/2.

## 3. What the PSP already counts (double-counting control)

V_PSP (GTH, largecore, LDA) contains: core Hartree (K2), core-valence
exchange-correlation at atomic LDA level (K3 + the static shadow of A/B/C),
and atomic scattering at the reference valence energy (norm conservation =>
correct LINEAR energy dependence of the static scattering). Therefore the
correction may keep only:
  (i) the omega-dependence of Sigma_cv beyond linear-static => QP weight z;
  (ii) the k-dependence of Sigma_cv ACROSS THE OCCUPIED BAND relative to the
       Fermi surface (a local static potential — LDA xc, local PSP — is exactly
       k-independent for plane waves; the nonlocal-projector k-dependence is
       atomic-reference static scattering, distinct from the 2nd-order
       correlation object computed here).
Operational QP equation (anchored at the Fermi surface, Luttinger-consistent):
  E_QP(nk) - E_F = z_nk [ (eps_nk - E_F) + DeltaSigma_nk ],
  DeltaSigma_nk = Sigma_cv(k, E_F) - Sigma_cv(k_F, E_F),
  z_nk = 1 / (1 + lambda_nk),  lambda_nk = -dSigma_cv/dw |_{E_F}.
Any constant static piece (and the PSP's content) cancels in DeltaSigma; at
k = k_F the correction vanishes identically (E_F preserved).

## 4. Open numerical questions (to be settled by Stage-1 numbers)

- magnitude of channels A vs C vs B on Li; target: Gamma-point depth ratio
  z_Gamma = E_QP/eps_KS depth = 2.60/3.48 = 0.747 (anchor).
- valence screening of the small-q Coulomb in channel A: Thomas-Fermi
  eps(q) = 1 + k_TF^2/q^2 (parameter-free; k_TF^2 = 4 k_F/pi). Derivable as the
  RPA static limit of the K1 sector dressing the K4 vertex; decide by control.
- whether the on-shell k-dependence (DeltaSigma) or the omega-slope (z) carries
  the narrowing. Both computed; no selection by anchor-fitting — the full
  QP equation of sec. 3 is used as derived.
