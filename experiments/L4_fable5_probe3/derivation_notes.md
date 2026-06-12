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

> [MAINTAINER COMMENT 2026-06-12 — ERROR SITE]: this dismissal, together with the
> later E0-closure substitution in the hybridization channel, is where the
> derivation went wrong. See scratch/MAINTAINER_NOTE.md before continuing.

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

## 5. GUIDED SESSION (2026-06-12) — deriving the hybridization vertex from the
##    two-body-induced transfer-amplitude class

STATUS DISCLOSURE: the identification of the correct amplitude CLASS (the
interaction-induced, spectator-mediated transfer G(p) = <chi p|v|chi chi>, with
orbital-average/zero-diagonal subtraction) was supplied by the maintainer note
(scratch/MAINTAINER_NOTE.md). Everything below — the collapse lemma made precise,
the algebra, the closed forms, the multi-shell generalization, the numerical
evaluation, and the agreement/disagreement verdict — is worked out here without
further guidance and without calibration.

### 5.1 No one-body vertex exists (the collapse lemma, made precise)

Let Phi be the reference determinant and Phi_f = a+_p a_{chi sigma} Phi the
core->valence transfer. Slater-Condon for the FULL Hamiltonian:
  <Phi_f|H|Phi> = <p|h1|chi> + sum_{q occ, q != chi sigma} [<pq|v|chi q> - delta_{s} <pq|v|q chi>]
The excluded self term q = chi sigma contributes <p chi|v|chi chi> - <p chi|v|chi chi> = 0,
so it can be re-included for free, and the sum closes to the full Fock operator:
  <Phi_f|H|Phi> = <p|F|chi> = eps_c <p|chi>      (Brillouin-type collapse).
If core and valence are eigenstates of a COMMON mean field, the transfer
amplitude is pure overlap x eps_c — exactly the static piece the reference
already counts. Hence: NO one-body operator can be the hybridization vertex.
Session-2/3 symptoms (cusp-dominated closure energy from -Z/r pieces, ~3x
reference sensitivity, Hermitian collapse) were fingerprints of searching in
that class.

### 5.2 What survives: the reference mismatch is the same-shell pair term

Our reference is NOT a common mean field: eps_c, chi_c come from the ATOMIC
LDA problem; the valence Bloch states from the PSEUDIZED crystal problem in
which the core appears only as a static density. The one term treated
differently in the two references is the interaction of the transferring core
electron with its OWN SHELL PARTNER(S):
  - in eps_c it enters as the static orbital average U_c = <chi chi|v|chi chi>
    (the atomic LDA eigenvalue contains the shell-partner Hartree);
  - in the crystal reference the same interaction sits inside the static
    frozen-core Hartree of V_PSP, counted once at the reference energy.
Neither reference contains the DYNAMICS of the pair: the amplitude for the
shell pair to break, |chi chi> -> |chi p>, the spectator staying behind. That
amplitude is generated by the W1/W3 quartic algebra (sec. 0) and is exactly the
channel-C-class object already written in sec. 1:
  G_c(p) = <chi p| v |chi chi> = <p| V_H[|chi_c|^2] |chi_c>.
Transfer mediated by the spectator's orbital Hartree potential. Properties
(contrast with the one-body class): V_H[chi^2](0) = finite (no nuclear cusp;
for Li V_H(0) = zeta), so G is cusp-free; it involves only shell-c objects, so
it is per-shell; and it contains no bare-Z piece, so it is weakly
reference-sensitive.

### 5.3 Zero-diagonal subtraction (q->0 / static bookkeeping)

The orbital-average (static) component of this transfer operator is its
diagonal U_c — the piece already counted once in eps_c (ledger item 3 pointer:
"an excitation-energy definition"). Subtract it by orthogonalizing the final
state to the core (the OPW discipline of sec. 1, T -> O(q)):
  M_c(p) = <OPW_p| V_H[rho_spec,c] |chi_c> = G_c(p) - U_c * chit_c(p),
  <chi_c| M_hat |chi_c> = 0  by construction (zero-diagonal rule).
The long-range monopole of V_H[rho_spec] (-> N_spec/r) is the slowly-varying
piece whose effect is static screening, already inside the PSP core Hartree;
the subtraction removes precisely its orbital-average shadow.

### 5.4 Spin / spectator counting; multi-shell generalization

The band electron hybridizes with each spatial core orbital m: (2l+1) channels
per shell — identical to the locked Q structure. For the transfer of one
electron out of shell c with N_c electrons, the spectator density is
  rho_spec,c(r) = (N_c - 1) * u_c(r)^2 / (4 pi r^2)   (spherical average),
i.e. exactly |chi|^2 (ONE spectator, opposite spin, no exchange term) for the
s^2 shells that dominate (Li-1s, Na-2s, Al-2s). For p^6 shells the direct
mediator is 5 spectators; same-spin exchange reduces the effective mediation
(2 of 5 spectators) — evaluated direct-only with the reduction bounded in the
ledger (p-shell contributions are subdominant; checked numerically).
Pair-state normalization: the singlet pair |chi chi> couples to the
(normalized) broken-pair singlet (1/sqrt2)(chi p + p chi) with amplitude
sqrt(2) G_orth; but the band electron enters ONE spin channel, and per spin
channel the Slater-Condon amplitude is G_orth(p) exactly once. No extra
factor. (Recorded; the sqrt2 variant is reported as a sensitivity, not used.)

### 5.5 Li closed form (Stage-1)

chi = sqrt(zeta^3/pi) e^{-zeta r}, zeta = 2.6875;
V_H[chi^2](r) = (1/r)[1 - e^{-2 zeta r}(1 + zeta r)];  V_H(0) = zeta.
chi V_H = sqrt(zeta^3/pi) [ e^{-zeta r}/r - e^{-3 zeta r}/r - zeta e^{-3 zeta r} ],
so the transform is elementary:
  G(p) = 4 pi sqrt(zeta^3/pi) [ 1/(p^2+zeta^2) - 1/(p^2+9 zeta^2) - 6 zeta^2/(p^2+9 zeta^2)^2 ].
  G(0) = (264/81) pi sqrt(zeta^3/pi) / zeta^2;  G(0)/chit(0) = (33/81) zeta = 1.095 Ha.
  U = (5/8) zeta = 1.6797 Ha.
  M(p) = G(p) - U chit(p);  effective scale  E_eff(p) = M(p)/chit(p):
    E_eff(0) = (33/81 - 5/8) zeta = -0.2176 zeta = -0.585 Ha,
    E_eff(p->inf) -> (1 - 5/8) zeta = 0.375 zeta = +1.008 Ha
  (G and chit share the p^-4 cusp...-free tail ratio zeta; subtraction leaves 3/8 zeta).
The derived per-shell vertex is therefore BOUNDED by ~|0.6|–|1.0| Ha on Li over
all p — to be compared against the gate-demanded 3.35 Ha (sec. 5.7).

### 5.6 Production form (drop-in replacement of E0)

  lambda_nk = sum_c (2l_c+1)/(4 pi Omega) * n_atoms * sum_G |c_G|^2 *
              M_c(|k+G|)^2 / (eps_nk - eps_c)^2,
identical to the locked model with  E0 * cf_c(kappa)  ->  M_c(kappa),
  M_c(kappa) = 4 pi int j_l(kappa r) [V_spec,c(r) - U_c] u_c(r) r dr,
  V_spec,c = Poisson[rho_spec,c],  U_c = int V_spec,c u_c^2 dr,
all from the same radial LDA solver that defines eps_c (reference-consistent).
Per-shell effective vertex for comparison with the gate-demanded numbers:
  E0_eff,c(nk) = sqrt[ sum_G |c_G|^2 M_c(|k+G|)^2 / sum_G |c_G|^2 cf_c(|k+G|)^2 ].

### 5.7 Numerical results (this session) — VERDICT: HONEST DISAGREEMENT

Convention audit (SETUP warning heeded): pushing the LOCKED vertex pi*cf_c
through this session's code path reproduces the production numbers exactly
(Li lam_Gamma = 0.2977; Na RMSE 0.094 eV; Al RMSE 0.219 eV). The comparison
below is therefore apples-to-apples; no contraction-convention factor hides.

Derived per-shell effective vertex E0_eff (Bloch-weighted, gate states) vs the
gate-demanded magnitudes:

| shell  | U_c (Ha) | derived E0_eff (M, subtracted) | unsubtracted G | demanded |
|--------|----------|-------------------------------|----------------|----------|
| Li 1s  | 1.613    | 0.620                         | 0.994          | 3.35     |
| Na 2s  | 1.151    | 0.314                         | 0.837          | 3.18     |
| Al 2s  | 1.440    | 0.371                         | 1.069          | 2.65     |

Gate metrics with the derived vertex (production formula, M_c drop-in):
  Li: lam_Gamma = 0.0116 -> z_Gamma = 0.989 vs anchor 0.747 (needed lam 0.339);
  Na: ARPES RMSE 0.355 eV (locked 0.094, baseline 0.078, bare KS 0.413);
  Al: ARPES RMSE 0.322 eV (locked 0.219, baseline 0.248, bare KS 0.414).
The derived channel moves every gate in the RIGHT DIRECTION but is a factor
5.4/10.1/7.1 too small in |M| per shell (30-300x in lambda).

Sensitivity variants (reported, not used — none rescues the gates):
  - sqrt(2) pair normalization: x1.41 on |M| -> Li E0_eff 0.88. Still 3.8x short.
  - full-shell mediator (N_c instead of N_c-1; both G and U double): exactly 2M
    -> Li 1.24, lam_Gamma 0.046, z 0.956. Still fails.
  - both combined (2*sqrt2): Li 1.75 / Na 0.89 / Al 1.05 — short AND ordering wrong.
  - cross-shell mediator (V_H[all core minus self], zero-diagonal subtracted):
    Li UNCHANGED 0.620 (no other shells exist); Na 2s 3.60, Al 2s 4.40
    -> Na near demanded but Al overshoots 66% (RMSE 0.257/0.654 — Al worse than
    bare KS). Not consistent across elements, and structurally unable to fix Li.

Analysis of the disagreement:
1. Li ceiling (closed form): within this class the vertex scale is bounded by
   the spectator orbital-Hartree potential, max V_H(0) = zeta = 2.69 Ha; the
   Bloch-weighted ratio G/chit at the momenta the gates probe (|k+G| <~ 3/bohr)
   is 0.41–0.6 zeta, and the zero-diagonal subtraction removes 0.625 zeta of it.
   The demanded 3.35 Ha = 1.25 zeta EXCEEDS even the unsubtracted mediator
   maximum. On Li — one spectator electron, no other shells — no bookkeeping
   within "transfer mediated by the spectator core electron's Hartree
   potential" can reach the demanded magnitude. This is a closed-form
   obstruction, not a numerical one.
2. The demanded set tracks no atomic same-shell scalar: 2*U_c(LDA) =
   3.23/2.30/2.88 Ha (Li/Na/Al) vs demanded 3.35/3.18/2.65 — Li coincides to
   4% (numerically 2U_1s^hyd = (5/4)zeta = 3.359), Na/Al do not, and the
   demanded Na > Al ordering INVERTS the atomic tightness ordering of the same
   2s shell type. The demanded vertex therefore carries environment (valence/
   lattice) dependence that no purely atomic same-shell object has. Notably
   r_s(Na)=3.93 > r_s(Al)=2.07: the ordering matches valence-screening
   strength, hinting the missing weight is the DYNAMICAL part of the
   q->0/valence-mediated piece of the same W1 class (the channel-B direct term
   v_q <chi|e^{iqr}|k> that was parked, whose static orbital-average we
   subtract; its plasmon-pole dynamics at the core pole is NOT static and was
   dropped with it). Untested here — would require the crystal RPA response;
   recorded as the leading hypothesis, not pursued (no time-boxed fit allowed
   or wanted).
3. Ratios demanded/derived = 5.4/10.1/7.1 are non-universal -> no missing
   constant prefactor (spin, pair normalization, mediator count) can close the
   gap.

CONSEQUENCE (per maintainer task item 3): the production formula is NOT
updated; E0 = pi remains CALIBRATED in solution/method.md's ledger. The
derived spectator-mediated vertex is the correct OBJECT CLASS by every
qualitative criterion (cusp-free: it is finite everywhere; weakly
reference-sensitive: hydrogenic vs LDA Li orbital moves E0_eff by 5.8%
[0.584 vs 0.620 Ha; computed], vs the ~3x of the one-body candidates;
per-shell resolved; right direction on all three gates) but its
same-shell-mediated magnitude is 5-10x too small to be the whole vertex.

## 4. Open numerical questions (to be settled by Stage-1 numbers)

- magnitude of channels A vs C vs B on Li; target: Gamma-point depth ratio
  z_Gamma = E_QP/eps_KS depth = 2.60/3.48 = 0.747 (anchor).
- valence screening of the small-q Coulomb in channel A: Thomas-Fermi
  eps(q) = 1 + k_TF^2/q^2 (parameter-free; k_TF^2 = 4 k_F/pi). Derivable as the
  RPA static limit of the K1 sector dressing the K4 vertex; decide by control.
- whether the on-shell k-dependence (DeltaSigma) or the omega-slope (z) carries
  the narrowing. Both computed; no selection by anchor-fitting — the full
  QP equation of sec. 3 is used as derived.
