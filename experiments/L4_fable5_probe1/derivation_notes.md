# Derivation notes (session 5, 2026-06-11)

## Enumeration of second-order contractions of the core-valence Coulomb term
Field split by orthogonal projection onto core orbitals (respects antisymmetry):
psi = psi_c + psi_v, P_c = sum_c |phi_c><phi_c| per site. Core subspace = the
OCCUPIED core orbitals only => core propagator G_c is pure hole propagator.
Vertices from (1/2) V psibar psibar psi psi, classified by core-field count:
 (i)   H_vvvv  : kept in valence EFT (LDA/eDMFT territory; NOT ours).
 (ii)  H_2c2v  : n_v V n_c (direct) and exchange-shaped B(r)Bbar(r').
        1st order => static core Hartree + static core-valence exchange => in V_PSP.
        2nd order needs one forward G_c => ZERO (core filled; no c->c' fluctuations
        within the occupied subspace).
 (iii) H_1c    : V psibar_v psibar_v psi_v psi_c + h.c.  (core electron promoted
        into the VALENCE-space continuum by valence density). ALL dynamical core
        polarization lives here under the projector split.
 (iv)  H_pair  : (1/2) V psibar_v psibar_v psi_c psi_c + h.c. (double core
        excitation). 2nd-order denominators ~ 2|eps_c| and tiny overlaps => negligible (ledger).
Second order in H_1c => 2p1h and 2h1p self-energy diagrams with intermediate
states |x k'; c^-1> where BOTH particles x,k' live in the valence space and the
pair is ANTISYMMETRIZED:
  dE1 = (1/2) sum_{c; x,k'} |<vc|V|k'x> - <vc|V|xk'>|^2 / (eps_v+eps_c-eps_x-eps_k')
  dE2 (2h1p, core-hole intermediate with large denominator |eps_c|) : negligible
      except the piece already in v1's f' kernel term.
|D|^2 part of dE1 + f'-term = EXACTLY v1's Fan kernel with the core fluctuation
propagator. NEW pieces: -2Re(D X*) + |X|^2 (exchange contractions), Pauli
blocking of x against the valence sea. Both estimated SMALL and partially
cancelling; to be quantified in the ledger, not the missing factor.

## Where the missing magnitude must live: the CONTRACTION, not the diagram
SETUP.md warns: correctly derived coupling appears wrong by a large factor if
Bloch-state contraction conventions are inconsistent. v1 contracts V_b with
PSEUDO (plane-wave) valence states. But the derived coupling V_b(r) is
concentrated in the core region (r <~ 1.5 a0), where the TRUE valence states are
the all-electron ones: orthogonal to the core, oscillating with local momentum
p_loc(r) = sqrt(2(eps - V_at(r))) -- the SAME local momentum as the core-excited
orbital u_x (eigenstates of the same potential). The radial coupling integrals
int u_eps'l'(r) w_bL(r) u_epsl(r) dr  with w_bL ~ u_x u_c r_</r_>^{L+1}
are therefore resonantly enhanced for distorted waves vs smooth Bessel j_l:
plane-wave contraction misses the oscillation matching.
=> v3: evaluate the SAME derived kernel entirely in the one-center all-electron
representation (renormalized-atom / distorted-wave): externals = AE radial
solutions at the band energy weighted by the Bloch state's local l-weights;
intermediates = the atomic KS spectrum (bound + box continuum), energies aligned
by eps' = eps_at + U0, U0 = EF - kF^2/2.
AUDIT IDENTITY (binding): with free externals AND free intermediates the
one-center evaluation must reproduce v1's plane-wave q-integral state by state.

## Closure rearrangement (numerical exactness of the counterterm)
K(e,e',D) -> 1/D as e'->inf, so Sum_k' |M|^2 K accumulates the closure sum
Sum_k'|<k'|V_b|ext>|^2/D_b = <ext|V_b^2|ext>/D_b over ALL intermediates.
Rearrange per channel:
  dE = C + Dyn,   C = sum_b <ext| (V_b)^2 |ext> / Delta_b   (EXACT, completeness;
      one-center radial: int [w_bL(r) u_ext(r)]^2 dr with angular reduction
      sum_lp (2lp+1) 3j(l,L,lp;000)^2 = 1)
  Dyn = sum_b sum_k' |M|^2 [ (1-f')/(e-e'-D) + f'/(e-e'+D) ]   (decays ~ |M|^2/e')
C is the (positive) removal of the adiabatic closure polarization; its state
dependence enters through the core-region amplitude of the external state.

## Angular reduction of the one-center matrix element (re-derivation, session 5)
ME = <u_lp,j' Y_lp,m' | V_b | u_l Y_l,mi>, V_b = sum_M (4pi/(2L+1)) w_L(r) g_b(M) Y_LM.
Sums over m_c,m_x (3j orthogonality, fixed M): sum |g_b(M)|^2 = (2lx+1)(2lc+1)(3j_xLc;000)^2/(4pi).
Gaunt for the valence side: sum over M,m' of (3j; -m' M mi)^2 = 1/(2l+1) [one m' per M].
|h|^2 carries (2lp+1)(2L+1)(2l+1)/(4pi) (3j_pLl;000)^2 => the (2l+1) CANCELS:
  sum_allm |ME|^2 = Rad^2 * W_b * (2lp+1) * 3j(l,L,lp;000)^2 / (2L+1)^2
=> pref_correct = W_b (2lp+1) tjp2 / (2L+1)^2.   onecenter.py had an extra
1/(2l_i+1): BUG (suppressed l>=1 external channels by 1/(2l+1)).
Cross-check (exact, kernel-free): closure sum over complete intermediates and
external plane wave must give S_direct = (n_at/Omega) sum_b (W_b/Delta_b)
(4pi/(2L+1)^2) int w_L(r)^2 r^2 dr; verified algebraically with sum_l (2l+1) j_l^2 = 1.

## Cross (exchange) contraction — explicit evaluation plan
Unrestricted (x,k') double sum already counts D^2 and X^2 (relabeling); the only
missing 2nd-order piece is
  dE_cross = -2 sum_{c,x,k'>F} D X / (eps_v + eps_c - eps_x - eps_k')  [same-spin core only => spin weight 1 vs 2 for D^2]
  D = <v c|V|k' x>  (multipole L between (v k') density at r and (c x) at r')
  X = <v c|V|x k'>  (multipole Lam between (v x) and (c k'))
Angular factor A(L,Lam; l_v l_c l_x l_k') = sum over all m of Gaunt products,
computed NUMERICALLY with a general Wigner-3j (Racah formula) — no 6j algebra
trusted by hand. Radial: R_L^D = int u_v u_k' w_L[(c,x)] dr, R_Lam^X = int u_v u_x w_Lam[(c,k')] dr.
Calibration: the D^2 term through the same code path must reproduce the
validated full-box value.

# Session 6 (Stage 1: Li)

## Target numbers (anchor)
LDA depth 3.48 eV, eDMFT 2.60 eV => narrowing 0.88 eV, z_G = 2.60/3.48 = 0.747,
lambda_G = 1/z - 1 = 0.338 (on-shell slope -dSigma/dw at band bottom ~ 0.25 if
linearized as dE = s*(eps-EF); self-consistent QP eq gives lambda = 0.338).

## 2nd-order contraction enumeration: 2h1p is exactly ZERO
Re-derived: with H_1c = sum <v1 v2|V|v3 c> a+_v1 a+_v2 a_c a_v3 + h.c., acting
on |gs + external particle>: the h.c. term requires a+_c on the FILLED core =>
0. The only 2nd-order self-energy topology is 2p1h with one core hole
(intermediates |x k'; c^-1>), both time orders included via the (1-f')/f'
kernel. Complete 2nd-order set: D^2, X^2 (both inside the unrestricted (x,k')
double sum), cross -2DX (same-spin core only), pair channel
<xx'|V|cc> (no external line => vacuum diagram, only Pauli-blocking residue,
bounded small), Pauli exclusion of core states from the k' sum.

## Li closed-form model (hydrogenic)
phi_1s = (Zs^3/pi)^{1/2} e^{-Zs r}, Zs = Z - 5/16 = 2.6875 (variational).
eps_1s(HF/Koopmans) = Zs^2/2 - Z Zs + (5/8)Zs = -2.7715 Ha  (HF exact -2.792).
U = F0(1s,1s) = (5/8) Zs = 1.6797 Ha.
E(1s^2) = -Zs^2 = -7.2227 Ha (exact -7.2799).
Excited x(nl): sees Z-1=2 => eps_x = -(Z-1)^2/(2n^2); Delta(1s->2p) =
-0.5 + 2.7715 = 2.2715 Ha = 61.8 eV (expt Li+ 1s2p ~61-62 eV). KS-LDA
eigenvalue difference instead gives ~1.83 Ha (49.8 eV) — KS underestimates Delta
=> overestimates the channel; ledger item.
alpha(Li+) ~ 0.192 a.u. (expt); dipole channel with alpha this small CANNOT
explain lambda ~ 0.34 by long-range polarization coupling — the question is the
SHORT-RANGE (in-core) part of the coupling with AE-contracted valence states.

## Immediate numeric scout (before closed forms): per-channel lambda_G(Li)
Use AtomChannels(3,1) + OneCenterSigma free-electron frame, l=0 AE local wave
at the band bottom; finite-difference -dSigma/dw. Channels for 1s core:
L = lx = 0 (1s->ns), 1 (1s->np), 2 (1s->nd).

## Step S6.1 — Closure/variance bound on the H_1c^2 class (THEOREM, Li)
For ANY contraction of the second-order core-hole self-energy (2p1h: D^2, X^2,
-2DX, with or without distorted-wave externals), the on-shell slope obeys
  lambda = -dSigma/dw |_{eps_G} <= sum_b <v| w_b^2 |v> / Delta_b^2   (no-recoil)
because every denominator is >= Delta_b (recoil eps_k'-eps_v >= 0 above EF) and
the k'-sum is bounded by completeness. Numerically on Li (AE-distorted l=0 wave
at the band bottom, WS-normalized, KS-LDA channels, Lmax=3):
  lambda_max(l=0) = 0.0112 ; lambda_max(l=1) = 0.0052 ;
  Sigma_ad(l=0) = 1.007 eV ; Sigma_ad(l=1) = 0.427 eV.
Anchor lambda(z_G=0.747) = 0.338. CONCLUSION (robust, factor >= 30): the omega-
dependence of the frozen-core 2p1h self-energy CANNOT produce z_G = 0.75. Any
contraction with core-scale denominators (2h1p crossed pieces have
|denominator| ~ |eps_c| as well) is bounded identically.

## Step S6.2 — Re-interpretation of the anchor (state- vs omega-dependence)
The anchor is a DEPTH ratio: E_G(eDMFT)/E_G(LDA) = 2.60/3.48 = 0.747. The
narrowing dE = [Sigma(k_F, EF) - Sigma(Gamma, eps_G)] = +0.88 eV mixes
STATE-dependence (k, channel weights, core-region amplitude) and omega-
dependence. S6.1 kills the pure-omega route only. v5 (Dyn-only, AE-contracted)
gave Na narrowing 0.28 of needed 0.46 (61%), NOT 8%: the per-channel lambda_G
scout compared the wrong quantity. Mandatory check now: full v5 on the Li band
-> predicted Gamma depth vs 2.60 eV (running).

## Step S6.3 — Li closed forms (hydrogenic/variational reference), exact rationals
Zs = Z - 5/16 = 43/16 = 2.6875 (variational 1s^2 screening).
u_1s(r) = 2 Zs^{3/2} r e^{-Zs r};  phi_1s = sqrt(Zs^3/pi) e^{-Zs r}.
eps_1s (Koopmans) = Zs^2/2 - Z Zs + (5/8) Zs = -1419/512 = -2.77148 Ha.
U = F0(1s,1s) = (5/8) Zs = 215/128 = 1.67969 Ha.
Excited core states: one electron in field Z2 = Z-1 = 2 (the other 1s screens):
u_2p(r) = (Z2^{5/2}/(2 sqrt6)) r^2 e^{-Z2 r/2};  eps_2p = -Z2^2/8 = -1/2 Ha.
Delta(1s->2p) = eps_2p - eps_1s = 1163/512 = 2.27148 Ha = 61.81 eV
   (expt Li+ 1s->2p 62.2 eV: 0.6% off — reference validated).
Dipole vertex (radial): d = <u_2p|r|u_1s> = 22544384 sqrt(129)/714924299
   = 0.358157 a.u.
alpha contribution of 2p alone = 4 d^2/(3 Delta) = 0.0753 a.u.; full np+continuum
sum gives alpha(Li+) ~ 0.19 a.u. (expt 0.192): polarization channels validated.

## Step S6.4 — Li dipole coupling vertex, closed form
w_1(r) = int u_2p(s) u_1s(s) r_<^1/r_>^2 ds  (potential of the 1s->2p transition density)
 = (688 sqrt(129)/714924299) * [32768(e^{ar} - 1) - 120832 r - 222784 r^2
   - 205379 r^3] e^{-a r} / r^2,  a = Zs + Z2/2 = 59/16.
Limits: w_1 -> d/r^2 (r->inf, d = 0.358157); w_1 -> (688 sqrt129/10443) r (r->0).
Values: w_1(0.5)=0.225, w_1(1)=0.199, w_1(1.5)=0.132, w_1(2)=0.085 Ha.
The vertex is O(0.2 Ha) over the core region but Delta = 2.27 Ha => per-channel
lambda ~ <w^2>/Delta^2 stays O(10^-2): consistent with the S6.1 bound.

## Step S6.5 — COMPLETE second-order contraction enumeration on Li (formal closure)
Vertex classes by core-field count (16 contractions of the quartic -> 5 classes):
 C0 vvvv (valence EFT; not part of the core integrate-out),
 C1 (1 core field) H_1c = V psibar_v psibar_v psi_v psi_c + h.c.,
 C2a n_c n_v (direct), C2b B-bar B (exchange-shaped), C2c pair cc->vv + h.c.,
 C3 (3 core fields), C4 core-core.
First order: C2a + C2b = static core Hartree + CV exchange => inside V_PSP
(reference). C4 = constant. C1, C2c, C3 vanish on the reference.
Second order on the added particle (Li 1s^2, closed shell):
 (1) C1xC1 2p1h (core hole c, valence pair x,k'): D^2 + X^2 (unrestricted
     (x,k') sum) and cross -2DX (same-spin core). Denominators
     -(Delta_b + recoil); BOTH time orders via the (1-f')/f' kernel.
 (2) C1xC1 2h1p orderings: need a+_c on the filled 1s^2 shell after normal
     ordering => only the f' kernel piece survives (already in (1)).
 (3) C2xC2 connected: requires a core PARTICLE propagator within the occupied
     projector subspace => ZERO identically.
 (4) odd core-field products (C1xC2, C3xC2, ...) => ZERO (unpaired core field).
 (5) C2c^2: double core excitation = vacuum fluctuation; enters Sigma_v only by
     Pauli blocking: |<v x|V|1s 1s>|^2 / (2|eps_1s|+...) — bounded below.
 (6) C3^2: needs an EMPTY core orbital => ZERO for the closed 1s^2 shell.
 (7) Pauli exclusion of core states from the k' sum (pauli=True term).
=> The COMPLETE dynamical set at O(V_cv^2) is {D^2, X^2, -2DX, pair-Pauli,
   core-state exclusion}: all carry core-scale denominators and obey the S6.1
   closure bound. CONSEQUENCE: the 25% Li narrowing CANNOT be an omega-slope
   (z) effect of any O(V^2) core contraction; it must (and does) come from the
   STATE-dependence of the SAME contractions — the resonant all-electron
   amplitude of the Bloch states in the core region, rising from the band
   bottom toward EF as eps approaches the atomic 2s resonance (Li: the AE 2s
   lies inside the occupied band; this is why Li narrows MOST despite
   alpha(Li+) = 0.19 being the SMALLEST core polarizability — penetration, not
   polarizability, sets the scale).

## Step S6.6 — The missing organizing principle: ADIABATIC integrate-out (Born-
## Oppenheimer in W/Delta), not bare 2nd-order recoil
Scale hierarchy on Li: occupied bandwidth W = 3.5 eV << Delta(1s->2p) = 61.8 eV.
The core fluctuations are FAST relative to the valence dispersion: integrating
the core out at FIXED valence configuration (adiabatic/BO limit of the core
path integral) is the correct leading order in W/Delta. For a valence electron
at position R the core's second-order response is
  V_pol(R) = - sum_b w_b(R)^2 / Delta_b      (exact closure form; for r->inf
  V_pol -> -alpha/2r^4, the CPP, since sum_b,dip w_b^2/Delta_b -> alpha/2r^4)
and the leading valence correction is the STATE-DEPENDENT expectation
  dE_pol(k) = <k| V_pol |k> = -C(k)   (AE-contracted Bloch state),
with the NONADIABATIC (retardation/recoil) correction R(k) = C(k)+Sigma_K(k)
>= 0 suppressed by W/Delta as a controlled ledger term. The strict O(V^2)
recoiled self-energy Sigma_K is NOT the leading answer in the W/Delta ordering:
the bare 2p1h kernel lets the intermediate valence electron recoil to
eps_k' ~ Zs^2/2 ~ 100 eV, which double-counts kinetic cost that the
(neglected at O(V^2)) ladder between the recoiled electron and the core hole
removes; the adiabatic limit resums exactly this physics (standard
strong-vs-weak-coupling polaron dichotomy, resolved here by W << Delta).
Bookkeeping: GTH/largecore is fit to the PLAIN LDA atom: contains NO V_pol
(LDA's local CV-correlation proxy ~ -0.05 eV on the atomic 2s, declared as
ledger error); so V_pol is omitted by the static reference and must be added
whole, with the ATOMIC-reference subtraction handled by the EF refill +
declared atomic-limit ledger entry.
Predicted sign/profile: C(k) grows from the band bottom toward EF (AE core-
region amplitude rises toward the 2s/3s resonance) => -C(k) lowers EF states
more than the bottom => NARROWING. Extracting C(k) on Li now via
dE_ct - dE_noct (fullbox aefree, pauli).
