# Derivation notes (micro-steps)

Atomic units (Ha, bohr) throughout.

## Step 0 — Setup and partition

All-electron action (given):
L = ∫ ψ̄ [∂_τ − ∇²/2 + V_Lat − μ] ψ + (1/2)∫ ψ̄ψ̄'ψ'ψ / |r−r'|.

Partition the one-particle Hilbert space per ion site R:
core subspace C_R = span{χ_{Ra}} = the N_c occupied orbitals of the free ion
(N_c = Z_nuclear − z_valence; Li: 1s↑1s↓), valence = orthogonal complement.
ψ(r) = Σ_{Ra} χ_{Ra}(r) d_{Ra} + ψ_v(r).

## Step 1 — Coulomb terms by leg count on the core

V splits by number of core legs: V_vvvv (valence-only; kept in valence theory,
treated by the pinned KS-LDA reference), V_cccc (core-only; defines ion),
V_vvcc (2+2), V_vccc, V_vvvc (3+1). Integrating out d at MEAN FIELD level
(frozen core ground state) gives the static Hartree + Fock exchange + Pauli
orthogonality of the frozen core acting on valence = what the norm-conserving
GTH pseudopotential V_PSP encodes. COUNTED — excluded from the correction.

The correction = effect of core FLUCTUATIONS: δV = V_cv − ⟨V_cv⟩_core-GS,
entering at second order (first order vanishes by construction).

## Step 2 — Second-order valence self-energy: complete contraction list

Reference: |0⟩ = |core GS (filled shells)⟩ ⊗ |KS valence Slater sea⟩.
Second-order self-energy of valence band state k (Goldstone, antisymmetrized
vertices ⟨pq||rs⟩ = ⟨pq|v|rs⟩ − ⟨pq|v|sr⟩). ALL second-order diagrams with at
least one core index, for a 1s² core (a, b ∈ {1s↑, 1s↓}):

(A) 2p1h, core hole a=1s:  Σ_A(k,ω) = Σ_{a} Σ_{p,k' unocc} |⟨k'p||k a⟩|² / (ω + ε_a − ε_p − ε_{k'})
    [p = excited orbital reached from 1s, k' = valence intermediate.
     Process: valence k scatters k→k', core excited 1s→p; OR exchange pairing.]
(B) 2h1p, one core hole + one valence hole: Σ_B = Σ_a Σ_{k₂ occ, p unocc} |⟨a k₂||k p⟩|² / (ω − ε_a − ε_{k₂} + ε_p)
(C) 2h1p, double core hole: Σ_C = Σ_{p unocc} |⟨1s↑1s↓||k p⟩|² / (ω − 2ε_1s + ε_p)
There are no other second-order contractions touching the core: 2p1h needs ≥1
hole (valence-only hole ⇒ pure valence diagram = valence correlation, assigned
to the KS-LDA reference), and a 1s² core admits exactly one or two core holes.
Each |⟨..||..⟩|² expands into direct², direct×exchange (×2), exchange² — the
complete second-order contraction set on the two-electron core is
{A-DD, A-DX, A-XX, B-DD, B-DX, B-XX, C} per core spin channel.

## Step 3 — Energy denominators and the fast-core limit

ε_1s ≈ −B (core binding, B ≈ 2.7 Ha for Li); valence energies |ε| ≲ 0.3 Ha.
A: ω + ε_1s − ε_p − ε_k' = ω − ε_k' − Δ_p, with Δ_p ≡ ε_p − ε_1s ≥ Δ_min ≈ 2 Ha.
B: ω − ε_1s − ε_k₂ + ε_p ≈ ω − ε_k₂ + Δ_p' (Δ_p' = ε_p − ε_1s again).
C: ω − 2ε_1s + ε_p ≈ ω + 2B + ε_p  (≈ 5.5 Ha for Li).

The static parts Σ(k, ω_ref) are representable by a static potential — counted
in V_PSP (the GTH fit reproduces the atomic valence levels, i.e. the static
core dressing AT the atomic reference energy). The part static pseudopotentials
OMIT is the ω-dependence. Linearize around the reference:
Σ(k,ω) ≈ Σ(k,ω_ref) + (ω−ω_ref) ∂Σ/∂ω,  λ_k ≡ −∂Σ/∂ω|_{ε_k} > 0
(every 2nd-order term contributes −∂/∂ω |M|²/(ω−X) = +|M|²/(ω−X)² to λ).

QP equation: E = ε_k^KS + Σ(E) − Σ(ω_ref) ⇒ E = z_k ε_k + (1−z_k) ω_ref,
z_k = 1/(1+λ_k). Same for the Fermi level ⇒ measured from E_F the reference
ω_ref drops out:
    E_QP − E_F = z_k (ε_k^KS − E_F),   z_k = 1/(1 + λ_k).
(For k-dependent λ: E_F^QP anchored at the Fermi surface, where ε=ε_F, so
E_QP−E_F^QP = z_k(ε_k−ε_F) + (z_k − z_{k_F})(ε_F − ω_ref); second term dropped —
ledger item; vanishes at the FS and is O(λ-variation × |ε_F−ω_ref|).)

## Step 4 — λ_k: closure over valence intermediates

In λ_A, denominators (ε_k − ε_{k'} − Δ_p)² ≈ Δ_p² [error O(2(ε_k−ε_k')/Δ_p),
controlled, ledger]. In λ_B same with +Δ_p. A sums k' over unocc, B over occ;
to leading order in 1/Δ they COMBINE into an unrestricted closure
Σ_{k' all} |k'⟩⟨k'| = 1 − P_core (per spin), with weight 1/Δ_p².

Vertex unpacking, fixed excited orbital p, core orbital a=1s, valence k:
 D(k')  = ⟨k'|V_p|k⟩,  V_p(r1) = ∫ p*(2) 1s(2) /|r1−r2| d2   (transition potential)
 X(k')  = ⟨k'| x_p k⟩,  [x_p k](r1) = 1s(r1) ∫ p*(2) k(2)/|r1−r2| d2
Spin sum over a (k fixed spin): |D−X|² + |D|² = 2|D|² − 2Re(D*X) + |X|².
Closure over k':
 λ_DD(k) = 2 Σ_p (1/Δ_p²) ⟨V_p k| (1−P_1s) |V_p k⟩
 λ_DX(k) = −2 Σ_p (1/Δ_p²) Re ⟨V_p k| (1−P_1s) |x_p k⟩
 λ_XX(k) =  Σ_p (1/Δ_p²) ⟨x_p k| (1−P_1s) |x_p k⟩
 λ_C(k)  =  Σ_p |⟨p ⊗ 1s | v | k ⊗ 1s ... ⟩|²/(2B+ε_p−...)² (double-core-hole, computed, small)
λ_k = Σ_sites [λ_DD + λ_DX + λ_XX + λ_C]   (multi-shell cores: also Σ over core shells a).

## Step 5 — p-sums via spectral kernel (per angular channel)

p runs over the complete unoccupied spectrum reached from the core orbital.
Kernel: G2_l(r,r') = Σ_{p∈l} R_p(r) R_p(r') / Δ_p², computed by exact
diagonalization of the radial Hamiltonian h_l in a large box (continuum
discretized; 1/Δ²-weighted sums converge).

Choice of h_l for excited orbitals (DECLARED): the excited electron moves in
the potential of the nucleus + the REMAINING core (core with one hole in the
shell being excited) — i.e. the physical (excitonic) ionic excitation:
Δ_p = B_a + ε_p[V_hole], B_a = binding energy of core orbital a (ionic HF /
exact two-electron closed form for Li). This resums the core-hole ladder; for
the He-like Li core it reproduces the exact Li⁺ spectrum (excited electron in
the field of nucleus + single 1s = an exact one-electron problem up to
exchange with the lone 1s electron, ledger).

## Step 6 — multipole decomposition (DD term)

1/|r1−r2| = Σ_L (4π/(2L+1)) (r_<^L/r_>^{L+1}) Σ_M Y_LM(r̂1)Y*_LM(r̂2).
1s → p=(n_p, l_p=L, M): V_p(r1) = (4π/(2L+1)) Y_LM(r̂1) (1/√4π) w_L^{(p)}(r1),
 w_L^{(p)}(r) = ∫ R_p(r') R_1s(r') (r_<^L/r_>^{L+1}) r'² dr'.
Σ_M |Y_LM|² = (2L+1)/4π ⇒ Σ_{p∈L,M} |V_p(r)|²/Δ_p² = Σ_{n_p} w_L²(r)/[(2L+1) Δ²]
F(r) ≡ Σ_L F_L(r),  F_L(r) = (1/(2L+1)) ⟨1s| u_L(r,·) G2_L u_L(r,·) |1s⟩,
 u_L(r,r') = r_<^L / r_>^{L+1}.
λ_DD(k) = 2 ∫ |ψ_k^{AE}(r)|² F(r) d³r  −  (P_1s correction, computed; small)
with ψ^{AE} the ALL-ELECTRON Bloch amplitude on the site (see Step 8).

## Step 7 — Li closed forms (Stage 1)

Core orbital: R_1s(r) = 2 Z_s^{3/2} e^{−Z_s r}, Z_s = Z − 5/16 = 43/16 = 2.6875
(variational optimum for the He-like ion; exact closed form).
Core total energy: E(1s²) = Z_s² − 2 Z Z_s + (5/8) Z_s = −7.2227 Ha.
Self-Coulomb integral: J_ss = ⟨1s1s|v|1s1s⟩ = (5/8) Z_s = 1.6797 Ha.
Core orbital removal energy: B = E(Li²⁺) − E(Li⁺) = −Z²/2 + 7.2227 = 2.7227 Ha
 (exact value 2.7798; ledger ±2% on B ⇒ ∓4% on λ — wait, enters via Δ², see ledger).
Hole potential for excited orbitals: V_h(r) = −Z/r + (1/r)[1 − e^{−2Z_s r}(1+Z_s r)]
 (nucleus + Hartree of the single remaining 1s electron; → −(Z−1)/r = −2/r).
Excitation energies: Δ_{nl} = B + ε_{nl}[V_h]; leading closed forms
 ε_{2p} ≈ −(Z−1)²/8 = −0.5 ⇒ Δ_{1s→2p} ≈ 2.22 Ha (exp. Li⁺ 1s2p ≈ 2.28);
 continuum Δ = B + k²/2.
Coupling vertex (dipole channel example, closed form):
 w_1^{(2p)}(r) and ⟨1s|r|2p⟩ analytic for hydrogenic orbitals.
z_Γ = 1/(1+λ_Γ) with λ_Γ from the complete set — numerical target ≈ 0.75.

## Step 8 — Bloch contraction: all-electron on-site amplitude

The vertex weight F(r) lives in/near the core; pseudo Bloch states lack the
AE core-region amplitude (orthogonality lobes + correct nucleus amplitude).
Consistent contraction: around each site expand the pseudo state exactly
(Rayleigh): ψ̃_kn(r+R_site) = Σ_lm A_lm(r) Y_lm(r̂),
 A_lm(r) = 4π i^l e^{i(k)·R}... Σ_G c_nk(G) e^{i(k+G)·R_site} j_l(|k+G| r) Y*_lm(k+G).
For r ≥ r_match (outside the core) pseudo = AE. Replace, inside the site
sphere, the radial profile of channel l by the AE frozen-core valence solution
φ_l^{AE}(r;ε̄) matched to A_lm on a shell at r_match:
 α_lm = ∫_shell A_lm(r) φ_l^{AE}(r) r²dr / ∫_shell [φ_l^{AE}]² r² dr,
 ψ^{AE}-on-site channel radial function = α_lm φ_l^{AE}(r).
Channel weights w_l(k,n) = Σ_m |α_lm|². Then (F spherical ⇒ lm-diagonal):
 λ_DD(k,n) = 2 Σ_l w_l(k,n) ∫_0^{r_sph} [φ_l^{AE}(r)]² F(r) r² dr + outside-sphere part
 (outside r_sph the state is the raw pseudo state; F there is the small
  multipole tail — evaluated with the plane-wave state, ledger).

## Step 9 — exchange terms angular reduction

[to be derived next: λ_DX, λ_XX reduced to radial integrals with 3j factors;
for the s-valence channel exact; for l_v ≥ 1 exact 3j algebra retained]
