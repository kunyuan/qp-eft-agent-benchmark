# First-principles z for Li with both 1s and 2s integrated out

**Ordering (the paradigm-correct one):** isolated-atom reference → exact
dual-fermion integration of all atomic configurations → effective action for
the inhomogeneous electron gas → KS from the paper's EFT (LDA-class
resummation of the *residual, inter-site* interaction) → fluctuation expansion
of the frequency-dependent self-energy with KS propagators.
In this ordering every diagram has a unique home (on-site vertices belong to
the atomic reference; diagrams with ≥1 inter-site interaction line belong to
the residual gas; static parts belong to KS; dynamic fluctuation parts to the
expansion below) — **no double counting by construction**, because the KS of
arXiv 2604.25199 is a diagrammatically defined resummation, not an empirical
functional.

Companion documents: `AUDIT-EFT_for_pseudopotential.md` (the 1s-channel
machinery, verified), `MEMO-zcore-normalization-issue.md` (the volume fix).

---

## Step 0 — Hamiltonian and site partition (exact)

Full action for Li metal (bcc, cell volume v = a³/2 = 145.8 a₀³):

    H = Σ_i [p_i²/2 − Σ_R Z/|r_i − R|] + Σ_{i<j} 1/r_ij + H_ion-ion.

Introduce a per-site orthonormalized atomic basis {u_1s^R, u_2s^R} ⊕
delocalized remainder (Löwdin; 1s–1s overlap corrections are O(e^{−Za}),
negligible; 2s–2s overlaps are absorbed into the hopping below). Split:

    H = Σ_R H_at^R + T + V_inter                                   (0.1)

- H_at^R: the complete Li atom at site R — its own ionic attraction and **all
  Coulomb matrix elements with all four indices on R** (so the full on-site
  U-physics of both shells lives in the reference).
- T: inter-site one-body terms (hopping + tails of neighboring ionic
  potentials).
- V_inter: Coulomb matrix elements with indices on ≥2 sites. By construction
  V_inter has **no on-site component**.

Status: exact rewriting.

## Step 1 — Reference: exact isolated atom (exact, ΔSCF inputs)

Atomic eigenstates define the reference correlation functions. The 2s-channel
one-body function follows from the Lehmann representation on the atomic Fock
space. Acting on the doublet ground state |1 1̄ 2σ⟩:

    c_{2σ}|11̄2σ⟩ = |11̄⟩            (removal cost  I = 5.392 eV)
    c†_{2σ̄}|11̄2σ⟩ = |11̄22̄⟩         (addition cost −A = −0.618 eV)

Paramagnetic average over the doublet (n_σ̄ = 1/2):

    g_2s(iω) = (1/2)/(iω + I) + (1/2)/(iω + A) + O(1% satellites)   (1.1)

The pole separation is the local charge-fluctuation energy

    U ≡ E(Li⁺) + E(Li⁻) − 2E(Li) = I − A = 4.774 eV                 (1.2)

— a ΔSCF definition, exactly parallel to the note's ΔE_c = |ε|−J for the
core. Equivalently Σ_at,2s(iω) = U/2 + (U/2)²/(iω − μ̃), μ̃ = −(I+A)/2.
The ~1% satellites of (1.1) (removal leaving core-excited 1s2s
configurations) are precisely the note's atomic Δ = 0.01 (its Eqs. 52/72) and
are carried by the 1s channel below. The atomic 4-point vertex γ⁽⁴⁾ is exactly
computable on the same Fock space; its static value in the 2s charge channel
is U.

Status: exact (inputs are atomic ΔSCF data).

## Step 2 — Exact dual-fermion integration (identity)

Hubbard–Stratonovich on the one-body difference T and exact integration of
Π_R H_at^R (the note's Eqs. 5–16, verified in the audit) gives the dual
action; V_inter rides along unchanged (it couples dual densities; no on-site
part). The exact lattice relation is

    Σ_latt(k,iω) = Σ_at(iω) + Σ_d(k,iω) [1 + g_at(iω) Σ_d(k,iω)]⁻¹  (2.1)

with Σ_d built from the atomic vertices γ⁽⁴⁾, γ⁽⁶⁾, … and dual propagators.
Setting Σ_d = 0 reproduces Hubbard-I — which at half filling predicts a Mott
gap of U = 4.8 eV (E±(k) = [ε_k ± √(ε_k²+U²)]/2), i.e. an insulator. Li is a
metal: the expansion **around the atomic limit** is the divergent
organization at U ≈ W. This forces the reorganization of Step 3 (resum
itinerancy first), and it is the sharp statement of *why the 2s cannot be a
static pseudopotential*: unlike the 1s (poles 70 eV away, V/ΔE ≈ 0.1), the 2s
fluctuation poles sit inside the band — no scale separation.

Status: identity; the observation about organization is exact (Li is a metal).

## Step 3 — KS resummation of the residual gas (the paper's one approved step)

Apply the EFT of arXiv 2604.25199 to the residual theory: resum the
LDA class of self-energy diagrams of **V_inter** (plus the static parts of the
local physics) into a static one-particle Hamiltonian:

    h_KS = ε_k(T) + Σ_at(static) + δV_pp^{1s,static}/v + v_LDA[n]    (3.1)

Bookkeeping at this point — diagrams consumed so far, each exactly once:
(i) all diagrams internal to the atom (inside g_at, γ);
(ii) the static/Hartree part of the on-site U (the U·n level shift — inside
     Σ_at(static), i.e. the band-center position);
(iii) the LDA class of V_inter (the KS potential; static by construction);
(iv) the static 1s pseudopotential with the volume-corrected projection.
Output for Li: a nearly-free-electron half-filled band. Model it by the KS
parabolic DOS ρ(ε) = (3/2)√ε/E_c^{3/2}, per-spin norm 1, E_F = 50.11/r_s² =
4.743 eV (r_s = 3.25), half filling ⇒ E_c = E_F/2^{−2/3} = 7.53 eV.

Status: the single approved approximation (KS); everything absorbed here is
static, hence disjoint from the dynamic diagrams below.

## Step 4 — Fluctuation expansion of the ω-dependence (second order, systematic)

What is left in Σ_latt(k,iω) − (static parts) at the Fermi surface, ordered
by diagram class:

**(a) 1s core channel** (the note's δV_pp, ω-linear part, volume-corrected):

    λ_1s = |f(0)|²-weighted sum / (v ΔE_c²) = 0.003                  (4.1)

(All machinery verified in the audit; f(0) = −1.88207, ΔSCF ΔE_c.)

**(b) 2s charge-fluctuation channel.** The fluctuation vertex is
U δn̂_↑ δn̂_↓ (first order vanishes: ⟨δn⟩ = 0, the mean field already sits in
h_KS — this is where "no double counting" is *visible*: the static moment of
this vertex is zero by construction). Leading (second-order) self-energy with
KS propagators:

    Σ_2s(iω) = U² ∫dε₁dε₂dε₃ ρ(ε₁)ρ(ε₂)ρ(ε₃)
               × [f₁f₂(1−f₃) + (1−f₁)(1−f₂)f₃] / (iω − ε₁ − ε₂ + ε₃) (4.2)

    λ_2s = −∂Σ_2s/∂(iω)|₀
         = U² [ ∫_{occ,occ,unocc} ρρρ/s² + ∫_{unocc,unocc,occ} ρρρ/s² ],
      s = ε₁+ε₂−ε₃ (measured from μ).                                (4.3)

No volume factor lurks here: with the local DOS normalized to one state per
site per spin, (4.3) is manifestly dimensionless (the 1/N's of the Bloch sums
cancel identically — checked).

Numerics (three independent methods — corrected-grid convolution, direct 3D
quadrature, median-of-means MC — agree to 3 digits):

    Li:  U = 4.774 eV, E_F = 4.743 eV → λ_2s = 0.474
         (kinematic split: occ-occ-unocc 0.226, unocc-unocc-occ 0.248)
    Na:  U = 4.591, E_F = 3.245 → λ_2s = 0.936
    K :  U = 3.840, E_F = 2.122 → λ_2s = 1.532

**(c) Not included at this order** (each a well-defined higher class, all
within the same EFT): dynamic inter-site diagrams (second order in V_inter —
the long-wavelength/collective class; requires the paper's EFT at its next
order because of the textbook infrared organization of the bare long-range
series), mixed U×V_inter diagrams (lattice screening of U), the full
frequency structure of γ⁽⁴⁾ beyond its static limit (exactly computable on
the atomic Fock space — a refinement, not a new approximation), k-resolved
vertex form factors, real KS bands instead of the parabolic model.

Status: systematic truncation; the expansion parameter is λ_2s itself.

## Step 5 — Assembly and error budget

    z⁻¹ = 1 + λ_1s + λ_2s

| element | λ_1s | λ_2s | **z (this work)** | published | expansion health |
|---|---|---|---|---|---|
| Li | 0.003 | 0.474 | **0.68** | 0.75 | λ = 0.47, marginal-controlled |
| Na | 0.004 | 0.936 | 0.52 | — | λ ≈ 0.9, critical |
| K  | 0.006 | 1.532 | 0.39 | 0.66 | λ > 1, truncation broken |

Error budget for Li: third-order fluctuation ~ λ² ≈ 0.2 → ±0.1 on z (the
dominant term, removable by computing the next order or by ED/QMC of the
downfolded model as an exact benchmark); static-γ⁽⁴⁾ limit; parabolic-DOS
model; omitted inter-site dynamic class (a candidate home for the remaining
0.68 → 0.75 residual).

## Conclusions

1. With the correct ordering the theory is **asymptotically exact**: Steps
   0–2 are identities, Step 3 is the single approved approximation, Step 4 is
   a systematic series whose error is measured by its own leading term. There
   is no undefined seam.
2. **z ≈ 0.7 for Li is a valence number, not a core number**: integrating out
   1s alone gives 0.997; adding the 2s charge-fluctuation channel gives 0.68
   with zero adjustable inputs (U and ΔE_c both ΔSCF; bands from KS). The
   published 0.75 is approximately *derived* rather than matched.
3. The structural hierarchy for the framework: **1s → static δV_pp** (scale
   separation, V/ΔE ≈ 0.1, effectively exact); **2s local fluctuations →
   dynamic second-order self-energy** (no scale separation — Hubbard-I /
   static projection is qualitatively wrong at half filling); **residual gas
   → KS**. For heavier alkalis the 2s-channel truncation must be pushed
   beyond second order (λ_K > 1) — a computation, not a new framework.
