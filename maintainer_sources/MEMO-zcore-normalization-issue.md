# Possible normalization bug in the z_core derivation (cell-volume factor)

**From:** Kunyuan (analysis carried out with Claude during a multi-day audit of the
qp-eft-agent-benchmark probing campaign)
**Date:** 2026-06-12
**Concerns:** arXiv 2604.25199 Eq. (3)/(30)–(31) and the eft-psp pipeline
**Status:** documented discrepancy + request for checking — we may well be
missing an argument; please correct us if so.

---

## TL;DR

Following the dual-fermion derivation step by step for Li, the dynamical
core correction comes out as λ_Γ ≈ 0.003 (z ≈ 0.997), **not** λ ≈ 0.34
(z ≈ 0.75). The factor between them is exactly the unit-cell volume
v ≈ 146 a.u. Tracing through the repo history, the factor appears to be
lost at one identifiable step in `legacy/notes/old-lithium.tex`: the
single-ion momentum-space kernel V_dyn(K,K′) is used directly as the
band-structure matrix element, although the same note's appendix derives
the lattice-sum prefactor (2π)³/v·δ(k−k′). The textbook statement of the
same rule: a crystal potential's Fourier *coefficient* is the single-ion
Fourier *transform* divided by the cell volume, V_G = Ṽ(G)/v.

## Sharpest formulation (EFT_for_pseudopotential.pdf, by equation number)

The full derivation document contains all three elements in one place:
- **Eq. (76)**: the crystal pseudopotential as the per-ion kernel projected and
  summed over sites — the correct starting point.
- **Eq. (88)–(89)**: the OPW basis metric, S_GG' = **v·δ_GG'** − ψ̃₁ₛψ̃₁ₛ' — the
  cell volume v is explicitly part of the basis normalization (and the second
  term is precisely the PK orthogonality correction).
- **Eq. (93)–(94)**: the spectrum rescaling ε_K → ε_K/(1+Δ) is read directly
  off the single-ion kernel V_dyn(K,K′) **without applying the metric** — the
  generalized eigenproblem with S would divide by v.
- **Eq. (96)**: δV_pp(Γ) = −11.409 a.u. − 0.461·iω — the static piece is
  −310 eV (unphysical as a valence PSP matrix element); dividing by v = 145.9
  gives −2.1 eV (normal) and slope 0.0032, matching the atomic estimate below.
- **Eq. (A2)**: the lattice-sum prefactor (2π)³/v·δ(k−k′) is derived in the
  appendix but not threaded into (94).

## The three exhibits (all in `eft-psp/legacy/notes/old-lithium.tex`)

1. **Line ~164 — the original atomic-basis estimate (correct):**
   Δ = (11|12)²/(E₁+J₁₁)² = (7.3 eV / 71.5 eV)² = **0.01**.
   This agrees with everything below.

2. **Appendix — the lattice prefactor is derived but then not applied:**
   "The full pseudopotential sums over the contributions from each ion …
   corresponds to an overall prefactor which forces crystal momentum
   conservation: **(2π)³/v · δ(k−k′)**."

3. **The momentum-space rescaling step — where the slip happens:**
   the note computes the single-ion kernel
   V_dyn(K,K′) = (64π/α)·A_K A_K′/(iω+ΔE) (dimension Ha·Bohr³ — a
   Fourier transform, not a matrix element) and reads the spectrum
   rescaling directly off it. Its own Γ-point evaluation is the smoking
   gun:

   δV_pp^reg(Γ; 0,0) = **−11.409 a.u.** − 0.461·iω

   A static valence-pseudopotential matrix element of −310 eV is not
   physical. Divided by v = 145.9: **−2.1 eV** (a normal PSP scale) and
   slope 0.461/145.9 = **0.0032** — simultaneously consistent with the
   atomic estimate (exhibit 1) and with the independent derivations below.

## Independent first-principles cross-checks (all give the 1/v answer)

- **Exact one-body lattice theorem:** N localized levels (one per cell,
  depth ΔE) hybridizing with a band through a localized vertex f(r):
  momentum conservation couples each Bloch state to exactly one coherent
  core combination with strength f/√v ⇒ λ = f²/(v·ΔE²). Exact, no
  approximation.
- **Crystal many-body Slater–Condon** (hole channel, normal-ordered,
  linked-cluster): M_R = f·e^{iθ}/√V per site, coherent sum ⇒ same.
- **Closure sum rule:** the second-order coupling weight into the
  core-hole sector is ⟨i|V P V|i⟩, independent of how the sector is
  diagonalized — intermediate-state effects (excitonic binding, sea
  relaxation, shake-up) redistribute weight within the sector but cannot
  increase the total. This closes the "collective enhancement" loophole
  at second order.
- **Atomic limit:** the same four-term kernel applied to the Li atom's 2s
  gives z ≈ 0.994–1.0000 (reference-consistent: ⟨2s|(V_H−J)|1s⟩ ≈
  −0.01 Ha after orthogonalization; ≤0.21 Ha without), consistent with
  atomic photoionization (no ~25% deep satellites).
- **Dimensional analysis:** λ = |F|²/ΔE² carries Bohr³; only the 1/v
  version is dimensionless.

We also verified the parts of the derivation that are *right* and, in our
view, genuinely valuable: the dual-fermion construction, the four-term
effective potential, the (V_H,c − J_c) vertex with the −J_c subtraction
(we re-derived f(0) = −1.88207 for Li in closed form, matching Eq. (38)),
and the ΔSCF excitation energies.

## Spectroscopic implications of the published value (for cross-checking)

z_Γ = 0.75 with a single pole at ΔE implies, literally, ~26% of each Li
valence QP's spectral weight in satellites ~3.8 Ha (≈104 eV) below ε_F,
equivalently an effective ~45 eV Coulomb mixing element between N−1
states differing by a core↔valence substitution (atomic Auger-type scale:
≤1 eV; our coherent lattice value: ≈3.7 eV). We are not aware of such
features in Li/Na valence-band XPS/EELS, but would be glad to be
corrected.

## Two questions we could not settle, where you may know more

1. Is there a derivation (beyond the materials in the paper, supplement,
   and eft-psp notes) that justifies using the single-ion kernel without
   the 1/v — e.g., a collective/q→0 enhancement argument? Our closure
   bound appears to exclude this at second order, but we may be missing
   the intended mechanism.
2. If the 1/v is indeed missing: the cross-element success (7 elements,
   esp. K with z = 0.66) cannot be mimicked by "1/v + one global
   constant" (Ω spans 4.3× across the validated set), so the agreement
   would need a different explanation — possibly that the alkali
   bandwidth narrowing originates in valence correlations with the
   f/ΔE trend co-varying. The eDMFT attribution would then need a
   second look as well.

## Pointers

- `eft-psp/legacy/notes/old-lithium.tex` — lines ~164 (atomic estimate),
  ~368–440 (momentum-space rescaling), appendix (lattice prefactor).
- `eft-psp/lithium/freq_correction.jl` — production code; header comment
  cites the old-lithium.tex rescaling equation.
- `eft-psp/lithium/dyson_exact.jl` — rank-1 postulate
  Σ(ω) = |f⟩⟨f|/(ω+ΔE) stated in the header.
- Paper: Eq. (30)–(31) (projection step), Fig. 2 caption ("dimensionless
  ratio f_K/ΔE_c"), Table IV.

Happy to walk through any step of this in detail — and genuinely hoping
there is a counter-argument we have missed.
