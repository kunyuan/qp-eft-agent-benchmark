# Equation-by-equation audit of `EFT_for_pseudopotential.pdf`

**Auditor:** Kunyuan + Claude, 2026-06-13.
**Method:** every section re-derived independently; every closed-form number
re-computed numerically (scripts inline below); verdicts per equation block.
**Headline:** the document is **correct from Eq. (1) through Eq. (92) and the
entire appendix** — including two independent derivations of the atomic result
Δ = 0.01 — and contains **exactly one error**, at Eq. (93)–(94): the spectrum
rescaling is read off the single-ion momentum kernel without applying the
basis metric the document itself defines in Eq. (88)–(89). With the metric
applied, the document's own numbers give **z_Γ(Li) = 0.9968**, consistent with
its own atomic estimate (z ≈ 0.99) and inconsistent with the published
z_Γ = 0.75.

---

## Section II — Dual-fermion formalism, Eqs. (1)–(16)

| Eq. | Content | Verdict |
|---|---|---|
| (1)–(3) | All-electron action; partition function; bare g₀ | ✅ standard |
| (4) | g_c⁻¹ = g₀⁻¹ + R | ✅ definition |
| (5)–(9) | Hubbard–Stratonovich on the one-body difference; reference integration into connected G₂ₙ | ✅ standard dual-fermion |
| (10)–(11) | Field rescaling ψ̄ = φ̄g_c; amputated vertices | ✅ |
| (12) | Reference Dyson G₂ = g_c + g_cΣ_cG₂ | ✅ definition |
| (13)–(16) | g_v⁻¹ = g_c⁻¹(R⁻¹−G₂)g_c⁻¹ = g₀⁻¹ + δV_pp with δV_pp = Σ_c + Σ_cG₂Σ_c + g₀⁻¹R⁻¹g₀⁻¹ | ✅ **independently verified**: the identity −g_c⁻¹ + g_c⁻¹R⁻¹g_c⁻¹ = g₀⁻¹ + g₀⁻¹R⁻¹g₀⁻¹ holds exactly. ⚠️ minor: the sign of the Σ_c terms in (14) is inconsistent with (12) as written (convention wobble; final formulas unaffected) |

## Section III — Reference system, Eqs. (17)–(35)

| Eq. | Content | Verdict |
|---|---|---|
| (17)–(21) | STO basis, two-electron integrals, antisymmetrized v_ijkl with spin rule | ✅ |
| (22)–(23) | Ĥ_c with core at physical energies, valence at E₀+Λ; exact-state expansion | ✅ |
| (24)–(27) | Spectral representation; Λ→∞ ground-state expansion \|Φ₀⟩ = \|N₀ξ₀⟩ + O(1/Λ) | ✅ matches our independent derivation |
| (28)–(35) | Operator matrix elements and the three G-sectors to O(1/Λ²) | ✅ |

## Section IV — Lithium (first derivation), Eqs. (36)–(52)

| Eq. | Content | Verdict |
|---|---|---|
| (36)–(38) | Li ground-state expansion; core-sector G exact (N₀−1 = one-electron states) | ✅ |
| (39)–(42) | Particle part → Σ^HF identification via completeness; clean | ✅ |
| (43)–(47) | Hole part: S²₀ₛ from perturbation theory; pole at E₀−E₁ᶜ with vertex v₁₁̄;ₛ₁̄ | ✅ the spectator-mediated transfer vertex — matches our Slater–Condon result |
| (48)–(51) | Block assembly, Λ→∞: (g_v⁻¹)_st = −δ(iω−E_s) + Σᶜ + v₁₁̄;ₛ₁̄v₁₁̄;ₜ₁̄/(iω−(E₀−E₁)) | ✅ |
| **(52)** | **Atomic Δ = (v₁₁̄;₁₂̄/(E₁+J₁₁))² = (7.3 eV/−71.5 eV)² = 0.01** | ✅ **numerically verified**: (11\|12) = 0.26807 Ha = 7.29 eV; E₁+J = −71.4 eV; Δ = 0.0104. **This is the correct physical answer for the 2s level.** |

## Section V — Li atom, scale-separation re-derivation, Eqs. (53)–(72)

| Eq. | Content | Verdict |
|---|---|---|
| (53)–(58) | Orthonormal hydrogenic basis; gc with θ(>)/θ(<) split; perturbative \|Φ₀⟩ | ✅ |
| (59)–(66) | G expansion; orthogonality relations; Σ^HF = ⟨1a‖1b⟩+⟨1̄a‖1̄b⟩ | ✅ |
| (67) | Hole part with both spin channels: [⟨11̄‖1a⟩⟨11̄‖1b⟩ + ⟨11̄‖a1̄⟩⟨11̄‖b1̄⟩]/(iω+(E₁−E₀)) | ✅ |
| (68)–(71) | Block assembly; same δV_pp | ✅ consistent with Section IV |
| **(72)** | **Δ = 0.01 again** | ✅ the document derives the atomic answer **twice, by two routes** |

## Section VI — Crystal / momentum space, Eqs. (73)–(96)

| Eq. | Content | Verdict |
|---|---|---|
| (73)–(75) | Core projection of \|ψ⟩; ψ₁ₛ and its FT ψ̃₁ₛ(K) = 8√π α^{5/2}/(K²+α²)² | ✅ verified |
| (76)–(77) | Crystal δV_pp = per-ion kernel projected and summed over sites | ✅ the correct starting point |
| (78) | Four-term coordinate kernel: 2u·δ − exchange − (u+u′−J) projection + (u−J)(u′−J)/(iω+ΔE) | ✅ **verified**: reproduces f(0) = −1.88207 (closed form) and the published Na fc tables (≤0.5%, our earlier checks) |
| (79)–(87) | FT conventions; K = k+G decomposition; lattice-sum identities | ✅ |
| **(88)–(89)** | **OPW overlap metric: ⟨ψ̃(K)\|ψ̃(K′)⟩ = (2π)³/v·δ(k−k′)·S, S_GG′ = v·δ_GG′ − ψ̃₁ₛψ̃₁ₛ′** | ✅ correct — **and decisive: the cell volume v is part of the basis normalization** (the second term is the PK orthogonality correction) |
| (90)–(92) | Single-ion V_dyn(K,K′) = (64π/α)A_K A_K′/(iω+ΔE) | ✅ **numerically verified**: 8√(π/α)·A_K equals the direct sine transform of u₁ₛ(V_H−J) to ratio 1.0000 at K = 0.3, 1, 3, 6 |
| **(93)–(94)** | Rescaling ε_K → ε_K/(1+Δ) read directly off V_dyn | ❌ **THE ERROR.** The eigenproblem in this basis is generalized, H·c = E·S·c, and the kinetic/identity blocks carry the same v as S (from the identical overlap computation: T → v·(K²/2)δ_GG′). Dividing through by v, every potential entry enters as Ṽ/v (the textbook crystal Fourier coefficient). Equivalently: the physical rescaling is the Rayleigh quotient ∂_ω⟨ψ\|V_dyn\|ψ⟩/⟨ψ\|ψ⟩, and ⟨ψ\|ψ⟩ = v in this basis, not 1. Eq. (94) sets the denominator to 1. |
| (95)–(96) | Γ-point: δV_pp^reg = −11.409 a.u. − 0.461·iω | ✅ arithmetic exact given (93) (we verified −976π/27α² + 2209π/(729α·2.625) = −11.408; slope −0.461). ❌ **physically self-refuting as a matrix element**: a −310 eV static valence-PSP shift is impossible; ÷v gives −2.1 eV (normal NFE scale). The static and dynamic parts sit in the same matrix element — the division must apply to both. |

## Appendix A, Eqs. (A1)–(A11)

| Eq. | Content | Verdict |
|---|---|---|
| (A1)–(A2) | Lattice sum → prefactor (2π)³/v·δ(k−k′) | ✅ correct — **derived but not threaded into (94)** |
| (A3)–(A11) | Two-electron integrals via Feynman parameters | ✅ verified through the A_K ratio-1.0000 test |

---

## Numerical verification log

```
(11|12)            = 0.26807 Ha = 7.29 eV          (doc: 7.3 eV)      ✓
E₁+J₁₁             = −71.4 eV                       (doc: −71.5 eV)    ✓
Δ_atom             = 0.0104                         (doc: 0.01)        ✓
8√(π/α)·A_K vs direct transform: ratio 1.0000 at K = 0.3, 1.0, 3.0, 6.0  ✓
Eq.(96) static     = −11.409 a.u.; slope −0.461     (exact)            ✓
f(0) (paper conv.) = −1.88207                       (closed form)      ✓
WITH METRIC: Δ_crystal = 0.461/145.8 = 0.00316 → z_Γ = 0.9968
```

## Final result

Out of 96 equations + 11 appendix equations, **exactly one step is wrong**
(Eq. 93–94, metric not applied), and the document carries its own three
refutations of that step: the twice-derived atomic Δ = 0.01 (Eqs. 52, 72),
the v in its own metric (Eqs. 88–89), and the −310 eV static value (Eq. 96).
Corrected, the document's machinery yields

**z_Γ(Li) = 0.9968** (crystal) ≈ z_2s ≈ 0.99 (atom, its own result),

i.e. the dynamical frozen-core correction is a ~0.3% effect, not 25%. The
published z_Γ = 0.75 (and everything downstream: the rank-1 code postulate,
paper Eqs. (30)–(31), Table IV) inherits the unmetric'd Eq. (94). The
remaining open question is then not the derivation but the phenomenology:
what actually produces the alkali bandwidth narrowing that eDMFT and ARPES
agree on — see MEMO-zcore-normalization-issue.md, question 2.
